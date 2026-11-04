"""Run the physical research workflow from a validated YAML/JSON configuration."""

import argparse
from dataclasses import fields
from pathlib import Path
import time

import numpy as np
import pandas as pd
import yaml

from enzymekinetics_prediction.paths import VENDOR_DIR
from .artifacts import artifact, content_key, digest_file, write_json
from .dynamics import MDSettings, run_md
from .electrostatics import pqr_observables, prepare_pqr, run_apbs
from .kinetics import CycleRates, eyring, tunneling
from .quantum import QuantumSettings, molecular_descriptors
from .structure import pocket_features, run_p2rank
from .tables import write_table

ALLOWED = {
    "schema_version",
    "protein_id",
    "protein_pdb",
    "substrate_smiles",
    "temperature_K",
    "ph",
    "ionic_strength_M",
    "conditions_assumed",
    "pocket",
    "quantum",
    "electrostatics",
    "dynamics",
    "kinetics",
}


def load_config(path):
    path = Path(path).resolve()
    config = yaml.safe_load(path.read_text())
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError("Expected schema_version: 1 configuration")
    unknown = set(config) - ALLOWED
    if unknown:
        raise ValueError(f"Unknown pipeline fields: {sorted(unknown)}")
    for field in [
        "protein_id",
        "protein_pdb",
        "substrate_smiles",
        "temperature_K",
        "ph",
    ]:
        if field not in config:
            raise ValueError(f"Missing pipeline field: {field}")
    if (
        not np.isfinite(config["temperature_K"])
        or config["temperature_K"] <= 0
        or not 0 <= config["ph"] <= 14
    ):
        raise ValueError("Invalid temperature or pH")
    pdb = Path(config["protein_pdb"]).expanduser()
    config["protein_pdb"] = str(
        (path.parent / pdb).resolve() if not pdb.is_absolute() else pdb.resolve()
    )
    if not Path(config["protein_pdb"]).is_file():
        raise FileNotFoundError(config["protein_pdb"])
    for section in ["pocket", "quantum", "electrostatics", "dynamics", "kinetics"]:
        if not isinstance(config.get(section, {}), dict):
            raise ValueError(f"{section} must be a mapping")
    for section, keys in {
        "pocket": {"enabled", "prediction_csv", "threads", "timeout", "config"},
        "quantum": {"enabled", *[f.name for f in fields(QuantumSettings)]},
        "electrostatics": {
            "enabled",
            "backend",
            "query_A",
            "dielectric",
            "timeout",
            "threads",
            "grid_points",
            "focus_length_A",
        },
        "dynamics": {
            "enabled",
            "system_xml",
            "restart_checkpoint",
            *[f.name for f in fields(MDSettings)],
        },
        "kinetics": {
            "enabled",
            "delta_g_kj_mol",
            "barrier_source",
            "kon_M_inv_s",
            "koff_s",
            "chemistry_reverse_s",
            "release_s",
            "imaginary_wavenumber_cm",
            "tunneling_method",
        },
    }.items():
        if set(config.get(section, {})) - keys:
            raise ValueError(f"Unknown {section} settings")
    for section, field in [
        ("pocket", "prediction_csv"),
        ("dynamics", "system_xml"),
        ("dynamics", "restart_checkpoint"),
    ]:
        value = config.get(section, {}).get(field)
        if value:
            p = Path(value).expanduser()
            config[section][field] = str(
                (path.parent / p).resolve() if not p.is_absolute() else p.resolve()
            )
    return config


def run_pipeline(config, output):
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(
            "Choose a new or empty output directory; results are immutable per run"
        )
    output.mkdir(parents=True, exist_ok=True)
    pdb = Path(config["protein_pdb"])
    manifest = {
        "schema_version": 1,
        "config": config,
        "pdb_sha256": digest_file(pdb),
        "run_key": content_key({"config": config, "pdb": digest_file(pdb)}),
        "stages": {},
        "scientific_status": "research; physical backend execution is not assay/mechanism validation",
    }
    write_json(output / "manifest.json", manifest)
    stages = manifest["stages"]
    uid = config["protein_id"]
    temperature, ph = config["temperature_K"], config["ph"]
    ionic = config.get("ionic_strength_M")

    def execute(name, function):
        started = time.monotonic()
        print(f"Running {name}", flush=True)
        try:
            result = function()
            write_json(output / f"{name}.json", result)
            stages[name] = {
                "status": "complete",
                "artifact": f"{name}.json",
                "seconds": time.monotonic() - started,
            }
            return result
        except Exception as exc:
            stages[name] = {
                "status": "failed",
                "error": str(exc),
                "error_type": type(exc).__name__,
                "seconds": time.monotonic() - started,
            }
            return None
        finally:
            write_json(output / "manifest.json", manifest)

    pocket = None
    settings = config.get("pocket", {})
    if settings.get("enabled", True):

        def calculate_pocket():
            path = settings.get("prediction_csv") or run_p2rank(
                pdb,
                output / "p2rank",
                VENDOR_DIR / "p2rank",
                settings.get("threads", 1),
                settings.get("timeout", 180),
                settings.get("config"),
            )
            return pocket_features(path, pdb)

        pocket = execute("pocket", calculate_pocket)
        if pocket:
            columns = {
                f"pocket_{k}": v
                for k, v in pocket["values"].items()
                if np.isscalar(v) or v is None
            }
            units = {
                f"pocket_{k}": pocket["units"][k]
                for k in pocket["values"]
                if f"pocket_{k}" in columns
            }
            write_table(
                output / "pocket.csv",
                pd.DataFrame([{"uniprot": uid, **columns}]),
                "pocket",
                ["uniprot"],
                units,
                pocket,
            )
    qconfig = config.get("quantum", {})
    if qconfig.get("enabled", True):
        quantum = execute(
            "quantum",
            lambda: molecular_descriptors(
                config["substrate_smiles"],
                QuantumSettings(**{k: v for k, v in qconfig.items() if k != "enabled"}),
                output / "cache/quantum",
            ),
        )
        if quantum:
            values = {
                f"qm_{k}": v
                for k, v in quantum["values"].items()
                if np.isscalar(v) or v is None
            }
            units = {
                f"qm_{k}": quantum["units"][k]
                for k in quantum["values"]
                if f"qm_{k}" in values
            }
            write_table(
                output / "quantum.csv",
                pd.DataFrame(
                    [{"canonical_smiles": quantum["inputs"]["smiles"], **values}]
                ),
                "qm",
                ["canonical_smiles"],
                units,
                quantum,
            )
    electro = config.get("electrostatics", {})
    if electro.get("enabled", False):

        def calculate_electrostatics():
            if ionic is None:
                raise ValueError(
                    "Explicit ionic_strength_M is required; do not infer assay ionic strength from pH"
                )
            query = electro.get("query_A")
            if query is None:
                if pocket is None or "center_A" not in pocket["values"]:
                    raise ValueError("Provide query_A or a successful pocket center")
                query = pocket["values"]["center_A"]
            pqr = prepare_pqr(
                pdb, output / "pqr/protein.pqr", ph, electro.get("timeout", 180)
            )
            backend = electro.get("backend", "apbs")
            if backend == "apbs":
                return run_apbs(
                    pqr,
                    output / "apbs",
                    query,
                    temperature_K=temperature,
                    ionic_strength_M=ionic,
                    dielectric=electro.get("dielectric", 78.5),
                    timeout=electro.get("timeout", 180),
                    threads=electro.get("threads", 1),
                    grid_points=electro.get("grid_points", 65),
                    focus_length_A=electro.get("focus_length_A"),
                )
            if backend == "screened-coulomb":
                return pqr_observables(
                    pqr,
                    query,
                    ph=ph,
                    temperature_K=temperature,
                    ionic_strength_M=ionic,
                    dielectric=electro.get("dielectric", 78.5),
                )
            raise ValueError("Choose apbs or explicitly approximated screened-coulomb")

        result = execute("electrostatics", calculate_electrostatics)
        if result:
            values = {
                "electro_potential_mean_V": float(
                    np.mean(result["values"]["potential_V"])
                ),
                "electro_field_norm_mean_V_per_A": float(
                    np.linalg.norm(result["values"]["field_V_per_A"], axis=1).mean()
                ),
                "electro_conditions_assumed": float(
                    config.get("conditions_assumed", False)
                ),
            }
            units = {
                "electro_potential_mean_V": "Volt",
                "electro_field_norm_mean_V_per_A": "Volt/Angstrom",
                "electro_conditions_assumed": "flag",
            }
            row = {
                "uniprot": uid,
                "ph": ph,
                "temperature": temperature - 273.15,
                "ionic_strength_M": ionic,
                **values,
            }
            write_table(
                output / "electrostatics.csv",
                pd.DataFrame([row]),
                "electro",
                ["uniprot", "ph", "temperature", "ionic_strength_M"],
                units,
                result,
            )
    mdconfig = config.get("dynamics", {})
    if mdconfig.get("enabled", False):

        def calculate_md():
            if ionic is None:
                raise ValueError(
                    "Explicit ionic_strength_M is required for this MD protocol"
                )
            parameters = {
                "temperature_K": temperature,
                "ph": ph,
                "ionic_strength_M": ionic,
            }
            parameters.update(
                {
                    k: v
                    for k, v in mdconfig.items()
                    if k not in {"enabled", "system_xml", "restart_checkpoint"}
                }
            )
            return run_md(
                pdb,
                output / "md",
                MDSettings(**parameters),
                mdconfig.get("system_xml"),
                restart_checkpoint=mdconfig.get("restart_checkpoint"),
            )

        result = execute("dynamics", calculate_md)
        if result:
            values = {f"md_{k}": v for k, v in result["values"].items()}
            units = {f"md_{k}": v for k, v in result["units"].items()}
            row = {
                "uniprot": uid,
                "ph": ph,
                "temperature": temperature - 273.15,
                "ionic_strength_M": ionic,
                **values,
            }
            write_table(
                output / "dynamics.csv",
                pd.DataFrame([row]),
                "md",
                ["uniprot", "ph", "temperature", "ionic_strength_M"],
                units,
                result,
            )
    kinetic = config.get("kinetics", {})
    if kinetic.get("enabled", False):

        def calculate_kinetics():
            if not kinetic.get("barrier_source"):
                raise ValueError(
                    "A justified barrier_source must accompany activation free energy"
                )
            factor = (
                tunneling(
                    kinetic["imaginary_wavenumber_cm"],
                    temperature,
                    kinetic.get("tunneling_method", "wigner"),
                )
                if "imaginary_wavenumber_cm" in kinetic
                else 1.0
            )
            chemical_rate = eyring(kinetic["delta_g_kj_mol"], temperature, factor)
            cycle = CycleRates(
                kinetic["kon_M_inv_s"],
                kinetic["koff_s"],
                chemical_rate,
                kinetic["chemistry_reverse_s"],
                kinetic["release_s"],
            )
            v = {
                **cycle.parameters(),
                "kchem_s_inv": chemical_rate,
                "transmission_factor": factor,
            }
            return artifact(
                "TST+explicit-cycle",
                {"barrier_source": kinetic["barrier_source"]},
                kinetic,
                v,
                {
                    "kcat_s_inv": "1/s",
                    "Km_M": "mol/L",
                    "kcat_over_Km_M_inv_s": "L/(mol*s)",
                    "kchem_s_inv": "1/s",
                    "transmission_factor": "dimensionless",
                },
                limitations=[
                    "Conditional on supplied activation free energy and microscopic rates; no barrier inferred from product energy or from kcat labels"
                ],
            )

        execute("kinetics", calculate_kinetics)
    manifest["status"] = (
        "failed"
        if any(s["status"] == "failed" for s in stages.values())
        else "complete"
    )
    if not stages:
        manifest["status"] = "failed"
        manifest["error"] = "No stages enabled"
    write_json(output / "manifest.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_pipeline(load_config(args.config), args.output)
    print(f"Pipeline {result['status']}: {args.output}")
    if result["status"] != "complete":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
