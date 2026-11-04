"""Label-independent features with explicit missingness and geometry provenance."""

from collections import Counter
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OneHotEncoder

AA = "ACDEFGHIKLMNPQRSTVWY"
FEATURE_VERSION = 3


@lru_cache(maxsize=40000)
def canonical_smiles(smiles):
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToSmiles(mol) if mol is not None else smiles


def enzyme_substrate_pairs(frame):
    column = "reactant_smiles" if "reactant_smiles" in frame else "substrate_smiles"
    return list(
        zip(
            frame["sequence"].astype(str),
            map(canonical_smiles, frame[column].fillna("").astype(str)),
            strict=True,
        )
    )


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def radial_quantiles(coords, count=32):
    """Rotation/translation/permutation invariant CA radial distribution, in Å."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or not len(coords):
        raise ValueError("Expected nonempty [N, 3] coordinates")
    if not np.isfinite(coords).all():
        raise ValueError("Nonfinite coordinates")
    radii = np.linalg.norm(coords - coords.mean(axis=0), axis=1)
    return np.quantile(radii, np.linspace(0, 1, count)).astype(np.float32)


def load_geometry(path):
    """Only the first PDB model; no substitution of ligand geometry for protein."""
    path = Path(path)
    if not path.is_file():
        return np.zeros(32, dtype=np.float32), False, None
    coords = []
    try:
        for line in path.read_text(errors="strict").splitlines():
            if line.startswith("ENDMDL"):
                break
            if (
                line.startswith("ATOM  ")
                and line[12:16].strip() == "CA"
                and line[16] in " A"
            ):
                coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
        if not coords:
            return np.zeros(32, dtype=np.float32), False, sha256(path)
        return radial_quantiles(coords), True, sha256(path)
    except (ValueError, UnicodeError):
        return np.zeros(32, dtype=np.float32), False, sha256(path)


@lru_cache(maxsize=40000)
def ligand_features(smiles):
    """Keep all dot-separated reactants, including water/co-substrates."""
    result = np.zeros(263, dtype=np.float32)
    # Parsing status is returned as a feature and reported in run coverage.
    # Avoid thousands of repeated RDKit messages about isolated protons.
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        result[-1] = 1  # explicit invalid/missing marker
        return result
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=256)
    result[:256] = generator.GetFingerprintAsNumPy(mol)
    result[256:262] = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
        len(Chem.GetMolFrags(mol)),
    ]
    if not np.isfinite(result).all():
        raise ValueError(f"Nonfinite ligand descriptors: {smiles}")
    return result


@lru_cache(maxsize=20000)
def sequence_features(sequence):
    counts = Counter(sequence)
    n = max(1, len(sequence))
    return np.asarray(
        [counts[a] / n for a in AA]
        + [np.log1p(len(sequence)), sum(counts[a] for a in AA) / n],
        dtype=np.float32,
    )


class FeatureBuilder(BaseEstimator, TransformerMixin):
    """Fit only EC vocabulary on training records; other features are deterministic.

    Products are deliberately excluded: the predictor requires enzyme, substrate
    and conditions, not a known reaction product or target-derived descriptors.
    """

    def __init__(self, data_dir, embeddings=None, physical_tables=None):
        self.data_dir = data_dir
        self.embeddings = embeddings
        self.physical_tables = physical_tables

    def fit(self, frame, y=None):
        self.environment_features_ = True
        self.encoder_ = OneHotEncoder(
            handle_unknown="ignore", sparse_output=False, dtype=np.float32
        )
        self.encoder_.fit(self._ec(frame))
        self.embedding_dim_ = (
            len(next(iter(self.embeddings.values()))) if self.embeddings else 0
        )
        self.geometry_ = {}
        self.geometry_provenance_ = {}
        from enzymekinetics_prediction.physics.tables import load_tables

        self.physical_lookup_ = load_tables(self.physical_tables)
        self.feature_blocks_ = {}
        offset = 0
        for name, width in [
            ("sequence", 22),
            ("substrate", 263),
            ("conditions", 4),
            ("embedding", self.embedding_dim_ + 1),
            ("geometry", 4),
            ("EC", len(self.encoder_.categories_[0])),
        ]:
            self.feature_blocks_[name] = (offset, offset + width)
            offset += width
        self.feature_blocks_["titration"] = (offset, offset + 3)
        offset += 3
        seen = set()
        for i, table in enumerate(self.physical_lookup_):
            columns = table["metadata"]["feature_columns"]
            if seen.intersection(columns):
                raise ValueError("Physical tables repeat feature columns")
            seen.update(columns)
            name = table["metadata"]["block"] + f"_{i}"
            self.feature_blocks_[name] = (offset, offset + len(columns) + 1)
            offset += len(columns) + 1
        return self

    @staticmethod
    def _ec(frame):
        return np.asarray(
            [".".join(str(ec).split(".")[:2]) for ec in frame["ec"].fillna("unknown")],
            dtype=object,
        ).reshape(-1, 1)

    def transform(self, frame):
        rows, geometry, flags = [], [], []
        for _, row in frame.iterrows():
            sequence = str(row["sequence"])
            uid = str(row.get("uniprot", "")).strip()
            smiles = row.get("reactant_smiles", row.get("substrate_smiles", ""))
            smiles = "" if pd.isna(smiles) else str(smiles)
            temperature = pd.to_numeric(row.get("temperature"), errors="coerce")
            ph = pd.to_numeric(row.get("ph"), errors="coerce")
            cond = [temperature, ph, float(pd.isna(temperature)), float(pd.isna(ph))]
            embed = None
            if self.embeddings:
                embed = self.embeddings.get(sequence)
                if embed is None:
                    embed = self.embeddings.get(uid)
                if embed is None:
                    embed = self.embeddings.get(f"AF-{uid}-F1-model_v4")
            missing_embed = embed is None
            if missing_embed:
                embed = np.zeros(self.embedding_dim_, dtype=np.float32)
            if uid not in self.geometry_:
                path = Path(self.data_dir) / "alphafold_pdb" / f"{uid}.pdb"
                q, valid, digest = load_geometry(path)
                self.geometry_[uid] = (q, valid)
                self.geometry_provenance_[uid] = {
                    "path": str(path),
                    "sha256": digest,
                    "valid": valid,
                }
            q, valid = self.geometry_[uid]
            ligand = ligand_features(smiles)
            # Fixed structural summaries are shared by all tabular comparators.
            struct = [q.mean() / 50, q.std() / 50, q[-1] / 50, float(not valid)]
            rows.append(
                np.concatenate(
                    [
                        sequence_features(sequence),
                        ligand,
                        cond,
                        embed,
                        [float(missing_embed)],
                        struct,
                    ]
                )
            )
            geometry.append(q / 50)
            flags.append([ligand[-1], float(missing_embed), float(not valid)])
        numeric = np.asarray(rows, dtype=np.float32)
        self.last_geometry_ = np.asarray(geometry, dtype=np.float32)
        self.last_flags_ = np.asarray(flags, dtype=np.float32)
        matrix = np.column_stack(
            [numeric, self.encoder_.transform(self._ec(frame))]
        ).astype(np.float32)
        if getattr(self, "environment_features_", False):
            from enzymekinetics_prediction.physics.conditions import sequence_titration

            environment = [
                sequence_titration(
                    str(row.sequence),
                    float(row.ph) if pd.notna(row.get("ph")) else float("nan"),
                )
                for _, row in frame.iterrows()
            ]
            matrix = np.column_stack(
                [matrix, np.asarray(environment, dtype=np.float32)]
            )
        self.last_physical_coverage_ = {}
        for table in getattr(self, "physical_lookup_", []):
            from enzymekinetics_prediction.physics.tables import key_for

            joins = table["metadata"]["join_on"]
            values = []
            for _, row in frame.iterrows():
                context = row.to_dict()
                if "canonical_smiles" in joins:
                    smiles = row.get("reactant_smiles", row.get("substrate_smiles", ""))
                    context["canonical_smiles"] = canonical_smiles(
                        "" if pd.isna(smiles) else str(smiles)
                    )
                key = key_for(context, joins)
                values.append(
                    table["lookup"].get(key, np.full(table["feature_count"], np.nan))
                )
            block = np.asarray(values, dtype=np.float32)
            missing = 1 - np.isfinite(block).mean(axis=1)
            self.last_physical_coverage_[table["path"]] = float(1 - missing.mean())
            matrix = np.column_stack([matrix, block, missing]).astype(np.float32)
        return matrix
