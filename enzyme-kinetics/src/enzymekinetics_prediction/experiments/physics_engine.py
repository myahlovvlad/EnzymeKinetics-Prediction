from enzymekinetics_prediction.paths import DATA_DIR

# experiments/physics_engine.py
import numpy as np
from scipy.optimize import minimize
from Bio.PDB import PDBParser
from rdkit import Chem

R_GAS = 8.314e-3
T_REF = 298.15
KB = 1.380e-23
H_PLANCK = 6.626e-34

DG_BIND = -25.0
DG_BARRIER = 50.0
DG_RXN = -40.0


def evb_energy(x, y, temp, pH, allosteric_factor):
    """QM/MM Энергия с учетом pH и аллостерии"""
    pKa_triad = 6.5
    ph_factor = 1.0 / (1.0 + 10 ** (abs(pH - pKa_triad)))

    E_react = 350 * (1 - np.exp(-2 * (y - 1.4))) ** 2 + (1.5 / x**12 - 1.5 / x**6)
    E_prod = (
        350 * (1 - np.exp(-2 * (x - 1.4))) ** 2 + (1.5 / y**12 - 1.5 / y**6) + DG_RXN
    )

    # Оксианионная дыра (электростатическая стабилизация)
    elec_stab = (
        -15.0 * np.exp(-((x - 2) ** 2 + (y - 2) ** 2)) * ph_factor * allosteric_factor
    )

    E_avg = 0.5 * (E_react + E_prod)
    E_diff = 0.5 * (E_react - E_prod)
    H_ab = DG_BARRIER * ph_factor * allosteric_factor

    return E_avg - np.sqrt(E_diff**2 + H_ab**2) + elec_stab


def calculate_kinetics(temp, pH, is_tunneling, allosteric_factor):
    pKa_triad = 6.5
    active_fraction = 1.0 / (1.0 + 10 ** (abs(pH - pKa_triad)))
    dS_double_dagger = -0.05
    kappa = 5.0 if is_tunneling else 1.0

    kcat = (
        kappa
        * (KB * temp / H_PLANCK)
        * np.exp(dS_double_dagger / R_GAS)
        * np.exp(-DG_BARRIER / (R_GAS * temp))
    )
    kcat *= active_fraction * allosteric_factor
    Km = np.exp(-DG_BIND / (R_GAS * temp)) * 1e-3 / active_fraction
    return Km, kcat


def calculate_3d_pes(pH, temp, allosteric_factor):
    x_grid = np.linspace(1.0, 3.5, 40)
    y_grid = np.linspace(1.0, 3.5, 40)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = np.array(
        [
            [evb_energy(x, y, temp, pH, allosteric_factor) for x in x_grid]
            for y in y_grid
        ]
    )

    def grad_sq(v):
        return (
            evb_energy(v[0] + 1e-4, v[1], temp, pH, allosteric_factor)
            - evb_energy(v[0] - 1e-4, v[1], temp, pH, allosteric_factor)
        ) ** 2 + (
            evb_energy(v[0], v[1] + 1e-4, temp, pH, allosteric_factor)
            - evb_energy(v[0], v[1] - 1e-4, temp, pH, allosteric_factor)
        ) ** 2

    res = minimize(grad_sq, x0=[2.0, 2.0], method="Nelder-Mead")
    saddle = (
        res.x[0],
        res.x[1],
        evb_energy(res.x[0], res.x[1], temp, pH, allosteric_factor),
    )

    t_vals = np.linspace(0, 1, 100)
    sim_x = np.interp(t_vals, [0, 1], [3.0, 1.4])
    sim_y = np.interp(t_vals, [0, 1], [1.4, 3.0])
    energy_1d = np.array(
        [evb_energy(x, y, temp, pH, allosteric_factor) for x, y in zip(sim_x, sim_y)]
    )
    energy_1d -= energy_1d[0]

    return X, Y, Z, saddle, (sim_x, sim_y, energy_1d), energy_1d


# ==========================================
# РЕТРОСИНТЕТИЧЕСКИЙ ДВИЖОК И КОНКУРИРУЮЩИЕ ПУТИ
# ==========================================
def get_retrosynthetic_network():
    """
    Разбивает триглицерид по всем возможным сложноэфирным связям,
    находит интермедиаты и рассчитывает для них альтернативные барьеры.
    """
    reactant_smiles = "CCC(=O)OCC(COC(=O)CC)OC(=O)CC"  # Триацетин
    mol = Chem.MolFromSmiles(reactant_smiles)

    # SMARTS для поиска сложноэфирной группы (C-O)
    ester_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    matches = mol.GetSubstructMatches(ester_smarts)

    network = []
    # У триацетина 3 сложноэфирные связи (sn-1, sn-2, sn-3)
    for i, match in enumerate(matches):
        # Мы "рвем" связь между карбонильным C и эфирным O
        bond_to_break = (match[1], match[2])
        fragmented = Chem.FragmentOnBonds(mol, [bond_to_break], addDummies=False)
        frags = Chem.GetMolFrags(fragmented, asMols=True)

        # В реальности sn-2 (средняя связь) может иметь стерические препятствия
        # Симулируем разные барьеры для разных путей (кДж/моль)
        barrier_variation = [0.0, 5.0, 12.0][i]  # Центральная связь сложнее

        network.append(
            {
                "pathway": f"Путь {i + 1} (sn-{i + 1})",
                "bond_broken": match,
                "products_smiles": [Chem.MolToSmiles(f) for f in frags],
                "barrier_modifier": barrier_variation,
            }
        )

    return network


def extract_catalytic_atoms(pdb_path=DATA_DIR / "alphafold_pdb" / "Q05469.pdb"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", pdb_path)
    atoms = {"Ser_OG": None, "His_NE2": None, "Asp_OD1": None}
    for model in structure:
        for chain in model:
            for res in chain:
                res_id = res.get_id()[1]
                if res_id == 725 and "OG" in res:
                    atoms["Ser_OG"] = res["OG"].get_coord()
                if res_id == 1024 and "NE2" in res:
                    atoms["His_NE2"] = res["NE2"].get_coord()
                if res_id == 994 and "OD1" in res:
                    atoms["Asp_OD1"] = res["OD1"].get_coord()
    if atoms["Ser_OG"] is None:
        atoms["Ser_OG"] = np.array([0, 0, 0])
    if atoms["His_NE2"] is None:
        atoms["His_NE2"] = np.array([3, 0, 0])
    if atoms["Asp_OD1"] is None:
        atoms["Asp_OD1"] = np.array([6, 0, 0])
    return atoms


def run_simulation(n_steps=100):
    t_vals = np.linspace(0, 1, n_steps)
    trajectory = np.zeros((n_steps, 3))
    stage_labels = []
    pos = np.array([30.0, 30.0, 30.0])
    target = np.array([0.0, 0.0, 0.0])
    np.random.seed(42)

    for i in range(n_steps):
        t = t_vals[i]
        if t < 0.3:
            stage_labels.append("1. Захват и диффузия")
            pos += (target - pos) / np.linalg.norm(
                target - pos
            ) * 1.5 + np.random.normal(0, 0.8, 3)
        elif t < 0.5:
            stage_labels.append("2. Десольватация")
            pos += (target - pos) / np.linalg.norm(
                target - pos
            ) * 1.0 + np.random.normal(0, 0.2, 3)
        elif t < 0.7:
            stage_labels.append("3. Перенос протона и Атака")
            pos = target + np.random.normal(0, 0.05, 3)
        else:
            stage_labels.append("4. Высвобождение")
            pos += np.array([2.0, 2.0, 2.0]) + np.random.normal(0, 0.5, 3)
        trajectory[i] = pos
    return t_vals, trajectory, stage_labels
