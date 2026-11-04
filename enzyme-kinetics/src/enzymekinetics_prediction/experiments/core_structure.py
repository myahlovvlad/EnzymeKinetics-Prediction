from enzymekinetics_prediction.paths import DATA_DIR
import logging
import numpy as np
from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)


def load_protein_and_lid(pdb_path=DATA_DIR / "alphafold_pdb" / "Q05469.pdb"):
    logger.info(f"Загрузка PDB: {pdb_path}")
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("p", pdb_path)
    except Exception as e:
        logger.error(f"Ошибка чтения PDB: {e}", exc_info=True)
        return np.array([]), [], np.array([]), []

    prot_coords, prot_colors = [], []
    lid_coords, lid_indices = [], []
    lid_range = range(700, 720)

    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    res_id = res.get_id()[1]
                    coord = res["CA"].get_coord()

                    if res_id in [725, 994, 1024]:
                        prot_colors.append("#FFD700")
                    elif res_id in lid_range:
                        prot_colors.append("#00FFFF")
                        lid_coords.append(coord)
                        lid_indices.append(len(prot_coords))
                    else:
                        prot_colors.append("#808080")
                    prot_coords.append(coord)

    logger.info(
        f"Белок загружен. Атомов: {len(prot_coords)}. Найдено атомов крышки: {len(lid_coords)}"
    )
    return np.array(prot_coords), prot_colors, np.array(lid_coords), lid_indices


def generate_lipid_droplet(center, radius=25.0, n_points=500):
    logger.info("Генерация модели липидной капли.")
    phi = np.random.uniform(0, 2 * np.pi, n_points)
    costheta = np.random.uniform(-1, 1, n_points)
    u = np.random.uniform(0, 1, n_points)
    theta = np.arccos(costheta)
    r = radius * np.cbrt(u)
    x = r * np.sin(theta) * np.cos(phi) + center[0]
    y = r * np.sin(theta) * np.sin(phi) + center[1]
    z = r * np.cos(theta) + center[2]
    return np.column_stack((x, y, z))
