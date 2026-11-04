import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)-7s] %(message)s',
    datefmt='%H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("DockingAnim")

# ==========================================
# 1. Загрузка координат
# ==========================================
def load_protein(pdb_path):
    logger.info(f"Шаг 1.1: Чтение PDB-файла: {pdb_path}")
    if not os.path.exists(pdb_path):
        logger.error(f"Файл не найден: {pdb_path}")
        sys.exit(1)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    prot_coords = []
    
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    prot_coords.append(res["CA"].get_coord())
                    
    prot_coords = np.array(prot_coords)
    if len(prot_coords) == 0:
        logger.error("Не найдено ни одного атома CA (альфа-углерода) в структуре!")
        sys.exit(1)
        
    prot_center = np.mean(prot_coords, axis=0)
    logger.info(f"Шаг 1.2: Белок успешно загружен. Найдено {len(prot_coords)} узлов (атомов CA).")
    logger.info(f"         Центр масс белка (X,Y,Z): [{prot_center[0]:.2f}, {prot_center[1]:.2f}, {prot_center[2]:.2f}]")
    return prot_coords, prot_center

def prepare_ligand(smiles):
    logger.info(f"Шаг 2.1: Парсинг SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error("RDKit не смог распарсить SMILES!")
        sys.exit(1)
        
    mol = Chem.AddHs(mol)
    logger.info(f"Шаг 2.2: Добавлены атомы водорода. Всего атомов: {mol.GetNumAtoms()}")
    
    logger.info("Шаг 2.3: Генерация 3D-конформера (ETKDG)...")
    rc = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if rc != 0:
        logger.warning("ETKDG не сошелся, пробуем случайные координаты...")
        rc = AllChem.EmbedMolecule(mol, randomSeed=42)
        if rc != 0:
            logger.error("Не удалось создать 3D-структуру лиганда!")
            sys.exit(1)
            
    conf = mol.GetConformer()
    lig_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    
    # Центрируем лиганд вокруг (0,0,0)
    lig_center = np.mean(lig_coords, axis=0)
    lig_coords -= lig_center
    logger.info(f"Шаг 2.4: Лиганд отцентрирован. Размер графа: {lig_coords.shape}")
    return lig_coords

# ==========================================
# 2. Настройка анимации
# ==========================================
def setup_animation(prot_coords, prot_center):
    logger.info("Шаг 3.1: Инициализация 3D-графики Matplotlib...")
    n_frames = 120  
    fig = plt.figure(figsize=(10, 8), facecolor='#111111')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#111111')
    ax.set_axis_off()
    
    # Начальная и конечная точка
    start_pos = prot_center + np.array([40.0, 40.0, 40.0])
    end_pos   = prot_center + np.array([2.0, 2.0, 2.0])
    
    # Отрисовка белка (статичная часть)
    logger.info("Шаг 3.2: Отрисовка облака атомов белка...")
    ax.scatter(prot_coords[:, 0], prot_coords[:, 1], prot_coords[:, 2], 
               c='#00FFFF', s=8, alpha=0.3, label='Фермент')
    
    # Инициализация лиганда (ВАЖНО: без запятой перед равно!)
    logger.info("Шаг 3.3: Инициализация объекта лиганда...")
    scatter_lig = ax.scatter([], [], [], c='#FFA500', s=50, depthshade=False)
    
    # HUD (Текст на экране)
    text_dist = ax.text2D(0.05, 0.95, "", transform=ax.transAxes, color='white', fontsize=12, fontweight='bold')
    text_status = ax.text2D(0.05, 0.90, "", transform=ax.transAxes, color='#FFA500', fontsize=12)
    
    # Границы осей
    margin = 60
    ax.set_xlim(prot_center[0] - margin, prot_center[0] + margin)
    ax.set_ylim(prot_center[1] - margin, prot_center[1] + margin)
    ax.set_zlim(prot_center[2] - margin, prot_center[2] + margin)
    
    return fig, ax, scatter_lig, text_dist, text_status, start_pos, end_pos, n_frames

# ==========================================
# 3. Функция обновления кадров
# ==========================================
def create_update_func(ax, scatter_lig, text_dist, text_status, lig_coords, prot_center, start_pos, end_pos, n_frames):
    
    def update(frame):
        # Логируем каждый 10-й кадр, чтобы не спамить
        if frame % 10 == 0 or frame == n_frames - 1:
            logger.info(f"  -> Рендеринг кадра {frame+1}/{n_frames} ({(frame/(n_frames-1))*100:.0f}%)")
            
        t = frame / (n_frames - 1)
        
        # Плавное замедление в конце (easing)
        ease_t = 1 - (1 - t)**3 
        
        # Интерполяция положения
        current_center = (1 - ease_t) * start_pos + ease_t * end_pos
        
        # Вращение лиганда (2 полных оборота)
        angle = t * np.pi * 4
        c, s = np.cos(angle), np.sin(angle)
        rot_matrix = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        
        # Применяем вращение и сдвиг
        rotated_lig = lig_coords @ rot_matrix.T
        final_lig_coords = rotated_lig + current_center
        
        # Обновляем данные на графике
        scatter_lig._offsets3d = (final_lig_coords[:, 0], final_lig_coords[:, 1], final_lig_coords[:, 2])
        
        # Обновляем текст
        dist = np.linalg.norm(current_center - prot_center)
        text_dist.set_text(f"Дистанция до центра: {dist:.2f} Å")
        
        if t < 0.85:
            text_status.set_text("Статус: Преодоление водной шубы...")
            text_status.set_color('#FFA500') # Оранжевый
            scatter_lig.set_color('#FFA500')
        else:
            text_status.set_text("Статус: УСПЕШНЫЙ ДОКИНГ (Посадка)")
            text_status.set_color('#00FF00') # Зеленый
            scatter_lig.set_color('#00FF00')
            
        # Вращение камеры
        ax.view_init(elev=20, azim=frame * 2.0)
        
        return scatter_lig, text_dist, text_status
        
    return update

# ==========================================
# 4. Главный пайплайн
# ==========================================
def main():
    start_time = time.time()
    logger.info("=== СТАРТ СИМУЛЯЦИИ МОЛЕКУЛЯРНОГО ДОКИНГА ===")
    
    pdb_path = "data/alphafold_pdb/Q05469.pdb"
    smiles = "CCC(=O)OCC(COC(=O)CC)OC(=O)CC.O"
    
    # 1. Подготовка данных
    prot_coords, prot_center = load_protein(pdb_path)
    lig_coords = prepare_ligand(smiles)
    
    # 2. Настройка графики
    fig, ax, scatter_lig, text_dist, text_status, start_pos, end_pos, n_frames = setup_animation(prot_coords, prot_center)
    
    # 3. Создание анимации
    update_func = create_update_func(ax, scatter_lig, text_dist, text_status, lig_coords, prot_center, start_pos, end_pos, n_frames)
    
    logger.info("Шаг 4.1: Запуск рендеринга анимации. Пожалуйста, подождите...")
    ani = FuncAnimation(fig, update_func, frames=n_frames, interval=50, blit=False)
    
    # 4. Сохранение
    output_file = "docking_simulation.gif"
    logger.info(f"Шаг 4.2: Сохранение GIF в файл '{output_file}'. Это может занять 10-30 секунд...")
    
    try:
        ani.save(output_file, writer=PillowWriter(fps=15))
        elapsed = time.time() - start_time
        logger.info(f"✅ ГОТОВО! Анимация успешно сохранена за {elapsed:.1f} секунд.")
        logger.info(f"Файл находится по пути: {os.path.abspath(output_file)}")
    except Exception as e:
        logger.error(f"Ошибка при сохранении файла: {e}")
        
if __name__ == "__main__":
    main()