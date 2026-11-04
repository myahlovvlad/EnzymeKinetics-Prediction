import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger("CinematicDocking")

# ==========================================
# 1. Подготовка данных
# ==========================================
logger.info("Загрузка структуры...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "data/alphafold_pdb/Q05469.pdb")

prot_coords, prot_colors = [], []
# Остатки каталитической триады (условные для Q05469)
catalytic_ids = [725, 994, 1024]

for model in structure:
    for chain in model:
        for res in chain:
            if "CA" in res:
                res_id = res.get_id()[1]
                prot_coords.append(res["CA"].get_coord())
                # Красим pH: красный (кислые), синий (основные), серый (нейтральные)
                resname = res.get_resname()
                if res_id in catalytic_ids:
                    prot_colors.append('#FFD700') # Золотой - активный центр
                elif resname in ['ASP', 'GLU']:
                    prot_colors.append('#FF4500') # Кислотные
                elif resname in ['LYS', 'ARG', 'HIS']:
                    prot_colors.append('#1E90FF') # Основные
                else:
                    prot_colors.append('#808080') # Гидрофобные

prot_coords = np.array(prot_coords)
pocket_center = np.mean(prot_coords[[i for i, c in enumerate(prot_colors) if c == '#FFD700']], axis=0) if any(c == '#FFD700' for c in prot_colors) else np.mean(prot_coords, axis=0)

logger.info("Генерация лиганда...")
mol = Chem.MolFromSmiles("CCC(=O)OCC(COC(=O)CC)OC(=O)CC.O")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()
lig_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
lig_coords -= np.mean(lig_coords, axis=0)

# ==========================================
# 2. Математика траектории и энергии
# ==========================================
n_frames = 180
t_vals = np.linspace(0, 1, n_frames)

# Профиль энергии (PES): 0 (старт) -> -20 (связывание) -> -5 (барьер) -> -40 (продукты)
energy = np.piecewise(t_vals, 
    [t_vals < 0.4, (t_vals >= 0.4) & (t_vals < 0.7), (t_vals >= 0.7) & (t_vals < 0.8), t_vals >= 0.8],
    [lambda x: -20 * (x / 0.4), 
     lambda x: -20, 
     lambda x: -20 + 15 * ((x - 0.7) / 0.1), 
     lambda x: -5 - 35 * ((x - 0.8) / 0.2)]
)

# Создаем статичную "шубу" для белка (точки вокруг белка)
protein_shell = prot_coords + np.random.normal(0, 4, prot_coords.shape)
protein_shell_colors = ['cyan'] * len(protein_shell)

# ==========================================
# 3. Генерация кадров
# ==========================================
frames = []
logger.info("Рендеринг сцен...")

for i, t in enumerate(t_vals):
    # --- Траектория лиганда ---
    start_pos = pocket_center + np.array([30, 30, 30])
    ease_t = 1 - (1 - t)**2
    curr_pos = (1 - ease_t) * start_pos + ease_t * pocket_center
    
    # Броуновское движение (затухает при связывании)
    noise = np.random.normal(0, 0.5, 3) * (1 - t) * 2 if t < 0.4 else np.zeros(3)
    curr_pos += noise
    
    # Вращение
    angle = t * np.pi * 6
    c, s = np.cos(angle), np.sin(angle)
    rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    curr_lig = lig_coords @ rot.T + curr_pos

    # --- Гидратная шуба лиганда (десольватация) ---
    w_rad = 10 * (1 - t) + 1  # Радиус уменьшается с 11 до 1
    theta, phi = np.linspace(0, 2*np.pi, 12), np.linspace(0, np.pi, 12)
    xw = w_rad * np.outer(np.cos(theta), np.sin(phi)).flatten() + curr_pos[0]
    yw = w_rad * np.outer(np.sin(theta), np.sin(phi)).flatten() + curr_pos[1]
    zw = w_rad * np.outer(np.ones(np.size(theta)), np.cos(phi)).flatten() + curr_pos[2]

    # --- Химический акт (разрыв связи) ---
    lig_color = 'orange'
    split = 0
    if t > 0.7: lig_color = '#FF4500' # Переходное состояние (красный)
    if t > 0.8:
        split = (t - 0.8) * 40
        curr_lig[:len(curr_lig)//2, 0] -= split
        curr_lig[len(curr_lig)//2:, 0] += split
        lig_color = '#00FF00' # Продукты (зеленые)

    # --- Межмолекулярные взаимодействия (пунктир) ---
    # Если лиганд близко к карману (t > 0.4), рисуем связи с остатками
    inter_x, inter_y, inter_z = [], [], []
    if 0.4 <= t < 0.8:
        # Соединяем центр лиганда с тремя каталитическими точками
        for cat_idx in [j for j, c in enumerate(prot_colors) if c == '#FFD700'][:3]:
            inter_x.extend([curr_pos[0], prot_coords[cat_idx][0], None])
            inter_y.extend([curr_pos[1], prot_coords[cat_idx][1], None])
            inter_z.extend([curr_pos[2], prot_coords[cat_idx][2], None])

    # --- Камера и Текст (Сцены) ---
    cam_eye = dict(x=2, y=2, z=1.5) # По умолчанию
    title = ""
    if t < 0.3:
        title = "<b>Сцена 1:</b> Диффузия и гидратация<br>Броуновское движение в вязкой среде"
        cam_eye = dict(x=2.5, y=2.5, z=2)
    elif t < 0.6:
        title = "<b>Сцена 2:</b> Адсорбция на поверхности<br>Десольватация лиганда, зум на карман"
        cam_eye = dict(x=1, y=1, z=0.5) # Приближаем
    elif t < 0.8:
        title = "<b>Сцена 3:</b> Индуцированное соответствие<br>Перестройка гидратной шубы, H-связи (пунктир)"
        cam_eye = dict(x=0.5, y=0.5, z=0.2)
    else:
        title = "<b>Сцена 4:</b> Катализ и высвобождение<br>Разрыв связи C-O, образование продуктов"
        cam_eye = dict(x=1.5, y=1.5, z=1.5)

    # Сборка кадра
    frame_data = [
        # 0. Белок (точки)
        go.Scatter3d(x=prot_coords[:,0], y=prot_coords[:,1], z=prot_coords[:,2],
                     mode='markers', marker=dict(size=3, color=prot_colors), name='Фермент (pH)', showlegend=False),
        # 1. Гидратная шуба белка
        go.Scatter3d(x=protein_shell[:,0], y=protein_shell[:,1], z=protein_shell[:,2],
                     mode='markers', marker=dict(size=2, color='cyan', opacity=0.15), name='Шуба фермента', showlegend=False),
        # 2. Гидратная шуба лиганда
        go.Scatter3d(x=xw, y=yw, z=zw,
                     mode='markers', marker=dict(size=2, color='blue', opacity=0.3), name='Шуба субстрата', showlegend=False),
        # 3. Лиганд
        go.Scatter3d(x=curr_lig[:,0], y=curr_lig[:,1], z=curr_lig[:,2],
                     mode='markers', marker=dict(size=7, color=lig_color), name='Субстрат'),
        # 4. Взаимодействия
        go.Scatter3d(x=inter_x, y=inter_y, z=inter_z,
                     mode='lines', line=dict(color='white', width=3, dash='dot'), name='H-связи', showlegend=False),
        
        # 5. График энергии (2D)
        go.Scatter(x=t_vals[:i+1], y=energy[:i+1], mode='lines', line=dict(color='lime', width=3), name='PES'),
        go.Scatter(x=[t], y=[energy[i]], mode='markers', marker=dict(size=10, color='red'), name='Текущее'),
        
        # 6. 2D Карта активного центра (Упрощенная проекция полей)
        go.Scatter(x=prot_coords[:,0], y=prot_coords[:,1], mode='markers', 
                   marker=dict(size=2, color=prot_colors, opacity=0.5), name='2D Проекция', showlegend=False)
    ]
    
    frame = go.Frame(
        data=frame_data, 
        name=str(i),
        layout=go.Layout(
            scene_camera=dict(eye=cam_eye),
            title_text=title
        )
    )
    frames.append(frame)

# ==========================================
# 4. Настройка макета и сборка HTML
# ==========================================
fig = make_subplots(
    rows=2, cols=2,
    specs=[[{'type': 'scene', 'rowspan': 2}, {'type': 'xy'}], [None, {'type': 'xy'}]],
    subplot_titles=("3D Молекулярный докинг", "Профиль энергии (PES)", "2D Карта полей активного центра"),
    column_widths=[0.6, 0.4], row_heights=[0.5, 0.5]
)

# Добавляем начальный кадр
for idx, trace in enumerate(frames[0].data):
    if idx < 5: fig.add_trace(trace, row=1, col=1)
    elif idx < 7: fig.add_trace(trace, row=1, col=2)
    else: fig.add_trace(trace, row=2, col=2)

fig.update_layout(
    template='plotly_dark', height=800, width=1200,
    scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False), aspectmode='data'),
    xaxis=dict(title='Координата реакции (Время)', showgrid=False),
    yaxis=dict(title='Энергия (кДж/моль)', zeroline=True, zerolinecolor='white'),
    xaxis2=dict(visible=False), yaxis2=dict(visible=False),
    updatemenus=[{
        "buttons": [
            {"args": [None, {"frame": {"duration": 80, "redraw": True}, "fromcurrent": True}], "label": "▶ Запуск", "method": "animate"},
            {"args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}], "label": "⏸ Пауза", "method": "animate"}
        ],
        "direction": "left", "pad": {"r": 10, "t": 87}, "showactive": False,
        "type": "buttons", "x": 0.05, "xanchor": "right", "y": 0, "yanchor": "top"
    }]
)

fig.frames = frames
output_file = "docking_cinematic.html"
fig.write_html(output_file, auto_play=False)
logger.info(f"✅ ГОТОВО! Открой {output_file} в браузере.")