import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger("WebDocking")

# ==========================================
# 1. Загрузка и подготовка данных
# ==========================================
logger.info("Чтение структуры белка...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "data/alphafold_pdb/Q05469.pdb")
prot_coords = []
prot_colors = []
# Цветовая схема по свойствам (pH-зависимость): красный (кислые), синий (основные), серый (гидрофобные)
acidic = ['ASP', 'GLU']
basic = ['LYS', 'ARG', 'HIS']

for model in structure:
    for chain in model:
        for res in chain:
            if "CA" in res:
                prot_coords.append(res["CA"].get_coord())
                resname = res.get_resname()
                if resname in acidic:
                    prot_colors.append('red')
                elif resname in basic:
                    prot_colors.append('blue')
                else:
                    prot_colors.append('gray')

prot_coords = np.array(prot_coords)
prot_center = np.mean(prot_coords, axis=0)

# Находим "карман" (упрощенно: берем точку чуть ниже центра масс)
pocket_center = prot_center + np.array([0, 0, -5]) 

logger.info("Генерация лиганда...")
mol = Chem.MolFromSmiles("CCC(=O)OCC(COC(=O)CC)OC(=O)CC.O")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()
lig_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
lig_coords -= np.mean(lig_coords, axis=0) # Центрируем

# ==========================================
# 2. Математическая модель (Траектория и Энергия)
# ==========================================
n_frames = 100
t_vals = np.linspace(0, 1, n_frames)

# Броуновское движение (вязкость воды)
np.random.seed(42)
brownian_noise = np.cumsum(np.random.normal(0, 0.1, (n_frames, 3)), axis=0)

# Траектория лиганда (издалека к карману)
start_pos = pocket_center + np.array([35, 35, 35])
trajectory = []
for t in t_vals:
    ease_t = 1 - (1 - t)**2.5 # Плавное замедление
    pos = (1 - ease_t) * start_pos + ease_t * pocket_center
    trajectory.append(pos)
trajectory = np.array(trajectory) + brownian_noise * (1 - t_vals[:, None]) # Шум затухает при связывании

# Модель Потенциальной Энергии (PES)
# Барьер активации (при t=0.7) и энергия связывания
energy = np.zeros(n_frames)
for i, t in enumerate(t_vals):
    if t < 0.7:
        # Приближение к ферменту (притяжение)
        energy[i] = -50 * np.sin(t * np.pi / 1.4)
    else:
        # Преодоление барьера и реакция (гидролиз)
        barrier = 20 * np.exp(-((t - 0.75) / 0.05)**2) # Гауссов барьер
        energy[i] = -50 + barrier - 40 * (t - 0.7)

# ==========================================
# 3. Построение кадров анимации
# ==========================================
frames = []
logger.info("Генерация кадров анимации (пожалуйста, подождите)...")

for i in range(n_frames):
    t = t_vals[i]
    
    # 1. Позиция лиганда и вращение
    angle = t * np.pi * 4
    c, s = np.cos(angle), np.sin(angle)
    rot_matrix = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    current_lig = lig_coords @ rot_matrix.T + trajectory[i]
    
    # 2. Растворение водной шубы (сферы)
    # В начале радиус шубы 12 Å, в конце (при связывании) 2 Å (десольватация)
    water_radius = 12 - 9 * t
    
    # Формируем 3D сферу для воды (упрощенно: облако точек)
    theta = np.linspace(0, 2*np.pi, 20)
    phi = np.linspace(0, np.pi, 20)
    x_w = water_radius * np.outer(np.cos(theta), np.sin(phi))
    y_w = water_radius * np.outer(np.sin(theta), np.sin(phi))
    z_w = water_radius * np.outer(np.ones(np.size(theta)), np.cos(phi))
    
    water_x = x_w.flatten() + trajectory[i][0]
    water_y = y_w.flatten() + trajectory[i][1]
    water_z = z_w.flatten() + trajectory[i][2]

    # 3. Химический акт (разрыв связи при t > 0.75)
    lig_alpha = 1.0
    if t > 0.75:
        # Визуально "разрываем" лиганд: отделяем часть атомов
        # В реальности здесь был бы гидролиз, мы просто раздвигаем атомы
        split_factor = (t - 0.75) * 20
        current_lig[:len(current_lig)//2, 0] -= split_factor
        current_lig[len(current_lig)//2:, 0] += split_factor
        lig_alpha = 1.0 - (t - 0.75) * 2 # Лиганд начинает исчезать (продукты уходят)

    # === Сборка кадра ===
    frame_data = [
        # Белок (Сцена 0)
        go.Scatter3d(x=prot_coords[:,0], y=prot_coords[:,1], z=prot_coords[:,2],
                     mode='markers', marker=dict(size=2, color=prot_colors, opacity=0.6),
                     name='Фермент (pH-поля)'),
        # Водная шуба (Сцена 0)
        go.Scatter3d(x=water_x, y=water_y, z=water_z,
                     mode='markers', marker=dict(size=2, color='cyan', opacity=0.1),
                     name='Гидратная шуба'),
        # Лиганд (Сцена 0)
        go.Scatter3d(x=current_lig[:,0], y=current_lig[:,1], z=current_lig[:,2],
                     mode='markers', marker=dict(size=8, color='orange', opacity=lig_alpha),
                     name='Субстрат'),
        # График PES (Сцена 1)
        go.Scatter(x=t_vals[:i+1], y=energy[:i+1], mode='lines', line=dict(color='lime', width=3),
                   name='Энергия'),
        # Точка на графике (Сцена 1)
        go.Scatter(x=[t], y=[energy[i]], mode='markers', marker=dict(size=12, color='red'),
                   name='Текущее состояние')
    ]
    
    # Заголовки с описанием физики
    annotations = [
        dict(x=0.25, y=0.95, xref="paper", yref="paper", text="<b>Сцена 1: Подход и десольватация</b><br>Броуновское движение, вязкость среды", showarrow=False, font=dict(color="white", size=14)),
        dict(x=0.75, y=0.95, xref="paper", yref="paper", text="<b>Сцена 2: Профиль энергии (PES)</b><br>Уравнение Леннард-Джонса + Барьер Аррениуса", showarrow=False, font=dict(color="white", size=14))
    ]

    frames.append(go.Frame(data=frame_data, name=str(i), layout=go.Layout(annotations=annotations)))

# ==========================================
# 4. Настройка макета и сборка HTML
# ==========================================
logger.info("Сборка интерактивного HTML-дашборда...")

# Начальный кадр
initial_data = frames[0].data

fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'xy'}]],
    subplot_titles=("3D Пространство (Вода, pH, Вязкость)", "Поверхность Потенциальной Энергии"),
    column_widths=[0.65, 0.35]
)

# Добавляем начальные данные
for trace in initial_data[:3]:
    fig.add_trace(trace, row=1, col=1)
for trace in initial_data[3:]:
    fig.add_trace(trace, row=1, col=2)

# Настройка осей и внешнего вида
fig.update_layout(
    template='plotly_dark',
    scene=dict(
        xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
        aspectmode='data'
    ),
    xaxis=dict(title='Координата реакции (Время)', color='white'),
    yaxis=dict(title='Энергия (кДж/моль)', color='white'),
    updatemenus=[{
        "buttons": [
            {"args": [None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}], "label": "Запустить симуляцию", "method": "animate"},
            {"args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}], "label": "Пауза", "method": "animate"}
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "type": "buttons",
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }]
)

# Внедряем кадры
fig.frames = frames

# Сохраняем в HTML
output_file = "docking_dashboard.html"
fig.write_html(output_file, auto_play=False)
logger.info(f"✅ ГОТОВО! Открой файл '{output_file}' в браузере.")