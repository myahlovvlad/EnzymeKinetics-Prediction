import logging
import streamlit as st
import numpy as np
import plotly.graph_objects as go
from enzymekinetics_prediction.experiments import core_chemistry
from enzymekinetics_prediction.experiments import core_structure
from enzymekinetics_prediction.experiments import core_physics

# Настройка логирования для Streamlit
logger = logging.getLogger(__name__)

st.set_page_config(layout="wide", page_title="Иллюстрация ферментативного катализа")
st.title("🧬 Иллюстрация ферментативного катализа")
st.caption(
    "Упрощённая анимация с заданными параметрами. Не выполняет DFT, QM/MM "
    "или молекулярную динамику и не использует обученную модель кинетики."
)

st.sidebar.header("Условия")
pH_val = st.sidebar.slider("pH", 4.0, 10.0, 7.0, 0.1)
temp_val = st.sidebar.slider("Температура (K)", 273.0, 340.0, 298.15, 1.0)
interfacial_act = st.sidebar.slider(
    "Активация липидной каплей (0=вода, 1=капля)", 0.0, 1.0, 0.0, 0.1
)

# --- Загрузка с try/except ---
try:
    prot_coords, prot_colors, lid_coords, lid_indices = (
        core_structure.load_protein_and_lid()
    )
    if not len(prot_coords):
        raise ValueError("Структура Q05469 не найдена или не содержит координат CA.")
    lipid_droplet = core_structure.generate_lipid_droplet(
        center=prot_coords.mean(axis=0) + [0, 0, -30]
    )
    reactions = core_chemistry.get_retrosynthetic_network(
        "CCC(=O)OCC(COC(=O)CC)OC(=O)CC", ec_class="3.1.1.1"
    )
except Exception as e:
    st.error(f"Критическая ошибка при загрузке данных: {e}")
    st.stop()

# --- Анимация ---
n_steps = 100
t_warped = core_physics.get_time_warping(n_steps)

frames = []
for i, t in enumerate(t_warped):
    display_prot = prot_coords.copy()
    lid_angle = min(1.0, max(0, (t - 0.3) / 0.3)) * np.pi / 2
    if len(lid_indices) > 0:
        for idx in lid_indices:
            display_prot[idx][0] += lid_angle * 10.0
            display_prot[idx][1] += lid_angle * 5.0

    curr_lig = np.array([[0, 0, 0]])
    if t < 0.2:
        curr_lig = np.random.normal(30, 2, (10, 3))
    elif t < 0.7:
        curr_lig = np.random.normal(0, 0.5, (10, 3))
    else:
        curr_lig = np.random.normal(10, 2, (10, 3))

    frame = go.Frame(
        data=[
            go.Scatter3d(
                x=display_prot[:, 0],
                y=display_prot[:, 1],
                z=display_prot[:, 2],
                mode="markers",
                marker=dict(size=2, color=prot_colors, opacity=0.6),
            ),
            go.Scatter3d(
                x=lipid_droplet[:, 0],
                y=lipid_droplet[:, 1],
                z=lipid_droplet[:, 2],
                mode="markers",
                marker=dict(size=2, color="yellow", opacity=0.1 * interfacial_act),
            ),
            go.Scatter3d(
                x=curr_lig[:, 0],
                y=curr_lig[:, 1],
                z=curr_lig[:, 2],
                mode="markers",
                marker=dict(size=5, color="red"),
            ),
        ],
        name=str(i),
    )
    frames.append(frame)

fig = go.Figure(frames=frames)
for trace in frames[0].data:
    fig.add_trace(trace)
fig.update_layout(
    template="plotly_dark",
    height=600,
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False),
        aspectmode="data",
    ),
    updatemenus=[
        {
            "buttons": [
                {
                    "args": [
                        None,
                        {
                            "frame": {"duration": 100, "redraw": True},
                            "fromcurrent": True,
                        },
                    ],
                    "label": "▶ Запуск",
                    "method": "animate",
                },
                {
                    "args": [
                        [None],
                        {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"},
                    ],
                    "label": "⏸ Пауза",
                    "method": "animate",
                },
            ],
            "type": "buttons",
            "x": 0.1,
            "y": 0,
        }
    ],
)
st.plotly_chart(fig, use_container_width=True)

# --- Безопасный вывод интермедиатов ---
st.subheader("♻️ Систематизация путей и интермедиатов (EC 3.1.1.1)")
if len(reactions) == 0:
    st.warning(
        "Ретросинтетический движок не нашел путей реакции для данного субстрата."
    )
else:
    cols = st.columns(len(reactions))  # Теперь безопасно, так как len > 0
    for i, rxn in enumerate(reactions):
        with cols[i]:
            st.info(
                f"**{rxn['site']}**\n\nТип: {rxn['intermediate_type']}\n\nБарьер: {50.0 + rxn['activation_energy_shift']} кДж/моль"
            )
            st.code(rxn["final_products"], language="text")
