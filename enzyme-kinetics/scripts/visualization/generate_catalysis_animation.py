from enzymekinetics_prediction.paths import OUTPUT_DIR
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def main():
    print("Генерация анимации ферментативного катализа...")

    # Инициализация графического окна
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-8, 8)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("#111111")
    ax.set_facecolor("#111111")

    # Константы геометрии активного центра (Фермент)
    theta = np.linspace(0.6 * np.pi, 2.4 * np.pi, 100)
    r_base = 5.0

    # Графические элементы
    (enzyme_line,) = ax.plot(
        [], [], color="#00FFCC", lw=4, label="Enzyme (Active Site)"
    )
    (substrate_line,) = ax.plot([], [], color="#FF3366", lw=6, label="Substrate")
    (product1_line,) = ax.plot([], [], color="#FFCC00", lw=5, label="Product 1")
    (product2_line,) = ax.plot([], [], color="#00CCFF", lw=5, label="Product 2")
    tunnel_effect = ax.scatter(
        [], [], color="#FFFFFF", s=0, alpha=0.8, edgecolors="#9933FF"
    )

    title_text = ax.text(
        0, 7, "", color="white", fontsize=14, ha="center", weight="bold"
    )
    state_text = ax.text(0, -7, "", color="#AAAAAA", fontsize=11, ha="center")

    total_frames = 160

    def init():
        enzyme_line.set_data([], [])
        substrate_line.set_data([], [])
        product1_line.set_data([], [])
        product2_line.set_data([], [])
        return enzyme_line, substrate_line, product1_line, product2_line

    def animate(frame):
        # 1. Фаза приближения субстрата (0-40 кадры)
        if frame < 40:
            t = frame / 40.0
            # Фермент в покое
            r = r_base + 0.3 * np.sin(5 * theta)
            x_enz = r * np.cos(theta)
            y_enz = r * np.sin(theta)

            # Субстрат приближается по оси Y сверху вниз
            y_sub_center = 7.0 - 7.0 * t
            x_sub = np.linspace(-1.5, 1.5, 20)
            y_sub = y_sub_center + 0.3 * np.sin(2 * x_sub)

            substrate_line.set_data(x_sub, y_sub)
            product1_line.set_data([], [])
            product2_line.set_data([], [])
            tunnel_effect.set_offsets(np.empty((0, 2)))

            title_text.set_text("Stage 1: Substrate Diffusion")
            state_text.set_text(
                f"Substrate approaching active site... Distance: {7.0 - 7.0 * t:.1f} Å"
            )

        # 2. Фаза Индуцированного соответствия (40-80 кадры)
        elif frame < 80:
            t = (frame - 40) / 40.0
            # Карман фермента сжимается вокруг субстрата (Induced Fit)
            deformation = 0.8 * t
            r = r_base + 0.3 * np.sin(5 * theta)
            # Сужаем "устье" кармана
            mask = (theta > 0.6 * np.pi) & (theta < 0.9 * np.pi)
            mask2 = (theta > 2.1 * np.pi) & (theta < 2.4 * np.pi)
            r[mask] -= deformation
            r[mask2] -= deformation

            x_enz = r * np.cos(theta)
            y_enz = r * np.sin(theta)

            # Субстрат внутри кармана слегка деформируется под действием электростатики
            x_sub = np.linspace(-1.5, 1.5, 20)
            y_sub = 0.0 + 0.5 * t * np.sin(2 * x_sub)

            substrate_line.set_data(x_sub, y_sub)
            product1_line.set_data([], [])
            product2_line.set_data([], [])

            # Вспышка квантового туннелирования в конце фазы (кадры 75-80)
            if frame > 72:
                # Генерируем "искры" квантового перехода протона
                spark_x = np.random.uniform(-0.5, 0.5, 10)
                spark_y = np.random.uniform(-0.2, 0.5, 10)
                tunnel_effect.set_offsets(np.c_[spark_x, spark_y])
                tunnel_effect.set_sizes([100] * 10)
                title_text.set_text("Stage 2: Quantum Tunneling Event")
                state_text.set_text(
                    "Proton/Hydride tunneling over activation barrier (KIE physics)"
                )
            else:
                tunnel_effect.set_offsets(np.empty((0, 2)))
                title_text.set_text("Stage 2: Induced Fit & Binding")
                state_text.set_text(
                    "Enzyme complex conformation optimizing. Electrostatic scaling active."
                )

        # 3. Фаза каталитического расщепления и выхода продуктов (80-160 кадры)
        else:
            t = (frame - 80) / 80.0
            # Фермент медленно релаксирует обратно
            r = (
                r_base
                + 0.3 * np.sin(5 * theta)
                - 0.8 * (1.0 - t) * (theta > 0.6 * np.pi) * (theta < 0.9 * np.pi)
            )  # упрощенная релаксация
            x_enz = r * np.cos(theta)
            y_enz = r * np.sin(theta)

            substrate_line.set_data([], [])
            tunnel_effect.set_offsets(np.empty((0, 2)))

            # Разделение субстрата на Product 1 (влево-вверх) и Product 2 (вправо-вверх)
            x_p1 = np.linspace(-3.0 - 5.0 * t, -1.0 - 3.0 * t, 15)
            y_p1 = 0.5 + 6.0 * t + 0.2 * np.cos(3 * x_p1)

            x_p2 = np.linspace(1.0 + 3.0 * t, 3.0 + 5.0 * t, 15)
            y_p2 = 0.5 + 6.0 * t + 0.2 * np.sin(3 * x_p2)

            product1_line.set_data(x_p1, y_p1)
            product2_line.set_data(x_p2, y_p2)

            title_text.set_text("Stage 3: Product Release")
            state_text.set_text(
                "Products leaving active center. Enzyme restoring structure. kcat logic simulated."
            )

        enzyme_line.set_data(x_enz, y_enz)
        return enzyme_line, substrate_line, product1_line, product2_line, tunnel_effect

    # Создание анимации
    ani = animation.FuncAnimation(
        fig, animate, init_func=init, frames=total_frames, interval=40, blit=True
    )

    # Сохранение анимации в файл
    output_path = OUTPUT_DIR / "animations" / "enzyme_catalysis.gif"
    os.makedirs("data", exist_ok=True)

    print(f"Сохранение анимации в {output_path}... Это может занять до минуты.")
    ani.save(output_path, writer="pillow", fps=25)
    print("Анимация успешно сохранена!")


if __name__ == "__main__":
    main()
