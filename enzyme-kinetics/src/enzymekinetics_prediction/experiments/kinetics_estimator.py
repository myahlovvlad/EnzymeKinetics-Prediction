# experiments/kinetics_estimator.py
import numpy as np

# Физические константы
R = 8.314e-3  # кДж/(моль*К)
T = 298.15  # К
kB = 1.380649e-23  # Дж/К
h = 6.62607015e-34  # Дж*с
NA = 6.022e23  # 1/моль


def estimate_kon(radius_enz_nm, radius_lig_nm, D_enz, D_lig, electrostatic_factor=1.0):
    """
    Оценка k_on (M^-1 s^-1) по уравнению Смолуховского.
    radius в нанометрах, D (коэфф. диффузии) в м^2/с.
    """
    # Перевод радиуса в метры
    r_total_m = (radius_enz_nm + radius_lig_nm) * 1e-9
    D_total = D_enz + D_lig

    # SI expression gives m^3/(mol*s); 1 m^3 = 1000 L.
    kon = 4 * np.pi * r_total_m * D_total * NA * 1e3

    # Электростатический фактор (белок притягивает лиганд)
    kon *= electrostatic_factor

    return kon  # Обычно ~10^7 - 10^9


def estimate_koff(delta_G_bind_kJ, kon):
    """
    Оценка k_off (s^-1) из энергии связывания.
    delta_G_bind должен быть отрицательным (например, -30 кДж/моль).
    """
    # Kd = exp(dG / RT)
    Kd = np.exp(delta_G_bind_kJ / (R * T))

    # k_off = Kd * k_on
    koff = Kd * kon
    return koff, Kd


def estimate_kcat(delta_G_dagger_kJ, temp=T, tunneling_factor=1.0):
    """
    Оценка k_cat (s^-1) через уравнение Эйринга.
    delta_G_dagger — высота барьера в кДж/моль (положительное число, например, 60).
    """
    # k_cat = (kB * T / h) * exp(-dG / RT)
    kcat = (kB * temp / h) * np.exp(-delta_G_dagger_kJ / (R * temp))

    # Учет квантового туннелирования
    kcat *= tunneling_factor

    return kcat  # Обычно ~10^0 - 10^5


def calculate_michaelis_menten_parameters():
    """
    Пример: как нейросеть KinPred предсказывает параметры и мы собираем их в Km.
    """
    # 1. Выход нейросети (Графовая архитектура):
    # Предсказала радиусы и электростатику
    predicted_radius_enz = 3.5  # нм
    predicted_radius_lig = 0.5  # нм
    D_enz = 5e-11  # м^2/с (для белка)
    D_lig = 5e-10  # м^2/с (для малой молекулы)
    electrostatics = 2.5  # Белок сильно притягивает лиганд

    kon = estimate_kon(
        predicted_radius_enz, predicted_radius_lig, D_enz, D_lig, electrostatics
    )

    # 2. Выход нейросети (Архитектура жесткости):
    # Предсказала энергию связывания
    predicted_dG_bind = -35.0  # кДж/моль

    koff, Kd = estimate_koff(predicted_dG_bind, kon)

    # 3. Выход нейросети (Path-GNN):
    # Предсказала барьер химической реакции
    predicted_dG_dagger = 55.0  # кДж/моль

    kcat = estimate_kcat(predicted_dG_dagger, tunneling_factor=3.0)

    # 4. Сборка Km по фундаментальной формуле
    Km = (koff + kcat) / kon

    print("=== Предсказанные кинетические параметры ===")
    print(f"k_on  (Ассоциация)    : {kon:.2e} M^-1 s^-1")
    print(f"k_off (Диссоциация)   : {koff:.2e} s^-1")
    print(f"k_cat (Катализ)       : {kcat:.2e} s^-1")
    print(f"Kd    (Сродство)      : {Kd:.2e} M")
    print("-" * 45)
    print(f"Km (Михаэлиса)        : {Km:.2e} M")
    print(f"Эффективность kcat/Km : {kcat / Km:.2e} M^-1 s^-1")

    # Анализ лимитирующей стадии
    if kcat < koff:
        print(
            "\n⚠️ Лимитирующая стадия: ХИМИЧЕСКАЯ РЕАКЦИЯ (kcat < koff). Km ~ Kd (истинное сродство)."
        )
    else:
        print(
            "\n⚠️ Лимитирующая стадия: ДИССОЦИАЦИЯ (koff < kcat). Km > Kd (сродство замаскировано)."
        )


if __name__ == "__main__":
    calculate_michaelis_menten_parameters()
