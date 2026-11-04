import numpy as np

R_GAS = 8.314e-3
KB = 1.380e-23
H_PLANCK = 6.626e-34

# Физические константы
DG_BARRIER = 50.0
ZPE_CORRECTION = -5.0  # Энергия нулевых колебаний (снижает барьер)


def evb_energy(x, y, temp, pH, allosteric_factor, interfacial_factor):
    pKa_triad = 6.5
    ph_factor = 1.0 / (1.0 + 10 ** (abs(pH - pKa_triad)))

    # Учитываем ZPE
    effective_barrier = (DG_BARRIER + ZPE_CORRECTION) * ph_factor * allosteric_factor

    # interfacial_factor отвечает за активацию на границе раздела (липидная капля)
    E_react = 350 * (1 - np.exp(-2 * (y - 1.4))) ** 2 + (1.5 / x**12 - 1.5 / x**6)
    E_prod = 350 * (1 - np.exp(-2 * (x - 1.4))) ** 2 + (1.5 / y**12 - 1.5 / y**6) - 40.0

    E_avg = 0.5 * (E_react + E_prod)
    E_diff = 0.5 * (E_react - E_prod)
    H_ab = effective_barrier * interfacial_factor

    return E_avg - np.sqrt(E_diff**2 + H_ab**2)


def get_time_warping(n_steps=200):
    """
    Деформирует время так, чтобы 4-я фаза (высвобождение продукта)
    шла медленнее, так как это лимитирующая стадия.
    """
    t = np.linspace(0, 1, n_steps)
    t_warped = np.zeros_like(t)
    for i, val in enumerate(t):
        if val < 0.7:
            t_warped[i] = val * (0.7 / 0.7)  # Первые фазы идут нормально
        else:
            # Высвобождение продукта растягивается
            t_warped[i] = 0.7 + ((val - 0.7) ** 2) * (0.3 / 0.09)
    return t_warped


def calculate_kinetics(temp, pH, is_tunneling, allosteric_factor, interfacial_factor):
    active_fraction = 1.0 / (1.0 + 10 ** (abs(pH - 6.5)))
    kappa = 5.0 if is_tunneling else 1.0
    eff_barrier = (DG_BARRIER + ZPE_CORRECTION) * (1.0 / interfacial_factor)

    kcat = (
        kappa
        * (KB * temp / H_PLANCK)
        * np.exp(-eff_barrier / (R_GAS * temp))
        * active_fraction
        * allosteric_factor
    )
    Km = np.exp(25.0 / (R_GAS * temp)) * 1e-3 / active_fraction
    return Km, kcat
