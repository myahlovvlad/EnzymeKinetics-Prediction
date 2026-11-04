# experiments/equations.py
import sympy as sp


def get_latex_equations():
    """Возвращает словарь с LaTeX-кодом уравнений для каждого этапа."""

    # 1. Стокс-Эйнштейн (Диффузия)
    D, kB, T, eta, r = sp.symbols("D k_B T eta r")
    eq_diff = sp.Eq(D, kB * T / (6 * sp.pi * eta * r))

    # 2. Леннард-Джонс (Адсорбция)
    eps, sigma, r_ij = sp.symbols("epsilon sigma r_ij")
    eq_lj = sp.Eq(
        sp.Symbol("V_{LJ}"), 4 * eps * ((sigma / r_ij) ** 12 - (sigma / r_ij) ** 6)
    )

    # 3. Уравнение Эйринга (Теория переходного состояния, вместо простого Аррениуса)
    k_rate, kappa, kB2, T2, h, dG = sp.symbols("k kappa k_B T h Delta_G")
    eq_eyring = sp.Eq(
        k_rate, kappa * (kB2 * T2 / h) * sp.exp(-dG / (sp.Symbol("R") * T2))
    )

    # 4. Михаэлис-Ментен с ингибированием субстратом (Конкурирующий механизм)
    v0, Vmax, S, Km, Ki = sp.symbols("v_0 V_{max} [S] K_M K_i")
    eq_mm = sp.Eq(v0, Vmax * S / (Km + S + (S**2) / Ki))

    # 5. Уравнение Хендерсона-Хассельбаха (Влияние pH на протонирование триады)
    pH, pKa = sp.symbols("pH pK_a")
    eq_ph = sp.Eq(pH, pKa + sp.log(sp.Symbol("[A^-]") / sp.Symbol("[HA]")))

    return {
        "diffusion": sp.latex(eq_diff),
        "adsorption": sp.latex(eq_lj),
        "chemical": sp.latex(eq_eyring),
        "release": sp.latex(eq_mm),
        "ph": sp.latex(eq_ph),
    }
