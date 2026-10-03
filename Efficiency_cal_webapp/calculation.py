import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def calculate_efficiency(filename):

    # ============================================================
    # Excel file
    # ============================================================

    data = pd.read_excel(filename)

    E = data.iloc[:, 0].to_numpy()
    k = data.iloc[:, 1].to_numpy() * 1e-2
    N = data.iloc[:, 2].to_numpy()
    N_rel_unc = data.iloc[:, 3].to_numpy() * 1e-2

    t_live = data.iloc[:, 4].to_numpy()
    t_ref = data.iloc[:, 7].to_numpy()

    A_0 = data.iloc[:, 8].to_numpy() * 1e3
    A_0_unc = data.iloc[:, 9].to_numpy() * 1e3

    T = data.iloc[:, 10].to_numpy()

    # ============================================================
    # Activity calculation
    # ============================================================

    A_0_rel_unc = A_0_unc / A_0

    A = A_0 * 2 ** (-t_ref / T)

    A_rel_unc = A_0_rel_unc

    A_abs_unc = A_rel_unc * A

    # ============================================================
    # Efficiency
    # ============================================================

    eta = N / (A * t_live * k)

    eta_rel_unc = np.sqrt(
        N_rel_unc ** 2 +
        A_rel_unc ** 2
    )

    eta_abs_unc = eta * eta_rel_unc

    # ============================================================
    # Weighted fitting
    # ============================================================

    E_eV = E * 1e3

    logE = np.log(E_eV)
    logEta = np.log(eta)

    sigma_logEta = eta_abs_unc / eta

    deg = 3

    A_matrix = np.zeros((len(logE), deg + 1))

    for k_poly in range(deg + 1):
        A_matrix[:, deg - k_poly] = logE ** k_poly

    W = np.diag(1 / sigma_logEta)

    Aw = W @ A_matrix
    bw = W @ logEta

    Q, R = np.linalg.qr(Aw, mode="reduced")

    p = np.linalg.solve(R, Q.T @ bw)

    cov = np.linalg.inv(R) @ np.linalg.inv(R).T

    sigma_p = np.sqrt(np.diag(cov))

    # ============================================================
    # Difference between measured and fitted values
    # ============================================================

    logEta_fit_data = np.polyval(p, logE)

    eta_fit_data = np.exp(logEta_fit_data)

    abs_diff = eta - eta_fit_data

    rel_diff = 100 * abs_diff / eta_fit_data

    results = pd.DataFrame({
        "Energy [keV]": E_eV * 1e-3,
        "Measured efficiency": eta,
        "Fitted efficiency": eta_fit_data,
        "Absolute difference": abs_diff,
        "Relative difference [%]": rel_diff
    })

    # ============================================================
    # Fitted curve and uncertainty band
    # ============================================================

    E_fit = np.linspace(
        0.8 * np.min(E_eV),
        1.2 * np.max(E_eV),
        1000
    )

    logE_fit = np.log(E_fit)

    logEta_fit = np.polyval(p, logE_fit)

    sigma_logEta_fit = np.zeros_like(logE_fit)

    for i, x in enumerate(logE_fit):

        J = np.array([
            x**3,
            x**2,
            x,
            1
        ])

        sigma_logEta_fit[i] = np.sqrt(
            J @ cov @ J.T
        )

    eta_fit = np.exp(logEta_fit)

    sigma_eta_fit = eta_fit * sigma_logEta_fit

    # ============================================================
    # Plot
    # ============================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    # Uncertainty band

    lower = np.maximum(
        eta_fit - sigma_eta_fit,
        1e-10
    )

    upper = eta_fit + sigma_eta_fit

    ax.fill_between(
        E_fit,
        lower,
        upper,
        alpha=0.15,
        label="Fit uncertainty"
    )

    # Fitted curve

    ax.plot(
        E_fit,
        eta_fit,
        linewidth=1.5,
        label="Fitted curve"
    )

    # Upper/lower uncertainty lines

    ax.plot(
        E_fit,
        upper,
        "--",
        linewidth=0.8
    )

    ax.plot(
        E_fit,
        lower,
        "--",
        linewidth=0.8
    )

    # Measurement points

    ax.errorbar(
        E_eV,
        eta,
        yerr=eta_abs_unc,
        fmt="o",
        markersize=5,
        linewidth=1,
        label="Measured data"
    )

    # Logarithmic axes

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Energy [eV]")
    ax.set_ylabel("Efficiency [-]")

    ax.set_title(
        r"$\ln(\epsilon) = a + b\ln(E) + c(\ln(E))^2 + d(\ln(E))^3$"
    )

    ax.grid(True, which="both")

    ax.legend()

    fig.tight_layout()

    # Save plot

    plot_filename = "static/efficiency_plot.png"

    fig.savefig(
       plot_filename,
       dpi=200
    )

    plt.close(fig)

    # ============================================================
    # Return results
    # ============================================================

    return {
        "E": E_eV,
        "eta": eta,
        "eta_abs_unc": eta_abs_unc,
        "eta_rel_unc": eta_rel_unc,
        "p": p,
        "cov": cov,
        "sigma_p": sigma_p,
        "results": results,
        "plot_filename": plot_filename
    }
def calculate_efficiency_at_energy(E_keV, p, cov):
    """
    Hatásfok és bizonytalanság számítása adott energián.

    E_keV: energia keV-ben
    p: illesztési paraméterek [d, c, b, a]
    cov: az illesztési paraméterek kovariancia-mátrixa
    """

    # Energia eV-ban
    E_eV = E_keV * 1e3

    # log(E)
    x = np.log(E_eV)

    # Jacobian
    J = np.array([
        x**3,
        x**2,
        x,
        1
    ])

    # Illesztett log(hatásfok)
    log_eta = J @ p

    # Hatásfok
    eta = np.exp(log_eta)

    # Bizonytalanság logaritmikus térben
    sigma_log_eta = np.sqrt(
        J @ cov @ J.T
    )

    # Hatásfok abszolút bizonytalansága
    sigma_eta = eta * sigma_log_eta

    # Relatív bizonytalanság
    relative_uncertainty = 100 * sigma_eta / eta

    return {
        "energy_keV": E_keV,
        "eta": eta,
        "sigma_eta": sigma_eta,
        "relative_uncertainty": relative_uncertainty
    }