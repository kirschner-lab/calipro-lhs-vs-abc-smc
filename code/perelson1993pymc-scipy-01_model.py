"""
Perelson AS, Kirschner DE, De Boer R. Dynamics of HIV infection of CD4+ T
cells. Math Biosci. 1993 Mar;114(1):81-125.
doi: 10.1016/0025-5564(93)90043-a. PMID: 8096155.
"""

import pickle

import numpy as np
import pandas as pd
import pymc as pm
from scipy.integrate import odeint


T_STEADY_STATE = 5  # years, Assumed time to CD4+ steady state.


def read_obs_hiv_cd4_timeseries():
    """Read the patient cd4 time series spreadsheet."""
    df = pd.read_csv("../data/pantaleo1995-figure1.csv")
    df = df.rename({"cd4_cells_per_mm3": "cd4"}, axis='columns')
    df.year += T_STEADY_STATE
    df["day"] = df.year * 365
    # Convert subject strings to float for PyMC's Distribution class that has a
    # convert_observed_data() function that assumes all data should be float.
    df["group"] = df["group"].str.extract(r"subject(\d+)")
    return df


# Sorting by time is necessary for merging columns later on.
DF_OBSERVED = read_obs_hiv_cd4_timeseries().sort_values("day")
# Collapse data into single representative sample.
DF_OBSERVED_1 = DF_OBSERVED[DF_OBSERVED["group"] == "1"]
OBSERVED = DF_OBSERVED_1.cd4.to_numpy()
# Only solve the ODE system at timepoints for which we have data.  T_VALS is in
# days because all the parameter units are in days.
T_VALS = DF_OBSERVED_1.day.to_numpy()
# Initial values.
Y0 = [1e3, 0, 0, 1e-3]


def hiv(y, t, s, r, T_max, mu_T, mu_b, mu_V, k_1, k_2, N):
    """Rates of change RHS (equations 5a-5d, page 87)."""
    return np.array([
        # T
        s - mu_T*y[0] + r*y[0]*(1 - (y[0] + y[1] + y[2])/T_max) -
        k_1*y[3]*y[0],
        # T_li
        k_1*y[3]*y[0] - mu_T*y[1] - k_2*y[1],
        # T_ai
        k_2*y[1] - mu_b*y[2],
        # V
        N*mu_b*y[2] - k_1*y[3]*y[0] - mu_V*y[3],
    ])


def simulate_hiv(rng, s, r, T_max, mu_T, mu_b, mu_V, k_1, k_2, N, size=None):
    # Pad zero value time to match the Y0 initial condition.
    times = np.concatenate((np.zeros((1,)), T_VALS))
    ret = odeint(hiv, Y0, times, rtol=0.1, mxstep=100,
                 args=(s, r, T_max, mu_T, mu_b, mu_V, k_1, k_2, N))
    # Apply the dimensional reduction here to make the model comparable to the
    # data.  Sum all the T cell counts.
    return np.sum(ret[1:, 0:3], axis=1)


if __name__ == "__main__":
    # Simulate from the model.
    with pm.Model() as model_hiv:
        # Table 1, page 88.
        # day^{-1}mm^{-3}, Rate of supply of CD4+ cells from precursors.
        s = pm.Gamma("s", alpha=1.985656, beta=5.681687)
        # day^{-1}, Rate of growth for the CD4+ cells.
        r = pm.Gamma("r", alpha=4.530347876, beta=0.006990707)
        # mm^{-3}, Maximum CD4+ cells.
        T_max = pm.NegativeBinomial("T_max", n=14.0126, p=0.02432633)
        # day^{-1}, Death rate of uninfected and latently CD4+ cells.
        mu_T = pm.Gamma("mu_T", alpha=2.10552523, beta=0.01068658)
        # day^{-1}, Death rate of actively infected CD4+ cells.
        mu_b = pm.Gamma("mu_b", alpha=1.9856561, beta=0.1363606)
        # day^{-1}, Death rate of free virus.
        mu_V = pm.Gamma("mu_V", alpha=1.985657, beta=1.363605)
        # mm^{3}day^{-1}, Rate constant for CD4+ becoming infected.
        k_1 = pm.Gamma("k_1", alpha=1.985657, beta=1.363605e-5)
        # day^{-1}, Rate latently to actively infected conversion.
        k_2 = pm.Gamma("k_2", alpha=1.594566171, beta=0.002008847)
        # Number of free virus produced by lysing a CD4+ cell.
        N = pm.NegativeBinomial("N", n=13.5, p=0.01477833)

        # Instead of specifying a likelihood function, simulate from the
        # model.
        sim = pm.Simulator("sim",
                           simulate_hiv,
                           params=(s, r, T_max, mu_T, mu_b, mu_V, k_1, k_2, N),
                           epsilon=10,
                           observed=OBSERVED)
        # Collect inference data.
        idata_hiv = pm.sample_smc(cores=8)

    with open("perelson1993-01_output-idata_hiv.pkl", "wb") as file_:
        pickle.dump(idata_hiv, file_)
