"""
Perelson AS, Kirschner DE, De Boer R. Dynamics of HIV infection of CD4+ T
cells. Math Biosci. 1993 Mar;114(1):81-125.
doi: 10.1016/0025-5564(93)90043-a. PMID: 8096155.
"""
import numpy as np
import pandas as pd
import pymc as pm
from pymc.ode import DifferentialEquation
import pytensor.tensor as pt
from scipy.spatial.distance import canberra
import sunode
import sunode.wrappers.as_pytensor


T_STEADY_STATE = 4  # years, Assumed time to CD4+ steady state.
# Convert T_VALS years to days because all the parameter units are in
# days.
T_VALS = np.arange(0, 10 + T_STEADY_STATE, 1e-4) * 365

def read_obs_hiv_cd4_cd8():
    """Select and rename cd4 and cd8 columns from spreadsheet."""
    df = pd.read_csv("../data/pantaleo1995-table1.csv")
    df = df.loc[:, df.columns.str.endswith("mm3")]
    df = df.rename(lambda x: x.split('_')[0], axis='columns')
    # return df.to_dict("list")
    return df


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
OBS_TS = read_obs_hiv_cd4_timeseries().sort_values("day")


def ode_hiv(t, y, p):
    """SUNDIALS CVODE equation wrapper."""
    # Rates of change RHS (equations 5a-5d, page 87).
    return {
        'T': p.s - p.mu_T*y.T +
        p.r*y.T*(1 - (y.T + y.T_li + y.T_ai)/p.T_max) -
        p.k_1*y.V*y.T,
        'T_li': p.k_1*y.V*y.T - p.mu_T*y.T_li - p.k_2*y.T_li,
        'T_ai': p.k_2*y.T_li - p.mu_b*y.T_ai,
        'V': p.N*p.mu_b*y.T_ai - p.k_1*y.V*y.T - p.mu_V*y.V,
    }


def distance(obs, sim):
    """Distance from simulated to observed data.

    Dimensions:

    1. T = time
    2. V = variables
    3. V' = aggregated variables

    Parameters
    ----------
    sim: numpy array [T x V]
    obs: pandas dataframe [T x V']

    """
    import pdb; pdb.set_trace()
    df = (
        pd.merge_asof(OBS_TS, sim["data"], on="day")
        .dropna(axis=0, how="any")
    )
    return canberra(df.cd4, df.cd4_sim)


def summary_stat(df):
    """."""
    import pdb; pdb.set_trace()
    return


def tribble(columns, *data):
    """Construct a pandas DataFrame using more legible row-wise assignment."""
    # https://stackoverflow.com/a/54368508
    return pd.DataFrame(
        data=list(zip(*[iter(data)] * len(columns))),
        columns=columns,
    )


# Enable debug logging for pyabc.
# logging.basicConfig(level=logging.DEBUG)
# Prior quantiles from
priors_quantiles = tribble(
    ["param", "median", "iqr_lb", "iqr_ub", "distr"],
    "s", 10, 5, 15, "gamma",
    "r_1", 0.03, 0.02, 0.04, "gamma",
    "T_max", 562, 485, 716, "nbinom",
    "mu_T", 0.02, 0.01, 0.03, "gamma",
    "mu_B", 0.24, 0.12, 0.36, "gamma",
    "mu_V", np.NaN, np.NaN, np.NaN, "gamma",
    "k_1", 2.4e-5, 1.2e-5, 3.6e-5, "gamma"
    "k_2", 3e-3, 2e-3, 4e-3, "gamma"
    "N", 900, np.NaN, np.NaN, "nbinom",
)

# FIXME: Automate fitting these distributions to the above quantiles.
# x = []
# q = [0.25, 0.5, 0.75]

# def gamma_fit_error(**kwargs):
#     x_est = stats.gamma.ppf(q, **kwargs)
#     mat = 

# optimize.minimize(
#     gamma_fit_error,
#     constraints=(optimize.Bounds([], [])),
# )

def simulate_ode_hiv(rng, unnamed_params, size):
    """Simulator function for ABC-SMC sampler."""
    # Solve the ODEs.
    params=dict(
        zip(
            ['s', 'r', 'T_max', 'mu_T', 'mu_b', 'mu_V', 'k_1', 'k_2', 'N'],
            unnamed_params)
    )
    params = {k: (v, ()) for k, v in params.items()}
    dy, _, problem, solver, _, _ = sunode.wrappers.as_pytensor.solve_ivp(
        y0={
            # Table 1, page 88.
            'T': (1e3, ()),   # mm^{-3}, Uninfected CD4+ cells (T).
            'T_li': (0, ()),  # mm^{-3}, Latently infected CD4+ cells (T_li).
            'T_ai': (0, ()),  # mm^{-3}, Actively infected CD4+ cells (T_ai).
            'V': (1e-3, ()),  # mm^{-3}, HIV cells (V).
        },
        params=params,
        rhs=ode_hiv,
        tvals=T_VALS,
        t0=T_VALS[0],
    )
    import pdb; pdb.set_trace()
    return dy


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

    # Instead of using a likelihood function, simulate using approximate
    # Bayesian computing (ABC).
    params=[
        pt.as_tensor_variable(param)
        for param in [s, r, T_max, mu_T, mu_b, mu_V, k_1, k_2, N]
    ],
    pm.Simulator("HIVModel",
                 fn=simulate_ode_hiv,
                 params=params,
                 distance=distance,
                 observed=OBS_TS)
    idata = pm.sample_smc(draws=2, chains=1, cores=1, progressbar=False)
