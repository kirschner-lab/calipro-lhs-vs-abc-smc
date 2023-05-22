"""
Perelson AS, Kirschner DE, De Boer R. Dynamics of HIV infection of CD4+ T
cells. Math Biosci. 1993 Mar;114(1):81-125.
doi: 10.1016/0025-5564(93)90043-a. PMID: 8096155.
"""
# import logging
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyabc
# from pyabc import sampler
# from pyabc.visualization import plot_data_default
from scikits.odes.odeint import odeint
from scipy.spatial.distance import canberra
# from scipy import optimize
# from scipy import stats


INITIAL_VALUES = np.array([
    # Table 1, page 88.
    1e3,   # mm^{-3}, Uninfected CD4+ cells (T).
    0,     # mm^{-3}, Latently infected CD4+ cells (T_li).
    0,     # mm^{-3}, Actively infected CD4+ cells (T_ai).
    1e-3,  # mm^{-3}, HIV cells (V).
])
T_STEADY_STATE = 4  # years, Assumed time to CD4+ steady state.
# Convert T_OUT years to days because all the parameter units are in
# days.
T_OUT = np.arange(0, 10 + T_STEADY_STATE, 1e-4) * 365


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
    return df


# Sorting by time is necessary for merging columns later on.
OBS_TS = read_obs_hiv_cd4_timeseries().sort_values("day")


def concretize_ode_hiv(pars):
    # Access the parameters.
    p = SimpleNamespace(**pars)

    def ode_hiv(t, y, dy):
        """SUNDIALS CVODE equation wrapper.

        This python interface does not support passing parameters via
        void* userdata like the SUNDIALS CVODES C API does, therefore
        embed the ODE equations as a function within model_hiv to
        access the parameters varied by pyABC.

        """
        # Alias the current solution values.
        T = y[0]
        T_li = y[1]
        T_ai = y[2]
        V = y[3]
        # Rates of change RHS (equations 5a-5d, page 87).
        dT = p.s - p.mu_T*T + p.r*T*(1 - (T + T_li + T_ai)/p.T_max) - p.k_1*V*T
        dT_li = p.k_1*V*T - p.mu_T*T_li - p.k_2*T_li
        dT_ai = p.k_2*T_li - p.mu_b*T_ai
        dV = p.N*p.mu_b*T_ai - p.k_1*V*T - p.mu_V*V
        # Alias the equation LHS.
        dy[0] = dT
        dy[1] = dT_li
        dy[2] = dT_ai
        dy[3] = dV

    return ode_hiv


def model_hiv(pars):
    """ODE model of HIV.

    The ODE solves 4 equations in Perelson 1993.

    Parameters
    ----------
    pars: dict
        Dictionary of the model parameters.

    Returns
    -------
    sample: any
        The sampled timeseries of 4 variables.

    """
    ode_hiv = concretize_ode_hiv(pars)
    sim = odeint(ode_hiv, T_OUT, INITIAL_VALUES)
    return {"data": sim}


def distance(sim, obs):
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
    df_sim = (
        pd.DataFrame(data=np.c_[sim["data"].values.t.T, sim["data"].values.y],
                     columns=["day", "dT", "dT_li", "dT_ai", "V"])
        .assign(cd4_sim=lambda x: x.dT + x.dT_li + x.dT_ai)
    )
    df = (
        pd.merge_asof(OBS_TS, df_sim, on="day")
        .dropna(axis=0, how="any")
    )
    return canberra(df.cd4, df.cd4_sim)


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

# Functions in scipy.stats are available for priors in pyabc.
priors = pyabc.Distribution(
    # Table 1, page 88.
    # day^{-1}mm^{-3}, Rate of supply of CD4+ cells from precursors.
    s=pyabc.RV("gamma", a=1.985656, scale=5.681687),
    # day^{-1}, Rate of growth for the CD4+ cells.
    r=pyabc.RV("gamma", a=4.530347876, scale=0.006990707),
    # mm^{-3}, Maximum CD4+ cells.
    T_max=pyabc.RV("nbinom", n=14.0126, p=0.02432633),
    # day^{-1}, Death rate of uninfected and latently CD4+ cells.
    mu_T=pyabc.RV("gamma", a=2.10552523, scale=0.01068658),
    # day^{-1}, Death rate of actively infected CD4+ cells.
    mu_b=pyabc.RV("gamma", a=1.9856561, scale=0.1363606),
    # day^{-1}, Death rate of free virus.
    mu_V=pyabc.RV("gamma", a=1.985657, scale=1.363605),
    # mm^{3}day^{-1}, Rate constant for CD4+ becoming infected.
    k_1=pyabc.RV("gamma", a=1.985657, scale=1.363605e-5),
    # day^{-1}, Rate latently to actively infected conversion.
    k_2=pyabc.RV("gamma", a=1.594566171, scale=0.002008847),
    # Number of free virus produced by lysing a CD4+ cell.
    N=pyabc.RV("nbinom", n=13.5, p=0.01477833),
)
model = pyabc.model.FunctionModel(model_hiv, name="Perelson1993HIV")
abc = pyabc.ABCSMC(
    models=model,
    parameter_priors=priors,
    distance_function=distance,
    population_size=500,
    # sampler=sampler.SingleCoreSampler(),
)
abc.new("sqlite:///../results/perelson1993abc.db")
abc.run(
    # minimum_epsilon=0.1,
    max_nr_populations=10,
)

# plot_data_default(abc)
