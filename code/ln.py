import pandas as pd
import pyabc as abc
from pyabc.visualization import plot_data_default


def read_obs():
    """Convert spreadsheet into scipy sparse matrix."""


def model_lymph_node(pars):
    """ODE model of lung lymph nodes during tuberculosis infection.

    The ODE solves 6 equations

    Parameters
    ----------
    pars: dict
        Dictionary of the model parameters.

    Returns
    -------
    sample: any
        The sampled timeseries of 6 variables.
    """
    
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
    obs: sparse scipy array [T x V']
    """
    return


priors = abc.Distribution(
    alpha4a=abc.RV("uniform", 0.25, 0.96),
    w1=abc.RV("uniform", 0.26, 0.36),
    w2=abc.RV("uniform", 0.94, 1.89),
    k2=abc.RV("uniform", 0.83, 2.31),
    k3=abc.RV("uniform", 0.034, 0.045),
    c8=abc.RV("uniform", 16.4, 40.9),
    c9=abc.RV("uniform", 1622, 7868),
    gamma4=abc.RV("uniform", 124, 255),
    muMR=abc.RV("uniform", 0.0044, 0.0057),
    k17=abc.RV("uniform", 0.088, 0.238),
    N=abc.RV("uniform", 5, 25),
    n=abc.RV("uniform", 1.74, 2.25),
    muMI=abc.RV("uniform", 0.003, 0.0038),
    muMA=abc.RV("uniform", 0.15, 0.19),
    alpha19=abc.RV("uniform", 0.81, 1.36),
    alpha20=abc.RV("uniform", 0.23, 0.43),
    Nfrac=abc.RV("uniform", 0.00086, 0.0011),
    k15=abc.RV("uniform", 0.029, 0.109),
    k18=abc.RV("uniform", 0.00029, 0.00078),
    muBE=abc.RV("uniform", 4.1e-9, 7.4e-9),
    k5=abc.RV("uniform", 0.27, 0.90),
    eta4=abc.RV("uniform", 2.93e4, 5.23e4),
    h5=abc.RV("uniform", 2.25, 4.72),
    xi3=abc.RV("uniform", 0.99, 3.60),
    k7=abc.RV("uniform", 0.14, 0.58),
)
model = abc.model.FunctionModel(model_lymph_node, name="LymphNode")
abc = abc.ABCSMC(model, priors, distance, populationsize=1e3)
