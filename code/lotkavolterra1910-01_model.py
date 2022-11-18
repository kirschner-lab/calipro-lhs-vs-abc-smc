import numpy as np
import pymc as pm
from scipy.integrate import odeint


def lotka_volterra(X, t, a, b, c, d):
    '''Return growth rate of prey and predator populations.'''
    return np.array([a*X[0] - b*X[0]*X[1],
                     -c*X[1] + d*b*X[0]*X[1]])


def simulate_lotka_volterra(rng, a, b, c, d, size=None):
    return odeint(lotka_volterra, X0, time, rtol=0.01, args=(a, b, c, d))


if __name__ == '__main__':
    # Ground truth parameters to create data before adding any noise.
    a = 1.0
    b = 0.1
    c = 1.5
    d = 0.75

    # Initial populations of prey and predator.
    X0 = [10.0, 5.0]
    # Time points.
    n_steps = 100
    time_end = 15
    time = np.linspace(0, time_end, n_steps)

    # Observations with noise.
    observed = (simulate_lotka_volterra(None, a, b, c, d) +
                np.random.normal(size=(n_steps, 2)) +
                np.random.normal(size=(n_steps, 2)))

    # Simulate from the model.
    with pm.Model() as model_lv:
        a = pm.HalfNormal('a', 1.0)
        b = pm.HalfNormal('b', 1.0)
        c = pm.ConstantData('c', 1.5)
        d = pm.ConstantData('d', 0.75)
        # Instead of specifying a likelihood function, simulate from the model.
        sim = pm.Simulator('sim', simulate_lotka_volterra,
                           params=(a, b, c, d), epsilon=10,
                           observed=observed)
        # Collect inference data.
        idata_lv = pm.sample_smc()

    idata_lv.to_netcdf("lotkavolterra1910-01_output-idata_lv.nc")
