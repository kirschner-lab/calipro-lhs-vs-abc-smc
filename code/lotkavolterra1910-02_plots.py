import arviz as az
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


def lotka_volterra(X, t, a, b, c, d):
    '''Return growth rate of prey and predator populations.'''
    return np.array([a*X[0] - b*X[0]*X[1],
                     -c*X[1] + d*b*X[0]*X[1]])


def simulate_lotka_volterra(rng, a, b, c, d, size=None):
    return odeint(lotka_volterra, X0, time, rtol=0.01, args=(a, b, c, d))


if __name__ == '__main__':
    # Initial populations of prey and predator.
    X0 = [10.0, 5.0]
    # Time points.
    n_steps = 100
    time_end = 15
    time = np.linspace(0, time_end, n_steps)

    # Load files saved by loktavolterra1910_model.py
    idata_lv = az.from_netcdf('lotkavolterra1910-01_output-idata_lv.nc')
    observed_lv = idata_lv.observed_data['sim'].to_numpy()
    c = float(idata_lv.constant_data['c'])
    d = float(idata_lv.constant_data['d'])

    # Plot model chains, posterior, and overlay ODEs on the noisy data.
    az.summary(idata_lv)
    az.plot_trace(idata_lv, kind='rank_vlines')
    plt.show()
    az.plot_posterior(idata_lv)
    plt.show()

    posterior = idata_lv.posterior.stack(samples=('draw', 'chain'))
    _, ax = plt.subplots(figsize=(14, 6))
    ax.plot(observed_lv[:, 0], 'o', label='prey', c='C0', mec='k')
    ax.plot(observed_lv[:, 1], 'o', label='predator', c='C1', mec='k')
    # ax.plot(simulate_lotka_volterra(None,
    #                                 posterior['a'].mean(),
    #                                 posterior['b'].mean(),
    #                                 c,
    #                                 d),
    #         linewidth=3)
    ax.set_xlabel('time')
    ax.set_ylabel('population')
    ax.legend()
    for i in np.random.randint(0, n_steps, 75):
        sim = simulate_lotka_volterra(None,
                                      posterior['a'][i],
                                      posterior['b'][i],
                                      c,
                                      d)
        ax.plot(sim[:, 0], alpha=0.1, c='C0')
        ax.plot(sim[:, 1], alpha=0.1, c='C1')
    plt.savefig('/Users/pnanda/Sync/lab-reports/2022-11-17/img/pymc.pdf')
    plt.show()
