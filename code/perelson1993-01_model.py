import pymc as pm


class HIVModel(pm.Model):
    def __init__(self, name=''):
        '''These weak priors fitted to CD4 counts in Figure 1 of Pantaleo et
        al. 1995 by varying only a single parameter at a time and using
        parametric bootstrap.  Model was re-ran until getting a million values
        for each parameter.

        '''
        super().__init__(name)

        # Parameters - continuous
        #
        # Notation for all parameter ranges to be fitted are "median
        # (interquartile range) units", followed by the fitted distribution
        # parameters to those 3 data points of the median and
        # interquartile range.
        #
        # 10 (5 - 15) day^{-1}mm^{-3}
        #    shape    scale
        # 1.985656 5.681687
        pm.Gamma('s', alpha=1.985656, sigma=5.681687)

        # 0.03 (0.02 - 0.04) day^{-1}
        #       shape       scale
        # 4.530347876 0.006990707
        pm.Gamma('r', alpha=4.530347876, sigma=0.006990707)

        # 0.02 (0.01 - 0.03) day^{-1}
        #      shape      scale
        # 2.10552523 0.01068658
        pm.Gamma('mu_T', alpha=2.10552523, sigma=0.01068658)

        # 0.24 (0.12 - 0.36) day^{-1}
        #     shape     scale
        # 1.9856561 0.1363606
        pm.Gamma('mu_b', alpha=1.9856561, sigma=0.1363606)

        # 2.4 (1.2 - 3.6) day^{-1}
        #    shape    scale
        # 1.985657 1.363605
        pm.Gamma('mu_V', alpha=1.985657, sigma=1.363605)

        # 2.4e-5 (1.2e-5 - 3.6e-5)
        unscaled = pm.Gamma('mu_V', alpha=1.985657, sigma=1.363605)
        pm.Deterministic('k_1', unscaled ** 1e-5)

        # 3e-3 (2e-3 - 4e-3)
        #       shape       scale
        # 2.253761026 0.001426931
        pm.Gamma('k_2', alpha=2.253761026, sigma=0.001426931)

        # Parameters - discrete

        # 562 (485–716), Table 1 in doi:10.1016/j.medj.2022.06.009
        tcells = pm.NegativeBinomial('T', mu=562, alpha=14.0126)
        pm.Truncated('T_max', tcells, lower=1454)  # 99.999% of tcells.
        pm.NegativeBinomial('N', mu=900, alpha=13.5)
