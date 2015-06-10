import numpy as np
import scipy.stats
import scipy.optimize
try:
    import matplotlib.pyplot as plt
except:
    pass
# S = SpotAnaliser.Spot_Analysis('data/3-GFP-MreB-488nm-10.tif')
# Y = S.time_pool[2][1]
# Y1 = S.Signal_Process(Y)
# X = range(len(Y))


class GaussianMixture:

    def multinomial_normal(self, x, params):
        par = np.reshape(params, (len(params) / 3, 3))
        # print x, params, par
        n = len(par)
        curve = np.zeros(len(x))
        # print curve
        for i in range(n):
            mu, sigma, k = par[i]
            curve += k * scipy.stats.norm.pdf(x, loc=mu, scale=sigma)
        return curve

    def multinomial_normal_fusion(self, x, params):
        par = np.reshape(params, (len(params) / 4, 4))
        # print x, params, par
        n = len(par)
        curve = np.zeros(len(x))
        # print curve
        for i in range(n):
            mu, sigma, k, spread = par[i]
            xi = np.searchsorted(x, mu)
            curve[0:xi] += k * scipy.stats.norm.pdf(x[0:xi], loc=mu, scale=sigma)
            xj = np.searchsorted(x, mu + spread)
            curve[xi:xj] += k * scipy.stats.norm.pdf([mu] * len(curve[xi:xj]), loc=mu, scale=sigma)
            curve[xj:] += k * scipy.stats.norm.pdf(x[xj:], loc=mu + spread, scale=sigma)
        return curve

    def residuals(self, coeffs, y, x):
        # print coeffs
        model = self.multinomial_normal(x, coeffs)
        # print model
        return y - model

    def residuals_fusion(self, coeffs, y, x):
        # print coeffs
        model = self.multinomial_normal_fusion(x, coeffs)
        # print model
        return y - model

    def Objective(self, coeffs, y, x):
        # print coeffs
        residuals = self.residuals(coeffs, y, x)
        # print model
        return ((residuals) ** 2).sum()

    def Objective_fusion(self, coeffs, y, x):
        # print coeffs
        residuals = self.residuals_fusion(coeffs, y, x)
        # print model
        return ((residuals) ** 2).sum()

    def BackFit(self, X, Y, peaks, min_component=1, fusion=True, backfit=True):
        "Takes peaks as starting points and reduce the number of parameters, compute AIC and BIC to extract the best model"
        optimized = []
        sigma0 = 3
        k0 = 1
        spread0 = 1
        # n = len(peaks)
        mus = peaks
        drop = []
        while len(drop) - len(mus) < 1:
            params = []
            bounds = []
            Map = {}
            j = 0
            for i in range(len(mus)):
                if i not in drop:
                    Map[j] = i
                    j += 1
                    params.append(mus[i])
                    bounds.append([0, None])
                    params.append(sigma0)
                    bounds.append([0, 20])
                    params.append(k0)
                    bounds.append([0.01, None])
                    if fusion:
                        params.append(spread0)
                        bounds.append([0, 20])
            # print Map
            # print "Starting parameters", params
            if not params:
                break
            # optimized.append(scipy.optimize.leastsq(self.residuals, params, args=(Y, X)))
            if fusion:
                optimized.append(scipy.optimize.minimize(self.Objective_fusion, params, args=(Y, X), method="SLSQP", bounds=bounds, options={'maxiter': 500}))
            else:
                optimized.append(scipy.optimize.minimize(self.Objective, params, args=(Y, X), method="SLSQP", bounds=bounds, options={'maxiter': 500}))
            # bic = self.__BIC(X, Y, optimized[-1].x, Lambda=4, fusion=fusion)
            # if bic < old_bic:
                # old_bic = bic
            # elif bic == old_bic:
                # break
            if fusion:
                learned = np.reshape(optimized[-1].x, (len(optimized[-1].x) / 4, 4)).T
            else:
                learned = np.reshape(optimized[-1].x, (len(optimized[-1].x) / 3, 3)).T
            drop.append(Map[np.argmin(learned[2])])
            if not backfit:
                break
            # print "Learned parameters"
            # print optimized[-1]
        return optimized

    def Plot(self, X, Y, P, fusion=True, makeplot=False, showplot=False, ax=None, color='r'):
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(X, Y)
        if fusion:
            ax.plot(X, self.multinomial_normal_fusion(X, P), '--', color=color, alpha=0.6)
        else:
            ax.plot(X, self.multinomial_normal(X, P), '--', color=color, alpha=0.6)
        if fusion:
            par = np.reshape(P, (len(P) / 4, 4))
        else:
            par = np.reshape(P, (len(P) / 3, 3))
        for i in range(len(par)):
            if fusion:
                mu, sigma, k, spread = par[i]
            else:
                mu, sigma, k = par[i]
            if fusion:
                curve = self.multinomial_normal_fusion(X, par[i])
            else:
                curve = k * scipy.stats.norm.pdf(X, loc=mu, scale=sigma)
            ax.plot(X, curve, '-.', color=color, alpha=0.4)

        if makeplot:
            if not fig:
                fig.savefig('GaussianMixture.png')
            else:
                return ax
        if showplot:
            plt.show()

    def __BIC(self, X, Y, P, Lambda=1, fusion=True):
        if fusion:
            RSS = (self.residuals_fusion(P, Y, X) ** 2).sum()
        else:
            RSS = (self.residuals(P, Y, X) ** 2).sum()
        n = len(X)
        k = len(P)
        BIC = n * np.log(RSS) + Lambda * k * np.log(n)
        return BIC

    def BIC(self, X, Y, models, Lambda=1, fusion=True):
        res = []
        for P in models:
            # params = np.reshape(P[0], (len(P[0]) / 3, 3))
            if fusion:
                RSS = (self.residuals_fusion(P[0], Y, X) ** 2).sum()
            else:
                RSS = (self.residuals(P[0], Y, X) ** 2).sum()
            n = len(X)
            k = len(P[0])
            BIC = n * np.log(RSS) + Lambda * k * np.log(n)
            res.append(BIC)
        return res

    def AIC(self, X, Y, models, Lambda=1):
        res = []
        for P in models:
            # params = np.reshape(P[0], (len(P[0]) / 3, 3))
            RSS = (self.residuals(P[0], Y, X) ** 2).sum()
            n = len(X)
            k = len(P[0])
            AIC = Lambda * 2 * k + n * np.log(RSS / n)
            res.append(AIC)
        return res
