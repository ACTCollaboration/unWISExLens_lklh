"""Taken from cobaya BoltzmanBase class
https://github.com/CobayaSampler/cobaya/blob/master/cobaya/theories/_cosmo/boltzmannbase.py"""

import numpy as np
from scipy.interpolate import RectBivariateSpline

class PowerSpectrumInterpolator(RectBivariateSpline):
    r"""2D spline interpolation object (scipy.interpolate.RectBivariateSpline) to evaluate matter power spectrum as function of z and k. *This class is adapted from CAMB's own P(k) interpolator, by Antony Lewis; it's mostly interface-compatible with the original.*

    :param z: values of z for which the power spectrum was evaluated.
    :param k: values of k for which the power spectrum was evaluated.
    :param P_or_logP: Values of the power spectrum (or log-values, if logP=True).
    :param logP: if True (default: False), log of power spectrum are given and used
        for the underlying interpolator.
    :param logsign: if logP is True, P_or_logP is log(logsign*Pk)
    :param extrap_kmax: if set, use power law extrapolation beyond kmax up to
        extrap_kmax; useful for tails of integrals.
    """

    def __init__(self, z, k, P_or_logP, extrap_kmax=None, extrap_kmin=None, logP=False, logsign=1):
        self.islog = logP
        #  Check order
        z, k = (np.atleast_1d(x) for x in [z, k])
        if len(z) < 4:
            raise ValueError('Require at least four redshifts for Pk interpolation.'
                             'Consider using Pk_grid if you just need a a small number'
                             'of specific redshifts (doing 1D splines in k yourself).')
        i_z = np.argsort(z)
        i_k = np.argsort(k)
        self.logsign = logsign
        self.z, self.k, P_or_logP = z[i_z], k[i_k], P_or_logP[i_z, :][:, i_k]
        self.zmin, self.zmax = self.z[0], self.z[-1]
        self.kmin, self.kmax = self.k[0], self.k[-1]
        logk = np.log(self.k)
        # Continue until extrap_kmax using a (log,log)-linear extrapolation
        if extrap_kmax and extrap_kmax > self.kmax:
            logk = np.hstack([logk, np.log(self.kmax) * 0.1 + np.log(extrap_kmax) * 0.9,np.log(extrap_kmax)])
            if logP:
                logPnew = np.empty((P_or_logP.shape[0], P_or_logP.shape[1] + 2))
                logPnew[:, :-2] = P_or_logP
                diff = (logPnew[:, -3] - logPnew[:, -4]) / (logk[-3] - logk[-4])
                if np.any(diff)<0:
                    raise ValueError("No log extrapolation possible! divergent behavior")
                delta = diff * (logk[-1] - logk[-3])
                logPnew[:, -1] = logPnew[:, -3] + delta
                logPnew[:, -2] = logPnew[:, -3] + delta * 0.9

                P_or_logP = logPnew
            else:
                Pnew = np.empty((P_or_logP.shape[0], P_or_logP.shape[1] + 2))
                Pnew[:, :-2] = P_or_logP
                if np.all(np.sign(P_or_logP[:, -5:])==1):
                    logFactor = 1
                elif np.all(np.sign(P_or_logP[:, -5:]) == -1):
                    logFactor = -1
                else:
                    raise ValueError("Can only log extrapolate if last few elements are strictly positive or negative")

                diff = (np.log(logFactor * P_or_logP[:, -3]) - np.log(logFactor * P_or_logP[:, -4])) / (logk[-3] - logk[-4])
                if np.any(diff)<0:
                    raise ValueError("No log extrapolation possible! divergent behavior")

                delta = diff * (logk[-1] - logk[-3])
                Pnew[:, -1] = logFactor * np.exp(np.log(logFactor * Pnew[:, -3]) + delta)
                Pnew[:, -2] = logFactor * np.exp(np.log(logFactor * Pnew[:, -3]) + delta * 0.9)

                P_or_logP = Pnew
            self.kmax = extrap_kmax  # Added for consistency with CAMB

        super().__init__(self.z, logk, P_or_logP)

    def P(self, z, k, grid=None):
        """
        Get the power spectrum at (z,k).
        """
        if grid is None:
            grid = not np.isscalar(z) and not np.isscalar(k)
        if self.islog:
            return self.logsign * np.exp(self(z, np.log(k), grid=grid))
        else:
            return self(z, np.log(k), grid=grid)

    def logP(self, z, k, grid=None):
        """
        Get the log power spectrum at (z,k). (or minus log power spectrum if
        islog and logsign=-1)
        """
        if grid is None:
            grid = not np.isscalar(z) and not np.isscalar(k)
        if self.islog:
            return self(z, np.log(k), grid=grid)
        else:
            return np.log(self(z, np.log(k), grid=grid))

    def piecewise_evaluate(self, z, k, pad=None, log=False):
        if pad is None:
            pad = (np.nan, np.nan)

        out = None
        if not log:
            out = self.P(z, k, grid=False)
            out[k < self.kmin] = pad[0]
            out[k > self.kmax] = pad[1]
        else:
            out = self.logP(z, k, grid=False)
            out[k < self.kmin] = pad[0]
            out[k > self.kmax] = pad[1]

        return out


class ScaledPowerSpectrumInterpolator(PowerSpectrumInterpolator):

    def __init__(self, scaling_function, pk_interp_lin, pk_interp_nonlin=None):
        self.__pk_interp_lin = pk_interp_lin
        self.__pk_interp_nonlin = pk_interp_nonlin
        self.__scaling_function = scaling_function

    def P(self, z, k, grid=None):
        if grid is None:
            grid = not np.isscalar(z) and not np.isscalar(k)
        if grid:
            scaling = self.__scaling_function(z)[...,None]
        else:
            scaling = self.__scaling_function(z)

        lin_power = self.__pk_interp_lin.P(z, k, grid=grid)
        if self.__pk_interp_nonlin is None:
            return scaling * lin_power
        else:
            return scaling * lin_power + (self.__pk_interp_nonlin.P(z, k, grid=grid) - lin_power)


    def logP(self, z, k, grid=None):
        if grid is None:
            grid = not np.isscalar(z) and not np.isscalar(k)
        if grid:
            scaling = self.__scaling_function(z)[...,None]
        else:
            scaling = self.__scaling_function(z)

        lin_power = self.__pk_interp_lin.P(z, k, grid=grid)
        if self.__pk_interp_nonlin is None:
            return np.log(scaling * lin_power)
        else:
            return np.log(scaling * lin_power + (self.__pk_interp_nonlin.P(z, k, grid=grid) - lin_power))