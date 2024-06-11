import numpy as np
from scipy import special
from ..auxiliary.pk_interpolator import PowerSpectrumInterpolator
from .model_helpers_unWISExLens import cosmo_from_camb
from scipy.interpolate import interp1d as interpolate


class CrossRedshiftCosmoCorrection(object):

    def __init__(self):
        pass

    def __compute_cosmo_factor(self, *args, **kwargs):
        raise NotImplementedError()

    def correction_factor(self, *args, **kwargs):
        raise NotImplementedError()

    def __call__(self, *args, **kwargs):
        raise NotImplementedError()


class CrossRedshiftCosmoCorrectionExact(CrossRedshiftCosmoCorrection):

    def __init__(self, zbins, smin=2.5, smax=10, nbins=3, log_bins=True, N_integration=1024, kmin=1e-4, kmax=10, fid_cosmo=None, fid_pk=None):
        super().__init__()
        self.__zbins = zbins
        self.__z_avg = (zbins[1:] + zbins[:-1]) / 2
        self.__dz = np.diff(zbins)

        if log_bins:
            self.__sbins = np.logspace(np.log10(smin), np.log10(smax), nbins + 1)
        else:
            self.__sbins = np.linspace(smin, smax, nbins + 1)

        self.__s_avg = 0.5 * (self.__sbins[1:] + self.__sbins[:-1])

        self.__gauss_x, self.__gauss_w = np.polynomial.legendre.leggauss(N_integration)

        self.__kmin = kmin
        self.__kmax = kmax
        self.__logk_vals = (np.log(kmax) - np.log(kmin)) / 2 * self.__gauss_x + (np.log(kmax) + np.log(kmin)) / 2
        self.__k_vals = np.exp(self.__logk_vals)#(kmax - kmin) / 2 * self.__gauss_x + (kmax + kmin) / 2

        if fid_cosmo is not None:
            self.set_fiducial_factor(fid_cosmo, fid_pk)
        else:
            self.__theta_bins = None
            self.__theta_avg = None
            self.__fid_factor = None

    def __geom_factor(self, k, r):
        return np.sum(np.diff(self.__theta_bins, axis=0)[:, :, None] * (self.__theta_bins[1:, :, None] * special.j1(k[None, None, :] * r[None, :, None] * self.__theta_bins[1:, :, None]) - self.__theta_bins[:-1, :, None] * special.j1(k[None, None, :] * r[None, :, None] * self.__theta_bins[:-1, :, None])) / (k[None, None, :] * r[None, :, None] * self.__theta_avg[:, :, None] * np.pi * (self.__theta_bins[1:, :, None]**2 - self.__theta_bins[:-1, :, None]**2)), axis=0)

    def __compute_cosmo_factor(self, cosmo, pk):
        """

        Parameters
        ----------
        pk : PowerSpectrumInterpolator
        cosmo : cosmo_from_camb
        """

        return (np.nansum(self.__k_vals[None, :] ** 2 * pk.P(self.__z_avg, self.__k_vals, grid=True) * self.__geom_factor(self.__k_vals, cosmo.comoving_angular_diameter_distance(cosmo.chi(self.__z_avg))) * self.__gauss_w[None, :], axis=-1) * (np.log(self.__kmax) - np.log(self.__kmin)) / 2 * self.__dz * cosmo.H(self.__z_avg))**(-1 / 2)

    def correction_factor(self, cosmo, pk, z_vals, normalise=False):
        if self.__fid_factor is None:
            raise ValueError("Correction factor can not be computed. No fiducial cosmology was provided.")

        corr = self.__compute_cosmo_factor(cosmo, pk)/self.__fid_factor
        if normalise:
            norm = np.sum(corr * np.diff(self.__zbins))/(np.max(self.__zbins) - np.min(self.__zbins))
        else:
            norm = 1.0
        out = np.ones_like(z_vals)
        out[(z_vals >= np.min(self.__zbins)) & (z_vals <= np.max(self.__zbins))] = np.interp(z_vals[(z_vals >= np.min(self.__zbins)) & (z_vals <= np.max(self.__zbins))], self.__z_avg, corr / norm)
        return out

    def set_fiducial_factor(self, fid_cosmo, fid_pk):
        self.__theta_bins = self.__sbins[:, None] / fid_cosmo.h / fid_cosmo.comoving_angular_diameter_distance(fid_cosmo.chi(self.__z_avg))[None, :]
        self.__theta_avg = self.__s_avg[:, None] / fid_cosmo.h / fid_cosmo.comoving_angular_diameter_distance(fid_cosmo.chi(self.__z_avg))[None, :]

        self.__fid_factor = self.__compute_cosmo_factor(fid_cosmo, fid_pk)

    def __call__(self, cosmo, pk, z_vals):
        return self.correction_factor(cosmo, pk, z_vals)


class CrossRedshiftCosmoCorrectionApprox(CrossRedshiftCosmoCorrection):

    def __init__(self, zmin=0.0, zmax=4.0, Nz=1000, fid_cosmo=None):
        super().__init__()
        self.__z_vals = np.linspace(zmin, zmax, Nz)
        if fid_cosmo is not None:
            self.set_fiducial_factor(fid_cosmo)
        else:
            self.__fid_factor = None

    @staticmethod
    def __compute_cosmo_factor(cosmo, z_vals):
        return cosmo.H(z_vals)**(1 / 2) / cosmo.chi(z_vals)

    def correction_factor(self, cosmo, pk, z_vals, normalise=False):
        return np.interp(z_vals, self.__z_vals, self.__compute_cosmo_factor(cosmo, self.__z_vals)/self.__fid_factor, left=1.0, right=1.0)

    def set_fiducial_factor(self, fid_cosmo):
        self.__fid_factor = self.__compute_cosmo_factor(fid_cosmo, self.__z_vals)

    def __call__(self, cosmo, pk, z_vals):
        return self.correction_factor(cosmo, pk, z_vals)


