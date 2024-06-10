import numpy as np
import os
from ..auxiliary.pk_interpolator import PowerSpectrumInterpolator


class cosmo_from_camb(object):
    r"""
    Class to store basic cosmology information from :doc:`camb:index`

    Parameters
    ----------
    camb_data : CAMBdata
        CAMB computation results of the type :class:`camb.results.CAMBdata`
    N_interpolation : int
        number of interpolation points for :math:`z(\chi)`
    zmin : float
        minimum redshift (default=0.0)
    zmax : float
        maximum redshift (default=4.0)
    include_nu_OmegaM : bool, optional
        include neutrinos in matter density (defaults to True)
    """

    def __init__(self, camb_data, include_nu_OmegaM=True):

        if camb_data.curv > 0:
            self.__comoving_angular_diameter_distance = lambda chi: camb_data.curvature_radius * np.sin(chi / camb_data.curvature_radius)
        elif camb_data.curv < 0:
            self.__comoving_angular_diameter_distance = lambda chi: camb_data.curvature_radius * np.sinh(chi / camb_data.curvature_radius)
        else:
            self.__comoving_angular_diameter_distance = lambda chi: chi

        self.__z_of_chi = camb_data.redshift_at_comoving_radial_distance
        self.__chi_of_z = camb_data.comoving_radial_distance
        self.__h_of_z = camb_data.h_of_z
        self.__curv = camb_data.curv

        self.__H0 = camb_data.h_of_z(0.0)
        self.__h = camb_data.Params.H0 / 100
        self.__Omega_m = (camb_data.Params.ombh2 + camb_data.Params.omch2 + camb_data.Params.omnuh2 * include_nu_OmegaM) / self.__h ** 2

        self.__chi_star = camb_data.conformal_time(0) - camb_data.tau_maxvis
        self.__z_star = camb_data.redshift_at_comoving_radial_distance(self.__chi_star)

        self.__camb_data = camb_data

    def H(self, z):
        r"""
        Hubble parameter as a function of redshift in units of :math:`\rm{Mpc}^{-1}`

        Parameters
        ----------
        z : array_like
            redshift
        Returns
        -------
        array_like
        """
        return self.__h_of_z(z)

    def chi(self, z):
        r"""
        Comoving distance as a function of redshift in units of :math:`\rm{Mpc}`

        Parameters
        ----------
        z : array_like
            redshift
        Returns
        -------
        array_like
        """
        return self.__chi_of_z(z)

    def z_of_chi(self, chi):
        r"""
        Redshift as a function of comoving distance.

        Parameters
        ----------
        chi : array_like
            Comoving distance in :math:`\rm{Mpc}`
        Returns
        -------
        array_like
        """
        return self.__z_of_chi(chi)

    def comoving_angular_diameter_distance(self, chi):
        r"""
        Angular diameter distance as a function of redshift in units of :math:`\rm{Mpc}`

        Parameters
        ----------
        chi : array_like
            comoving distance in :math:`\rm{Mpc}`
        Returns
        -------
        array_like
        """
        return self.__comoving_angular_diameter_distance(chi)

    @property
    def H0(self, unit="1/Mpc"):
        r"""
        Hubble rate

        Parameters
        ----------
        unit : str
            Available units are "1/Mpc" or "km/s/Mpc"
        Returns
        -------
        float
        """
        if unit == "1/Mpc":
            return self.__H0
        elif unit == "km/s/Mpc":
            return self.__h * 100
        else:
            raise Exception(f"Unit {unit} for H0 unknown.")

    @property
    def h(self):
        """
        Hubble parameter

        Returns
        -------
        float
        """
        return self.__h

    @property
    def Omega_m(self):
        """
        Matter density in units of the critical density

        Returns
        -------
        float
        """
        return self.__Omega_m

    @property
    def chi_star(self):
        """
        Comoving distance to surface of last scattering

        Returns
        -------
        float
        """
        return self.__chi_star

    @property
    def z_star(self):
        """
        Redshift at the surface of last scattering

        Returns
        -------
        float
        """
        return self.__z_star

    @property
    def curvature(self):
        """
        Curvature of the universe

        Returns
        -------
        float
        """
        return self.__curv


class dNdz(object):
    """
    Class for storing :math:`dN/dz` information

    Parameters
    ----------
    xmatch : callable
        interpolation of cross-match :math:`dN/dz`
    xcorr : callable
         interpolation of cross-correlation :math:`dN/dz`
    delta_xcorr_dndz_pcs : callable, optional
         interpolation of principal components of :math:`\Delta dN/dz`. First component must be mean :math:`\Delta dN/dz` then principal components in order of decreasing significance.

    """

    def __init__(self, xmatch, xcorr, delta_xcorr_dndz_pcs=None):
        self.__xmatch = xmatch
        self.__xcorr = xcorr
        self.__delta_xcorr_dndz_pcs = delta_xcorr_dndz_pcs

    def __call__(self, z, cross_correlation=False):
        """
        Evaluate :math:`dN/dz`. By default will return cross-match :math:`dN/dz`. Set `cross-correlation` to True to get cross-correlation redshifts. If principal components are provided and `cross-correlation=True` the function will return cross-correlation redshifts stacked with the principal components for :math:`dN/dz` marginalisation.

        Parameters
        ----------
        z : array_like
            redshifts
        cross_correlation : bool
            whether to use cross-correlation redshifts (will default to using cross-match redshifts)
        Returns
        -------
        array_like
        """
        if cross_correlation:
            return self.bdNdz(z)
        else:
            return self.dNdz(z)

    def dNdz(self, z):
        return self.__xmatch(z)

    def bdNdz(self, z, pcs=False):
        if not pcs:
            return self.__xcorr(z)
        else:
            if self.__delta_xcorr_dndz_pcs is not None:
                return np.hstack([self.__xcorr(z)[:, None], self.__delta_xcorr_dndz_pcs(z)])
            else:
                return self.bdNdz(z, pcs=False)[:, None]

    @property
    def n_pcs(self):
        if self.__delta_xcorr_dndz_pcs is not None:
            return len(self.__delta_xcorr_dndz_pcs(0.0))
        else:
            return 0


class CleftInterpolationHelper(object):
    def __init__(self,
                 bias_coevolution,
                 use_second_order=True,
                 use_shear=True,
                 use_third_order=False,
                 nk=100,
                 kmin=1.0e-3,
                 kmax=0.3,
                 extrap_kmax=None,
                 threads=None,
                 cutoff=10,
                 jn=5,
                 N=2700):
        self._bias_coevolution = bias_coevolution
        self._use_second_order = use_second_order
        self._use_shear = use_shear
        self._use_third_order = use_third_order

        self._N = N
        self._jn = jn
        self._cutoff = cutoff
        self._nk = nk
        self._kmin = kmin
        self._kmax = kmax
        self._extrap_kmax = extrap_kmax

        self._cleft = None
        if threads is None:
            try:
                self._threads = len(os.sched_getaffinity(0))
            except AttributeError:
                self._threads = None
        else:
            self._threads = threads

    def compute_cleft_spectra(self, k, z, Pk, get_table=False):
        try:
            from velocileptors.LPT.cleft_fftw import CLEFT
        except ImportError:
            raise ImportError("velocileptors not installed. Please install velocileptors to use this functionality")

        cleft_table = None
        if self._cleft is None:
            self._cleft = CLEFT(k, Pk[0, :], shear=self._use_shear, third_order=self._use_third_order, threads=self._threads, N=self._N, cutoff=self._cutoff, jn=self._jn)

        for i in range(len(z)):
            self._cleft.update_power_spectrum(k, Pk[i, :])
            self._cleft.make_ptable(kmin=self._kmin, kmax=self._kmax, nk=self._nk, )
            if i == 0:
                cleft_table = np.zeros((len(z), *self._cleft.pktable.shape))
            cleft_table[i] = self._cleft.pktable

        if not get_table:
            return self.interpolate_cleft_spectra(z, cleft_table[0, :, 0], cleft_table[..., 1:])
        else:
            return cleft_table

    def interpolate_cleft_spectra(self, z, k_vals, cleft_table, offset=3):
        cleft_interpolations = []
        for i in range(offset, offset + 3 * self._use_second_order + 3 * self._use_shear + 2 * self._use_third_order + self._use_second_order * self._use_shear):
            sel = np.all(np.isfinite(cleft_table[:, :, i]), axis=0)
            if np.all(cleft_table[:, sel, i] < 0):
                cleft_interpolations.append(PowerSpectrumInterpolator(z, k_vals[sel], np.log(-cleft_table[:, sel, i]), extrap_kmax=max([self._extrap_kmax if self._extrap_kmax is not None else self._kmax, self._kmax]), logP=True, logsign=-1))
            elif np.all(cleft_table[:, sel, i] > 0):
                cleft_interpolations.append(PowerSpectrumInterpolator(z, k_vals[sel], np.log(cleft_table[:, sel, i]), extrap_kmax=max([self._extrap_kmax if self._extrap_kmax is not None else self._kmax, self._kmax]), logP=True, logsign=1))
            else:
                cleft_interpolations.append(PowerSpectrumInterpolator(z, k_vals[sel], cleft_table[:, sel, i], logP=False, logsign=1))

        return cleft_interpolations

    def multiply_eft_biases(self, cleft_pk, z_vals, fid_bias_evol, fid_bias_evol_2=None):
        bias_matrix_mg = np.zeros((cleft_pk.shape[-1], len(z_vals)))  # b2/2, 0.0, 0.0, bs/2, 0.0, 0.0, 0.0, b3/2, 0.0
        # the terms in bias_matrix_gg_b are multiplied by b f(z)dN/dz hence it includes only terms that are proportional to b1 (i.e. P_b1b2, P_b1bs, etc.)
        # since b is the eulerian bias rather than the lagrangian ( b f(z) =b1(z)+1 ) for each P_b1bX we have to a) subtract out bX*P_b1bX in the nob1 term
        # and b) divide bX by two in the b1 term.
        bias_matrix_gg_b = np.zeros((cleft_pk.shape[-1], len(z_vals)))  # 0, b*b2, 0, 0, b*bs, 0, 0, 0, b*b3
        bias_matrix_gg_nob = np.zeros((cleft_pk.shape[-1], len(z_vals)))  # b2, -b2, b2**2, bs, -bs, b2*bs, bs**2, b3, -b3

        if fid_bias_evol_2 is None:
            fid_bias_evol_2 = fid_bias_evol

        if cleft_pk.shape[-1] in [3, 7, 9]:
            b1 = fid_bias_evol(z_vals) - 1
            b1_2 = fid_bias_evol_2(z_vals) - 1

            b2 = self._bias_coevolution.get_b2(b1)
            b2_2 = self._bias_coevolution.get_b2(b1_2)

            bias_matrix_mg[0] = b2 / 2

            bias_matrix_gg_nob[0] = b2
            bias_matrix_gg_b[1] = b2
            bias_matrix_gg_nob[1] = -b2
            bias_matrix_gg_nob[2] = b2 * b2_2

        if cleft_pk.shape[-1] in [7, 9]:
            bs = self._bias_coevolution.get_bs(b1)
            bs_2 = self._bias_coevolution.get_bs(b1_2)

            bias_matrix_mg[3] = bs / 2

            bias_matrix_gg_nob[3] = bs
            bias_matrix_gg_b[4] = bs
            bias_matrix_gg_nob[4] = -bs
            bias_matrix_gg_nob[5] = b2 * bs_2
            bias_matrix_gg_nob[6] = bs * bs_2

        if cleft_pk.shape[-1] == 9:
            b3 = self._bias_coevolution.get_b3(b1)

            bias_matrix_mg[7] = b3 / 2

            bias_matrix_gg_nob[7] = b3
            bias_matrix_gg_b[8] = b3
            bias_matrix_gg_nob[8] = -b3

        if cleft_pk.shape[-1] not in [3, 7, 9]:
            raise Exception("Length of CLEFT components unknown")

        return cleft_pk * bias_matrix_gg_b.T[:, None, :], cleft_pk * bias_matrix_gg_nob.T[:, None, :], cleft_pk * bias_matrix_mg.T[:, None, :]

    def assemble_cleft_coeff(self, cb2=1.0, cbs=1.0, cb3=1.0, cb2_2=None, cbs_2=None, cb3_2=None):
        if cb2_2 is not None or cbs_2 is not None or cb3_2 is not None:
            assert (cb2_2 is not None and (cbs_2 is not None or not self._use_shear) and (cb3_2 is not None or not self._use_third_order))
        else:
            assert (cb2_2 is None and cbs_2 is None and cb3_2 is None)
            cb2_2 = cb2
            cbs_2 = cbs
            cb3_2 = cb3

        if self._use_second_order and self._use_shear and self._use_third_order:
            return np.array([cb2, cb2, cb2 ** 2, cbs, cbs, cb2 * cbs, cbs ** 2, cb3, cb3])
        elif self._use_second_order and self._use_shear:
            return np.array([cb2, cb2, cb2 ** 2, cbs, cbs, cb2 * cbs, cbs ** 2])
        elif self._use_second_order:
            return np.array([cb2, cb2, cb2 ** 2])
        else:
            return None

    def __reduce__(self):
        self._cleft = None
        return super().__reduce__()


class CleftInterpolationHelperFreeCleft(CleftInterpolationHelper):

    def multiply_eft_biases(self, cleft_pk, z_vals, fid_bias_evol, fid_bias_evol_2=None):
        use_b2 = cleft_pk.shape[-1] >= 3
        use_bs = cleft_pk.shape[-1] >= 7
        use_b3 = cleft_pk.shape[-1] == 9

        if fid_bias_evol_2 is None:
            fid_bias_evol_2 = fid_bias_evol

        if cleft_pk.shape[-1] not in [3, 7, 9]:
            raise Exception("Length of CLEFT components unknown")

        bias_matrix_1 = np.zeros((len(z_vals), use_b2 + use_bs + use_b3, 2))  # b2, bs, b3
        bias_matrix_2 = np.zeros((len(z_vals), use_b2 + 2 * use_bs, 4))  # b2**2, b2*bs, bs**2

        if use_b2 and self._use_second_order:
            b1 = fid_bias_evol(z_vals) - 1
            b1_2 = fid_bias_evol_2(z_vals) - 1

            b2 = self._bias_coevolution.get_b2(b1)
            b2_2 = self._bias_coevolution.get_b2(b1_2)

            bias_matrix_1[:, 0] = np.vstack([b2, np.ones(len(z_vals))]).T
            bias_matrix_2[:, 0] = np.vstack([b2 * b2_2, b2, b2_2, np.ones(len(z_vals))]).T

        if use_bs and self._use_shear:
            bs = self._bias_coevolution.get_bs(b1)
            bs_2 = self._bias_coevolution.get_bs(b1_2)

            bias_matrix_1[:, 1] = np.vstack([bs, np.ones(len(z_vals))]).T
            bias_matrix_2[:, 1] = np.vstack([b2 * bs_2, b2, bs_2, np.ones(len(z_vals))]).T
            bias_matrix_2[:, 2] = np.vstack([bs * bs_2, bs, bs_2, np.ones(len(z_vals))]).T

        if use_b3 and self._use_third_order:
            b3 = self._bias_coevolution.get_b3(b1)

            bias_matrix_1[:, 2] = np.vstack([b3, np.ones(len(z_vals))]).T

        sel_1 = [0, 3, 7][:use_b2 + use_bs + use_b3]
        sel_2 = [1, 4, 8][:use_b2 + use_bs + use_b3]
        sel_3 = [2, 5, 6][:use_b2 + 2 * use_bs]

        return cleft_pk[:, :, sel_2, None] * bias_matrix_1[:, None, :, :], (cleft_pk[:, :, sel_1, None] - cleft_pk[:, :, sel_2, None]) * bias_matrix_1[:, None, :, :], cleft_pk[:, :, sel_3, None] * bias_matrix_2[:, None, :, :], cleft_pk[:, :, sel_1, None] * bias_matrix_1[:, None, :, :] / 2

    def assemble_cleft_coeff(self, cb2=(1.0, 0.0), cbs=(1.0, 0.0), cb3=(1.0, 0.0), cb2_2=None, cbs_2=None, cb3_2=None):
        cb2 = np.array(cb2)
        cbs = np.array(cbs)
        cb3 = np.array(cb3)

        if cb2_2 is not None or cbs_2 is not None or cb3_2 is not None:
            assert (cb2_2 is not None and (cbs_2 is not None or not self._use_shear) and (cb3_2 is not None or not self._use_third_order))
        else:
            assert (cb2_2 is None and cbs_2 is None and cb3_2 is None)
            cb2_2 = cb2
            cbs_2 = cbs
            cb3_2 = cb3

        if self._use_second_order and self._use_shear and self._use_third_order:
            return np.array([cb2 + cb2_2, cbs + cbs_2, cb3 + cb3_2]) / 2, np.array([np.outer(cb2, cb2_2).flatten(), np.outer(cb2 + cb2_2, cbs + cbs_2).flatten() / 4, np.outer(cbs, cbs_2).flatten()])
        elif self._use_second_order and self._use_shear:
            return np.array([cb2 + cb2_2, cbs + cbs_2]) / 2, np.array([np.outer(cb2, cb2_2).flatten(), np.outer(cb2 + cb2_2, cbs + cbs_2).flatten() / 4, np.outer(cbs, cbs_2).flatten()])
        elif self._use_second_order:
            return np.array([cb2 + cb2_2]) / 2, np.array([np.outer(cb2, cb2_2).flatten()])
        else:
            return None, None
