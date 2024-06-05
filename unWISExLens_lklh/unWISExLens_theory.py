from cobaya.theory import Theory
from cobaya.log import LoggedError
import itertools
import numpy as np
import inspect
from scipy.interpolate import pchip as interpolate, interp1d, UnivariateSpline
import cobaya
import os

from .auxiliary.dN_dz_aux import dN_dz_Helper
from .auxiliary.pk_interpolator import PowerSpectrumInterpolator, ScaledPowerSpectrumInterpolator
from .auxiliary.auxiliary_functions import combine_inputs

from .theory_modules.unWISExkappa_model import unWISExACT_theory_model
from .theory_modules.unWISExkappa_model_freeCLEFT import unWISExACT_theory_model as unWISExACT_theory_model_new

from .theory_modules.model_helpers_unWISExACT import cosmo_from_camb, dNdz, CleftInterpolationHelper, CleftInterpolationHelperFreeCleft
from .theory_modules.cross_dNdz_cosmo_correction import CrossRedshiftCosmoCorrectionExact, CrossRedshiftCosmoCorrectionApprox

from .auxiliary.bias_coevolution import BiasEvolution


class unWISExLens_theory(Theory):
    stop_at_error = True

    data_base_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")

    zmax = 3.0
    zmin = 0.0
    Nz = 256
    k_max = 10.0

    N_integration = 512
    use_linear_theory = False
    use_cleft = True
    use_free_cleft_model = True
    use_fiducial_cleft = True
    use_Az_parametrisation = False
    Az_parametrisation = {'type': 'bins', 'bin_edges': [1.09, 1.75, np.inf], 'params': [f'A{i}' for i in range(3)]}
    compute_gg_cross_spectra = False

    cleft_interp_helper = None
    fiducial_cleft_interpolations_nonu_nonu = None
    fiducial_cleft_interpolations_tot_nonu = None

    correct_clustering_redshift_cosmo = True
    correct_clustering_redshift_cosmo_approx = False
    clustering_redshift_fid_cosmo = {'H0': 67.66, 'ombh2': 0.02242, 'omch2': 0.11935351837638222, 'mnu': 0.06, 'nnu': 3.046, 'ns': 0.9665, 'As': 2.105209331337507e-09}
    _clustering_redshift_cosmo_correction = None

    use_mean_field_correction = True

    bxdN_dz_paths = {'Blue': 'aux_data/dndz/unWISE_blue_xcorr_bdndz.txt', 'Green': 'aux_data/dndz/unWISE_green_xcorr_bdndz.txt'}
    xmatch_dN_dz_paths = {'Blue': 'aux_data/dndz/unWISE_blue_xmatch_dndz.txt', 'Green': 'aux_data/dndz/unWISE_green_xmatch_dndz.txt'}
    delta_dndz_pcs_paths = {'Blue': 'aux_data/dndz/unWISE_blue_delta_bdndz_pcs.dat', 'Green': 'aux_data/dndz/unWISE_green_delta_bdndz_pcs.dat'}

    noise_bias_corr_paths = {'Blue': 'aux_data/noise_bias_correction/unWISE_blue_noise_bias_corr_free_cleft.dat', 'Green': 'aux_data/noise_bias_correction/unWISE_green_noise_bias_corr_free_cleft.dat', 'BluexGreen': 'aux_data/noise_bias_correction/unWISE_blueXgreen_noise_bias_corr_free_cleft.dat'}

    fid_bias_evol = {'Blue': 'lambda z:0.8 + 1.2*z', 'Green': 'lambda z: np.clip(1.6*z**2, 1, np.inf)'}
    _fid_bias_evol_list = []

    noise_bias_corr = []
    _cross_spec_noise_bias_corr = {}
    _samples = []
    _dNdz_list = []
    _theory_code = None
    _Az_function = None

    _want_clkk = False
    _want_kg_gg = False

    _ell_vals = np.array([])
    _ell_vals_clkk = np.array([])



    fiducial_cleft_interpolations = None  # deprecated (maintained for compatibility)
    use_pk_nu_cross = None  # deprecated (maintained for compatibility)
    N_interpolation = 1024  # deprecated (maintained for compatibility)

    def initialize(self):
        """called from __init__ to initialize"""

        self.log.info("Getting CLEFT helpers and interpolations ...")
        if self.use_cleft and self.cleft_interp_helper is None:
            assert (self.data_base_path is not None)
            if self.use_free_cleft_model:
                interp_helper_class = CleftInterpolationHelperFreeCleft
            else:
                interp_helper_class = CleftInterpolationHelper

            coevl = BiasEvolution(path=os.path.join(self.data_base_path, "aux_data/bias_coevol_coeff.txt"))

            self.cleft_interp_helper = interp_helper_class(coevl, kmax=3.0)

        if self.use_cleft and self.use_fiducial_cleft:
            if self.fiducial_cleft_interpolations_nonu_nonu is None or self.fiducial_cleft_interpolations_tot_nonu is None:
                k, z, fid_cleft_table = np.load(os.path.join(self.data_base_path, "aux_data/fiducial_cleft_spectra.npy"), allow_pickle=True)
                self.fiducial_cleft_interpolations_nonu_nonu = self.cleft_interp_helper.interpolate_cleft_spectra(z, fid_cleft_table[0, :, 0], fid_cleft_table[..., 1:])
                self.fiducial_cleft_interpolations_tot_nonu = self.fiducial_cleft_interpolations_nonu_nonu

        if self.use_Az_parametrisation:
            if self.Az_parametrisation['type'] == 'bins' or self.Az_parametrisation['type'] == 'bins+slope':

                if self.Az_parametrisation['bin_edges'][0] == 0.0:  # remove 0.0 bin edge
                    self.Az_parametrisation['bin_edges'] = self.Az_parametrisation['bin_edges'][1:]

                if self.Az_parametrisation['bin_edges'][-1] != np.inf:  # add infinity bin edge
                    self.Az_parametrisation['bin_edges'] = np.concatenate([self.Az_parametrisation['bin_edges'], [np.inf]])

                if self.Az_parametrisation['type'] == 'bins':
                    def func(z, *params):
                        params = np.array(params)
                        return params[np.digitize(z, self.Az_parametrisation['bin_edges'])]

                    self._Az_function = func

                elif self.Az_parametrisation['type'] == 'bins+slope':
                    def func(z, *params):
                        params = np.array(params)
                        return np.clip(params[:len(self.Az_parametrisation['bin_edges'])][np.digitize(z, self.Az_parametrisation['bin_edges'])] + params[len(self.Az_parametrisation['bin_edges']):][np.digitize(z, self.Az_parametrisation['bin_edges'])] * (z - np.concatenate([[0.0], self.Az_parametrisation['bin_edges']])[np.digitize(z, self.Az_parametrisation['bin_edges'])]), 0, np.inf)

                    self._Az_function = func
            else:
                raise LoggedError(self.log, f"Unknown Az parametrisation type {self.Az_parametrisation['type']}.")

    def initialize_with_provider(self, provider):
        """
        Initialization after other components initialized, using Provider class
        instance which is used to return any dependencies (see calculate below).
        """
        self.provider = provider

        camb = None
        for c in provider.model.components:
            if isinstance(c, cobaya.theories.camb.camb.CAMB):
                camb = c.camb

        if camb is None:
            import camb

        if self.correct_clustering_redshift_cosmo or self.correct_clustering_redshift_cosmo_approx:
            pars = camb.CAMBparams()
            pars.set_cosmology(**{key: self.clustering_redshift_fid_cosmo[key] for key in set(inspect.signature(pars.set_cosmology).parameters.keys()).intersection(set(self.clustering_redshift_fid_cosmo.keys()))})
            pars.InitPower.set_params(**{key: self.clustering_redshift_fid_cosmo[key] for key in set(inspect.signature(pars.InitPower.set_params).parameters.keys()).intersection(set(self.clustering_redshift_fid_cosmo.keys()))})
            pars.set_matter_power(redshifts=np.linspace(0, 4.0, 41), kmax=10.0)

            # Non-Linear spectra (Halofit)
            pars.NonLinear = camb.model.NonLinear_both
            pars.NonLinearModel.set_params(halofit_version='mead')
            results = camb.get_results(pars)
            k_nonlin, z_nonlin, pk_nonlin = results.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False, var1='delta_nonu', var2='delta_nonu')
            fid_pk_interp = PowerSpectrumInterpolator(z_nonlin, k_nonlin, pk_nonlin)
            fid_cosmo = cosmo_from_camb(results)

            zbins = np.concatenate([np.arange(0, 0.8, 0.05), np.arange(0.8, 4.0, 0.2), [4.0]])

            if self.correct_clustering_redshift_cosmo_approx and not self.correct_clustering_redshift_cosmo:
                self._clustering_redshift_cosmo_correction = CrossRedshiftCosmoCorrectionApprox(fid_cosmo)
            else:
                self._clustering_redshift_cosmo_correction = CrossRedshiftCosmoCorrectionExact(fid_cosmo, fid_pk_interp, zbins, smin=2.5, smax=10, nbins=3, log_bins=True, N_integration=1024, kmin=1e-4, kmax=10)

    def get_requirements(self):
        """
        Return dictionary of derived parameters or other quantities that are needed
        by this component and should be calculated by another theory class.
        """
        if self.use_Az_parametrisation:
            return dict(zip(self.Az_parametrisation['params'], itertools.repeat(None)))
        else:
            return {}

    def must_provide(self, **requirements):
        requires = {}
        zmin = self.zmin
        zmax = self.zmax
        z_vals = np.logspace(np.log10(self.zmin + 1), np.log10(self.zmax + 1), self.Nz) - 1

        pk_vars_pairs = set()

        if 'clkk' in requirements:
            self._want_clkk = True
            self._ell_vals_clkk = combine_inputs(requirements['clkk'].pop('ell_vals', []), self._ell_vals_clkk)
            if zmin > 0.0 or zmax < 1100:
                raise LoggedError(self.log, f"Clkk calculation requires zmin=0.0 and zmax>=1100.0, but zmin={zmin} and zmax={zmax} were provided.")
            pk_vars_pairs.add(("Weyl", "Weyl"))

        if 'cls' in requirements:
            self._want_kg_gg = True
            if 'samples' in requirements['cls']:
                self._samples = list(requirements['cls']['samples'])

            self._ell_vals = combine_inputs(requirements['cls'].pop('ell_vals', []), self._ell_vals)

            pk_vars_pairs.add(("Weyl", "Weyl"))
            pk_vars_pairs.add(("delta_nonu", "delta_nonu"))
            pk_vars_pairs.add(("Weyl", "delta_nonu"))

            assert (not self.use_cleft or not self.use_linear_theory), "Cannot use CLEFT when linear theory is activated."

            if self.use_cleft and not self.use_fiducial_cleft:
                requires['Pk_grid'] = {'z': z_vals[z_vals<=3.0], 'k_max': self.k_max, 'nonlinear': False, 'vars_pairs': [["delta_nonu", "delta_nonu"], ["delta_tot", "delta_nonu"]]}

        if self._want_clkk or self._want_kg_gg:
            # Nz = self.Nz if not self._want_clkk else int(np.round((zmax - zmin) / (self.zmax - self.zmin) * self.Nz))

            requires['Pk_interpolator'] = {'z': z_vals, 'k_max': self.k_max, 'nonlinear': not self.use_linear_theory, 'vars_pairs': pk_vars_pairs}
            if self.use_Az_parametrisation and not self.use_linear_theory:
                requires['Pk_interpolator'] = {'z': z_vals, 'k_max': self.k_max, 'nonlinear': (False, True), 'vars_pairs': pk_vars_pairs}
            requires['CAMBdata'] = None

            theory_params = {'zmax': self.zmax,
                             'zmin': self.zmin,
                             'k_max': self.k_max,
                             'N_integration': self.N_integration,
                             'cross_correlation_redshift_correction': self._clustering_redshift_cosmo_correction,
                             'cleft_interp_helper': self.cleft_interp_helper,
                             'ell_vals': self._ell_vals,
                             'want_gg_cross': self.compute_gg_cross_spectra,
                             'ell_vals_clkk': self._ell_vals_clkk}

            if self.use_free_cleft_model:
                self._theory_code = unWISExACT_theory_model_new(**theory_params)
            else:
                self._theory_code = unWISExACT_theory_model(**theory_params)

            self.load_sample_data()

        return requires

    def load_sample_data(self):
        self.noise_bias_corr = []

        self._dNdz_list = []
        for i, sample in enumerate(self._samples):
            # print(f"Loading dn/dz from path {self.bxdN_dz_paths[sample]}.")
            bxdN_dz, dump, dump = dN_dz_Helper.get_dn_dz(os.path.join(self.data_base_path, self.bxdN_dz_paths[sample]))

            dN_dz_xmatch, zmin, zmax = dN_dz_Helper.get_dn_dz(os.path.join(self.data_base_path, self.xmatch_dN_dz_paths[sample]))

            delta_dndz_pcs_data = np.loadtxt(os.path.join(self.data_base_path, self.delta_dndz_pcs_paths[sample]))
            delta_dndz_pcs_interp = interpolate(delta_dndz_pcs_data[:, 0], delta_dndz_pcs_data[:, 1:])

            self._dNdz_list.append(dNdz(dN_dz_xmatch, bxdN_dz, delta_dndz_pcs_interp))

            if self.use_mean_field_correction:
                noise_bias_corr = np.load(os.path.join(self.data_base_path, self.noise_bias_corr_paths[sample]), allow_pickle=True).item()

                noise_bias_ell_selection, dump = np.where(noise_bias_corr['ell'][:, np.newaxis] == self._theory_code.ell_vals[np.newaxis, :])
                assert (len(noise_bias_ell_selection) == len(self._theory_code.ell_vals))
                if self.cleft_interp_helper is not None:
                    self.noise_bias_corr.append({'gg': {'gg_bsq': noise_bias_corr['gg']['gg_bsq'][noise_bias_ell_selection],
                                                                'gg_b': noise_bias_corr['gg']['gg_b'][noise_bias_ell_selection],
                                                                'gmu_b': noise_bias_corr['gg']['gmu_b'][noise_bias_ell_selection]},
                                                         'kg': {'kg_b': noise_bias_corr['kg']['kg_b'][noise_bias_ell_selection]}})
                else:
                    self.noise_bias_corr.append({'gg': {'gg_bsq': noise_bias_corr['gg']['gg_bsq'][noise_bias_ell_selection],
                                                                'gg_b': np.zeros((len(self._theory_code.ell_vals), 1)),
                                                                'gmu_b': noise_bias_corr['gg']['gmu_b'][noise_bias_ell_selection]},
                                                         'kg': {'kg_b': noise_bias_corr['kg']['kg_b'][noise_bias_ell_selection]}})

                if self.compute_gg_cross_spectra:
                    for j, sample2 in enumerate(self._samples[i + 1:]):
                        j += i + 1

                        noise_bias_corr = np.load(os.path.join(self.data_base_path, self.noise_bias_corr_paths[f"{sample}x{sample2}"]), allow_pickle=True).item()
                        noise_bias_ell_selection, dump = np.where(noise_bias_corr['ell'][:, np.newaxis] == self._theory_code.ell_vals[np.newaxis, :])
                        assert (len(noise_bias_ell_selection) == len(self._theory_code.ell_vals))
                        if self.cleft_interp_helper is not None:
                            self._cross_spec_mean_field_corrections[(i, j)] = {'g1g2': {'g1b_g2b': noise_bias_corr['g1g2']['g1b_g2b'][noise_bias_ell_selection],
                                                                                        'g1b_g2': noise_bias_corr['g1g2']['g1b_g2'][noise_bias_ell_selection],
                                                                                        'g2b_g1': noise_bias_corr['g1g2']['g2b_g1'][noise_bias_ell_selection]},
                                                                               'g1mu2': {'gmu_b': noise_bias_corr['g1mu2']['gmu_b'][noise_bias_ell_selection]},
                                                                               'g2mu1': {'gmu_b': noise_bias_corr['g2mu1']['gmu_b'][noise_bias_ell_selection]}}
                        else:
                            self._cross_spec_mean_field_corrections[(i, j)] = {'g1g2': {'g1b_g2b': noise_bias_corr['g1g2']['g1b_g2b'][noise_bias_ell_selection],
                                                                                        'g1b_g2': np.zeros((len(self._theory_code.ell_vals), 1)),
                                                                                        'g2b_g1': np.zeros((len(self._theory_code.ell_vals), 1))},
                                                                               'g1mu2': {'gmu_b': noise_bias_corr['g1mu2']['gmu_b'][noise_bias_ell_selection]},
                                                                               'g2mu1': {'gmu_b': noise_bias_corr['g2mu1']['gmu_b'][noise_bias_ell_selection]}}
            if not callable(self.fid_bias_evol[sample]):
                try:
                    self._fid_bias_evol_list.append(eval(self.fid_bias_evol[sample]))
                except SyntaxError as ex:
                    raise LoggedError(self.log, f"The fiducial bias evolution {self.fid_bias_evol[sample]} can not be interpreted as function.")
            else:
                self._fid_bias_evol_list.append(self.fid_bias_evol[sample])

    def get_can_provide_params(self):
        return []

    def get_pk_interp(self, var_pair=("delta_tot", "delta_tot")):
        if not self.use_Az_parametrisation:
            return self.provider.get_Pk_interpolator(var_pair=var_pair, nonlinear=not self.use_linear_theory, extrap_kmax=self.k_max)
        else:
            lin_power_interp = self.provider.get_Pk_interpolator(var_pair=var_pair, nonlinear=False, extrap_kmax=self.k_max)
            if not self.use_linear_theory:
                nonlin_power_interp = self.provider.get_Pk_interpolator(var_pair=var_pair, nonlinear=True, extrap_kmax=self.k_max)
            else:
                nonlin_power_interp = None

            args = [self.provider.get_param(p) for p in self.Az_parametrisation['params']]
            return ScaledPowerSpectrumInterpolator(lambda z: self._Az_function(z, *args), lin_power_interp, nonlin_power_interp)

    def calculate(self, state, want_derived=True, **params_values_dict):
        cosmo = cosmo_from_camb(self.provider.get_CAMBdata())
        pk_interp_nonlin_weyl_weyl = self.get_pk_interp(var_pair=("Weyl", "Weyl"))
        if self._want_kg_gg:
            pk_interp_nonlin_nonu_nonu = self.get_pk_interp(var_pair=("delta_nonu", "delta_nonu"))
            pk_interp_nonlin_weyl_nonu = self.get_pk_interp(var_pair=("Weyl", "delta_nonu"))

            if self.use_cleft and not self.use_fiducial_cleft:
                k, z, Pk = self.provider.get_Pk_grid(var_pair=("delta_nonu", "delta_nonu"), nonlinear=False)
                cleft_interpolations = self.cleft_interp_helper.compute_cleft_spectra(k / cosmo.h, z, Pk * cosmo.h ** 3)

                k, z, Pk_tot_nonu = self.provider.get_Pk_grid(var_pair=("delta_tot", "delta_nonu"), nonlinear=False)
                cleft_interpolations_tot_nonu = self.cleft_interp_helper.compute_cleft_spectra(k / cosmo.h, z, Pk_tot_nonu * cosmo.h ** 3)

            elif self.use_cleft and self.use_fiducial_cleft:
                cleft_interpolations = self.fiducial_cleft_interpolations_nonu_nonu
                cleft_interpolations_tot_nonu = self.fiducial_cleft_interpolations_tot_nonu
            else:
                cleft_interpolations = None
                cleft_interpolations_tot_nonu = None

            state['cls'] = self._theory_code.compute_raw_spectra(cosmo, self._dNdz_list, pk_interp_nonlin_weyl_weyl, pk_interp_nonlin_weyl_nonu, pk_interp_nonlin_nonu_nonu, cleft_interpolations_tot_nonu, cleft_interpolations, self._fid_bias_evol_list if cleft_interpolations is not None else None)

        if self._want_clkk:
            state['clkk'] = self._theory_code.compute_clkk(cosmo, pk_interp_nonlin_weyl_weyl)

    def get_kg(self, samples=None, ell=None, **kwargs):
        if ell is not None:
            ell_selection, dump = np.where(self._theory_code.ell_vals[:, np.newaxis] == ell)
        else:
            ell_selection = np.full(len(self._theory_code.ell_vals), True)

        if samples is None:
            indices = list(range(len(self._samples)))
        else:
            indices = [self._samples.index(s) for s in samples]

        if not self.compute_gg_cross_spectra:
            raw_spectra = [self._current_state['cls'][i] for i in indices]
        else:
            raw_spectra = [self._current_state['cls'][0][i] for i in indices], dict([(k, self._current_state['cls'][1][k]) for k in itertools.combinations(indices, 2)])
        mean_field = [self.noise_bias_corr[i] for i in indices] if self.use_mean_field_correction else None

        return [out[ell_selection] for out in self._theory_code.evaluate_kg(raw_spectra, mean_field=mean_field, **kwargs)]

    def get_gg(self, ell=None, samples=None, want_gg_cross=False, **kwargs):
        if ell is not None:
            ell_selection, dump = np.where(self._theory_code.ell_vals[:, np.newaxis] == ell)
        else:
            ell_selection = np.full(len(self._theory_code.ell_vals), True)

        if samples is None:
            indices = list(range(len(self._samples)))
        else:
            indices = [self._samples.index(s) for s in samples]

        mean_field = [self.noise_bias_corr[i] for i in indices] if self.use_mean_field_correction else None

        if not self.compute_gg_cross_spectra or not want_gg_cross:
            raw_spectra = [self._current_state['cls'][i] for i in indices]
            return [out[ell_selection] for out in self._theory_code.evaluate_gg(raw_spectra, mean_field=mean_field, **kwargs)]
        else:
            raw_spectra = [self._current_state['cls'][0][i] for i in indices], dict([(k, self._current_state['cls'][1][k]) for k in itertools.combinations(indices, 2)])
            gg_auto, gg_cross = self._theory_code.evaluate_gg(raw_spectra, mean_field=mean_field, want_gg_cross=True, **kwargs)
            return [out[ell_selection] for out in gg_auto], gg_cross[:, :, ell_selection]

    def get_cls(self, ell=None, samples=None, want_gg_cross=False, **kwargs):
        if ell is not None:
            ell_selection, dump = np.where(self._theory_code.ell_vals[:, np.newaxis] == ell)
        else:
            ell_selection = np.full(len(self._theory_code.ell_vals), True)

        if samples is None:
            indices = list(range(len(self._samples)))
        else:
            indices = [self._samples.index(s) for s in samples]

        mean_field = [self.noise_bias_corr[i] for i in indices] if self.use_mean_field_correction else None

        if not self.compute_gg_cross_spectra or not want_gg_cross:
            raw_spectra = [self._current_state['cls'][i] for i in indices]
            out_gg, out_kg = self._theory_code.evaluate(raw_spectra, mean_field=mean_field, **kwargs)
            return [out[ell_selection] for out in out_gg], [out[ell_selection] for out in out_kg]
        else:
            raw_spectra = [self._current_state['cls'][0][i] for i in indices], dict([(k, self._current_state['cls'][1][k]) for k in itertools.combinations(indices, 2)])
            out_gg, gg_cross, out_kg = self._theory_code.evaluate_gg(raw_spectra, mean_field=mean_field, want_gg_cross=True, **kwargs)
            return [out[ell_selection] for out in out_gg], gg_cross[:, :, ell_selection], [out[ell_selection] for out in out_kg]

    def get_clkk(self, ell=None):
        if ell is not None:
            ell_selection, dump = np.where(self._theory_code.ell_vals_clkk[:, np.newaxis] == ell)
            return self._current_state['clkk'][ell_selection]
        else:
            return self._current_state['clkk']

    def get_theory(self):
        return self._theory_code
