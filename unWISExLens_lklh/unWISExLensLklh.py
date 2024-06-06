import warnings

import numpy as np
import healpy as hp
from cobaya.likelihood import Likelihood
from cobaya.yaml import yaml_load_file
import os

from .auxiliary.auxiliary_functions import density_space, select_from_matrix
from .auxiliary.binning_helpers import PowerSpectrumBinning, MatrixPowerSpectrumBinning, NaMasterPowerSpectrumBinning


class unWISExLensLklh(Likelihood):

    _binning_function_gg = []  # will default to simply binning the spectrum by taking the unweighted average in every bin
    _binning_function_kg = []
    _binning_function_kk = []
    _ell_vals = None
    _ell_vals_kk = None
    _ell_selections = []
    _ell_selections_kk = []
    _ell_conds = []
    _ell_conds_kk = []
    _lmax_kk = 0
    _lmin_kk = 0
    _lmax_TEB = 0
    _data_gg = []
    _data_kg = []
    _data_kk = []
    _cov = None
    _inv_cov = None
    _use_free_cleft_model = False
    _logp_const = 0.0

    _pixwin_correction_gg = None
    _pixwin_correction_kg = None
    _pixwin_correction_nside = 2048
    _n_pca = {'Blue_ACT': 3, 'Green_ACT': 5, 'Blue_Planck': 3, 'Green_Planck': 5}

    def initialize(self):
        """
         Prepare any computation, importing any necessary code, files, etc.

         e.g. here we load some data file
        """

        if self.path is None:
            self.data_base_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../data")
        else:
            self.data_base_path = self.path

        # clear
        self._binning_function_gg = []
        self._binning_function_kg = []
        self._binning_function_kk = []

        self._ell_vals = None
        self._ell_vals_kk = None
        self._ell_selections = []
        self._ell_selections_kk = []
        self._ell_conds = []
        self._ell_conds_kk = []
        self._lmax_kk = 0
        self._lmax_TEB = 0
        self._data_gg = []
        self._data_kg = []
        self._data_kk = []
        self._cov = None
        self._inv_cov = None
        self._use_free_cleft_model = False

        self._pixwin_correction_gg = None
        self._pixwin_correction_kg = None

        data_filename_dict = yaml_load_file(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config_files/data_filenames.yaml"))
        if self.want_lensing_norm_correction:
            covmat_filename_dict = yaml_load_file(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config_files/covmat_filenames_baseline.yaml"))
        else:
            covmat_filename_dict = yaml_load_file(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config_files/covmat_filenames_cmbmarg.yaml"))

        binning_instructions = yaml_load_file(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config_files/binning_setup.yaml"))

        if len(self.samples) > 0:

            self.log.info("Getting binning functions ...")
            for s in self.samples:
                assert (s in binning_instructions.keys()), f"No instructions provided for how to bin sample {s}"
                instructions = binning_instructions[s]

                if 'transfer_path' in instructions.keys():
                    transfer_function = np.loadtxt(os.path.join(self.data_base_path, f"aux_data/transfer_functions/{instructions['transfer_path']}"))
                else:
                    transfer_function = None

                if 'bandwindow_matrix_path' in instructions.keys():
                    bandwindow_matrices = np.load(os.path.join(self.data_base_path, f"aux_data/bandwindow_matrices/{instructions['bandwindow_matrix_path']}"), allow_pickle=True).item()

                    ell_bin_edges = np.array(instructions['ell_bin_edges'][:len(transfer_function[:, 1]) + 1] if transfer_function is not None else instructions['ell_bin_edges'])
                    self._binning_function_gg.append(NaMasterPowerSpectrumBinning(bandwindow_matrices['gg']['coupling'], bandwindow_matrices['gg']['bandwindow'], ell_bin_edges, transfer_function=transfer_function[:, 1] if transfer_function is not None else None))
                    self._binning_function_kg.append(NaMasterPowerSpectrumBinning(bandwindow_matrices['kg']['coupling'], bandwindow_matrices['kg']['bandwindow'], ell_bin_edges, transfer_function=transfer_function[:, 2] if transfer_function is not None else None))
                elif 'binning_matrix_path' in instructions.keys():
                    binning_matrix = np.loadtxt(os.path.join(self.data_base_path, f"aux_data/bandwindow_matrices/{instructions['binning_matrix_path']}"))
                    self._binning_function_gg.append(MatrixPowerSpectrumBinning(binning_matrix, transfer_function=transfer_function[:, 1] if transfer_function is not None else None))
                    self._binning_function_kg.append(MatrixPowerSpectrumBinning(binning_matrix, transfer_function=transfer_function[:, 2] if transfer_function is not None else None))
                else:
                    ell_bin_edges = np.array(instructions['ell_bin_edges'][:len(transfer_function[:, 1]) + 1] if transfer_function is not None else instructions['ell_bin_edges'])
                    self._binning_function_gg.append(PowerSpectrumBinning(ell_bin_edges, transfer_function=transfer_function[:, 1] if transfer_function is not None else None))
                    self._binning_function_kg.append(PowerSpectrumBinning(ell_bin_edges, transfer_function=transfer_function[:, 2] if transfer_function is not None else None))

            self.log.info("Setting ell-spacing ...")
            all_ells = np.sort(np.unique(np.concatenate([func.get_input_ells() for func in self._binning_function_kg] + [func.get_input_ells() for func in self._binning_function_gg]).astype('int')))

            if not self.interpolate_cls:
                self._ell_vals = all_ells
            else:
                if isinstance(self._binning_function_kg[0], NaMasterPowerSpectrumBinning):
                    assert (isinstance(self._binning_function_gg[0], NaMasterPowerSpectrumBinning))
                    ell_space_weights = np.zeros(len(all_ells))

                    max_ell = np.max([self.lranges_gg[s][1] for s in self.samples] + [self.lranges_kg[s][1] for s in self.samples])

                    ell_space_weights[all_ells <= max_ell] = 1.0
                    ell_space_weights[max_ell < all_ells] = all_ells[max_ell < all_ells][0] / all_ells[max_ell < all_ells]

                    self._ell_vals = np.concatenate([all_ells[all_ells <= 30], np.round(density_space(all_ells[30 < all_ells], ell_space_weights[30 < all_ells], self.cl_interp_N - np.sum(all_ells <= 30), endpoint=True), decimals=0).astype('int')])
                else:
                    self._ell_vals = all_ells

            self._pixwin_correction_gg = hp.pixwin(self._pixwin_correction_nside) ** 2
            self._pixwin_correction_kg = hp.pixwin(self._pixwin_correction_nside)

        self.log.info("Loading data ...")
        for i, s in enumerate(self.samples):
            cond_gg = (self.lranges_gg[s][0] <= self._binning_function_gg[i].get_ell_bin_edges()[:-1]) & (self._binning_function_gg[i].get_ell_bin_edges()[1:] < self.lranges_gg[s][1])
            cond_kg = (self.lranges_kg[s][0] <= self._binning_function_kg[i].get_ell_bin_edges()[:-1]) & (self._binning_function_kg[i].get_ell_bin_edges()[1:] < self.lranges_kg[s][1])

            self._ell_conds.append((cond_gg, cond_kg))

            data = np.loadtxt(os.path.join(self.data_base_path, f"bandpowers/{data_filename_dict[s]}"))
            data_ells = np.round(data[:, 0], 10)
            ell_selection_gg = [data_ells[j] in np.round(self._binning_function_gg[i].get_binned_ell_vals()[cond_gg], 10) for j in range(len(data_ells))]
            ell_selection_kg = [data_ells[j] in np.round(self._binning_function_kg[i].get_binned_ell_vals()[cond_kg], 10) for j in range(len(data_ells))]

            self._ell_selections.append((ell_selection_gg, ell_selection_kg))

            self._data_gg.append(data[ell_selection_gg, 1])
            self._data_kg.append(data[ell_selection_kg, 3])

        if self.include_lensing_auto_spectrum:
            for i, s in enumerate(self.lensing_auto_spectrum_samples):
                binning_matrix = np.loadtxt(os.path.join(self.data_base_path, f"aux_data/bandwindow_matrices/{binning_instructions[s]['binning_matrix_path']}"))
                self._binning_function_kk.append(MatrixPowerSpectrumBinning(binning_matrix))

                self._ell_conds_kk.append((self.lranges_kk[s][0] < self._binning_function_kk[i].get_binned_ell_vals()) & (self._binning_function_kk[i].get_binned_ell_vals() < self.lranges_kk[s][1]))
                data = np.loadtxt(os.path.join(self.data_base_path, f"bandpowers/{data_filename_dict[s]}"))
                self._ell_selections_kk.append(self._ell_conds_kk[i])  # this could be changed in the future for the case that the data provided doesn't have the same dimensions as the output of the theory code
                self._data_kk.append(data[self._ell_selections_kk[i]])

            self._lmax_kk = np.max([self.lranges_kk[s][1] for s in self.lensing_auto_spectrum_samples])
            self._lmin_kk = np.min([self.lranges_kk[s][0] for s in self.lensing_auto_spectrum_samples])

            if self.want_lensing_norm_correction:
                self._lmax_kk = 3000
                self._lmin_kk = 0

            if not self.lensing_auto_spectrum_from_camb:
                all_ells_kk = np.linspace(self._lmin_kk, self._lmax_kk, self._lmax_kk - self._lmin_kk + 1).astype('int')
                ell_space_weights_kk = np.zeros_like(all_ells_kk, dtype='float')

                ell_space_weights_kk[(max(self._lmin_kk, 30) <= all_ells_kk) & (all_ells_kk <= self._lmax_kk)] = all_ells_kk[(max(self._lmin_kk, 30) <= all_ells_kk) & (all_ells_kk <= self._lmax_kk)][0] / all_ells_kk[(max(self._lmin_kk, 30) <= all_ells_kk) & (all_ells_kk <= self._lmax_kk)]

                self._ell_vals_kk = np.concatenate([all_ells_kk[all_ells_kk <= 30], np.round(density_space(all_ells_kk[30 < all_ells_kk], ell_space_weights_kk[30 < all_ells_kk], self.cl_interp_N - np.sum(all_ells_kk <= 30), endpoint=True), decimals=0).astype('int')])

        self.log.info("Loading covariances ...")
        n_sample_bpw_list = [np.sum(self._ell_selections[i][0]) + np.sum(self._ell_selections[i][1]) for i in range(len(self.samples))]
        if self.include_lensing_auto_spectrum:
            n_sample_bpw_list += [np.sum(self._ell_selections_kk[i]) for i in range(len(self.lensing_auto_spectrum_samples))]
        n_sample_bpw = np.array(n_sample_bpw_list, dtype=int)
        covmat = np.zeros((np.sum(n_sample_bpw), np.sum(n_sample_bpw)))

        for i, s in enumerate(self.samples):

            cov = np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[s]))
            ell_selection = np.concatenate(self._ell_selections[i])
            if cov.shape != (len(ell_selection), len(ell_selection)):
                warnings.warn(f"Covmat and input shape are in disagreement! Using first {len(self._ell_selections[i][0])} and {len(self._ell_selections[i][1])} bins for Clgg and Clkg respectively")

                self._ell_selections[i] = np.concatenate([self._ell_selections[i][0], np.full(cov.shape[0] // 2 - len(self._ell_selections[i][0]), False)]), np.concatenate([self._ell_selections[i][1], np.full(cov.shape[0] // 2 - len(self._ell_selections[i][1]), False)])
                ell_selection = np.concatenate(self._ell_selections[i])

            covmat[np.sum(n_sample_bpw[:i]):np.sum(n_sample_bpw[:i + 1]), np.sum(n_sample_bpw[:i]):np.sum(n_sample_bpw[:i + 1])] = select_from_matrix(cov, ell_selection)  # /h_factor

            for j, s2 in enumerate(self.samples[i + 1:]):
                j += i + 1

                if f"{s}_X_{s2}" in covmat_filename_dict.keys():
                    cross_cov = select_from_matrix(np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s}_X_{s2}"])), np.concatenate(self._ell_selections[i]), np.concatenate(self._ell_selections[j]))
                elif f"{s2}_X_{s}" in covmat_filename_dict.keys():
                    cross_cov = select_from_matrix(np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s2}_X_{s}"])), np.concatenate(self._ell_selections[j]), np.concatenate(self._ell_selections[i])).T
                else:
                    warnings.warn(f"Cross covariance between {s} and {s2} not provided.")
                    continue

                covmat[np.sum(n_sample_bpw[:i]):np.sum(n_sample_bpw[:i + 1]), np.sum(n_sample_bpw[:j]):np.sum(n_sample_bpw[:j + 1])] = cross_cov
                covmat[np.sum(n_sample_bpw[:j]):np.sum(n_sample_bpw[:j + 1]), np.sum(n_sample_bpw[:i]):np.sum(n_sample_bpw[:i + 1])] = cross_cov.T

        if self.include_lensing_auto_spectrum:
            for i, s in enumerate(self.lensing_auto_spectrum_samples):
                cov = np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s}_Clkk"]))

                covmat[np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)]), np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)])] = select_from_matrix(cov, self._ell_selections_kk[i])  # / h_factor

                for j, s2 in enumerate(self.samples):

                    if f"{s}_Clkk_X_{s2}" in covmat_filename_dict.keys():
                        cross_cov = select_from_matrix(np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s}_Clkk_X_{s2}"])), self._ell_selections_kk[i], np.concatenate(self._ell_selections[j]))
                    else:
                        warnings.warn(f"Cross covariance between clkk_{s} and {s2} not provided.")
                        continue

                    covmat[np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)]), np.sum(n_sample_bpw[:j]):np.sum(n_sample_bpw[:j + 1])] = cross_cov  # / h_factor
                    covmat[np.sum(n_sample_bpw[:j]):np.sum(n_sample_bpw[:j + 1]), np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)])] = cross_cov.T  # / h_factor

                for j, s2 in enumerate(self.lensing_auto_spectrum_samples[i + 1:]):
                    j += i + 1

                    if f"{s}_Clkk_X_{s2}_Clkk" in covmat_filename_dict.keys():
                        cross_cov = select_from_matrix(np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s}_Clkk_X_{s2}_Clkk"])), self._ell_selections_kk[i], self._ell_selections_kk[j])
                    elif f"{s2}_Clkk_X_{s}_Clkk" in covmat_filename_dict.keys():
                        cross_cov = select_from_matrix(np.loadtxt(os.path.join(self.data_base_path, covmat_filename_dict[f"{s2}_Clkk_X_{s}_Clkk"])), self._ell_selections_kk[j], self._ell_selections_kk[i]).T
                    else:
                        warnings.warn(f"Cross covariance between {s} and {s2} not provided.")
                        continue

                    covmat[np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)]), np.sum(n_sample_bpw[:j + len(self.samples)]):np.sum(n_sample_bpw[:j + 1 + len(self.samples)])] = cross_cov  # / h_factor
                    covmat[np.sum(n_sample_bpw[:j + len(self.samples)]):np.sum(n_sample_bpw[:j + 1 + len(self.samples)]), np.sum(n_sample_bpw[:i + len(self.samples)]):np.sum(n_sample_bpw[:i + 1 + len(self.samples)])] = cross_cov.T  # / h_factor

        if self.apply_hartlap_correction:
            if 'Blue_ACT' in self.samples or 'Green_ACT' in self.samples:
                nsims_min = 400
            elif 'Blue_Planck' in self.samples or 'Green_Planck' in self.samples or ('Planck' in self.lensing_auto_spectrum_samples and self.include_lensing_auto_spectrum):
                nsims_min = 480
            else:
                nsims_min = 798

            h_factor = (nsims_min - np.sum(n_sample_bpw) - 2) / (nsims_min - 1)
            covmat /= h_factor
            self.log.info(f"Inflating covariance by hartlap correction: {1 / h_factor:.3f}")

        self._cov = covmat
        self._inv_cov = np.linalg.inv(covmat)

        # self._logp_const = -(np.log(2.0 * np.pi) * np.sum(n_sample_bpw) + np.linalg.slogdet(self._cov)[1]) / 2.0

        if self.theory_eval_kwargs is None:
            self.theory_eval_kwargs = {}

        self._use_free_cleft_model = np.any([self.scale_cleft_b2, self.shift_cleft_b2, self.scale_cleft_bs, self.shift_cleft_bs, self.scale_cleft_b3, self.shift_cleft_b3])

    def get_requirements(self):
        """
         return dictionary specifying quantities calculated by a theory code are needed
        """
        if len(self.samples)>0:
            reqs = {'cls':{'ell_vals':self._ell_vals, 'samples': np.unique([s.split("_")[0] for s in self.samples])}}
        else:
            reqs = {}

        if self.include_lensing_auto_spectrum:
            if self.lensing_auto_spectrum_from_camb:
                reqs['Cl'] = {'pp': self._lmax_kk}
            else:
                reqs['clkk'] = {'ell_vals': self._ell_vals_kk}

        if self.want_lensing_norm_correction:
            reqs['lensing_norm_correction'] = None

        return reqs

    def logp(self, **params_values):
        """
        Taking a dictionary of (sampled) nuisance parameter values params_values
        and return a log-likelihood.

        e.g. here we calculate chi^2
        """

        bias_list = [params_values['b_{}'.format(s)] for s in self.samples]
        n_shot_list = [params_values['n_shot_{}'.format(s)] for s in self.samples]
        s_mag_list = [params_values['s_{}'.format(s)] for s in self.samples]

        pca_coeff_list = None
        if self.do_pca_dndz_marg:
            pca_coeff_list = [np.array([params_values['{}_pca_{}'.format(s,i)] for i in range(self._n_pca[s])]) for s in self.samples]

        if self._use_free_cleft_model:
            cleft_coeff = []
            for i,s in enumerate(self.samples):
                coeff_tmp = {}
                for bX in ['b2', 'bs', 'b3']:
                    scale, shift = 1.0, 0.0
                    if getattr(self, f'scale_cleft_{bX}'):
                        scale = params_values[f'{s}_scale_cleft_c{bX}']
                    if getattr(self, f'shift_cleft_{bX}'):
                        shift = params_values[f'{s}_shift_cleft_c{bX}']

                    coeff_tmp[f"c{bX}"] = np.concatenate((np.atleast_1d(scale), np.atleast_1d(shift)))
                cleft_coeff.append(coeff_tmp)
        else:
            cleft_coeff = None

        residuals = []
        if len(self.samples) > 0:
            gg, kg = self.provider.get_cls(ell=self._ell_vals, bias=bias_list, s_mag=s_mag_list, pca_coeff=pca_coeff_list, cleft_coeff=cleft_coeff, samples=[s.split("_")[0] for s in self.samples], do_dndz_pca=self.do_pca_dndz_marg, **self.theory_eval_kwargs)

            for i, s in enumerate(self.samples):
                binned_gg = self._binning_function_gg[i](np.interp(self._binning_function_gg[i].get_input_ells(), self._ell_vals, gg[i]) * self._pixwin_correction_gg[self._binning_function_gg[i].get_input_ells()], white_noise=n_shot_list[i])[self._ell_conds[i][0]]

                all_ell_kg = np.interp(self._binning_function_kg[i].get_input_ells(), self._ell_vals, kg[i])
                if self.want_lensing_norm_correction:
                    all_ell_kg = self.provider.get_lensing_norm_correction(all_ell_kg, s, cross_spectrum=True)

                binned_kg = self._binning_function_kg[i](all_ell_kg * self.pixwin_correction_kg[self._binning_function_kg[i].get_input_ells()])[self._ell_conds[i][1]]

                residual_gg = self._data_gg[i] - binned_gg
                residual_kg = self._data_kg[i] - binned_kg
                residuals.append(np.concatenate([residual_gg, residual_kg]))

        if self.include_lensing_auto_spectrum:
            if self.lensing_auto_spectrum_from_camb:
                cl = self.provider.get_Cl(ell_factor=False, units='FIRASmuK2')
                clkk = cl['pp'] * (cl['ell'] * (cl['ell'] + 1.)) ** 2. / 4.
                all_ells = cl['ell']
            else:
                clkk = self.provider.get_clkk(ell=self._ell_vals_kk)
                all_ells = np.linspace(self._lmin_kk, self._lmax_kk, self._lmax_kk - self._lmin_kk + 1).astype('int')
                clkk = np.interp(all_ells, self._ell_vals_kk, clkk, left=0.0, right=0.0)

            for i,s in enumerate(self.lensing_auto_spectrum_samples):

                if self.want_lensing_norm_correction:
                    clkk_corrected = np.copy(clkk)
                    clkk_corrected = self.provider.get_lensing_norm_correction(clkk_corrected, s, cross_spectrum=False)
                else:
                    clkk_corrected = clkk

                clkk2bin = np.zeros_like(self._binning_function_kk[i].get_input_ells())
                sel_ells = np.where(np.isclose(self._binning_function_kk[i].get_input_ells()[:, None], all_ells[None, :]))
                clkk2bin[sel_ells[0]] = clkk_corrected[sel_ells[1]]
                binned_kk = self._binning_function_kk[i](clkk2bin)[self._ell_conds_kk[i]]
                residuals.append(self._data_kk[i] - binned_kk)

        chi2 = np.concatenate(residuals) @ self._inv_cov @ np.concatenate(residuals)

        return -chi2 / 2.0 + self._logp_const
