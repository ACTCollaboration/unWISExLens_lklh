import numpy as np
from .model_helpers_unWISExLens import cosmo_from_camb, dNdz
from ..auxiliary.pk_interpolator import PowerSpectrumInterpolator
from .unWISExkappa_model import unWISExLens_theory_model as base_model


class unWISExLens_theory_model(base_model):

    def __init__(self, want_poly_marg=False, *args, **kwargs):

        self._poly_marg = want_poly_marg
        super().__init__(*args, **kwargs)
    def compute_raw_spectra(self, cosmo, dNdz_list, pk_weyl_weyl, pk_weyl_dnonu, pk_dnonu_dnonu, cleft_interpolations_dtot_dnonu=None, cleft_interpolations_dnonu_dnonu=None, fid_bias_evol_list=None):
        """
        Function to compute raw spectral components for :math:`C_\ell^{gg}` and :math:`C_\ell^{\kappa g}`. Expects to be given input power spectra and cosmology and performs projection integrals.

        Parameters
        ----------
        cosmo : cosmo_from_camb
            object to hold cosmology objects
        dNdz_list : list [dNdz]
            list of dNdz objects for galaxy kernels
        pk_weyl_weyl : PowerSpectrumInterpolator
            non-linear Weyl potential power spectrum (including neutrinos). Used for lensing contributions (:math:`C_\ell^{\kappa \mu}` and :math:`C_\ell^{\mu \mu}`).
        pk_weyl_dnonu : PowerSpectrumInterpolator
            non-linear cross power spectrum of Weyl potential and matter density without neutrinos. Used for galaxy - lensing contributions (:math:`C_\ell^{\kappa g}` and :math:`C_\ell^{g \mu}`).
        pk_dnonu_dnonu : PowerSpectrumInterpolator
            non-linear power spectrum of the matter density without neutrinos. Used for galaxy - galaxy contribution (:math:`C_\ell^{g g}`).
        cleft_interpolations_dtot_dnonu : list [PowerSpectrumInterpolator]
            interpolations of higher order contributions to cross power spectrum of total matter density and matter density without neutrinos. If None, no higher order contributions are used
        cleft_interpolations_dnonu_dnonu :  list [PowerSpectrumInterpolator]
            interpolations of higher order contributions to power spectrum of matter density without neutrinos. If None, no higher order contributions are used
        fid_bias_evol_list : list [callable], optional
            List of fiducial bias evolutions

        Returns
        ----------
        list [dict]
            The function returns a list of dictionaries of the raw components of :math:`C_\ell^{\kappa g}` and :math:`C_\ell^{gg}`. Pass this list to :func:`~unWISExACT_theory_model.evaluate` along with the appropriate nuisance parameters to evaluate the final spectra
        """

        if cleft_interpolations_dtot_dnonu is not None:
            assert (self._cleft_interp_helper is not None and fid_bias_evol_list is not None and cleft_interpolations_dnonu_dnonu is not None)

        chi_min, chi_max = cosmo.chi(self._zmin), cosmo.chi(self._zmax)
        chi_vals = (chi_max - chi_min) / 2 * self._gauss_x + (chi_max + chi_min) / 2
        z_vals = cosmo.z_of_chi(chi_vals)

        hubble_vals = cosmo.H(z_vals)

        muKernel_evals = np.zeros((len(dNdz_list), len(chi_vals)))
        for i, dndz in enumerate(dNdz_list):
            muKernel_evals[i] = self._lensing_magnification_weights(chi_vals, cosmo, dndz.dNdz)

        kappaKernel_eval = self._kappa_kernel(chi_vals, z_vals, cosmo)

        kGrid = (self._ell_vals[None, :] + 0.5) / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None]

        matter2weyl_factor = 3/2 * cosmo.Omega_m * cosmo.H0**2 * (1+z_vals)[:,None] * kGrid**2/(3*cosmo.curvature - kGrid**2)

        halofit_pk_evals_weyl_weyl = pk_weyl_weyl.P(z_vals[:, None], kGrid, grid=False) / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None] ** 2
        halofit_pk_evals_weyl_dnonu = pk_weyl_dnonu.P(z_vals[:, None], kGrid, grid=False) / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None] ** 2
        halofit_pk_evals_dnonu_dnonu = pk_dnonu_dnonu.P(z_vals[:, None], kGrid, grid=False) / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None] ** 2
        if cleft_interpolations_dtot_dnonu is not None:
            cleft_pk_evals_weyl_dnonu = np.zeros((*halofit_pk_evals_weyl_dnonu.shape, len(cleft_interpolations_dnonu_dnonu)))
            cleft_pk_evals_dnonu_dnonu = np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, len(cleft_interpolations_dnonu_dnonu)))
            for i in range(len(cleft_interpolations_dtot_dnonu)):
                cleft_pk_evals_weyl_dnonu[:, :, i] = matter2weyl_factor * cleft_interpolations_dtot_dnonu[i].piecewise_evaluate(z_vals[:, None], kGrid / cosmo.h, pad=(0.0, 0.0)) / cosmo.h ** 3 / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None] ** 2
                cleft_pk_evals_dnonu_dnonu[:, :, i] = cleft_interpolations_dnonu_dnonu[i].piecewise_evaluate(z_vals[:, None], kGrid / cosmo.h, pad=(0.0, 0.0)) / cosmo.h ** 3 / cosmo.comoving_angular_diameter_distance(chi_vals)[:, None] ** 2

        if self._cross_correlation_redshift_correction is not None:
            cross_dndz_correction = self._cross_correlation_redshift_correction(cosmo, pk_dnonu_dnonu, z_vals)
        else:
            cross_dndz_correction = np.ones_like(z_vals)

        outputs = []
        outputs_cross = {}
        # the indexing below is as follows (z, l, pcs, cleft_pk)
        for i in range(len(dNdz_list)):

            bdN_dz_H = dNdz_list[i].bdNdz(z_vals, pcs=True) * cross_dndz_correction[:, None] * hubble_vals[:, None]
            dN_dz_H = dNdz_list[i].dNdz(z_vals) * hubble_vals

            #recompute norm after correction
            bdndz_norm = np.sum(bdN_dz_H * self._gauss_w[:, None], axis=0) * (chi_max - chi_min) / 2

            if cleft_interpolations_dtot_dnonu is not None:
                cleft_pk_gg_b, cleft_pk_gg_nob1, cleft_pk_gg_nob2, dump = self._cleft_interp_helper.multiply_eft_biases(cleft_pk_evals_dnonu_dnonu, z_vals, fid_bias_evol_list[i])
                dump, dump, dump, cleft_pk_mg = self._cleft_interp_helper.multiply_eft_biases(cleft_pk_evals_weyl_dnonu, z_vals, fid_bias_evol_list[i])
            else:
                cleft_pk_gg_b, cleft_pk_gg_nob1, cleft_pk_gg_nob2, cleft_pk_mg = np.zeros((*halofit_pk_evals_weyl_weyl.shape, 1, 1)), np.zeros((*halofit_pk_evals_weyl_weyl.shape, 1, 1)), np.zeros((*halofit_pk_evals_weyl_weyl.shape, 1, 1)), np.zeros((*halofit_pk_evals_weyl_weyl.shape, 1, 1))

            p_mg_b = bdN_dz_H[:, None, :] * halofit_pk_evals_weyl_dnonu[:, :, None]
            p_mg_nob = dN_dz_H[:, None, None, None] * cleft_pk_mg

            kg_spectra_b = np.nansum(p_mg_b * kappaKernel_eval[:, None, None] * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
            kg_spectra_nob = np.nansum(p_mg_nob * kappaKernel_eval[:, None, None, None] * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2
            kmu_spectra = np.nansum(muKernel_evals[i][:, None] * halofit_pk_evals_weyl_weyl * kappaKernel_eval[:, None] * self._gauss_w[:, None], axis=0) * (chi_max - chi_min) / 2

            p_gg_bsq = (bdN_dz_H[:, None, :, None] * bdN_dz_H[:, None, None, :] * halofit_pk_evals_dnonu_dnonu[:, :, None, None]).reshape(len(z_vals), len(self._ell_vals), (dNdz_list[i].n_pcs + 1)**2)
            p_gg_b = bdN_dz_H[:, None, :, None, None] * dN_dz_H[:, None, None, None, None] * cleft_pk_gg_b[:, :, None, :, :]
            p_gg_nob1 = (dN_dz_H**2)[:, None, None, None] * cleft_pk_gg_nob1
            p_gg_nob2 = (dN_dz_H**2)[:, None, None, None] * cleft_pk_gg_nob2

            gg_spectra_bsq = np.nansum(p_gg_bsq * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
            gg_spectra_b = np.nansum(p_gg_b * self._gauss_w[:, None, None, None, None], axis=0) * (chi_max - chi_min) / 2
            gg_spectra_nob1 = np.nansum(p_gg_nob1 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2
            gg_spectra_nob2 = np.nansum(p_gg_nob2 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2

            mumu_spectra = np.nansum((muKernel_evals[i]**2)[:, None] * halofit_pk_evals_weyl_weyl * self._gauss_w[:, None], axis=0) * (chi_max - chi_min) / 2
            gmu_spectra_b = np.nansum(muKernel_evals[i][:, None, None] * p_mg_b * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
            gmu_spectra_nob = np.nansum(muKernel_evals[i][:, None, None, None] * p_mg_nob * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2

            outputs.append({'kg': {"kg_b": kg_spectra_b, "kg_nob": kg_spectra_nob, "kmu": kmu_spectra},
                            'gg': {"gg_bsq": gg_spectra_bsq, "gg_b": gg_spectra_b, "gg_nob1": gg_spectra_nob1, "gg_nob2": gg_spectra_nob2, "gmu_b": gmu_spectra_b, "gmu_nob": gmu_spectra_nob, "mumu": mumu_spectra},
                            'bdndz_norm': bdndz_norm})

            if self._want_gg_cross:
                for j in range(i+1,len(dNdz_list)):
                    bdN_dz_H_j = dNdz_list[j].bdNdz(z_vals, pcs=True) * cross_dndz_correction[:, None] * hubble_vals[:, None]
                    dN_dz_H_j = dNdz_list[j].dNdz(z_vals) * hubble_vals

                    if cleft_interpolations_dtot_dnonu is not None:
                        cleft_pk_g1g2_b_cross, cleft_pk_g1g2_nob1_cross, cleft_pk_g1g2_nob2_cross, _ = self._cleft_interp_helper.multiply_eft_biases(cleft_pk_evals_dnonu_dnonu, z_vals, fid_bias_evol_list[i], fid_bias_evol_list[j])
                        cleft_pk_g2g1_b_cross, cleft_pk_g2g1_nob1_cross, cleft_pk_g2g1_nob2_cross, _ = self._cleft_interp_helper.multiply_eft_biases(cleft_pk_evals_dnonu_dnonu, z_vals, fid_bias_evol_list[j], fid_bias_evol_list[i])
                        _, _, _, cleft_pk_mg2 = self._cleft_interp_helper.multiply_eft_biases(cleft_pk_evals_weyl_dnonu, z_vals, fid_bias_evol_list[j], fid_bias_evol_list[i])
                    else:
                        cleft_pk_g1g2_b_cross, cleft_pk_g1g2_nob1_cross, cleft_pk_g1g2_nob2_cross = np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1)), np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1)), np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1))
                        cleft_pk_g2g1_b_cross, cleft_pk_g2g1_nob1_cross, cleft_pk_g2g1_nob2_cross, cleft_pk_mg2 = np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1)), np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1)), np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1)), np.zeros((*halofit_pk_evals_dnonu_dnonu.shape, 1, 1))

                    p_g1b_g2b = (bdN_dz_H_j[:, None, :, None] * p_mg_b[:, :, None, :]).reshape(len(z_vals), len(self._ell_vals), (dNdz_list[i].n_pcs + 1) * (dNdz_list[j].n_pcs + 1))
                    p_g1b_g2 = bdN_dz_H[:, None, :, None, None] * dN_dz_H_j[:, None, None, None, None] * cleft_pk_g2g1_b_cross[:, :, None, :, :]
                    p_g2b_g1 = bdN_dz_H_j[:, None, :, None, None] * dN_dz_H[:, None, None, None, None] * cleft_pk_g1g2_b_cross[:, :, None, :, :]
                    p_g1_g2_1 = (dN_dz_H * dN_dz_H_j)[:, None, None, None] * cleft_pk_g1g2_nob1_cross
                    p_g1_g2_2 = (dN_dz_H * dN_dz_H_j)[:, None, None, None] * cleft_pk_g1g2_nob2_cross
                    p_g2_g1_1 = (dN_dz_H * dN_dz_H_j)[:, None, None, None] * cleft_pk_g2g1_nob1_cross
                    p_g2_g1_2 = (dN_dz_H * dN_dz_H_j)[:, None, None, None] * cleft_pk_g2g1_nob2_cross

                    g1b_g2b_spectra = np.nansum(p_g1b_g2b * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
                    g1b_g2_spectra = np.nansum(p_g1b_g2 * self._gauss_w[:, None, None, None, None], axis=0) * (chi_max - chi_min) / 2
                    g2b_g1_spectra = np.nansum(p_g2b_g1 * self._gauss_w[:, None, None, None, None], axis=0) * (chi_max - chi_min) / 2
                    g1_g2_1_spectra = np.nansum(p_g1_g2_1 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2
                    g1_g2_2_spectra = np.nansum(p_g1_g2_2 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2
                    g2_g1_1_spectra = np.nansum(p_g2_g1_1 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2
                    g2_g1_2_spectra = np.nansum(p_g2_g1_2 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2

                    mu1mu2_spectra = np.nansum((muKernel_evals[i] * muKernel_evals[j])[:, None] * halofit_pk_evals_weyl_weyl * self._gauss_w[:, None], axis=0) * (chi_max - chi_min) / 2
                    g1mu2_spectra = np.nansum(muKernel_evals[j][:, None, None] * p_mg_b * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
                    g1mu2_spectra_nob = np.nansum(muKernel_evals[j][:, None, None, None] * p_mg_nob * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2

                    g2mu1_spectra = np.nansum(muKernel_evals[i][:, None, None] * bdN_dz_H_j[:, None, :] * halofit_pk_evals_weyl_dnonu[:, :, None] * self._gauss_w[:, None, None], axis=0) * (chi_max - chi_min) / 2
                    g2mu1_spectra_nob = np.nansum(muKernel_evals[i][:, None, None, None] * dN_dz_H_j[:, None, None, None] * cleft_pk_mg2 * self._gauss_w[:, None, None, None], axis=0) * (chi_max - chi_min) / 2

                    outputs_cross[(i, j)] = {'g1g2': {'g1b_g2b': g1b_g2b_spectra,
                                                      'g1b_g2': g1b_g2_spectra,
                                                      'g2b_g1': g2b_g1_spectra,
                                                      'g1_g2_1': g1_g2_1_spectra,
                                                      'g1_g2_2': g1_g2_2_spectra,
                                                      'g2_g1_1': g2_g1_1_spectra,
                                                      'g2_g1_2': g2_g1_2_spectra},
                                             'g1mu2': {'gmu_b': g1mu2_spectra, 'gmu_nob': g1mu2_spectra_nob},
                                             'g2mu1': {'gmu_b': g2mu1_spectra, 'gmu_nob': g2mu1_spectra_nob},
                                             'mu1mu2': mu1mu2_spectra}

        if self._want_gg_cross:
            return outputs, outputs_cross
        else:
            return outputs

    def evaluate(self, raw_spectra, get='all', want_gg_cross=False, noise_bias=None, cross_noise_bias=None, bias=None, s_mag=None, pca_coeff=None, do_dndz_pca=True, cleft_coeff=None, **kwargs):

        kg_kwargs = {}
        gg_kwargs = {}
        for key in kwargs.keys():
            if 'cl_kg' in key or 'cl_kmu' in key:
                kg_kwargs[key] = kwargs[key]
            elif 'cl_gg' in key or 'cl_gmu' in key or 'cl_mumu' in key:
                gg_kwargs[key] = kwargs[key]

        if self._want_gg_cross:
            raw_spectra, raw_spectra_cross = raw_spectra

        if noise_bias is None:
            noise_bias = [{'kg': {"kg_b": np.zeros_like(self._ell_vals)},
                           'gg': {"gg_bsq": np.zeros_like(self._ell_vals),
                                  "gg_b": np.zeros((len(self._ell_vals), 1, 1)) if self._cleft_interp_helper is None else np.zeros((len(self._ell_vals), *self._cleft_interp_helper.assemble_cleft_coeff()[0].shape)),
                                  "gmu_b": np.zeros_like(self._ell_vals)}} for i in range(len(raw_spectra))]

        out_gg = np.empty((len(raw_spectra), len(self._ell_vals)))
        out_kg = np.empty((len(raw_spectra), len(self._ell_vals)))

        if want_gg_cross:
            assert self._want_gg_cross
            out_gg_cross = np.empty((len(raw_spectra), len(raw_spectra), len(self._ell_vals)))

        for i in range(len(raw_spectra)):
            b = bias[i] if bias is not None else 1.0
            s = s_mag[i] if s_mag is not None else 0.4
            gg_kwargs['b_gmu'] = kwargs['bias_gmu'] if 'bias_gmu' in kwargs.keys() else None

            n_pcs = raw_spectra[i]['kg']['kg_b'].shape[-1]-1
            #first is fiducial dndz, next mean delta_dndz and then pcs
            if pca_coeff is None or not do_dndz_pca or n_pcs == 0:
                pca_coeff_final = np.concatenate([[1.0], np.zeros(n_pcs)])
            else:
                assert(len(pca_coeff[i])==n_pcs - 1)
                pca_coeff_final = np.concatenate([[1.0, 1.0], pca_coeff[i]])

            if self._cleft_interp_helper is not None:
                cleft_coeff_final1, cleft_coeff_final2 = self._cleft_interp_helper.assemble_cleft_coeff() if cleft_coeff is None else self._cleft_interp_helper.assemble_cleft_coeff(**cleft_coeff[i])
            else:
                cleft_coeff_final1, cleft_coeff_final2 = np.zeros((1,1)), np.zeros((1,1))

            b /= np.dot(raw_spectra[i]['bdndz_norm'], pca_coeff_final)

            if get=='gg' or get=='all':
                out_gg[i] = self.__gg(raw_spectra[i], b, s, cleft_coeff_final1, cleft_coeff_final2, pca_coeff_final, noise_bias[i], **gg_kwargs)

                if want_gg_cross:
                    for j in range(i+1, len(raw_spectra)):
                        b2 = bias[j] if bias is not None else 1.0
                        s2 = s_mag[j] if s_mag is not None else 0.4

                        n_pcs2 = raw_spectra[j]['kg']['kg_b'].shape[-1] - 1
                        # first is fiducial dndz, next mean delta_dndz and then pcs
                        if pca_coeff is None or not do_dndz_pca or n_pcs2 == 0:
                            pca_coeff_final2 = np.concatenate([[1.0], np.zeros(n_pcs2)])
                        else:
                            assert (len(pca_coeff[j]) == n_pcs2 - 1)
                            pca_coeff_final2 = np.concatenate([[1.0, 1.0], pca_coeff[j]])

                        if self._cleft_interp_helper is not None:
                            cleft_coeff_final12_1, cleft_coeff_final12_2 = self._cleft_interp_helper.assemble_cleft_coeff() if cleft_coeff is None else self._cleft_interp_helper.assemble_cleft_coeff(**cleft_coeff[i], **dict([(f"{bX}_2", cleft_coeff[j][bX]) for bX in  cleft_coeff[j].keys()]))
                            cleft_coeff_final21_1, cleft_coeff_final21_2 = self._cleft_interp_helper.assemble_cleft_coeff() if cleft_coeff is None else self._cleft_interp_helper.assemble_cleft_coeff(**cleft_coeff[j], **dict([(f"{bX}_2", cleft_coeff[i][bX]) for bX in cleft_coeff[i].keys()]))
                        else:
                            cleft_coeff_final12_1, cleft_coeff_final12_2 = np.zeros((1, 1)), np.zeros((1, 1))
                            cleft_coeff_final21_1, cleft_coeff_final21_2 = np.zeros((1, 1)), np.zeros((1, 1))

                        b2 /= np.dot(raw_spectra[j]['bdndz_norm'], pca_coeff_final2)

                        mf_cross = cross_noise_bias[(i, j)] if cross_noise_bias is not None else {'g1g2': {"g1b_g2b": np.zeros_like(self._ell_vals),
                                                                                                           "g1b_g2": np.zeros((len(self._ell_vals), 1, 1)) if self._cleft_interp_helper is None else np.zeros((len(self._ell_vals), *self._cleft_interp_helper.assemble_cleft_coeff()[0].shape)),
                                                                                                           "g2b_g1": np.zeros((len(self._ell_vals), 1, 1)) if self._cleft_interp_helper is None else np.zeros((len(self._ell_vals), *self._cleft_interp_helper.assemble_cleft_coeff()[0].shape))},
                                                                                                  "g1mu2": {'gmu_b': np.zeros_like(self._ell_vals)},
                                                                                                  "g2mu1": {'gmu_b': np.zeros_like(self._ell_vals)}}

                        out_gg_cross[i,j] = self.__gg_cross(raw_spectra_cross[(i,j)], b, b2, s, s2, pca_coeff_final, pca_coeff_final2, cleft_coeff_final12_1, cleft_coeff_final12_2, cleft_coeff_final21_1, cleft_coeff_final21_2, mf_cross)
                        out_gg_cross[j,i] = out_gg_cross[i,j]

            if get=='kg' or get=='all':
                out_kg[i] = self.__kg(raw_spectra[i], b, s, cleft_coeff_final1, pca_coeff_final, noise_bias[i], **kg_kwargs)

        if get=='all':
            if want_gg_cross:
                return out_gg, out_gg_cross, out_kg
            else:
                return out_gg, out_kg
        elif get=='gg':
            if want_gg_cross:
                return out_gg, out_gg_cross
            else:
                return out_gg
        elif get=='kg':
            return out_kg
        else:
            raise ValueError(f"Get options are 'all', 'gg' or 'kg', got unsupported option {get}.")

    @staticmethod
    def __kg(raw_spectra, b, s, cleft_coeff, pca_coeff, noise_bias, cl_kg=True, cl_kg_HF=None, cl_kg_CLEFT=None, cl_kmu=True):
        if cl_kg_HF is None:
            cl_kg_HF = cl_kg
        if cl_kg_CLEFT is None:
            cl_kg_CLEFT = cl_kg

        cleft_terms = np.einsum('lij, ij->l', raw_spectra['kg']['kg_nob'], cleft_coeff)

        return cl_kg_HF * (np.dot(raw_spectra['kg']['kg_b'], pca_coeff) - noise_bias['kg']['kg_b']) * b \
               + cl_kg_CLEFT * cleft_terms \
               + cl_kmu * raw_spectra['kg']['kmu'] * (5 * s - 2)

    @staticmethod
    def __gg(raw_spectra, b, s, cleft_coeff1, cleft_coeff2, pca_coeff, noise_bias, b_gmu=None, cl_gg=True, cl_gmu=True, cl_mumu=True, cl_gg_bsq=None, cl_gg_b=None, cl_gg_nob=None, cl_gmu_b=None, cl_gmu_nob=None):
        if cl_gg_bsq is None:
            cl_gg_bsq = cl_gg
        if cl_gg_b is None:
            cl_gg_b = cl_gg
        if cl_gg_nob is None:
            cl_gg_nob = cl_gg

        if cl_gmu_b is None:
            cl_gmu_b = cl_gmu
        if cl_gmu_nob is None:
            cl_gmu_nob = cl_gmu

        if b_gmu is None:
            b_gmu = b

        pca_coeff_final_sq = np.outer(pca_coeff, pca_coeff).flatten()

        cleft_terms_gg_b = np.einsum('lpij, ij->lp', raw_spectra['gg']['gg_b'], cleft_coeff1)
        cleft_terms_gg_nob = np.einsum('lij, ij->l', raw_spectra['gg']['gg_nob1'], cleft_coeff1) + np.einsum('lij, ij->l', raw_spectra['gg']['gg_nob2'], cleft_coeff2)
        cleft_terms_gmu_nob = np.einsum('lij, ij->l', raw_spectra['gg']['gmu_nob'], cleft_coeff1)

        return cl_gg_bsq * (np.dot(raw_spectra['gg']['gg_bsq'], pca_coeff_final_sq) - noise_bias['gg']['gg_bsq']) * b**2 \
               + cl_gg_b * (np.dot(cleft_terms_gg_b, pca_coeff) - np.einsum('lij, ij->l', noise_bias['gg']['gg_b'], cleft_coeff1)) * b \
               + cl_gg_nob * cleft_terms_gg_nob \
               + cl_gmu_b * 2 * (np.dot(raw_spectra['gg']['gmu_b'], pca_coeff) - noise_bias['gg']['gmu_b']) * b_gmu * (5 * s - 2) \
               + cl_gmu_nob * 2 * cleft_terms_gmu_nob * (5 * s - 2) \
               + cl_mumu * raw_spectra['gg']['mumu'] * (5 * s - 2)**2

    @staticmethod
    def __gg_cross(raw_spectra, b1, b2, s1, s2, pca_coeff1, pca_coeff2, cleft_coeff12_1, cleft_coeff12_2, cleft_coeff21_1, cleft_coeff21_2, noise_bias, cl_gg=True, cl_gmu=True, cl_mumu=True):

        pca_coeff_final_sq = np.outer(pca_coeff1, pca_coeff2).flatten()

        cleft_terms_g1b_g2 = np.einsum('lpij, ij->lp', raw_spectra['g1g2']['g1b_g2'], cleft_coeff21_1)
        cleft_terms_g2b_g1 = np.einsum('lpij, ij->lp', raw_spectra['g1g2']['g2b_g1'], cleft_coeff12_1)
        cleft_terms_g1_g2 = np.einsum('lij, ij->l', raw_spectra['g1g2']['g1_g2_1'], cleft_coeff12_1) + np.einsum('lij, ij->l', raw_spectra['g1g2']['g1_g2_2'], cleft_coeff12_2)
        cleft_terms_g2_g1 = np.einsum('lij, ij->l', raw_spectra['g1g2']['g2_g1_1'], cleft_coeff21_1) + np.einsum('lij, ij->l', raw_spectra['g1g2']['g2_g1_2'], cleft_coeff21_2)
        cleft_terms_g1_mu2 = np.einsum('lij, ij->l', raw_spectra['g1mu2']['gmu_nob'], cleft_coeff12_1)
        cleft_terms_g2_mu1 = np.einsum('lij, ij->l', raw_spectra['g2mu1']['gmu_nob'], cleft_coeff21_1)

        return cl_gg * ((np.dot(raw_spectra['g1g2']['g1b_g2b'], pca_coeff_final_sq) - noise_bias['g1g2']['g1b_g2b']) * b1 * b2
                        + (np.dot(cleft_terms_g1b_g2, pca_coeff1) - np.einsum('lij, ij->l', noise_bias['g1g2']['g1b_g2'], cleft_coeff21_1)) * b1 / 2
                        + (np.dot(cleft_terms_g2b_g1, pca_coeff2) - np.einsum('lij, ij->l', noise_bias['g1g2']['g2b_g1'], cleft_coeff12_1)) * b2 / 2
                        + cleft_terms_g1_g2 / 2 + cleft_terms_g2_g1 / 2) \
            + cl_gmu * ((np.dot(raw_spectra['g1mu2']['gmu_b'], pca_coeff1) - noise_bias['g1mu2']['gmu_b']) * b1 * (5 * s2 - 2)
                        + (np.dot(raw_spectra['g2mu1']['gmu_b'], pca_coeff2) - noise_bias['g2mu1']['gmu_b']) * b2 * (5 * s1 - 2)
                        + cleft_terms_g1_mu2 * (5 * s2 - 2) + cleft_terms_g2_mu1 * (5 * s1 - 2)) \
            + cl_mumu * raw_spectra['mu1mu2'] * (5 * s1 - 2) * (5 * s2 - 2)

