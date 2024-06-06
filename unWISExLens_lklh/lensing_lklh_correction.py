from cobaya.theory import Theory
import numpy as np
import os


class LensingLklhCorrectionHelper(object):

    def __init__(self, clTEB_fid, clkk_fid, norm_fid, dNorm_dCl, dN1_dCl, dN1_dClkk):

        self.__clTEB_fid = clTEB_fid
        self.__clkk_fid = clkk_fid
        self.__norm_fid = norm_fid
        self.__dNorm_dCl = dNorm_dCl
        self.__dN1_dCl = dN1_dCl
        self.__dN1_dClkk = dN1_dClkk

        self.__lmax = self.__dNorm_dCl.shape[-1] - 1
        self.__Lmax = self.__dNorm_dCl.shape[1] - 1
        self.__Lmax_N1corr = self.__dN1_dClkk.shape[0] - 1

        assert self.__Lmax >= self.__dN1_dClkk.shape[0] - 1, "Lmax of dN1_dClkk is greater than the maximum Lmax of the fiducial Cls."

    def compute_corrections(self, clTEB):
        clTEB_diff = clTEB[:,:self.__lmax+1] - self.__clTEB_fid[:,:self.__lmax+1]
        with np.errstate(divide='ignore', invalid='ignore'):
            norm_correction = - np.nan_to_num(np.einsum('tLl,tl->L', self.__dNorm_dCl, clTEB_diff) / self.__norm_fid, nan=0.0)
        N1_TEB_corr = np.einsum('tLl,tl->L', self.__dN1_dCl, clTEB_diff[:, :self.__dN1_dCl.shape[-1]])

        return norm_correction, N1_TEB_corr

    def apply_corrections(self, clkk_kg, norm_correction, N1_TEB_corr, is_cross=False):
        clkk_kg = np.copy(clkk_kg)

        if is_cross:
            clkk_kg[:len(norm_correction)] *= (1 + norm_correction[:len(clkk_kg)])
        else:
            assert (len(clkk_kg) >= self.__Lmax_N1corr + 1), "Input Clkk has Lmax less than the required Lmax of the N1 correction."
            N1_kk_corr = self.__dN1_dClkk[:, :len(clkk_kg)] @ (clkk_kg[:self.__dN1_dCl.shape[-1]] - self.__clkk_fid[:min([self.__dN1_dCl.shape[-1], len(clkk_kg)])])

            clkk_kg[:min([len(norm_correction), len(self.__clkk_fid)])] += 2*norm_correction[:min([len(clkk_kg), len(self.__clkk_fid)])] * self.__clkk_fid[:min([len(norm_correction), len(clkk_kg)])]
            clkk_kg[:len(N1_kk_corr)] += N1_kk_corr[:len(clkk_kg)]
            clkk_kg[:len(N1_TEB_corr)] += N1_TEB_corr[:len(clkk_kg)]

        return clkk_kg

    def correct_cross_spectrum(self, clkg, norm_correction, N1_TEB_corr):
        return self.apply_corrections(clkg, norm_correction, N1_TEB_corr, is_cross=True)

    def correct_auto_spectrum(self, clkk, norm_correction, N1_TEB_corr):
        return self.apply_corrections(clkk, norm_correction, N1_TEB_corr)

    def get_lmax(self):
        return self.__lmax

    def get_Lmax(self):
        return self.__Lmax


class LensingLklhCorrection(Theory):

    lklh_corr_base_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/aux_data/lklh_corr")

    lklh_correction_paths = {'ACT': {'fid_cls': 'cosmo2017_10K_acc3_lensedCls.dat',
                                                        'phi_cls': 'cosmo2017_10K_acc3_lenspotentialCls.dat',
                                                        'dNorm_dCl': 'norm_kk_correction_matrix_Lmin0_Lmax4000_new.npy',
                                                        'fAL': 'n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt',
                                                        'dN1_dkk': 'N1der_KK_lmin600_lmax3000_full.txt',
                                                        'dN1_dCl': {'tt': 'N1der_TT_lmin600_lmax3000_full.txt',
                                                                    'ee': 'N1der_EE_lmin600_lmax3000_full.txt',
                                                                    'bb': 'N1der_BB_lmin600_lmax3000_full.txt',
                                                                    'te': 'N1der_TE_lmin600_lmax3000_full.txt'},
                                                        },
                             'Planck': {'fid_cls': 'cosmo2017_10K_acc3_lensedCls.dat',
                                                                 'phi_cls': 'cosmo2017_10K_acc3_lenspotentialCls.dat',
                                                                 'dNorm_dCl': 'P18_norm_kk_correction_matrix_Lmin0_Lmax3000_new.npy',
                                                                 'fAL': 'PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt',
                                                                 'dN1_dkk': 'N1_planck_der_KK_lmin100_lmax2048.txt',
                                                                 'dN1_dCl': {'tt': 'N1_planck_der_TT_lmin100_lmax2048.txt',
                                                                             'ee': 'N1_planck_der_EE_lmin100_lmax2048.txt',
                                                                             'bb': 'N1_planck_der_BB_lmin100_lmax2048.txt',
                                                                             'te': 'N1_planck_der_TE_lmin100_lmax2048.txt'},
                                                                 },
                             }

    _lklh_correction_modules = {}
    _lmax_TEB = 0
    _Lmax_kk = 0
    _lklh_correction_sample_module_mapping = {}

    def initialize(self):

        self.log.info("Loading  ...")
        

    def initialize_with_provider(self, provider):
        """
        Initialization after other components initialized, using Provider class
        instance which is used to return any dependencies (see calculate below).
        """
        self.provider = provider

    def must_provide(self, **requirements):
        requires = {}

        if 'lensing_lklh_correction' in requirements:
            self.load_lklh_corr_data(requirements['lensing_lklh_correction']['samples'])
            requires['Cl'] = {'tt': self._lmax_TEB, 'te': self._lmax_TEB, 'ee': self._lmax_TEB, 'bb': self._lmax_TEB}

        return requires

    def load_lklh_corr_data(self, samples):

        lensing_recons = set([])
        for s in samples:
            if 'Planck' in s:
                lensing_recons.add('Planck')
                self._lklh_correction_sample_module_mapping[s] = 'Planck'
            if 'ACT' in s:
                lensing_recons.add('ACT')
                self._lklh_correction_sample_module_mapping[s] = 'ACT'

        for key in lensing_recons:
            f_ls, f_tt, f_ee, f_bb, f_te = np.loadtxt(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['fid_cls']), unpack=True)
            f_tt = f_tt / (f_ls * (f_ls + 1.)) * 2. * np.pi
            f_ee = f_ee / (f_ls * (f_ls + 1.)) * 2. * np.pi
            f_bb = f_bb / (f_ls * (f_ls + 1.)) * 2. * np.pi
            f_te = f_te / (f_ls * (f_ls + 1.)) * 2. * np.pi

            fd_ls, f_dd = np.loadtxt(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['phi_cls']), unpack=True, usecols=[0, 5])
            f_kk = f_dd * 2. * np.pi / 4.

            dNorm_dCl = np.load(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['dNorm_dCl']))
            fAL_ls, fAL = np.loadtxt(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['fAL']))
            dN1_dkk = np.loadtxt(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['dN1_dkk']))
            dN1_dCl = []
            for spec in ['tt', 'ee', 'bb', 'te']:
                dN1_dCl.append(np.loadtxt(os.path.join(self.lklh_corr_base_path, self.lklh_correction_paths[key]['dN1_dCl'][spec])))

            self._lklh_correction_modules[key] = LensingLklhCorrectionHelper(np.array([f_tt, f_ee, f_bb, f_te]), f_kk, fAL, dNorm_dCl, np.array(dN1_dCl), dN1_dkk)
            self._lmax_TEB = max([self._lmax_TEB, self._lklh_correction_modules[key].get_lmax()])
            self._Lmax_kk = max([self._Lmax_kk, self._lklh_correction_modules[key].get_Lmax()])
    
    def calculate(self, state, want_derived=True, **params_values_dict):

        cl = self.provider.get_Cl(ell_factor=False, units='FIRASmuK2')

        state['lensing_lklh_correction'] = {}
        cl_TEB = np.array([cl['tt'], cl['ee'], cl['bb'], cl['te']])
        for key in self._lklh_correction_modules.keys():
            state['lensing_lklh_correction'][key] = self._lklh_correction_modules[key].compute_corrections(cl_TEB)

    def get_lensing_lklh_correction(self, cls, sample, cross_spectrum=True):
        correction = self._current_state['lensing_lklh_correction'][self._lklh_correction_sample_module_mapping[sample]]
        return self._lklh_correction_modules[self._lklh_correction_sample_module_mapping[sample]].apply_corrections(cls, *correction, is_cross=cross_spectrum)

    def get_Lmax_kk(self, sample):
        return self._Lmax_kk

    def get_lmax_TEB(self, sample):
        return self._lmax_TEB