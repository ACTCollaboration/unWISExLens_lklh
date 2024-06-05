import numpy as np
from .auxiliary_functions import bin_spectrum

class PowerSpectrumBinning(object):
    def __init__(self, ell_bin_edges, transfer_function=None):
        self._ell_bin_edges = ell_bin_edges
        self.__transfer = transfer_function if transfer_function is not None else np.ones(len(ell_bin_edges)-1)

    def __call__(self, cells, white_noise=None):
        return self.bin(cells, white_noise=white_noise) * self.__transfer

    def bin(self, signal_cells, white_noise=None):
        assert (len(signal_cells) == len(self.get_input_ells()))
        if white_noise is None:
            return bin_spectrum(self.get_input_ells(), signal_cells, self._ell_bin_edges)
        else:
            assert(isinstance(white_noise, float))
            return bin_spectrum(self.get_input_ells(), signal_cells + white_noise, self._ell_bin_edges)

    def get_input_ells(self):
        return np.arange(max([np.min(self._ell_bin_edges) - 1, 0]), np.max(self._ell_bin_edges) + 1, 1).astype('int')

    def get_ell_bin_edges(self):
        return self._ell_bin_edges

    def get_binned_ell_vals(self):
        return (self._ell_bin_edges[:-1]+self._ell_bin_edges[1:])/2


class MatrixPowerSpectrumBinning(PowerSpectrumBinning):

    def __init__(self, binning_matrix, transfer_function=None, lmin=0.0):
        self.__binning_matrix = binning_matrix
        self._lmin = lmin
        self._lmax = lmin + self.__binning_matrix.shape[-1] - 1

        super().__init__(np.full(self.__binning_matrix.shape[0]+1, np.nan), transfer_function)

    def bin(self, signal_cells, white_noise=None):
        assert (len(signal_cells) == len(self.get_input_ells()))
        if white_noise is None:
            white_noise = 0.0
        assert (isinstance(white_noise, float))
        return self.__binning_matrix @ (signal_cells + white_noise)

    def get_input_ells(self):
        return np.arange(self._lmin, self._lmax + 1, 1)

    def get_binned_ell_vals(self):
        return self.__binning_matrix @ self.get_input_ells()


class NaMasterPowerSpectrumBinning(PowerSpectrumBinning):
    def __init__(self, coupling_matrix, bandpower_windows, ell_bin_edges, lmax=None, transfer_function=None):
        self.__coupling_matrix = coupling_matrix
        self.__decoupling_matrix = bandpower_windows@np.linalg.inv(self.__coupling_matrix)
        self.__w2 = np.sum(self.__coupling_matrix[0, :])
        if lmax is not None:
            assert(lmax <= self.__coupling_matrix.shape[-1]-1)
            self._lmin, self._lmax = 0, lmax
        else:
            self._lmin, self._lmax = 0, self.__coupling_matrix.shape[-1]-1

        super().__init__(ell_bin_edges, transfer_function)

    def bin(self, signal_cells, white_noise=None):
        if white_noise is None:
            white_noise = 0.0
        assert(len(signal_cells) == self._lmax - self._lmin + 1)
        assert(isinstance(white_noise, float))
        padded_cells = np.zeros((self.__coupling_matrix.shape[-1], *signal_cells.shape[1:]))
        padded_cells[self._lmin:self._lmax+1] = signal_cells

        return self.__decoupling_matrix@(self.__coupling_matrix@padded_cells + white_noise*self.__w2)

    def get_input_ells(self):
        return np.arange(self._lmin, self._lmax+1, 1)

    def get_bpw(self):
        return self.__decoupling_matrix @ self.__coupling_matrix
