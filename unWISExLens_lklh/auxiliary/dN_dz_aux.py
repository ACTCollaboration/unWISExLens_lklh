from scipy.interpolate import interp1d
import numpy as np
from scipy.integrate import quad as integrate, IntegrationWarning
import warnings

class dN_dz_Helper:
    @staticmethod
    def get_dn_dz(file_path, pad=0.0, n_norm_integral=1000):
        data = np.loadtxt(file_path)
        return dN_dz_Helper.__get_normalized_interp(data[:, 1], data[:, 0], pad=pad, n_norm_integral=n_norm_integral), np.min(data[:, 0]), np.max(data[:, 0])

    @staticmethod
    def __get_normalized_interp(dN_dz, z_vals, pad=0.0, n_norm_integral=1000):
        func = interp1d(z_vals, dN_dz, fill_value=pad, bounds_error=False)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=IntegrationWarning)
            norm = integrate(func, np.min(z_vals), np.max(z_vals), limit=n_norm_integral)[0]
        return lambda z: func(z)/norm