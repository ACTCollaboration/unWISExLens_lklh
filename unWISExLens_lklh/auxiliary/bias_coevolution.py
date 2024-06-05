import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

class BiasEvolution(object):
    def __init__(self, path):
        try:
            data = np.loadtxt(path, delimiter=",")
            b1_data = data[:,[0,1]]
            b2_data = data[:,[0,2]]
            b3_data = np.zeros_like(b1_data)
            bs_data = data[:,[0,3]]
        except:
            b1_data = self._sort(np.loadtxt(f"{path}b1_L.csv", skiprows=1, delimiter=","))
            b2_data = self._sort(np.loadtxt(f"{path}b2_L.csv", skiprows=1, delimiter=","))
            b3_data = self._sort(np.loadtxt(f"{path}b3_L.csv", skiprows=1, delimiter=","))
            bs_data = self._sort(np.loadtxt(f"{path}bs_L.csv", skiprows=1, delimiter=","))

        self.__splinesb2 = self.__spline(b1_data, b2_data)
        self.__splinesb3 = self.__spline(b1_data, b3_data)
        self.__splinesbs = self.__spline(b1_data, bs_data)

    @staticmethod
    def _sort(data):
        return data[np.argsort(data[:, 0])]

    @staticmethod
    def __spline(datax, datay):
        b_spl = InterpolatedUnivariateSpline(datax[:,1], datay[:,1])
        b_spl_const = InterpolatedUnivariateSpline(datax[:,1], datay[:,1], ext='const')

        return b_spl, b_spl_const, np.min(datax[:,1])

    def get_b2(self, b1):
        min_b1 = self.__splinesb2[-1]
        return np.piecewise(b1, [b1<min_b1, min_b1<=b1], [self.__splinesb2[1](b1[b1<min_b1]), self.__splinesb2[0](b1[min_b1<=b1])])

    def get_b3(self, b1):
        min_b1 = self.__splinesb3[-1]
        return np.piecewise(b1, [b1 < min_b1, min_b1 <= b1], [self.__splinesb3[1](b1[b1<min_b1]), self.__splinesb3[0](b1[min_b1<=b1])])

    def get_bs(self, b1):
        min_b1 = self.__splinesbs[-1]
        return np.piecewise(b1, [b1 < min_b1, min_b1 <= b1], [self.__splinesbs[1](b1[b1<min_b1]), self.__splinesbs[0](b1[min_b1<=b1])])
