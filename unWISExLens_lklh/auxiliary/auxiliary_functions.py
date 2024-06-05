import numpy as np
from scipy import stats, interpolate, integrate

def bin_spectrum(ells, cells, bin_edges, ell_weighted=False):
    if ell_weighted:
        sums = stats.binned_statistic(ells, ells, statistic='sum', bins=bin_edges)
        cl = stats.binned_statistic(ells, ells * cells, statistic='sum', bins=bin_edges)
        cl = cl[0] / sums[0]
    else:
        cl = stats.binned_statistic(ells, cells, statistic='mean', bins=bin_edges)[0]

    return cl


def select_from_matrix(mat, selectionA, selectionB=None):
    if selectionB is None:
        selectionB = selectionA
    assert(mat.shape == (len(selectionA), len(selectionB)))
    selection_matrix = np.outer(selectionA, selectionB)

    return mat[np.where(selection_matrix)].reshape(np.sum(selectionA), np.sum(selectionB))

def density_space(xs, ps, n, endpoint=False, order=1, random=False):
    """Draw samples with spacing specified by a density function.

    Copyright Han-Kwang Nienhuys (2020).
    License: any of: CC-BY, CC-BY-SA, BSD, LGPL, GPL.
    Reference: https://stackoverflow.com/a/62740029/6228891

    Parameters:

    - xs: array, ordered by increasing values.
    - ps: array, corresponding densities (not normalized).
    - n: number of output values.
    - endpoint: whether to include x[-1] in the output.
    - order: interpolation order (1 or 2). Order 2 will
      require dense sampling and a smoothly varying density
      to work correctly.
    - random: whether to return random samples, ignoring endpoint).
      in this case, n can be a shape tuple.

    Return:

    - array, shape (n,), with values from xs[0] to xs[-1]
    """

    cps = integrate.cumtrapz(ps, xs, initial=0)
    cps *= (1 / cps[-1])
    intfunc = interpolate.interp1d(cps, xs, kind=order, fill_value='extrapolate')
    if random:
        return intfunc(np.random.uniform(size=n))
    else:
        return intfunc(np.linspace(0, 1, n, endpoint=endpoint))


def combine_inputs(*arrays, sorted=True):
    if sorted:
        return np.sort(np.unique(np.concatenate(arrays)))
    else:
        return np.unique(np.concatenate(arrays))

