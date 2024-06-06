# unWISE x CMB lensing likelihood

This repository provides the public likelihood for the cross-correlation analysis using unWISE galaxies and CMB lensign reconstructions from ACT and Planck.

If you use this software and/or the associated data, please cite both of the following papers:
- [Farren, Krolewski, MacCrann, Ferraro et al ACT Collaboration (2023), arxiv:2309.05659](https://arxiv.org/abs/2309.05659)
- Farren, Krolewski, Qu, Ferraro et al ACT Collaboration (2024; in prep.), arxiv:2407.XXXXX

Furthermore, for the unWISE data cite:
- [Krolewski, Ferraro, Schlafly and White, arxiv:1909.07412](https://arxiv.org/abs/1909.07412)
- [Schlafly, Meisner and Green, arxiv:1901.03337](https://arxiv.org/abs/1901.03337)

For the ACT DR6 Lensing reconstructions and when using the lenesing auto-spectrum in the joined 3x2pt analysis cite:
- [Madhavacheril, Qu, Sherwin, MacCrann, Li et al ACT Collaboration (2023), arxiv:2304.05203](https://arxiv.org/abs/2304.05203)
- [Qu, Sherwin, Madhavacheril, Han, Crowley et al ACT Collaboration (2023), arxiv:2304.05202](https://arxiv.org/abs/2304.05202)

When using cross-correlations with the Planck PR4 lensing reconstruction and/or the Planck PR4 lensing auto-spectrum also cite:
- [Carron, Mirmelstein, Lewis (2022), arxiv:2206.07773, JCAP09(2022)039](https://arxiv.org/abs/2206.07773)


## Chains

Chains from Farren et al. 2023 are available for download on NERSC [here](https://portal.nersc.gov/project/act/act_x_unWISE_xcorr+3x2pt/). Chains from Farren et al. 2024 will also be made available under this address upon journal publication of the paper.

## Installation
### Option 1: Install from PyPI
*NOTE: not yet available*
You can install the likelihood directly with:

    pip install unWISExLens_lklk

### Option 2: Install from Github
If you wish to be able to make changes to the likelihood for development, first clone this repository. Then install with symbolic links:

    pip install -e . --user

## Data
The bandpowers, covariances and auxiliary data for this likliehood is available for download [here](https://portal.nersc.gov/project/act/act_x_unWISE_xcorr+3x2pt/data.zip). Download the data archive and extract it inside the cloned directory such that `unWISExLens_lklh/data/` contains three directories `bandpowers`, `covariances`, and `aux_data`.

## Use with Cobaya

This likelihood provides several versions of the cross-correlation and 3x2pt analysis using two redshift samples of unWISE data and CMB lensing reconstructions from ACT DR6 and *Planck* PR4. The analysis requires a dedicted theory module `unWISExLens_lklh.unWISExLensTheory` which has to be included in the theory block of the Cobaya-`.yaml`-file. The likelihood itself comes in several versions enabling analyses using only the cross-correlations (`XCorr`) or the full 3x2pt dataset (`ThreeXTwo`). You can choose to use both ACT DR6 and *Planck* PR4 (`XCorrACTPlanck` or `ThreeXTwoACTPlanck`) or ACT and *Planck* alone (`XCorr(ACT|Planck)` or `ThreeXTwo(ACT|Planck)`). The `XCorr` likelihood includes the galaxy-CMB lensing cross-correlations ($C_\ell^{\kappa g}$) along with the galaxy auto-correlations ($C_\ell^{gg}$) of the two samples while the `ThreeXTwo` likelihood additionally includes the CMB lensing auto-correlation ($C_\ell^{\kappa \kappa}$).

To use for example the 3x2pt dataset from ACT and *Planck* include the following in your `theory` and `likelihood` blocks.

```
theory:
  camb: ...
  unWISExLens_lklh.unWISExLensTheory: null
likelihood:
  unWISExLens_lklh.ThreeXTwoACTPlanck: null
```

Note that by default the likelihood includes marginalisation over the primary CMB power spectrum (see Farren et al. 2023 and Qu et al. 2023 for details). To combine with primary CMB data set the `want_lensing_lklh_correction` attribute of the likelihood to `True`. Furthermore, this requires the `LensingLklhCorrection` module to be loaded as a theory class. Separating this module into its own component ensures that the corrections are only evaluated for once for each set of cosmological parameters enabling one to take advantage of the parameter speed hierarchy to more efficiently marginalise over the galaxy nuisance parameters which can be evaluated faster than the cosmological parameters. The Cobaya `.yaml`-file should then contain the following

```
theory:
  camb: ...
  unWISExLens_lklh.unWISExLensTheory: null
  unWISExLens_lklh.unWISExLensTheory: null
likelihood:
  primaty_CMB_likelihoods: ...
  unWISExLens_lklh.ThreeXTwoACTPlanck:
    want_lensing_lklh_correction: True
```



