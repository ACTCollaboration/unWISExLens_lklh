# unWISE x CMB lensing likelihood

This repository provides the public likelihood for the cross-correlation analysis using unWISE galaxies and CMB lensing reconstructions from ACT and Planck.

If you use this software and/or the associated data, please cite both of the following papers:
- [Farren, Krolewski, MacCrann, Ferraro et al ACT Collaboration (2023), arxiv:2309.05659](https://arxiv.org/abs/2309.05659)
- [Farren, Krolewski, Qu, Ferraro et al ACT Collaboration (2024), arxiv:2409.02109](https://arxiv.org/abs/2409.02109)

Furthermore, for the unWISE data cite:
- [Krolewski, Ferraro, Schlafly and White, arxiv:1909.07412](https://arxiv.org/abs/1909.07412)
- [Schlafly, Meisner and Green, arxiv:1901.03337](https://arxiv.org/abs/1901.03337)

For the ACT DR6 Lensing reconstructions and when using the lensing auto-spectrum in the joined 3x2pt analysis cite:
- [Madhavacheril, Qu, Sherwin, MacCrann, Li et al ACT Collaboration (2023), arxiv:2304.05203](https://arxiv.org/abs/2304.05203)
- [Qu, Sherwin, Madhavacheril, Han, Crowley et al ACT Collaboration (2023), arxiv:2304.05202](https://arxiv.org/abs/2304.05202)

When using cross-correlations with the Planck PR4 lensing reconstruction and/or the Planck PR4 lensing auto-spectrum also cite:
- [Carron, Mirmelstein, Lewis (2022), arxiv:2206.07773, JCAP09(2022)039](https://arxiv.org/abs/2206.07773)


## Chains

Chains from Farren et al. 2023 and Farren et al. 2024 are available for download on NERSC [here](https://portal.nersc.gov/project/act/act_x_unWISE_xcorr+3x2pt/).

## Installation
### Option 1: Install from PyPI
You can install the likelihood directly with:

    pip install unWISExLens_lklk

### Option 2: Install from Github
If you wish to be able to make changes to the likelihood for development, first clone this repository. Then install with symbolic links:

    pip install -e . --user

*Note*: Up to *Cobaya* version 3.5.1 a minor bug prevents the defaults from being initialised correctly. It has been fixed [here](https://github.com/CobayaSampler/cobaya/pull/360), but you may have to update your sampler and/or install it from source.

## Data
The bandpowers, covariances and auxiliary data for this likliehood is available for download [here](https://portal.nersc.gov/project/act/act_x_unWISE_xcorr+3x2pt/data_unWISExLens.tar.gz). Download the data archive and extract it inside the cloned directory such that `unWISExLens_lklh/data/` contains three directories `bandpowers`, `covariances`, and `aux_data`. You can simply run the `get_unWISExLens_data.sh` script to achieve this automatically.

## Use with *Cobaya*

This likelihood provides several versions of the cross-correlation and 3x2pt analysis using two redshift samples of unWISE data and CMB lensing reconstructions from ACT DR6 and *Planck* PR4. The analysis requires a dedicted theory module `unWISExLens_lklh.unWISExLensTheory` which has to be included in the theory block of the *Cobaya* `.yaml`-file. The likelihood itself comes in several versions enabling analyses using only the cross-correlations (`XCorr`) or the full 3x2pt dataset (`ThreeXTwo`). You can choose to use both ACT DR6 and *Planck* PR4 (`XCorrACTPlanck` or `ThreeXTwoACTPlanck`) or ACT and *Planck* alone (`XCorr(ACT|Planck)` or `ThreeXTwo(ACT|Planck)`). The `XCorr` likelihood includes the galaxy-CMB lensing cross-correlations ($C_\ell^{\kappa g}$) along with the galaxy auto-correlations ($C_\ell^{gg}$) of the two samples while the `ThreeXTwo` likelihood additionally includes the CMB lensing auto-correlation ($C_\ell^{\kappa \kappa}$).

To use for example the 3x2pt dataset from ACT and *Planck* include the following in your `theory` and `likelihood` blocks.

```
theory:
  camb: ...
  unWISExLens_lklh.unWISExLensTheory: null
likelihood:
  unWISExLens_lklh.ThreeXTwoACTPlanck: null
```

Note that by default the likelihood includes marginalisation over the primary CMB power spectrum (see Farren et al. 2023 and Qu et al. 2023 for details). To combine with primary CMB data set the `want_lensing_lklh_correction` attribute of the likelihood to `True`. Furthermore, this requires the `LensingLklhCorrection` module to be loaded as a theory class. This module provides the likelihood corrections discussed in Appendix A of Farren et al. 2024. Separating this module ensures that the corrections are only evaluated once for each set of cosmological parameters enabling one to take advantage of the parameter speed hierarchy to more efficiently marginalise over the galaxy nuisance parameters which can be evaluated faster than the cosmological parameters. The *Cobaya* `.yaml`-file should then contain the following

```
theory:
  camb: ...
  unWISExLens_lklh.unWISExLensTheory: null
  unWISExLens_lklh.LensingLklhCorrection: null
likelihood:
  primary_CMB_likelihoods: ...
  unWISExLens_lklh.ThreeXTwoACTPlanck:
    want_lensing_lklh_correction: True
```

### Other important parameters

The theory code and likelihood have several options. Most of these should not be altered, but some may be of interest to perform alternative analyses.

Options for the likelihood include the following:

```
unWISExLens_lklh.(ThreeXTwo|XCorr)(ACTPlanck|ACT|Planck):
    samples: # which cross-correlations to use (Note that you will also have to adjust the galaxy nuisance parameters and priors which are automatically selected when using the predefined likelihoods)
        - Blue_ACT # names are in the form (Blue|Green)_(ACT|Planck)
        - ...

    lranges_(gg|kg|kk): # bandpower selection to use in analysis (We caution that the modelling has not been verfied outside the default range, so this should only be used to decrease the scale range, e.g. to restrict to (quasi-)linear scales)
        (Blue|Green)_(ACT|Planck): # for Clgg and Clkg
            - lmin
            - lmax
        (ACT|Planck): # for Clkk
            - lmin
            - lmax

    lensing_auto_spectrum_from_camb: (true|false) # whether to use CAMB to compute the lensing auto-power spectrum or compute it internally using the limber approximation (useful when reconstructing growth of perturbations, sigma8(z))

    #parameters for varying the LPT contributions to Clgg and Clkg (see Farren et al. 2023 for model details; should all be set to false if use_linear_theory = true in theory module or use_cleft = false)
    scale_cleft_b2: false
    shift_cleft_b2: true
    scale_cleft_bs: false
    shift_cleft_bs: true
    
```
The theory module provides the following options:
```
unWISExLens_lklh.unWISExLensTheory:
    use_linear_theory: (false|true) # use linear theory P(k) (will deactivate LPT corrections and HALOFIT)
    use_cleft: (true|false) # whether to use LPT terms
    use_free_cleft_model: (true|false) # use free parameters in LPT expansion rather than bias-coevolutoon relations (see options to likelihood above)
    use_fiducial_cleft: (true|false) # use LPT terms evaluated at fiducial cosmology (requires velocileptor package to deactivate)
    use_Az_parametrisation: (false|true) # use free scaling of the linear power spectrum in bins to reconstruct growth of structure in model agnostic way
    Az_parametrisation:
      type: bins
      bin_edges:
      - 1.09
      - 1.75
      - .inf
      params:
      - A0
      - A1
      - A2
```

### Recommended theory accuracy

For CAMB calls, we recommend the following (or higher accuracy):
- `lmax`: 4000
- `lens_margin`:1250
- `lens_potential_accuracy`: 4
- `AccuracyBoost`:1
- `lSampleBoost`:1
- `lAccuracyBoost`:1
- `halofit_version`:`mead2016`

## Example `.yaml`-files and starting covmat

Along with the likelihood we provide an example `.yaml`-file to perform the likelihood analysis with *Cobaya*. After installation the command

    cobaya-run example_LCDM_unWISExLens_3x2pt+CMB2pt.yaml

will perform a joined analysis of the unWISE 3x2pt data using both ACT and *Planck* lensing along with primary CMB observations from Planck. The file should be adapted for the desired analysis. We also provide a starting covariance matrix for this run (`LCDM_unWISExLens_3x2pt+CMB2pt.covmat`) which can be used as a starting pont for other analyses.

As a simple check on the installation the user may wish to run

    cobaya-run test_unWISExLens_lklh.yaml

which will evaluate the likelihood at a single point. If all components are installed correctly the test should yield a total log-posterior of $-62.1652$.

## Notes

At present the liklihood is compatible only with `camb`. Compatiblity with `class` and emulators replacing the boltzman solver are under development.

