The latest release of the MSEtool package is available on [CRAN](https://CRAN.R-project.org/package=MSEtool).

## MSEtool 1.3.1 (ongoing)
- New plots for `SRA_scope` markdown report, and draft vignette is now available.
- Since FMSY is the estimated rate parameter in `SP` (surplus production model), a log-Jacobian transform is needed for the r prior and has been added. It is assumed n is fixed in the model.

## MSEtool 1.3.0

### Updates to SRA_scope
- Defaults for `SRA_scope` are now more robust (set maximum F in model, higher std. dev. for likelihood of mean lengths).
- Users can choose to use `SRA_scope` conditioned on either observed catch or observed effort. 
- `SRA_scope` returns an S4 object of class `SRA` with a `plot()` method that generates a markdown report of model fits.

### Assessment models
- A prior for r is now possible with `SP` and `SP_SS` using life history information (priors in natural mortality and steepness, as well as maturity/weight at age). To use this feature, set argument `use_r_prior = TRUE`.
- Default process error standard deviation for `SP_SS` is reduced to 0.1.
- `cDD` and `cDD_SS` are more robust when catch is very, very small (F is set to 0). This is important for management procedures that shut down fishing.
- Minor updates to simplify TMB code.
- Minor revision to `make_MP`.

### Other edits
- Vignette links are now available through the MSEtool help page. Type `?MSEtool` into the console.
- Minor fixes to `multiMSE`.
- `SS2OM` now has an option for selecting male or female life history parameters.

## MSEtool 1.2.1
- Fixed error in Solaris build.

## MSEtool 1.2.0
For the new features described below, DLMtool version 5.3.1 is recommended.

### multiMSE
- The initial release for multi-fleet and multi-stock operating models and MSEs are released in this version, with `multiMSE` being the core function. The multiMSE vignette will be quite useful and can be accessed at `browseVignettes("MSEtool")`.

### Assessment models

Quite a few additions and changes have been made to the Assessment models. See the help manual and vignettes for descriptions of these new Assessment functions.

- The continuous delay-differential model with deterministic and stochastic recruitment (`cDD` and `cDD_SS`, respectively) have been added as new Assessment models to the package. The continuous formulation should be more stable in high F situations.
- A virtual population analysis `VPA` model has also been added to the package.
- The surplus production model `SP` assumes continuous production and estimates continuous F's, similar to ASPIC. This formulation will be more stable in high F situations. The Fox model can be implemented by setting the production function exponent `n = 1`.
- A wrapper function for `spict` (state-space surplus production model) has been written and is available in the `DLMextra` package (located on Github). While reporting functions are available in MSEtool, the output of the wrapper function can still be used with the diagnostic functions in the spict package.
- `SCA` and `SCA2` estimate annual F's and include a likelihood function for the catch. In previous versions, `SCA` matched the predicted catch to observed catch. This feature has been transfered over to the `SCA_Pope` function.
- Summary of assessment results can be obtained by the `plot` function which now generates a markdown report. This will be useful for diagnosing model fits and evaluating parameter estimates.

### Other edits and additions, including:
- A scoping function `SRA_scope` fits an assessment model to catch, indices, and age/length comps to inform historical effort, recruitment deviations, and depletion for data-moderate operating models. Multiple fits are done based on the different life history parameters assumed in the operating model. This function is intended to be an alternative to `DLMtool::StochasticSRA`.
- `profile` and `retrospective` functions for profiling the likelihood and retrospective analyses, respectively, of assessment models are now improved.
- The `compare_models` function has been added to compare time series estimates, e.g. B/BMSY and F/FMSY, among different assessment models.
- Manual starting values (the `start` argument) for parameters of assessment models can be expressions and subsequently evaluated in the assessment function. This can be very helpful when passing starting values in the `make_MP` function.
- The `CASAL2OM` function can be used to generate an operating model from CASAL assessments.
- The `SS2OM` and `SS2Data` functions are updated for the latest versions of r4ss on Github.
- More options are provided to increase flexibility for ramped harvest control rules (`HCR_ramp`).


## MSEtool 1.1.0

### New additions, including:
- Re-parameterization of dome selectivity in SCA so that estimated parameters are age-based after transformation.
- Additional argument in SCA for lognormal distribution for age comps.
- A more efficient method is used to report convergence diagnostics of assessment models when running in closed-loop simulation.

### Minor edits, including:
- By default, steepness is now fixed in the SCA and SCA2 assessment functions.
- By default, nine data-rich MPs are now included in the package. See the help documentation: ?`Data-rich-MP`
- A generic function for ramped harvest control rules (`HCR_ramp`) is now included. Users can input the desired limit and target reference points.
- `make_MP` adds dependencies to the MP so that `DLMtool::Required` returns the appropriate dependencies. Dependencies are dynamic based on the configuration of the assessment model. For example, `Data@steep` is a dependency for a SCA-based model only if steepness is fixed.


## MSEtool 1.0.0

- Initial CRAN release.
- Assessment models: Delay-Difference (DD_TMB, DD_SS); Surplus Production (SP, SP_SS); Statistical Catch-at-Age (SCA, SCA2)
- Harvest control rules: HCR_MSY, HCR40_10, HCR60_20
- `simmov` function for multiple-area movement models (age-independent)
- Functions for converting Stock Synthesis and iSCAM assessments to OM and Data objects (classes inherited from DLMtool)
