The latest release of the MSEtool package is available on [CRAN](https://CRAN.R-project.org/package=MSEtool).

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

### Other minor edits and additions, including:
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
