### Epidemiological modeling of the circadian immunity project

#### Dependencies

To make sure the `DEoptim` package used is bug-free for parallelized optimization, install the package included in this repo by

``` R CMD INSTALL DEoptim_2.2-5.tar.gz ```

#### Scripts

* `R0_mods.R`: __Family of $R_0$-covariate models__

_Flu analysis_

* `states_flu_analysis.tsv`: __List of states included in the flu analysis__

* `SIRS_model.R`: __Core functions of the SIRS model__

* `fit_postest.R`: __Main model fitting script__

[//]: # (Usage: `$ nohup ./fit_all.R [state_ls.tsv] [R0_model] [lambda] &`)

[//]: #   (`[state_code]`: 2-letter state code)
[//]: # 
[//]: #   (`[R0_model]`: functional form of R0, see `SIRS_model.R` code for options)
[//]: # 
[//]: #   (`[lambda]`: penalty for neg-log-likelihood when model prediction exceeds cap)

* `hessian_flufit.R`: __Curvature analysis of parameter variance using the Hessian__

* `fit_postest_boot.R`: __Fit resampled parameter values through bootstrap__

* `sim_SIRS.R`: __Forward simulations of the flu SIRS models__

_COVID hospitalization analysis_

* `states_COVID_hosp_analysis.tsv`: __List of states included in the COVID hospitalization analysis__

* `SEIH_mod.R`: __Core functions of the SIR-derived hospitalization model__

* `fit_hosp.R`: __Main model fitting script__

* `fit_hosp_iter.R`: __Iterative optimization script__

* `hessian_hospfit.R`: __Curvature analysis of parameter variance using the Hessian__

* `fit_hosp_boot.R`: __Fit resampled parameter values through bootstrap__

* `sim_SEIH.R`: __Forward simulations of the COVID hospitalization models__