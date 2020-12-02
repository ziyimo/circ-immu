## circ-immu

#### Scripts

* `SIRS_dev.R`: __Sandbox__

* `SIRS_model.R`: __Core functions of the SIRS model__

  `R0_*`: A library of R0 models, see code for details
  
  `SIRS_R0`: Differention equations of SIRS model given pre-computed R0 values
  
  `SIRSvar_pred`: Run SIRS model given i) R0 model, ii) R0 parameters, iii) a list of variables, and iv) population size, returns raw model predictions (p)
  
  `p2q`: Scaling and capping of model predictions given data, returns binomial prob (q)
  
  `binom_Lp`: Wrapper for single state binomial likelihood
  
  `load_state_epi`: Function for loading state epidemiological data

* `fit_all.R`: __Fit model jointly to many states__

  Usage: `$ nohup ./fit_all.R [state_ls.tsv] [R0_model] [share_incpt] [lambda] &`

    `[state_code]`: 2-letter state code

    `[R0_model]`: functional form of R0, options: `exp`, `cdexp`, `bell`, `linEE`, `linGE`, `mixGE`

    `[share_incpt]`: `0` or `1`, whether states share the same "intercept" (peak in sunrise model) or not

    `[lambda]`: penalty for neg-log-likelihood when model prediction exceeds cap

* `fit_state.R`: __Fit model to one state__

  Usage: `$ nohup ./fit_state.R [state_code] [R0_model] [lambda] &`

    `[state_code]`: 2-letter state code

* `plot_SIRS.R`: __Some plotting functions__ (deprecated, to be updated)

* `tally.R`: __Aggregate state level fitting results__ (deprecated, to be updated)
