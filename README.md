## circ-immu

#### Scripts

* `SIRS_dev.R`: __Sandbox__

* `SIRS_model.R`: __Core functions of the SIRS model__

  `R0_sig`: Single variable sigmoid function for R0
  
  `R0_sig2`: Two variable sigmoid function for R0
  
  `SIRS_R0sig`: Differention equations of SIRS model given pre-computed R0 values
  
  `SIRS1var_pred`: Run SIRS model with single variable given sigmoid parameters, returns raw model predictions (p)
  
  `SIRS2var_pred`:  Run SIRS model with two variables given sigmoid parameters, returns raw model predictions (p)
  
  `p2q`: Scaling and capping of model predictions given data, returns binomial prob (q)
  
  `binom2_L`: Wrapper for 2 variable binomial likelihood
  
  `binom1_L`: Wrapper for 1 variable binomial likelihood

* `fit_SIRS.R`: __Fit sigmoid model__

    Usage: `$ nohup ./fit_SIRS.R [state_code] [variable] [optim_arg] [scaler] &`

    		`[state_code]`: 2-letter state code

    		`[variable]`: one of `sun`, `cli` or `both`
    		
    		`[optim_arg]`: see script

        `[scaler]`: a constant to down scale model predicted `p`, (0, 1]

* `plot_SIRS.R`: __Some plotting functions__ (deprecated, to be updated)

* `tally.R`: __Aggregate state level fitting results__ (deprecated, to be updated)
