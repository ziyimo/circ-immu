## circ-immu

#### Scripts

* `SIRS_dev.R`: __Sandbox__

* `SIRS_model.R`: __Core functions of the SIRS model__ (sigmoid, single- and two- variable)

* `fit_SIRS.R`: __Fit sigmoid model__

    Usage: `$ nohup ./fit_SIRS.R [state_code] [variable] [optim_method] &`
    		`[state_code]`: 2-letter state code
    		`[variable]`: one of `sun`, `cli` or `both`
    		`[optim_method]`: one of `DE` (genetic) or `NM` (Nelder-Mead