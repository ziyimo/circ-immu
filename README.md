## circ-immu covid deaths

#### Dependencies

The scripts in this repo rely on `pso`, `DEoptim`, `deSolve`, `parallel` and `binom`.

#### Scripts

* `Covid_state_joint.R`: __Core functions of the SIRS model__

* `fit_states_death.R`: __Fit model jointly to many states__

  Usage: `$ Rscript ./fit_states_death.R [state_ls.tsv] [R0_model] [params] [threads] &> fit.log`

    `[state_ls.tsv]`: list of states, one per line

    `[R0_model]`: functional form of R0, e.g. `sd`

    `[params]`: model specific colon-separated list of initialization params, e.g. `0.001:0.8:1:-45:0.5:-13:0.32`

    `[threads]`: Number of threads to run in parallel, e.g. `24`

* `fit_cities_death_global_clean.R`: __Fit model to many cities with unique s0__

  Usage: `$ nohup ./fit_cities_death_global_clean.R  [city_ls.tsv] [R0_model] [params] [threads] &> fit.log`

    `[params]`: model specific colon-separated list of initialization params, e.g. `0.001:1:1.2:-289:0.57:-104:0.379`

