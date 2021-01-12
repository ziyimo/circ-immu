### Covid-19 modelling 

NOTE: Currently testing has been restricted to fitting a single variable (daylength). The script may not work as expected with other variables.

To fit a city using the daylength variable, use the below command.

``` Rscript fit_city.R "New York" day ```

The fit model can be visualised together with the R0 estimates by providing an additional parameter to the script.

``` Rscript fit_city.R "New York" day 108207.362674,1637736.177180,-17.320636,0.206693 ```
