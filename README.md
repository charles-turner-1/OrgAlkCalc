# OrgAlkCalc

This is a package which you can use to compute organic alkalinity from titrations.

## Usage
This toolbox contains two basic classes which may be used to perform organic 
alkalinity calculations: `OrgAlkTitration` and `OrgAlkTitrationBatch`.
- `OrgAlkTitrationBatch` is intended as a 'batch mode', which will allow the user
to perform calculations with no additional input.
- `OrgAlkTitration` is more granular, allowing the user to specify in detail how
they wish the calculation to be performed. 

### OrgAlkTitrationBatch

`OrgAlkTitrationBatch` allows the user to take a master spreadsheet and 
automatically perform all organic alkalinity calculations for all titrations 
contained in the master spreadsheet of interest. It is invoked as follows

`titr = OrgAlkCalc.OrgAlkTitrationBatch(master_spreadsheet_path,master_spreadsheet_filename, master_results_path,master_results_filename)`

A sample call is shown below:

`titr = OrgAlkCalc.OrgAlkTitrationBatch("~/Python/OrgAlkCalculations/","Master_Titration_file.xlsx"
,"~/Python/OrgAlkCalculations/","Master_Results_File.xlsx")`

This initialises the batch calculation object as `titr`. This will load all data 
contained in `/Python/OrgAlkCalculations/Master_Titration_file.xlsx`

It is then called using
`titr.batch_calculate()`
This will perform all calculations and write results to the master results file, 
in this case `~/Python/OrgAlkCalculations/Master_ResultsFile.xlsx`.

Alternatively, you may call `batch_calculate` with plotting enabled:
`titr.batch_calculate(plot_results=True)`
This will perform all the same calculations, but additionally plot titration 
curves measured and calculated results.

Each argument of the initialisation is now explained in turn: 

-  `master_spreadsheet_path`  (string)
    - The absolute path of the master spreadsheet. This tells the program 
       where to look for the master spreadsheet which informs the individual 
       calculations.
-  `master_spreadsheet_filename` (string)
    - The name of the master spreadsheet, eg. master_titration.xlsx 
-   `master_results_path`  (string)
    - This function will write results out to a master results file. As with
       master_spreadsheet_path, this argument tells the toolbox which directory
       to look for a master results file to write to.
-   `master_results_filename`  (string)
    - The name of the master results spreadsheet, eg. master_results.xlsx 


### OrgAlkTitration

`OrgAlkTitration` allows the user to take a master spreadsheet and automatically
perform all organic alkalinity calculations for a titration contained in the 
master spreadsheet of interest. It may be invoked as follows:

1. `titr = OrgAlkTitration()` initialises the `OrgAlkTitration` object.
2. `titr.read_master_spreadsheet(master_speadsheet_path,master_speadsheet_filename,
titration_name)` reads the master spreadsheet specified to find the titrations
associated with `titration_name`.
3. `titr.pipeline()` performs all necessary data processing before minimisation.
For additional granular control, inspect `titr.pipeline()`: functions invoked by 
it have additional parameters which may be altered by the user. 
4. `titr.repeat_minimise(minimiser_no,SSR_frac_change_limit,plot_results)` performs
the repeated minimisation in order to calculate output parameters. This must be 
run in order (ie. run `minimiser_no = 1`, `= 2`, `=3`, `=4`).
`SSR_frac_change_limit` specifies the fractional change at which the minimiser 
will stop running.
`plot_results` may be `true` or `false`: if `true`, once the repeated minimisation
has reached its fractional change limit, the data points and calculated titration
curve will be plotted.
5. `titr.select_output_params(row_to_select)`: this selects which output parameters
version we manually from `titr.df_minimiser_outputs` (this DataFrame can be 
inspected manually by the user). Alternatively, we may call 
`titr.select_output_params(batch_mode=True)`, in order to allow the minimiser to 
automatically select the 'best' output parameters, based on automated reliability 
checks.
6. `titr.write_results(master_results_path,master_results_filename)`: this writes
results to a master spreadsheet, specified by `master_results_path` and 
`master_results_filename`.


## Bugs and Errors

If you do find any bug, error, unexpected behavioural quirk or just anything 
which seems odd, please do get in contact: this toolbox is still very much a work
in progress.