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

`titr = OrgAlkCalc.OrgAlkTitrationBatch("~/Python/OrgAlkCalculations/"
                                         ,"Master_Titration_file.xlsx"
                                         ,"~/Python/OrgAlkCalculations/"
                                         ,"Master_Results_File.xlsx")`

This initialises the batch calculation object. This will load all data contained
in `/Python/OrgAlkCalculations/Master_Titration_file.xlsx`

It is then called using
`titr.batch_calculate()`
This will perform all calculations and write results to the master results file, 
in this case `~/Python/OrgAlkCalculations/Master_ResultsFile.xlsx`.

Alternatively, you may call `batch_calculate` with plotting enabled:
`titr.batch_calculate(plot_results=True)`

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