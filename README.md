[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9471829&assignment_repo_type=AssignmentRepo)

# IPDA Analyzer

This module will perform the entire IPDA analysis for you. 

What you need to do in order to use this program is:
* copy the .csv file into the input_files folder
* right click on the .csv file and copy the relative pathname
* In the IPDA_analyzer_client.py replace the pathname that is between the '' in input_IPDA_file = '' with your pathname
* it should look like this `input_IPDA_file = 'input_files/YOUR-FILE-NAME.csv'`

in addition to the .csv file from the ddPCR analyzer, you will need to tell the program at which concentration you ran the Rpp30 and HIV reactions
To do so:
* paste you sample names in the Excel DNA concentration template. 
* Make sure that your sample names are idential to your sample names on the ddPCR analyzer (and thus the .csv file)
* paste the concentrations as numbers (no units).
* upload this .xlsx file into the input_files folder
* right click on the .xlsx file and copy the relative pathname
* replace the pathname that is between the '' in input_DNAconcentration = '' with your pathname
* it should look like this `input_DNAconcentration = 'input_files/YOUR-FILE-NAME.xlsx'`


The output is a completely analyzed new .xlsx output_file that contains analyzed IPDA data.
Don't worry, this module won't alter the original raw data.
In addition, it will save a new .png output_file which lists the sample names on the x-axis and plots the intact HIV copies/10^6 cells on the y-axis. 
Note, these values are present in the .xlsx output_file on the Summary table, so you can also copy those into GraphPad Prism if you prefer.


# Additional information and limitations of the IPDA Analyzer

The analysis works for any .csv file exported from the analzer that contains single data, regardless of whether you exported as "Single" or "Both". 
In either case, the single values will be used for analysis as those result in more precise results. 

There are a few limitations of the analyzer but issue  can easily be avoided when the user is aware of how the program functions
* don't export "Merged" data only from the ddPCR analyzer. 
* Always include ***vic*** in the RPP30 vic Target name. Capitalization does not matter, spelling does.
* Always include ***Rpp30 Shear*** or ***Rpp30 Fam*** in the RPP30 Fam Target name. Capitalization does not matter, spelling does.
* Always include ***gag*** or ***Psi*** in the HIV Gag/Psi Target name. Capitalization does not matter, spelling does.
* If you run the same sample with different gag primers and env primers on the same plate, change the name of the sample for each new primer combination
* e.g. imagine you run sample WWHB031 with Jones-Gag/Jones-Env primers as well as Siliciano-Gag/Siliciano-Env primers
* in this case, name the sample WWHB031_withJones in the wells that received the Jones primers. And WWHB031_withSiliciano in the other wells.
* otherwise, the IPDA_analyzer will pool all those samples since he pools anything that's the same Sample and where the Target contains "gag".
* the output from the IPDA Analyzer is always called "Analyzed_IPDA_data.xlsx". So you must rename this file to something else before analyzing the next .csv file
* the IPDA Analyzer uses the reccomended minimum_required_droplets = 10000. However, if your IPDA did not run well and you cannot re-run it, you can adjust this threshold, e.g. to 8000. This might allow you to have more samples pass the quality control step.