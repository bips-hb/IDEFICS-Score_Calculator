# Example_Calculator.R  - Readme

The script shows the exemplary application of the function ScoreCalc() using an example data set. The following files are required for Example_Calculator.R  :

* Example_Calculator.R 
* IDEFICS_Score_Calculator.R
* all_para_tables.RData
* Example_Data.csv


To run the script the following steps must be conducted:

1.	Open the file Example_Calculator.R in RStudio; make sure that your environment is empty.
2.	Change the path in the second step to the directory containing all the above mentioned files (line 21).
3.	If files are not in the same folder, you need to change the file path of Example_Data.csv (line 29), of all_para_tables.RData (line 33), of IDEFICS_Score_Calculator.R (line 36) individually. Furthermore, the path of the output file z_Scored_Data.csv can be changed (line 40).
4.	Run the whole script of Example_Calculator.R. As a result, the data sets data_input and z_scored_data should be in your environment. The latter contains all newly calculated values and is also exported as csv-file into the working directory folder.
