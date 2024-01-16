###########################################################################
## Cattaneo, Jansson, and Nagasawa (2024)
## Bootstrap-assisted inference for generalized Grenander-type estimators
## Replication files
###########################################################################

-collect_tables_isoreg.R
This file constructs the table from the output of main_isoreg.R

-isoreg.cpp
This file is called within main_function_isoreg. 
It provides functions that implement our proposed boostrap-assisted inference procedure.
It is called within main_function_isoreg.R.

-main_function_isoreg.R
This file is called by main_isoreg.R and it defines functions used to compute various estimates. 

-main_isoreg.R
This is the main R file that implements simulations. The code was run in batch mode.