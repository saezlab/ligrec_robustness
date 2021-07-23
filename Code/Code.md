# Code:

1. **Functions**: This script defines all the functions each other script needs,
sorted by script name.

2. **GetData**: This script will retrieve all the data currently in use in this 
project and drops it in the Data folder for you.

3. **OmniPath Dilution Robustness**: This script dilutes OP with genes active in
the test data set that are not recorded signalling genes. It then investigates 
how much overlap between the top n highest ranked CC-interactions is affected by
the dilution rate.

4. **Snippets**: This script just represents some code snippets I used to test 
the outputs of OmniPath Dilution Robustness. None of the content in it is 
intended to be a long term part of the project. If it were relevant it would be 
moved into the right script.

5. **Tutorials**: This folder has two scripts in it, both of which run through 
reproductions of the official tutorials of dependencies of this project. The 
LIANA tutorial is useful for installing LIANA on windows, the Seurat one formats
the PBMC data set to that it is ready to be used with LIANA.
