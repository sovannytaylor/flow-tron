## flow-tron
This repository uses scripts made by Sovanny Taylor and Borna Novak to streamline flow analysis after gating done in flowjo.

## Workflow 
* Use LSRFortessa (T103) to analyze cells and collect flourescence and cell morphology measurements. 
* Use FLOWJO software to gate the populations and export each population as an excel file. 
* Use soph-tron to take the raw values for each cell from the corresponding flourescence channel and concatanate into one excel.
* Use the next scripts to make a dataframe and plot the values using either violinplot, barplot, or other various plots using seaborn or matplotlib packages. 

### Reproducing workflow
* See the associated template README.md for setup instructions - look over the flow_analysis_template readme file. Run the scripts in the order indicated by the script prefix.
* Make sure that you use the yml file for the environment if a new environment creation is required.


