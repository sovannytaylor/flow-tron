# getting necessary libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from loguru import logger
import os
from statannotations.Annotator import Annotator
logger.info('Import ok')

## ----------- File path -----------------------------

input_folder = 'python_results/plotting/'
output_folder = 'python_results/plotting/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

## ---------Set theme ----------------
# sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

## ------- Initialize the excel that you want to work with ------------

# List all Excel files in the folder (assuming only one Excel file)
input_files = [f for f in os.listdir(input_folder) if f.endswith('.xlsx')]
input_file = input_files[0]  # Get the only Excel file name
file_path = os.path.join(input_folder, input_file)  # Construct the full path
highpeptide_comp_df = pd.read_excel(file_path, engine='openpyxl')

logger.info('file initialization complete :)')

melted_df = highpeptide_comp_df 

####---------Plotting all of the conditions -------------------------------------------

fig, ax = plt.subplots(figsize=(10,10))

# Unique column names from melted_df
columns = melted_df["Condition"].unique()
colors = sns.color_palette("husl", len(columns))  # Generate distinct colors

# Loop through each column and plot normalized histogram
for col, color in zip(columns, colors):
    sns.histplot(
        data=melted_df[melted_df["Condition"] == col], 
        x="median_ldl_normalization", 
        label=col, 
        kde=True, 
        color=color, 
        alpha=0.6,
        stat="density"  # Normalizes the histogram
    )

# Fix the legend transparency
legend = ax.legend(title="Condition", fontsize=12, frameon=False)  # Remove background

plt.xlabel("Normalized Uptake")
plt.ylabel("Density")  # Updated to reflect normalization
plt.title("cPPs compete for LDL internalization")


plt.tight_layout()

# save the plot
plt.savefig(os.path.join(output_folder, 'LNP_round1_peptide_only.png'), dpi=300)

plt.show()


