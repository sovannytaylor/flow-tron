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

input_folder = 'all_rounds_data'
output_folder = 'python_results/summary_plots/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)



## ------- Initialize the excel that you want to work with ------------

# List all Excel files in the folder (assuming only one Excel file)
input_files = [f for f in os.listdir(input_folder) if f.endswith('.xlsx')]
input_file = input_files[0]  # Get the only Excel file name
file_path = os.path.join(input_folder, input_file)  # Construct the full path
lnp_comp_df = pd.read_excel(file_path, engine='openpyxl')

logger.info('file initialization complete :)')

## ---- grouping what are conditions ----------
conditions = ['unstained','lnp_only', 'lnp_hdl', 'lnp_ldl', 'lnp_vldl']


## ------------------creating dataframe------------------- 

# Your data
df = pd.DataFrame({
    'replicate': ['r1', 'r2', 'r4'],
    'unstained': [0.00, 0.00, 0.01],
    'lnp_only': [1.74, 1.28, 1.46],
    'lnp_hdl': [1.87, 1.83, 2.03],
    'lnp_vldl': [0.97, 1.36, 1.09],
    'lnp_ldl': [1.42, 1.21, 0.81]
})

# Melt to long format
melted_df = df.melt(id_vars='replicate', var_name='condition', value_name='value')

file_name = 'melted_data.xlsx'

# Combine the folder path and file name
output_path = f'{output_folder}{file_name}'

# Save the DataFrame to an Excel file
melted_df.to_excel(output_path, index=False)

print(f"Data saved to {output_path}")

### --- Bar plot ---------------------------

# Define a custom blue color palette with lighter, more aesthetically pleasing shades
blue_palette = sns.color_palette(["#A9CBE5", "#7FB9D4", "#5D8FAC", "#42728E", "#295D6D"])

# Get the number of unique replicates and generate shades of gray
replicates = melted_df['replicate'].unique()
gray_palette = sns.color_palette("gray", n_colors=len(replicates))

# Map the replicates to the gray palette
replicate_color_map = {replicate: gray_palette[i] for i, replicate in enumerate(replicates)}

# Set up plot
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Create a barplot with statistical error bars (use custom blue palette here instead of pastel)
ax = sns.barplot(data=melted_df, x='condition', y='value', ci='sd', palette=blue_palette)

# Overlay the strip plot with replicate points in the custom blue palette
sns.stripplot(data=melted_df, x='condition', y='value', hue='replicate',
              dodge=False, palette=blue_palette, size=8, edgecolor='gray', linewidth=1, alpha=0.7)

# Define control and pairs for statistical tests
control = 'lnp_only'
all_conditions = melted_df['condition'].unique()
pairs = [(control, cond) for cond in all_conditions if cond != control]

# Apply the statistical annotations (t-test)
annotator = Annotator(ax, pairs, data=melted_df, x='condition', y='value', order=all_conditions)
annotator.configure(test='t-test_ind', text_format='star', loc='outside')
annotator.apply_test()
annotator.annotate()

# Final touches
ax.set_title("LNP Uptake Across Conditions with Stats", fontsize=14)
ax.set_ylabel("Uptake Value")
ax.set_xlabel("Condition")
ax.legend_.remove()  # Remove the legend, or you can customize it if needed
sns.despine()

# Ensure layout is tight
plt.tight_layout()

# Save the plot
output_folder = './'  # Specify your output folder path
plt.savefig(f'{output_folder}lnp_competition_allreplicates_barplot.png', 
            bbox_inches='tight', pad_inches=0.1, dpi=300)

# Show the plot
plt.show()
## -----Trying more things ----------------

# Order of conditions
order = ['unstained', 'lnp_only', 'lnp_hdl', 'lnp_vldl', 'lnp_ldl']
control = 'lnp_only'
pairs = [(control, cond) for cond in order if cond != control]


# Define a custom blue color palette with lighter, more aesthetically pleasing shades
blue_palette = sns.color_palette(["#A9CBE5", "#7FB9D4", "#5D8FAC", "#42728E", "#295D6D"])

# Get the number of unique replicates and generate shades of gray
replicates = melted_df['replicate'].unique()
gray_palette = sns.color_palette("gray", n_colors=len(replicates))

# Map the replicates to the gray palette
replicate_color_map = {replicate: gray_palette[i] for i, replicate in enumerate(replicates)}

# Plot
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Strip plot with custom blue shades
ax = sns.stripplot(data=melted_df, x='condition', y='value', order=order,
                   jitter=True, size=8, palette=blue_palette, edgecolor='gray', linewidth=0.5)

# Point plot with error bars
sns.pointplot(data=melted_df, x='condition', y='value', order=order,
              ci='sd', join=False, color='black', markers='D', scale=1.2, errwidth=1)

# Stats (with t-test comparisons)
annotator = Annotator(ax, pairs, data=melted_df, x='condition', y='value', order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='outside')
annotator.apply_test()
annotator.annotate()

# Style
ax.set_title("LNP Uptake Across Conditions")
ax.set_ylabel("Uptake Value")
ax.set_xlabel("Condition")
sns.despine()

# Make sure everything fits nicely
plt.tight_layout()

# Save plot
plt.savefig(f'{output_folder}lnp_competition_allreplicates_dotplot.png', 
            bbox_inches='tight', pad_inches=0.1, dpi=300)

plt.show()

## ----- Violin plot ---------------------------

# Filter the melted_df to include only the specified conditions

# Your data and conditions
conditions = ['unstained', 'lnp_only', 'lnp_vldl', 'lnp_ldl', 'lnp_hdl']
control = 'lnp_only'
pairs = [(control, cond) for cond in conditions if cond != control]

x = 'condition'
y = 'value'
order = conditions

# Define a custom blue color palette with lighter, more aesthetically pleasing shades
blue_palette = sns.color_palette(["#A9CBE5", "#7FB9D4", "#5D8FAC", "#42728E", "#295D6D"])

# Get the number of unique replicates and generate shades of gray
replicates = melted_df['replicate'].unique()
gray_palette = sns.color_palette("gray", n_colors=len(replicates))

# Map the replicates to the gray palette
replicate_color_map = {replicate: gray_palette[i] for i, replicate in enumerate(replicates)}

# Create the figure with violin and scatter plot
plt.figure(figsize=(10,10))
ax = sns.violinplot(x=x, y=y, data=melted_df, width=0.6, palette=blue_palette, order=order)

# Add the scatterplot with edge color black and shades of gray for each replicate
sns.scatterplot(data=melted_df, x=x, y=y, edgecolor='k', ax=ax, hue='replicate', palette=replicate_color_map, legend=False, zorder=2)

# Apply the statistical annotations
annotator = Annotator(ax, pairs, data=melted_df, x=x, y=y, order=conditions)
annotator.configure(test='t-test_ind', text_format='star', loc='outside')
annotator.apply_test()
annotator.annotate()

# Formatting
plt.xticks(rotation=45, ha='right')
plt.ylabel('Intensity')
plt.title('LNP competition')
sns.despine()

# Make sure everything fits nicely
plt.tight_layout()
# plt.show()

# Save plot
plt.savefig(f'{output_folder}lnp_competition_allreplicates_violinplot.png', 
            bbox_inches='tight', pad_inches=0.1, dpi=300)
