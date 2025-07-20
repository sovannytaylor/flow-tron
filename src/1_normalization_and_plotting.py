import os
import matplotlib.pyplot as plt
from itertools import combinations
import pandas as pd
import seaborn as sns
from loguru import logger
import numpy as np
from statannotations.Annotator import Annotator
from scipy.stats import levene
from scipy import ndimage, stats
from scipy.stats import mannwhitneyu, kruskal, shapiro
pd.DataFrame.iteritems = pd.DataFrame.items


logger.info('Import ok')

input_folder = 'results/excel_clean_up/'
output_folder = 'results/plotting/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# -------------- defining functions --------------

def remove_outliers_iqr(df, col):
    Q1 = df[col].quantile(0.25)
    Q3 = df[col].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[col] >= lower_bound) & (df[col] <= upper_bound)]



## ------- Initialize the excel that you want to work with ------------
input_files = [f for f in os.listdir(input_folder) if f.endswith(('.xlsx', '.csv'))]
file_path = os.path.join(input_folder, input_files[0])

if file_path.endswith('.xlsx'):
    df = pd.read_excel(file_path, engine='openpyxl')
elif file_path.endswith('.csv'):
    df = pd.read_csv(file_path)
else:
    raise ValueError("Unsupported file format. Expected .xlsx or .csv.")
logger.info('file initialization complete :)')


## ----make the conditions the peptide---

#condense the labels 
df.columns = df.columns.str.extract(r'export_(.*?)_SSC singlets', expand=False)



#define conditions 
filename = ['Round1_1hour_LNP_CROT', 'Round1_1hour_LNP_GP30',
       'Round1_1hour_LNP_GR30', 'Round1_1hour_LNP_HDL', 'Round1_1hour_LNP_LDL',
       'Round1_1hour_LNP_LL37', 'Round1_1hour_LNP_MOLLUSC',
       'Round1_1hour_LNP_VLDL', 'Round1_1hour_LNP_only',
       'Round1_1hour_unstained', 'Round2_1hour_LNP_CROT',
       'Round2_1hour_LNP_GP30', 'Round2_1hour_LNP_GR30',
       'Round2_1hour_LNP_HDL', 'Round2_1hour_LNP_LDL', 'Round2_1hour_LNP_LL37',
       'Round2_1hour_LNP_MOLLUSC', 'Round2_1hour_LNP_VLDL',
       'Round2_1hour_LNP_only', 'Round2_1hour_unstained',
       'Round3_1hour_LNP_CROT', 'Round3_1hour_LNP_GP30',
       'Round3_1hour_LNP_GR30', 'Round3_1hour_LNP_HDL', 'Round3_1hour_LNP_LDL',
       'Round3_1hour_LNP_LL37', 'Round3_1hour_LNP_MOLLUSC',
       'Round3_1hour_LNP_VLDL', 'Round3_1hour_LNP_only',
       'Round3_1hour_unstained']



# Convert DataFrame to long format for plotting
melted_df = df.melt(var_name='filename', value_name='value')


# # #drop NA
melted_df = melted_df.dropna()

# Create Replicate column
melted_df['rep'] = melted_df['filename'].str.split('_').str[0]
melted_df['condition'] = melted_df['filename'].str.split('_').str[-1]
                                                             

# grab avg per replicate
avg_per_condition = melted_df.groupby(['condition', 'rep'])['value'].mean().reset_index()

#----------------Normalizing the data and plotting that ------------------

# outlier removal ?? 

melted_df = remove_outliers_iqr(melted_df, "value")

# Step 1: Compute mean of 'unstained' per replicate
unstained_means = (
    melted_df[melted_df["condition"] == "unstained"]
    .groupby("rep")["value"]
    .mean()
    .reset_index()
    .rename(columns={"value": "unstained_mean"})
)

# Step 2: Merge rep-specific unstained means back to main DataFrame
melted_df = melted_df.merge(unstained_means, on="rep", how="left")

# Step 3: Subtract unstained mean from all values, clip at 0
melted_df["background_correction"] = (
    melted_df["value"] - melted_df["unstained_mean"]
).clip(lower=0)

# Optional: Drop rows with missing values (e.g., if rep had no unstained)
melted_df = melted_df.dropna(subset=["background_correction"])

# Step 4: Compute per-rep mean for 'only' condition after subtracting
lnp_means = (
    melted_df[melted_df["condition"] == "only"]
    .groupby("rep")["background_correction"]
    .mean()
    .reset_index()
    .rename(columns={"background_correction": "lnp_mean"})
)

# Step 5: Merge those means back to the original DataFrame
melted_df = melted_df.merge(lnp_means, on="rep", how="left")

# Step 6: Normalize
melted_df["mean_lnp_normalization"] = (
    melted_df["background_correction"] / melted_df["lnp_mean"])



# --------- compiling the dataframe for the save ---------

avg_per_condition_normalized = melted_df.groupby(['condition','rep'])['mean_lnp_normalization'].mean().reset_index()

#save the normalized file 
output_file = os.path.join(output_folder, 'ANA-2556_LDLR_melted-df-lnp647-1hr.xlsx')
melted_df.to_excel(output_file, index=False)


#save the normalized file 
output_file = os.path.join(output_folder, 'summary_LNP647_1hr.xlsx')
avg_per_condition_normalized.to_excel(output_file, index=False)





## ----- Violin plot for normalized ---------------------------
blue_palette = sns.color_palette(["#A9CBE5", "#7FB9D4", "#5D8FAC", "#42728E", "#295D6D"])

# Define x and y axes
x = 'condition'
y = 'mean_lnp_normalization'

#define order 
order  = ['only', 'VLDL','LDL','HDL']

# Pairs for statistical testing
pairs = [('only', condition) for condition in order if condition != 'only']

# Create plot
plt.figure(figsize=(10, 10), constrained_layout=True)
ax = sns.boxplot(x=x, y=y, data=avg_per_condition_normalized, width=0.6, palette=blue_palette, order=order)

# Overlay replicate mean dots
sns.stripplot(
    x=x, y=y, data=avg_per_condition_normalized,
    order=order,
    color='black',
    size=6,
    jitter=True,
    ax=ax
)

# Add mean value labels on top of each violin
grouped_means = avg_per_condition_normalized.groupby(x)[y].mean()
for tick, label in enumerate(order):
    if label in grouped_means:
        mean_val = grouped_means[label]
        ax.text(
            tick, mean_val + 0.1, f'{mean_val:.2f}',  # adjust vertical offset if needed
            ha='center', va='bottom', fontsize=9, color='black', weight='bold'
        )

# Add statistical annotation
annotator = Annotator(ax, pairs, data=avg_per_condition_normalized, x=x, y=y, order=order)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='outside',
    line_height=0.02,     # Controls spacing between annotation lines
    verbose=2)
annotator.apply_test()
annotator.annotate()

# Formatting
plt.xticks(rotation=45, ha='right')
# plt.ylim(0, 3.5)
# ymax = melted_df[y].max()

plt.ylabel('Normalized Mean Uptake')
plt.title('Normalized Mean Uptake - 60 min - MannWhitney')
sns.despine()
# plt.tight_layout()

# Save plot
plt.savefig(f'{output_folder}lnp_60min_mean_violin_normalized_mannwhit.png', 
            bbox_inches='tight', pad_inches=0.1, dpi=300)

plt.savefig(f'{output_folder}lnp_60min_mean_violin_normalized_mannwhit.svg', format = 'svg')





### --- additional distrubtion stats --------
### trying to test distribution using kolmogorov smir

from scipy.stats import ks_2samp


# LNP VS LDL
group1 = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
group2 = melted_df[melted_df['condition'] == 'lnp_ldl']['mean_lnp_normalization']
ks_stat, ks_p = ks_2samp(group1, group2)
print(f"KS test statistic: {ks_stat}, p-value: {ks_p}")

# Subset your data correctly
lnp_only_data = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
lnp_ldl_data = melted_df[melted_df['condition'] == 'lnp_ldl']['mean_lnp_normalization']

# Plot KDEs
sns.kdeplot(lnp_only_data, label='LNP only')
sns.kdeplot(lnp_ldl_data, label='LNP with LDL')

# Annotate
plt.title(f'KS p-value = {ks_p:.4g}, KS STATISTIC: 0.358')
plt.xlabel('mean_lnp_normalization')
plt.ylabel('Density')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig(f'{output_folder}ldl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.svg', format='svg')
plt.savefig(f'{output_folder}ldl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.png', bbox_inches='tight', pad_inches=0.1, dpi=300)



# #LNP VS VLDL

# group1 = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
# group2 = melted_df[melted_df['condition'] == 'lnp_vldl']['mean_lnp_normalization']
# ks_stat, ks_p = ks_2samp(group1, group2)
# print(f"KS test statistic: {ks_stat}, p-value: {ks_p}")


# # Subset your data correctly
# lnp_only_data = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
# lnp_vldl_data = melted_df[melted_df['condition'] == 'lnp_vldl']['mean_lnp_normalization']

# # Plot KDEs
# sns.kdeplot(lnp_only_data, label='LNP only')
# sns.kdeplot(lnp_vldl_data, label='LNP with VLDL')

# # Annotate
# plt.title(f'KS p-value = {ks_p:.4g}, KS STATISTIC: 0.244')
# plt.xlabel('mean_lnp_normalization')
# plt.ylabel('Density')
# plt.legend()
# plt.tight_layout()
# # plt.show()
# plt.savefig(f'{output_folder}vldl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.svg', format='svg')
# plt.savefig(f'{output_folder}vldl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.png', bbox_inches='tight', pad_inches=0.1, dpi=300)


# #LNP VS HDL

# group1 = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
# group2 = melted_df[melted_df['condition'] == 'lnp_hdl']['mean_lnp_normalization']
# ks_stat, ks_p = ks_2samp(group1, group2)
# print(f"KS test statistic: {ks_stat}, p-value: {ks_p}")


# # Subset your data correctly
# lnp_only_data = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
# lnp_vldl_data = melted_df[melted_df['condition'] == 'lnp_hdl']['mean_lnp_normalization']

# # Plot KDEs
# sns.kdeplot(lnp_only_data, label='LNP only')
# sns.kdeplot(lnp_vldl_data, label='LNP with HDL')

# # Annotate
# plt.title(f'KS p-value = {ks_p:.4g}, KS STATISTIC: 0.151')
# plt.xlabel('mean_lnp_normalization')
# plt.ylabel('Density')
# plt.legend()
# plt.tight_layout()
# # plt.show()
# plt.savefig(f'{output_folder}hdl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.svg', format='svg')
# plt.savefig(f'{output_folder}hdl_lnp_competition_24hrs_KSTEST_MEANLNPNORM.png', bbox_inches='tight', pad_inches=0.1, dpi=300)


### ------ Plotting them together on one plot 


group1 = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
group2 = melted_df[melted_df['condition'] == 'lnp_ldl']['mean_lnp_normalization']
group2 = melted_df[melted_df['condition'] == 'lnp_vldl']['mean_lnp_normalization']
group3 = melted_df[melted_df['condition'] == 'lnp_hdl']['mean_lnp_normalization']
# ks_stat, ks_p = ks_2samp(group1, group2)
# print(f"KS test statistic: {ks_stat}, p-value: {ks_p}")

# Subset your data correctly
lnp_only_data = melted_df[melted_df['condition'] == 'lnp_only']['mean_lnp_normalization']
lnp_ldl_data = melted_df[melted_df['condition'] == 'lnp_ldl']['mean_lnp_normalization']
lnp_vldl_data = melted_df[melted_df['condition'] == 'lnp_vldl']['mean_lnp_normalization']
lnp_hdl_data = melted_df[melted_df['condition'] == 'lnp_hdl']['mean_lnp_normalization']



# Plot KDEs
# orange_palette = sns.color_palette(["#FDCB9E", "#FDB57B", "#FD9E58", "#FC8740", "#F66B1E"])
sns.kdeplot(lnp_only_data, label='LNP only', color="#0A0A0A" )
sns.kdeplot(lnp_ldl_data, label='LNP with LDL',color="#2BCF15")
sns.kdeplot(lnp_vldl_data, label='LNP with VLDL', color="#37BDBA")
sns.kdeplot(lnp_hdl_data, label='LNP with HDL',color="#DC1ABF")

# Add vertical line at control mean
control_mean = lnp_only_data.mean()
plt.axvline(control_mean, color='black', linestyle='--', linewidth=1.5, label='LNP only mean')

# Annotate
plt.title('LNP uptake is affected by lipoproteins')
plt.xlabel('norm EGFP expression')
plt.ylabel('Density')
plt.legend()
plt.tight_layout()
# plt.show()
# plt.show()
plt.savefig(f'{output_folder}full_lnp_competition_24hrs_KSTEST_MEANLNPNORM.svg', format='svg')
plt.savefig(f'{output_folder}full_lnp_competition_24hrs_KSTEST_MEANLNPNORM.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
