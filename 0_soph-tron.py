import glob
import os

import numpy as np
import pandas as pd

path = 'raw_data'
column = 'R670-A'
output = 'results/excel_clean_up'

if not os.path.exists(output):
    os.mkdir(output)

def read_csv(file_path):
    """
    Read a CSV file and return a DataFrame.
    """
    df = pd.read_csv(file_path)
    return df


def extract_columns(df, columns):
    """
    Extract specified columns from a DataFrame.
    """
    return df[columns]


def main():
    csv_files = sorted(glob.glob(f"{path}/*.csv"))

    # Create an empty DataFrame to store our results
    result_df = pd.DataFrame()

    # Process each file
    for file in csv_files:
        df = read_csv(file)
        # Get just the filename without the path for cleaner column names
        filename = os.path.basename(file)
        # Add the column from this file to our result DataFrame
        result_df[filename] = df[column]

    # Save the result to a new CSV file
    result_df.to_csv(f"{output}/{column}_combined.csv", index=False)
    print(f"Combined data saved to {column}_combined.csv")
    print(f"Number of columns in output: {len(result_df.columns)}")


if __name__ == "__main__":
    main()