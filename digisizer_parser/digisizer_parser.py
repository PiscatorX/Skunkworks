import argparse
from pathlib import Path

import pandas as pd


def extract_test3_data(xls_file):

    data = pd.read_excel(xls_file, engine = "xlrd", header = None)

    sample_id = str(data.loc[data.iloc[:, 0].eq("Sample:"),1].iat[0]).strip()

    header_row = data.index[data.iloc[:, 0].eq("Particle Diameter (µm)")][0]

    blank = data.iloc[header_row + 1:, 0].isna()

    last_row = blank.idxmax() if blank.any() else len(data)

    data = data.iloc[header_row:last_row, :4].copy()

    data.columns = data.iloc[0]

    data = data.iloc[1:].reset_index(drop = True)

    data.insert(0, "sample_id", sample_id)

    return data


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("xls_file", type = Path)
    parser.add_argument("-o", "--output", type = Path)
    args = parser.parse_args()
    extract_test3_data(args.xls_file).to_csv(args.output or args.xls_file.with_suffix(".csv"),index = False)
