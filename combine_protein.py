from glob import glob
import pandas as pd
import os

if __name__ == "__main__":
    df_all = pd.DataFrame()
    for f in glob("./data/*.csv"):
        filename = os.path.split(f)[1]
        df = pd.read_csv(f)
        for i, r in df.iterrows():
            if not pd.isnull(r["log2FC"]):
                df_all.at[r["Protein"], r["Label"]+ "_" + filename] = r["log2FC"]

    df_all = df_all.fillna(0)
    df_all.to_csv("combined.csv", index=True)