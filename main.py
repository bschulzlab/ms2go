#! python3
#import os

#os.environ["R_HOME"] = r""
#os.environ["path"] = r"C:\Users\localadmin\Anaconda3;C:\Users\localadmin\Anaconda3\Scripts;C:\Users\localadmin\Anaconda3\Library\bin;C:\Users\localadmin\Anaconda3\Library\mingw-w64\lib;C:\Users\localadmin\Anaconda3\Library\mingw-w64\bin;" + os.environ["path"]
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-3.6.3"

import argparse
import msstats
import gostats
import pandas as pd
from io import StringIO
from get_uniprot import UniprotParser, UniprotSequence
import csv


parser = argparse.ArgumentParser(description="Automated workflow for processing PeakView data through MSstats and GOstats")
parser.add_argument("-i",
                    "-input_file",
                    type=str,
                    help="Filepath to experiment description file where each row has 5 columns, ion, fdr, out, treatment, "
                         "control", dest="i")


msstats_pvalue_cutoff = 0.00001
gostats_pvalue_cutoff = 0.01
gostats_check = True


def get_uniprot_data(msstats):
    for i, r in msstats.iterrows():
        seq = UniprotSequence(r["Protein"], True)
        msstats.at[i, "Accession"] = str(seq)
    accessions = msstats["Accession"].unique()
    parser = UniprotParser(accessions, True)

    data = []
    for i in parser.parse("tab"):
        frame = pd.read_csv(StringIO(i), sep="\t")
        frame = frame.rename(columns={frame.columns[-1]: "Accession"})
        data.append(frame)
    data = pd.concat(data, ignore_index=True)
    unmatched = []
    for a in accessions:
        if a not in data["Accession"].values:
            unmatched.append(a)
    if unmatched:
        print("Non-Uniprot ID found:", unmatched)
    return data


def create_go_association_file(data, output_file):
    with open(output_file, "wt", newline="\n") as outfile:
        writer = csv.writer(outfile, dialect="excel", delimiter="\t")
        writer.writerow(["Geneontology IDs", "Gocode", "Entry"])
        for i, r in data.iterrows():
            if pd.notnull(r["Gene ontology IDs"]):
                go_ids = r["Gene ontology IDs"].split("; ")
                if go_ids:
                    for g in go_ids:
                        writer.writerow([g, "IEA", r["Entry name"]])


def perform_gostats(association, universe, study):
    gos = gostats.GOStats(association, universe, study, gostats_pvalue_cutoff)
    gos.initiate_gostats()
    result = gos.process()
    return result

if __name__ == "__main__":
    args = parser.parse_args()
    for i, ms in msstats.process_msstats(args.i):
        uniprot_data = get_uniprot_data(ms)
        uniprot_data.to_csv(i[2]+"_uniprot.txt", "\t", index=False)
        if gostats_check:
            for ind, g in ms.groupby(["Label"]):
                combined_msstats_uniprot = pd.merge(ms, uniprot_data, how="left", on=["Accession"])

                pvalue_cut_ms = g[g["adj.pvalue"] <= msstats_pvalue_cutoff]
                increase_set = pd.merge(pvalue_cut_ms[pvalue_cut_ms["log2FC"] > 0], uniprot_data, how="left", on=["Accession"])
                increase_set = increase_set[pd.notnull(increase_set["Gene ontology IDs"])]
                increase_set.to_csv(i[2]+ind+"increase.txt", sep="\t", index=False)
                decrease_set = pd.merge(pvalue_cut_ms[pvalue_cut_ms["log2FC"] < 0], uniprot_data, how="left", on=["Accession"])
                decrease_set = decrease_set[pd.notnull(decrease_set["Gene ontology IDs"])]
                decrease_set.to_csv(i[2]+ind+"decrease.txt", sep="\t", index=False)

                create_go_association_file(combined_msstats_uniprot, i[2]+ind+"gostats_association.txt")
                result_increase = perform_gostats(
                    i[2]+ind+"gostats_association.txt",
                    i[2]+"_uniprot.txt",
                    i[2]+ind+"increase.txt"
                )
                result_increase.to_csv(i[2]+ind+"gostats_increase.txt", sep="\t", index=False)
                result_decrease = perform_gostats(
                    i[2] + ind + "gostats_association.txt",
                    i[2] + "_uniprot.txt",
                    i[2] + ind + "decrease.txt"
                )
                result_decrease.to_csv(i[2] + ind + "gostats_decrease.txt", sep="\t", index=False)

