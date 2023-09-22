#! python3
#import os

#os.environ["R_HOME"] = r""
#os.environ["path"] = r"C:\Users\localadmin\Anaconda3;C:\Users\localadmin\Anaconda3\Scripts;C:\Users\localadmin\Anaconda3\Library\bin;C:\Users\localadmin\Anaconda3\Library\mingw-w64\lib;C:\Users\localadmin\Anaconda3\Library\mingw-w64\bin;" + os.environ["path"]
import os
import settings
from multiprocessing import Pool
if settings.R_HOME:
    # check if os is linux
    if os.name != "posix":
        os.environ["R_HOME"] = settings.R_HOME

import argparse
import msstats
import gostats
import pandas as pd
from io import StringIO
import csv
from uniprotparser.betaparser import UniprotParser, UniprotSequence

parser = argparse.ArgumentParser(description="Automated workflow for processing PeakView data through MSstats and GOstats")
parser.add_argument("-i",
                    "-input_file",
                    type=str,
                    help="Filepath to experiment description file where each row has 5 columns, ion, fdr, out, treatment, "
                         "control", dest="i")


msstats_pvalue_cutoff = settings.msstats_cutoff
gostats_pvalue_cutoff = settings.gostats_cutoff
gostats_check = settings.gostats_check

def split_base(work):
    work = pd.read_csv(work, sep="\t")
    new_work = []
    for i, r in work.iterrows():
        r["out"] = r["out"].rstrip("/")
        #print(r["out"])
        #if ";" in r["control"] or ";" in r["treatment"]:
        control = r["control"].split(";")
        treatment = r["treatment"].split(";")
        control_dict = {}
        treatment_dict = {}
        ion = pd.read_csv(r["ion"])
        # ind = ""
        # for c in ion.columns:
        #    if "MANND_3" in c:
        #        ind = c
        #ion = ion.drop(ind, axis=1)
        fdr = pd.read_csv(r["fdr"], dtype={"Decoy": str})
        #fdr = fdr.drop(ind, axis=1)
        # keep = []
        for c in ion.columns[9:]:
            for con in control:
                if con in c:
                    # keep.append(c)
                    if con not in control_dict:
                        control_dict[con] = []
                    control_dict[con].append(c)
            for tre in treatment:
                if tre in c:
                    # keep.append(c)
                    if tre not in treatment_dict:
                        treatment_dict[tre] = []
                    treatment_dict[tre].append(c)
        # print(keep)
        # keep_ion = list(ion.columns[:9]) + keep
        # keep_fdr = list(fdr.columns[:7]) + keep
        for t in treatment_dict:
            for c in control_dict:

                ion_new_name = r["ion"] + "{}vs{}.csv".format(t, c)
                fdr_new_name = r["fdr"] + "{}vs{}.csv".format(t, c)
                keep = list(ion.columns[:9]) + treatment_dict[t] + control_dict[c]
                #print(keep)
                ion_new = ion[keep]
                # ion["Protein"] = ion["Protein"].str.replace("CSL", "P00740")
                ion_new.to_csv(ion_new_name, index=False)
                keep = list(fdr.columns[:7]) + treatment_dict[t] + control_dict[c]
                fdr_new = fdr[keep]
                # fdr["Protein"] = fdr["Protein"].str.replace("CSL", "P00740")
                fdr_new.to_csv(fdr_new_name, index=False)
                new_work.append([ion_new_name, fdr_new_name, r["out"] + "/{}vs{}".format(t, c) + "/", t, c])
        #else:
            #new_work.append(r)
    df = pd.DataFrame(new_work, columns=["ion", "fdr", "out", "treatment", "control"])
    df.to_csv("workextt.txt", index=False, sep="\t")

def get_uniprot_data(msstats):

    accessions = msstats["Accession"].unique()
    parser = UniprotParser()
    data = []
    for i in parser.parse(accessions):
        frame = pd.read_csv(StringIO(i), sep="\t")
        frame = frame.rename(columns={"From": "Accession"})
        data.append(frame)
    if len(data) == 1:
        data = data[0]
    else:
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
            if pd.notnull(r["Gene Ontology IDs"]):
                go_ids = r["Gene Ontology IDs"].split("; ")
                if go_ids:
                    for g in go_ids:
                        writer.writerow([g, "IEA", r["Entry Name"]])


def perform_gostats(association, universe, study):
    gos = gostats.GOStats(association, universe, study, gostats_pvalue_cutoff)
    gos.initiate_gostats()
    result = gos.process()
    return result

if __name__ == "__main__":
    args = parser.parse_args()
    workfile_path = args.i
    if not workfile_path:
        workfile_path = settings.work
    if settings.split:
        split_base(workfile_path)
        workfile_path = "workextt.txt"
    for i, ms in msstats.process_msstats(workfile_path):
        ms = ms[ms["log2FC"] != "NA"]
        for i, r in ms.iterrows():
            seq = UniprotSequence(r["Protein"], True)
            ms.at[i, "Accession"] = str(seq)
        if not os.path.exists(i[2]+"_uniprot.txt"):
            uniprot_data = get_uniprot_data(ms)
            uniprot_data.to_csv(i[2]+"_uniprot.txt", "\t", index=False)
        else:
            uniprot_data = pd.read_csv(i[2]+"_uniprot.txt", sep="\t")
        uniprot_data = uniprot_data[pd.notnull(uniprot_data["Gene Ontology IDs"])]
        uniprot_data.to_csv(i[2]+"_universe.txt", "\t", index=False)
        if gostats_check and uniprot_data.shape[0] > 0:
            for ind, g in ms.groupby(["Label"]):
                print("Running GOstats for " + ind[0])

                combined_msstats_uniprot = g.merge(uniprot_data, how="left", left_on=["Protein"], right_on=[uniprot_data.columns[0]])

                pvalue_cut_ms = g[g["adj.pvalue"] <= msstats_pvalue_cutoff]
                increase_set = pd.merge(pvalue_cut_ms[pvalue_cut_ms["log2FC"] > 0], uniprot_data, how="left", left_on=["Protein"], right_on=[uniprot_data.columns[0]])
                increase_set = increase_set[pd.notnull(increase_set["Gene Ontology IDs"])]
                increase_set.to_csv(i[2]+ind[0]+"increase.txt", sep="\t", index=False)
                decrease_set = pd.merge(pvalue_cut_ms[pvalue_cut_ms["log2FC"] < 0], uniprot_data, how="left", left_on=["Protein"], right_on=[uniprot_data.columns[0]])
                decrease_set = decrease_set[pd.notnull(decrease_set["Gene Ontology IDs"])]
                decrease_set.to_csv(i[2]+ind[0]+"decrease.txt", sep="\t", index=False)

                create_go_association_file(combined_msstats_uniprot, i[2]+ind[0]+"gostats_association.txt")
                if increase_set[pd.notnull(increase_set["Gene Ontology IDs"])].shape[0] > 0:
                    result_increase = perform_gostats(
                        i[2]+ind[0]+"gostats_association.txt",
                        i[2]+"_universe.txt",
                        i[2]+ind[0]+"increase.txt"
                    )
                    result_increase.to_csv(i[2]+ind[0]+"gostats_increase.txt", sep="\t", index=False)
                else:
                    print("No GOIDs found for increase set")
                if decrease_set[pd.notnull(decrease_set["Gene Ontology IDs"])].shape[0] > 0:
                    result_decrease = perform_gostats(
                        i[2] + ind[0] + "gostats_association.txt",
                        i[2] + "_universe.txt",
                        i[2] + ind[0] + "decrease.txt"
                    )
                    result_decrease.to_csv(i[2] + ind[0] + "gostats_decrease.txt", sep="\t", index=False)
                else:
                    print("No GOIDs found for decrease set")

