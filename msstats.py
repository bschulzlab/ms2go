from collections import OrderedDict
import pathlib
from rpy2.robjects.vectors import DataFrame, FloatVector, ListVector, IntVector, StrVector
import numpy as np
import subprocess
import pandas as pd
import rpy2.robjects as ro
import os
import settings

def recurList(data):
    rDictTypes = [DataFrame, ListVector]
    rArrayTypes = [FloatVector, IntVector]
    rListTypes = [StrVector]
    if type(data) in rDictTypes:
        return OrderedDict(zip(data.names, [recurList(elt) for elt in data]))
    elif type(data) in rListTypes:
        return [recurList(elt) for elt in data]
    elif type(data) in rArrayTypes:
        return np.array(data)
    else:
        if hasattr(data, "rclass"):
            raise KeyError('Could not proceed, type {} is not defined'.format(type(data)))
        else:
            return data


class MSstats:
    def __init__(self, ion, fdr, out, binary_path, cutoff=0.01):
        self.ion = ion.replace("\\", "/")
        self.fdr = fdr.replace("\\", "/")
        self.out = os.path.join(out, "reformated.csv").replace("\\", "/")
        self.binary_path = binary_path
        self.cutoff = cutoff
        self.levels = []
        self.code_source = []
    
    def run_r_code(self, r_code):
        print(r_code)
        self.code_source.append(r_code)
        return ro.r(r_code)

    def write_r_source(self):
        with open(self.out + ".R", "wt", newline="\n") as R_file:
            for i in self.code_source:
                R_file.write(i+"\n")
    
    def reformat_ms(self):
        print("Using {} for reformating file: {}".format(self.binary_path, self.ion))
        if "linux" in self.binary_path:
            reformat = subprocess.run(
                [self.binary_path, "-ion", self.ion, "-fdr", self.fdr, "-out", self.out, "-t", str(self.cutoff)]
            )
        else:
            reformat = subprocess.run(
                [self.binary_path, "-ion", self.ion, "-fdr", self.fdr, "-out", self.out, "-t", str(self.cutoff)],
                shell=True
            )
        print(reformat.args)
        if reformat.returncode != 0:
            raise subprocess.CalledProcessError("Error reformatting input ion for MSstats.")
        else:
            print("Finished reformating " + self.ion)


    def get_levels(self):

        self.run_r_code('''library(MSstats)''')
        
        self.run_r_code('''ms<-read.csv("{}", sep=",")'''.format(self.out))
        self.run_r_code('''quant<-dataProcess(ms)''')
        self.levels = self.run_r_code('''levels(quant$GROUP_ORIGINAL)''')

    def generate_comparisons(self, comparisons):
        code_gen_list = []
        comparison_names = []
        comparison_list = []
        for i, c in enumerate(comparisons):
            comp = []
            treatment = []
            control = []
            for l in self.levels:
                if l in c["treatment"]:
                    comp.append("1")
                    treatment.append("(1){}".format(l))
                elif l in c["control"]:
                    comp.append("-1")
                    control.append("(-1){}".format(l))
                else:
                    comp.append("0")
            comparison = "comparison{}".format(str(i))
            comparison_list.append(comparison)
            code_gen_list.append('''{}<-matrix(c({}),nrow=1)'''.format(comparison,",".join(comp)))
            comparison_names.append('''"{}vs{}"'''.format(";".join(treatment), ";".join(control)))
        return pd.DataFrame({"comparison_id": comparison_list, "comparison_name": comparison_names, "comparison_matrix": code_gen_list})

    def process_comparisons(self, comparison_dataframe):
        assert len(comparison_dataframe.index) > 0
        for i, r in comparison_dataframe.iterrows():
            self.run_r_code(r["comparison_matrix"])
        if len(comparison_dataframe.index) > 1:
            self.run_r_code('''comparison<-rbind({})'''.format(",".join(comparison_dataframe["comparison_id"])))
        else:
            self.run_r_code('''comparison<-comparison0''')
        self.run_r_code('''row.names(comparison)<-c({})'''.format(",".join(comparison_dataframe["comparison_name"])))
        self.run_r_code('''testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quant,labeled=FALSE,interference=FALSE,featureVar=TRUE)''')
        self.run_r_code('''write.csv(testResultMultiComparisons$ComparisonResult, "{}")'''.format(self.out+"_msstats.csv"))
        return pd.read_csv(self.out+"_msstats.csv", index_col=0)


def process_msstats(filename):
    work = pd.read_csv(filename, sep="\t")

    # process input data in parallel using multiprocessing pool of 6
    from multiprocessing import Pool

    with Pool(settings.msstats_parallel) as p:
        data = p.map(process_msstats_worker, work.groupby(["ion", "fdr", "out"]))

        for i, r in data:
            if r is not None:
                yield i, r



def process_msstats_worker(group_data):
    i = group_data[0]
    d = group_data[1]
    pathlib.Path(i[2]).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(i[2] + "/reformated.csv"):
        ms = MSstats(i[0], i[1], i[2], settings.Reformat_MSstats)
        comparisons = []
        for i2, r in d.iterrows():
            comp = {}
            comp["treatment"] = r["treatment"].split(";")
            comp["control"] = r["control"].split(";")
            comparisons.append(comp)
        ms.reformat_ms()
        ms.get_levels()
        if not ms.levels:
            print("No samples can be used for comparison.")
            return i, None
        else:
            com_df = ms.generate_comparisons(comparisons)
            result = ms.process_comparisons(com_df)
            ms.write_r_source()
            return i, result
    else:
        print("Already analyzed: {}".format(i[2]))
        return i, pd.read_csv(i[2]+"reformated.csv_msstats.csv", index_col=0)
    # ro.r('''rm(list = ls(all.names = TRUE))''')


if __name__ == "__main__":
    in_file = "work2.txt"
    work = pd.read_csv(in_file, sep="\t")
    for i, d in work.groupby(["ion", "fdr", "out"]):
        ms = MSstats(i[0], i[1], i[2], "reformatMS_windows_amd64.exe")
        comparisons = []
        for i2, r in d.iterrows():
            comp = {}
            comp["treatment"] = r["treatment"].split(";")
            comp["control"] = r["control"].split(";")
            comparisons.append(comp)
        ms.reformat_ms()
        ms.get_levels()
        if not ms.levels:
            raise ValueError("No samples can be used for comparison.")
        else:
            com_df = ms.generate_comparisons(comparisons)
            result = ms.process_comparisons(com_df)


