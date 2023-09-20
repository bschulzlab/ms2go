from rpy2.robjects.vectors import DataFrame, FloatVector, ListVector, IntVector, StrVector
import numpy as np
import rpy2.robjects as ro
from collections import OrderedDict
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests


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

class GOStats:
    _ontology = ("BP", "CC", "MF")
    _methods = ["bonferroni", "fdr_bh"]

    def __init__(self, association_file, universe_file, study_file, organism, association_file_delimiter="\t", universe_delimiter="\t", study_delimiter="\t", pvalue_cutoff=0.001):
        self.pvalue_cutoff = pvalue_cutoff
        self.association_file = association_file.replace("\\", "/")
        self.association_file_delimiter = association_file_delimiter
        self.universe_file = universe_file.replace("\\", "/")
        self.universe_delimiter = universe_delimiter
        self.study_file = study_file.replace("\\", "/")
        self.study_delimiter = study_delimiter
        self.organism = organism
        self.code_source = []

    def run_r_code(self, r_code):
        self.code_source.append(r_code)
        return ro.r(r_code)

    def write_r_source(self):
        with open(self.study_file + ".R", "wt", newline="\n") as R_file:
            for i in self.code_source:
                R_file.write(i + "\n")

    def initiate_gostats(self):
        self.run_r_code('''library(GOstats)''')
        self.run_r_code('''library(GSEABase)''')
        build_GO = '''GODB<-read.table("{}", sep="{}", header=TRUE, fill=TRUE, comment.char = "")'''.format(self.association_file, self.association_file_delimiter)
        print(build_GO)
        self.run_r_code(build_GO)
        universe = '''universe<-read.table("{}", sep="{}", header=TRUE, fill=TRUE, comment.char = "")'''.format(self.universe_file, self.universe_delimiter)
        self.run_r_code(universe)
        study = '''study<-read.table("{}", sep="{}", header=TRUE, fill=TRUE, comment.char = "")'''.format(self.study_file, self.study_delimiter)
        self.run_r_code(study)

    def process(self):
        result = []
        gf = '''gf<-GOFrame(GODB, organism="{}")'''.format(self.organism)
        self.run_r_code(gf)
        gaf = '''gaf<-GOAllFrame(gf)'''
        self.run_r_code(gaf)
        gsc = '''gsc<-GeneSetCollection(gaf,setType=GOCollection())'''
        self.run_r_code(gsc)
        for o in self._ontology:
            para = '''params<-GSEAGOHyperGParams(
            name="My params",
            geneSetCollection=gsc,
            geneIds=study$Entry.Name,
            universeGeneIds=universe$Entry.Name,
            pvalueCutoff=1,
            testDirection="over",
            ontology="{}",
            conditional=FALSE)'''.format(o)
            print(para)
            self.run_r_code(para)
            proc = '''Over<-hyperGTest(params)'''
            self.run_r_code(proc)
            go_genes= '''geneIdsByCategory(Over)'''
            g = self.run_r_code(go_genes)
            g_py = recurList(g)
            summary_res = '''summary(Over)'''
            s = self.run_r_code(summary_res)
            s_py = pd.DataFrame(recurList(s))
            s_py["Ontology"] = pd.Series([o]*len(s_py.index), index=s_py.index)
            s_py = s_py.rename(columns={"GO"+o+"ID": "GOID"})
            for i, r in s_py.iterrows():
                if r["GOID"] in g_py:
                    s_py.at[i, "Proteins"] = ";".join(g_py[r["GOID"]])
            s_py = s_py.sort_values("Pvalue")
            for m in self._methods:
                r, p_corrected, als, alb = multipletests(s_py["Pvalue"].values, self.pvalue_cutoff, is_sorted=True,
                                                         method=m)
                s_py["Pvalue_" + m] = p_corrected
            result.append(s_py)
        if not result:
            print("No matches were found for GOstats analysis")
        else:
            self.write_r_source()
            return pd.concat(result, ignore_index=True)



if __name__ == "main":
    g = GOStats("C:/Users/localadmin/PycharmProjects/gostats/gostats_uniprot_parsed.txt",
                "C:/Users/localadmin/PycharmProjects/gostats/Universe.txt",
                "C:/Users/localadmin/PycharmProjects/gostats/DA downreg.txt", "Bos taurus")
    g.initiate_gostats()
    result = g.process()
