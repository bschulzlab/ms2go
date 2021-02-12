# MS2GO

This is a Python, Go and R 3.x pipelines for analysis of PeakView data.

The workflow would start by reading a tabulated `txt` file that define the locations of various files needed for the analysis as well as the label of treatment and control conditions. reformating PeakView ion output with fdr quality control into a format accepted by `MSstats`. Then executed `MSstats` analysis in R with the set treatment and control conditions. Then the program would use the internet to request information from the UniProt database about the proteins analyzed by `MSstats`.

If `GOstats` analysis is allowed to be performed, the `MSstats` dataset would first be cutoff at the set P-value then split into those that were found to be upregulated and downregulated. For each of these datasets, we performed `GOstats` on them separately for 3 pathways, `Molecullar functions`, `Cellular components`, and `Biological processes`. The output would automatically have its P-value adjusted by FDR BH and Bonferronni methods. 

## Note
`rpy2`, a Python package necessary for R and Python communication, can be tricky to install on Windows.
In order to get it installed, you first have to make sure that the R binary that you want to run is in windows environmental variable `PATH`.

If `rpy2` could not find your R installation, you will have to manually set R path using the the command below.
`set PATH=%PATH%;path-to-Rbinary` replace `path-to-Rbinary` with the actual path on your computer.

Other Python packages that the pipeline requires are `statsmodels`, `numpy`, `pandas`, `requests`

R packages that the pipeline requires are `MSstats` 2.x and `GOstats`. If an `MSstats` 3.x were converted to be converted to `MSstats` 2.x, `MSnbase` and `reshape` packages are needed for R too.

Within the folder, `settings.py` is used to for setting p-value cutoff for MSstats and GOstats analysis, as well as location of folder with R binary and the `reformatMSstats` program binary and whether or not the program should automatically perform GOstats analysis.

## Usage

The program accepts one argument, "-i" which is a tabulated `.txt` file with 5 columns, `ion`, `fdr`, `out`, `treatment`, `control`.

- ion: ion file in `csv` format from PeakView ion output.
- fdr: fdr file in `fdr` format from PeakView fdr output.
- out: output folder location. A `/` should be added at the end of the folder path if it does not have it. A new folder would be created in the folder path if the folder does not exist.
- treatment: treatment condition id. Assuming the sample labels are `condtion_replicate`.
- control: control condition id. Assuming the sample labels are `condtion_replicate`.

Operation example:

The tabulated file is `work.txt`

`python main.py -i work.txt`
