# MS2GO

This is a Python, Go and R pipelines for analysis of PeakView data.

# Note
`rpy2`, a Python package necessary for R and Python communication, can be tricky to install on Windows.
In order to get it installed, you first have to make sure that the R binary that you want to run is in windows environmental variable `PATH`.

If `rpy2` could not find your R installation, you will have to manually set R path using the the command below.
`set PATH=%PATH%;path-to-Rbinary` replace `path-to-Rbinary` with the actual path on your computer.

Other Python packages that the pipeline requires are `statsmodels`, `numpy`, `pandas`, `requests`

R packages that the pipeline requires are `MSstats` 2.x and `GOstats`. `MSstats` 2.x requires `MSnbase` and `reshape` which newer `MSstats` does not need.
