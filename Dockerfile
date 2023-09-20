FROM rocker/tidyverse:3.6.3


WORKDIR /app
RUN apt-get update
RUN apt-get install libnetcdf-dev -y
RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install("MSstats")'
ADD ./MSstats /app/MSstats
RUN cp -r /app/MSstats /usr/local/lib/R/site-library
RUN R -e 'BiocManager::install("mzR")'
RUN R -e 'BiocManager::install("reshaper")'
RUN R -e 'BiocManager::install("MSnbase")'
RUN R -e 'BiocManager::install("GOstats")'
