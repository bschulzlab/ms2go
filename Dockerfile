# docker build -t ms2go .
# docker run --rm -v ".\data:/data" ms2go -i "/data/GO Stats Input 2021_new.txt"
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /app

RUN apt-get update
RUN apt-get install libnetcdf-dev wget r-base-dev=3.6.3-2 libxml2-dev cmake -y

RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install(version = "3.9", ask = FALSE)'
RUN R -e 'BiocManager::install("MSstats")'
RUN apt-get install
RUN R -e 'BiocManager::install("mzR")'
RUN R -e 'BiocManager::install("reshape")'
RUN R -e 'BiocManager::install("XML")'
RUN R -e 'BiocManager::install("mzID")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/MALDIquant/MALDIquant_1.16.tar.gz")'
RUN R -e 'BiocManager::install("MSnbase")'
RUN R -e 'BiocManager::install("GOstats")'
RUN apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get install python3.9 python3.9-dev python3-pip -y
COPY . /app
RUN cp -r /app/MSstats /usr/local/lib/R/site-library
RUN sed -i 's+C:\\Program Files\\R\\R-3.6.3+/usr/lib/R+g' /app/settings.py
RUN sed -i 's/docker = False/docker = True/g' /app/settings.py
RUN sed -i 's+Reformat_MSstats = "reformatMS_windows_amd64.exe"+Reformat_MSstats = "/app/reformatMS_linux_amd64"+g' /app/settings.py

RUN python3.9 -m pip install -r requirements.txt
ENTRYPOINT ["python3.9", "main.py"]
#ENTRYPOINT ["tail", "-f", "/dev/null"]