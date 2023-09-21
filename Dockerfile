FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /app

RUN apt-get update
RUN apt-get install libnetcdf-dev wget r-base-dev=3.6.3-2 -y
RUN apt-get install cmake -y

RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install("MSstats")'

RUN R -e 'BiocManager::install("mzR")'
RUN R -e 'BiocManager::install("reshaper")'
RUN R -e 'BiocManager::install("XML")'
RUN R -e 'BiocManager::install("MSnbase")'
RUN R -e 'BiocManager::install("GOstats")'
COPY . /app
RUN cp -r /app/MSstats /usr/local/lib/R/site-library
RUN sed -i 's/docker = False/docker = True/g' /app/settings.py
RUN apt-get install python3-pip -y
RUN pip3 install -r requirements.txt
ENTRYPOINT ["bash"]