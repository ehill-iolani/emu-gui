# Set the base image to Ubuntu 20.04 LTS
FROM --platform=linux/amd64 rocker/shiny
# My authorship
LABEL maintainer="ehill@iolani.org"
LABEL version="1.0.0"
LABEL description="Emu GUI for Iolani School"

# Disable prompts during package installation
ENV DEBIAN_FRONTEND noninteractive

# Convenience packages
RUN apt update
RUN apt install -y curl git nano

# Install R packages
RUN R -e "install.packages(c('stringr', 'dplyr', 'ggplot2', 'plotly'))"

# Copy the app to the image
RUN rm -r /srv/shiny-server/*
RUN git clone https://github.com/ehill-iolani/emu-gui.git
RUN cp -r emu-gui/* /srv/shiny-server/
RUN rm -r emu-gui

# Miniconda installation
RUN cd tmp
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh --output miniconda.sh
RUN bash miniconda.sh -bu
ENV PATH="/root/miniconda3/bin:$PATH"
RUN conda update -y conda

# Setting up the emu file system
RUN mkdir /home/database && \
    pip install osfclient && \
    export EMU_DATABASE_DIR=/home/database && \
    cd ${EMU_DATABASE_DIR} && \
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar && \
    tar -xvf emu.tar
    
RUN conda create -y -n emu python=3.10
SHELL ["conda", "run", "-n", "emu", "/bin/bash", "-c"]

# Install emu
RUN conda init && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y emu && \
    conda install -y seqtk && \
    conda install -y -c bioconda nanofilt && \
    echo "conda activate emu" >> ~/.bashrc && \
    mkdir /home/data && \
    # Cleanup
    rm miniconda.sh && \
    rm /home/database/emu.tar

RUN echo "export EMU_DATABASE_DIR=/home/database" >> ~/.bashrc

# Set working directory
WORKDIR /home

# Make /home/ writeable to all "users"
RUN chmod -R 777 /home/

# Make conda executable to all "users"
RUN chmod -R 777 /miniconda3/bin/*
