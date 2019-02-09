FROM continuumio/miniconda3:4.5.12

# Configuration required for using Binder
ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

# create the notebook user
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Install some linux tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    apt-transport-https \
    gnupg \
    graphviz \
    software-properties-common \
    unzip

# Install wine
RUN dpkg --add-architecture i386 \
    && apt-add-repository 'deb https://dl.winehq.org/wine-builds/debian/ stretch main' \
    && wget -nc https://dl.winehq.org/wine-builds/winehq.key \
    && apt-key add winehq.key \
    && apt-get update \
    && apt-get install -y --install-recommends winehq-stable

# Install LAStools
RUN wget http://www.cs.unc.edu/~isenburg/lastools/download/LAStools.zip \
    -O lastools.zip \
    && unzip -q lastools.zip \
    -x "LAStools/*toolbox/*" "LAStools/example*/*" "LAStools/src/*" \
    "LAStools/data/*" \
    -d ${HOME} \
    && rm lastools.zip

# Install FUSION
RUN wget http://forsys.sefs.uw.edu/Software/FUSION/fusionlatest.zip \
    -O fusion.zip \
    && unzip -q fusion.zip -x "APScripts/*" -d ${HOME} \
    && rm fusion.zip


# Get the contents of our repo added to ${HOME}
COPY . ${HOME}

# install conda requirements
RUN conda config --add channels conda-forge \
    && conda env create -n pyFIRS -f ${HOME}/environment.yml \
    && rm -rf /opt/conda/pkgs/*

RUN echo "source activate pyFIRS" > ~/.bashrc
ENV PATH /opt/conda/envs/pyFIRS/bin:$PATH

# install pyFIRS in development mode
RUN pip install --no-cache-dir -e ${HOME}

# pangeo-related configuration for jupyter
RUN jupyter serverextension enable --py nbserverproxy --sys-prefix \
    && jupyter labextension install @jupyter-widgets/jupyterlab-manager \
    @jupyterlab/hub-extension@0.12 \
    @pyviz/jupyterlab_pyviz \
    jupyter-leaflet \
    dask-labextension

# set ownership of HOME and conda env to NB_USER
USER root
RUN chown -R ${NB_UID} ${HOME} /opt/conda/envs/pyFIRS
USER ${NB_USER}

# JUPYTERHUB_USER provided as an argument when the docker image is run
ARG JUPYTERHUB_USER
ENTRYPOINT[${HOME}/binder/start]
