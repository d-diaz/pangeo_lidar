FROM buildpack-deps:bionic

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && apt-get install -qq --yes --no-install-recommends \
    locales \
    && apt-get purge -qq \
    && apt-get clean -qq \
    && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

ENV SHELL /bin/bash

# Configuration required for using Binder
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV NB_USER $NB_USER
ENV HOME /home/${NB_USER}

# create the notebook user
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

WORKDIR ${HOME}

RUN wget --quiet -O - https://deb.nodesource.com/gpgkey/nodesource.gpg.key | \
    apt-key add - && DISTRO="bionic" \
    && echo "deb https://deb.nodesource.com/node_10.x $DISTRO main" \
    >> /etc/apt/sources.list.d/nodesource.list \
    && echo "deb-src https://deb.nodesource.com/node_10.x $DISTRO main" \
    >> /etc/apt/sources.list.d/nodesource.list

# Install some linux tools
RUN apt-get update -qq \
    && apt-get install -qq --yes --no-install-recommends \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    graphviz \
    less \
    nodejs \
    software-properties-common \
    && apt-get purge -qq && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

# Install wine
RUN dpkg --add-architecture i386 \
    && wget -nc https://dl.winehq.org/wine-builds/winehq.key \
    && apt-key add winehq.key \
    && apt-add-repository 'deb https://dl.winehq.org/wine-builds/ubuntu/ bionic main' \
    && apt-get update \
    && apt-get install --yes --install-recommends winehq-stable \
    && apt-get purge -qq && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

# Install LAStools
RUN wget http://www.cs.unc.edu/~isenburg/lastools/download/LAStools.zip \
    -O lastools.zip --quiet \
    && unzip -q lastools.zip -x "LAStools/*toolbox/*" "LAStools/example*/*" \
    "LAStools/src/*" "LAStools/data/*" -d ${HOME} \
    && rm lastools.zip

# Install FUSION
RUN wget http://forsys.sefs.uw.edu/Software/FUSION/fusionlatest.zip \
    -O fusion.zip --quiet \
    && unzip -q fusion.zip -x "APScripts/*" -d ${HOME}/fusion \
    && rm fusion.zip

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.12-Linux-x86_64.sh -O ~/miniconda.sh \
    && /bin/bash ~/miniconda.sh -b -p /opt/conda \
    && rm ~/miniconda.sh \
    && /opt/conda/bin/conda clean -tipsy \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate base" >> ~/.bashrc

EXPOSE 8888
ENV APP_BASE /srv
ENV CONDA_DIR ${APP_BASE}/conda
ENV NB_PYTHON_PREFIX ${CONDA_DIR}
ENV KERNEL_PYTHON_PREFIX ${NB_PYTHON_PREFIX}
ENV PATH ${CONDA_DIR}/bin:$HOME/.local/bin:${PATH}

# Get the contents of our repo added to ${HOME}
COPY . ${HOME}
# COPY ./binder ${HOME}/binder

# install conda requirements
RUN conda config --system --add channels conda-forge \
    && conda config --system --set auto_update_conda false \
    && conda config --system --set show_channel_urls true \
    && conda env update -n base -f ${HOME}/binder/pangeo_env.py-3.6.frozen.yml \
    && conda env update -n base -f ${HOME}/environment.yml \
    && conda clean -tipsy && pip install --no-cache-dir -e ${HOME}

# pangeo-related configuration for jupyter
RUN jupyter serverextension enable nteract_on_jupyter --sys-prefix \
    && jupyter serverextension enable --py nbserverproxy --sys-prefix \
    && jupyter labextension install @jupyter-widgets/jupyterlab-manager \
    @jupyterlab/hub-extension@0.12 \
    @pyviz/jupyterlab_pyviz \
    jupyter-leaflet \
    dask-labextension

USER root
RUN chown -R ${NB_USER}:${NB_USER} ${HOME} /srv/conda/envs/
USER ${NB_USER}

# specify jupyterhub configuration
RUN ls -lash ${HOME}/binder/start.sh
ENTRYPOINT /bin/bash ${HOME}/binder/start.sh
