FROM continuumio/miniconda3:4.5.12
ENV DEBIAN_FRONTEND=noninteractive

# Install some linux tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    apt-transport-https \
    gnupg \
    graphviz \
    locales \
    less \
    nodejs \
    software-properties-common \
    unzip \
    && apt-get purge -qq && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

# Install wine
RUN dpkg --add-architecture i386 \
    && apt-add-repository 'deb https://dl.winehq.org/wine-builds/debian/ stretch main' \
    && wget -nc https://dl.winehq.org/wine-builds/winehq.key \
    && apt-key add winehq.key \
    && apt-get update \
    && apt-get install -y --install-recommends winehq-stable \
    && apt-get purge && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV SHELL /bin/bash

# Configuration required for using Binder
ARG NB_USER
ARG NB_UID
ENV USER {NB_USER}
ENV HOME /home/${NB_USER}

# create the notebook user
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

WORKDIR ${HOME}
EXPOSE 8888
ENV APP_BASE /srv
ENV CONDA_DIR ${APP_BASE}/conda
ENV NB_PYTHON_PREFIX ${CONDA_DIR}
ENV KERNEL_PYTHON_PREFIX ${NB_PYTHON_PREFIX}
ENV PATH ${CONDA_DIR}/bin:$HOME/.local/bin:${PATH}

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

# Get the contents of our repo added to ${HOME}
COPY . ${HOME}

# install conda requirements
RUN conda config --add channels conda-forge \
    && conda env update -n base -f ${HOME}/environment.yml \
    && conda clean -tipsy && pip install --no-cache-dir -e ${HOME}

RUN echo "source activate" > ~/.bashrc

USER root
RUN chown -R ${NB_USER}:${NB_USER} ${HOME} /opt/conda/envs/
USER ${NB_USER}

# specify jupyterhub configuration
ENTRYPOINT ['~/binder/start.sh']
