FROM continuumio/miniconda3:4.5.12

# Configuration required for using Binder
ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

# create the notebook user jovyan
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Install some linux tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    unzip \
    software-properties-common \
    gnupg \
    apt-transport-https \
    graphviz

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
RUN jupyter serverextension enable --py nbserverproxy --sys-prefix
RUN labextension install @jupyter-widgets/jupyterlab-manager \
    @jupyterlab/hub-extension@0.12 \
    @pyviz/jupyterlab_pyviz \
    jupyter-leaflet \
    dask-labextension

# Replace DASK_DASHBOARD_URL with the proxy location
RUN sed -i -e "s|DASK_DASHBOARD_URL|/user/${JUPYTERHUB_USER}/proxy/8787|g" \
    ${HOME}/binder/jupyterlab-workspace.json
# Get the right workspace ID
RUN sed -i -e "s|WORKSPACE_ID|/user/${JUPYTERHUB_USER}/lab|g" \
    ${HOME}/binder/jupyterlab-workspace.json

# Import the workspace into JupyterLab
RUN jupyter lab workspaces import ${HOME}/binder/jupyterlab-workspace.json \
  --NotebookApp.base_url=user/${JUPYTERHUB_USER}

RUN exec "$@"

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
