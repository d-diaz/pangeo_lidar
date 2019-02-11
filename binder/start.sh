# #!/bin/bash
echo ${JUPYTERHUB_USER}
# Replace DASK_DASHBOARD_URL with the proxy location
sed -i -e "s|DASK_DASHBOARD_URL|/user/${JUPYTERHUB_USER}/proxy/8787|g" \
  ~/binder/jupyterlab-workspace.json
# Get the right workspace ID
sed -i -e "s|WORKSPACE_ID|/user/${JUPYTERHUB_USER}/lab|g" \
  ~/binder/jupyterlab-workspace.json

# Import the workspace into JupyterLab
jupyter lab workspaces import ~/binder/jupyterlab-workspace.json \
  --NotebookApp.base_url=user/${JUPYTERHUB_USER}

exec "$@"
