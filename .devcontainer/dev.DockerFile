FROM ghcr.io/mamba-org/micromamba-devcontainer:git-d175103

# Ensure that all users have read-write access to all files created in the subsequent commands.
ARG DOCKERFILE_UMASK=0000

# Install the Conda packages.
COPY --chown=$MAMBA_USER:$MAMBA_USER ${CONDA_ENV_FILE} /tmp/${CONDA_ENV_FILE}
RUN : \
    # Configure Conda to use the conda-forge channel
    && micromamba config append channels conda-forge \
    # Install and clean up
    && micromamba install --yes --name base --file /tmp/env.yml \
    && micromamba clean --all --yes \
    ;

# Activate the conda environment for the Dockerfile.
# <https://github.com/mamba-org/micromamba-docker#running-commands-in-dockerfile-within-the-conda-environment>
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Create and set the workspace folder
ARG CONTAINER_WORKSPACE_FOLDER=/workspaces/default-workspace-folder
RUN mkdir -p "${CONTAINER_WORKSPACE_FOLDER}"
WORKDIR "${CONTAINER_WORKSPACE_FOLDER}"
