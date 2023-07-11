FROM ghcr.io/mamba-org/micromamba:jammy

# Take root for the system-wide setup
USER root

# Install system-wide packages with `apt`
RUN export DEBIAN_FRONTEND=noninteractive && apt update --fix-missing && \
    apt install --no-install-recommends -y \
    ca-certificates \
    moreutils \
    bash-completion \
    less \
    watch \
    git \
    patch \
    git-lfs \
    nano \
    jq \
    gnupg2 \
    openssh-client \
    file \
    htop \
    zip \
    unzip \
    p7zip-full \
    curl \
    wget \
    lsof \
    iputils-ping \
    iproute2 \
    net-tools \
    dnsutils \
    socat \
    telnet \
    && apt clean && rm -rf /var/lib/apt/lists/*

# Install bash completions
RUN mkdir -p /etc/bash_completion.d \
    && sh -c "curl -L https://raw.githubusercontent.com/docker/compose/1.29.2/contrib/completion/bash/docker-compose > /etc/bash_completion.d/docker-compose" \
    && sh -c "curl -L https://raw.githubusercontent.com/docker/cli/v20.10.13/contrib/completion/bash/docker > /etc/bash_completion.d/docker" \
    && sh -c "curl -L https://raw.githubusercontent.com/git/git/v2.35.1/contrib/completion/git-completion.bash > /etc/bash_completion.d/git"

# Go back to regular user
USER $MAMBA_USER

# Install git prompt
# NOTE(hadim): likely a bit hacky
COPY --chown=$MAMBA_USER:$MAMBA_USER .devcontainer/bashrc /tmp/bashrc
RUN cat /tmp/bashrc >> ~/.bashrc

# Install the Conda packages
ARG CONDA_ENV_FILE
COPY --chown=$MAMBA_USER:$MAMBA_USER ${CONDA_ENV_FILE} /tmp/${CONDA_ENV_FILE}
RUN micromamba install --yes --name base --file /tmp/${CONDA_ENV_FILE} \
    && micromamba clean --all --yes

# Activate the conda environment for the Dockerfile.
# <https://github.com/mamba-org/micromamba-docker#running-commands-in-dockerfile-within-the-conda-environment>
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Create and set the workspace folder
ARG CONTAINER_WORKSPACE_FOLDER=/workspaces/default-workspace-folder
RUN mkdir -p "${CONTAINER_WORKSPACE_FOLDER}"
WORKDIR "${CONTAINER_WORKSPACE_FOLDER}"
