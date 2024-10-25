# Use the official C++ development container image
FROM mcr.microsoft.com/devcontainers/cpp:ubuntu-22.04

# Install required dependencies
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y --no-install-recommends \
        libarmadillo-dev \
        libsuperlu-dev \
        libhdf5-dev \
        libopenblas-dev \
        liblapack-dev \
        libarpack++2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# RUN mkdir /workspaces && chown -R vscode:vscode /workspaces

# USER vscode

# COPY .devcontainer/entrypoint.sh /entrypoint.sh
# RUN chmod +x /entrypoint.sh
# ENTRYPOINT ["/entrypoint.sh"]
