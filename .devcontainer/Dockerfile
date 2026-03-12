FROM julia:1.11

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Typst
RUN curl -fsSL https://github.com/typst/typst/releases/latest/download/typst-x86_64-unknown-linux-musl.tar.xz \
    | tar -xJ --strip-components=1 -C /usr/local/bin/ typst-x86_64-unknown-linux-musl/typst
