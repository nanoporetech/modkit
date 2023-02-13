FROM rust@sha256:02a53e734724bef4a58d856c694f826aa9e7ea84353516b76d9a6d241e9da60e

RUN apt-get update && apt-get install -y cmake

WORKDIR /app
COPY src ./src
COPY Cargo.toml .
COPY rust-toolchain.toml .
COPY .rustfmt.toml .
COPY tests ./tests
