FROM ubuntu:latest

# Install the dependencies.
RUN apt-get update -qq \
 && apt-get install -qq git vim wget python3 build-essential numactl libjemalloc-dev

# Copy the current directory to a location within the container & move there
COPY . /root/cpam
WORKDIR /root/cpam
