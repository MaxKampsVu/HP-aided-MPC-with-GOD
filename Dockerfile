FROM ubuntu:focal

RUN apt-get update

# Ensure tzdata installation (a dependency) is non-interactive.
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get install -y sudo vim software-properties-common wget build-essential cmake git python3-dev xxd libgmp-dev libssl-dev iproute2
RUN apt-get install -y --no-install-recommends libntl-dev

# Install nlohmann_json
WORKDIR /home
RUN git clone https://github.com/nlohmann/json.git \
  && cd json \
  && mkdir build && cd build \
  && cmake -DCMAKE_BUILD_TYPE=Release .. \
  && make \
  && make install

# Install boost
WORKDIR /home
RUN wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
RUN tar -xzvf boost_1_86_0.tar.gz && rm boost_1_86_0.tar.gz
RUN cd /home/boost_1_86_0 && ./bootstrap.sh --with-python=/usr/bin/python3
RUN cd /home/boost_1_86_0 && ./b2 install

# Install emp-tool and emp-ot
WORKDIR /home
RUN wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py \
  && python3 install.py --deps --tool \
  && python3 install.py --install --tool --ot

CMD ["bash"]
