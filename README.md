# Asynchronous HP-aided MPC with GOD 

This directory contains the implementation of the Synchronous/Asynchronous HP-aided MPC with GOD dishonest majority protocols Alhena and Wasat.
The protocol is implemented in C++17 and [CMake](https://cmake.org/) is used as the build system.

**NOTICE:** This is an academic proof-of-concept prototype and has not received careful code review. This implementation is NOT ready for production use.

## External Dependencies
The following libraries need to be installed separately and should be available to the build system and compiler.

- [GMP](https://gmplib.org/)
- [NTL](https://www.shoup.net/ntl/) (11.0.0 or later)
- [Boost](https://www.boost.org/) (1.72.0 or later)
- [Nlohmann JSON](https://github.com/nlohmann/json)
- [EMP Tool](https://github.com/emp-toolkit/emp-tool)
- [EMP OT](https://github.com/emp-toolkit/emp-ot)

### Docker
All required dependencies to compile and run the project are available through the docker image.
To build and run the docker image, execute the following commands from the root directory of the repository:

```sh
# Build the Asynchronous HP-aided MPC Docker image.
#
# Building the Docker image requires at least 4GB RAM. This needs to be set 
# explicitly in case of Windows and MacOS.
docker build -t async-hp-aided-mpc .

# Create and run a container.
#
# This should start the shell from within the container.
docker run --cap-add=NET_ADMIN -it -v $PWD:/code async-hp-aided-mpc

# The following command changes the working directory to the one containing the 
# source code and should be run on the shell started using the previous command.
cd /code
```

## Compilation
The project uses [CMake](https://cmake.org/) for building the source code. 
To compile, run the following commands from the root directory of the repository:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..

# The two main targets are 'benchmarks' and 'tests' corresponding to
# binaries used to run benchmarks and unit tests respectively.
make <target>
```

## Usage 


- `benchmarks/alhena_mpc`: Benchmark the performance of synchronous HP-aided MPC protocol with god (both offline and shared online phases) Alhena by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/alhena_offline`: Benchmark the performance of synchronous HP-aided MPC protocol (only offline phase) Alhena by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/alhena_online`: Benchmark the performance of the shared online phase synchronous/asynchronous HP-aided MPC protocol Alhena/Wasat by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/wasat_mpc`: Benchmark the performance of asynchronous HP-aided MPC protocol (both offline and shared online phases) Wasat by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/wasat_offline`: Benchmark the performance of asynchronous HP-aided MPC protocol (only offline phase) Wasat by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/wasat_online`: Benchmark the performance of the shared online phase synchronous/asynchronous HP-aided MPC protocol Alhena/Wasat by evaluating a circuit with a given depth and number of multiplication gates at each depth.

- `benchmarks/alhena_online_cheaters`: Benchmark the performance of the shared online phase synchronous/asynchronous HP-aided MPC protocol Alhena/Wasat by evaluating a circuit with a given depth and number of multiplication gates at each depth in the presence of cheaters. A cheater sends invalid data to the HP in every round of the protocol.

- `tests/*`: These programs contain unit tests for various parts of the codebase. 


Execute the following commands from the `build` directory created during compilation to run the programs:
```sh
# Benchmark Asynchronous HP-aided MPC MPC for honest-majority setting.
#
# The command below should be run on four different terminals with $PID set to
# 0, 1, 2, and 3 i.e., one instance corresponding to each party.
#
# The number of threads can be set using the '-t' option. '-g' denotes the 
# number of multiplication gates at each level, '-d' denotes the multiplicative depth 
# of the circuit and '-n' the number of parties participating in the protocol.
#
# The program can be run on different machines by replacing the `--localhost`
# option with '--net-config <net_config.json>' where 'net_config.json' is a
# JSON file containing the IPs of the parties. A template is given in the
# repository root.¨¨
./benchmarks/mpc -p $PID --localhost -g 100 -d 10 -n 5

# The `mpc` script in the repository root can be used to run the
# programs for all parties from the same terminal.
# For example, the previous benchmark can be run using the script as shown below.
# Here the last argument stands for either Alhena (0) or
# Wasat (1) and the preceding arguments stand for: 
#   - gates per circuit layer: 100
#   - number of circuit layers: 10
#   - number of parties: 5
#   - number of delayed parties: 0
#   - protocol: Alhena (0) 
./../mpc.sh 100 10 5 0 0

# Additionaly, we provide scripts with the same arguments for just the offline phase (`offline.sh`) 
# and just the online phase (`online.sh`). 

# For the `online` script the user can additionally specify the number of corrupt parties.
# A corrupt party sends invalid data to the HP in every round.
# Once it is caught cheating by the HP, the HP closes the communication channel to that party.
# The number of corrupt parties can be specified in the last argument, for example 3
./../online.sh 100 10 5 0 0 3

```
