# Asynchronous HP-aided MPC with GOD 

The asynchronous version of our new Protocol can be executed with the script:
```sh
async_asterisk_god_mpc.sh
```
The script takes the same arguments as the original protocols with abort. 

## Location of the implementation: 
- dmgod_offline_evaluator.cpp
- dmgod_online_evaluator.cpp

Note, that the online phase remains the same for both the synchronous and asynchronous version of our protocol. Only the offline phase changes in (`dmgod_offline_evaluator.cpp`).
To run the synchronous version of the protocol set the member variable `run_async_` in `dmgod_offline_evaluator.cpp` to false and recompile. (I will provide a separate script once I know how to set up a synchronous network environment). 

The new sharing type (which macs only the parties share) is located in `dmgod_sharing`. 

Moreover, for the new protocol the `ot_provider.cpp` been extended in `ot_provider_ha.cpp` (ot provider with honest abort) to return the hash of the messages sent during the OLE. Note, that some messages are not included in the hash yet, due to me still figuring out how to instantiate the Ferret OT instances such that they produce the same intermediate messages for all Parties. 

# Asynchronous HP-aided MPC (README of forked repo)

This directory contains the implementation of the Asyncronous HP-aided MPC (accepted in IEEE S&P 2026) fair protocols for both honest-majority and dishonest-majority settings.
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
c
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
A short description of the compiled programs is given below.
All of them provide detailed usage description on using the `--help` option.

- `benchmarks/hm_async_asterisk_mpc`: Benchmark the performance of Asynchronous HP-aided MPC protocol (both offline and online phases) in honest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/hm_async_asterisk_offline`: Benchmark the performance of Asynchronous HP-aided MPC protocol (only offline phase) in honest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/hm_async_asterisk_online`: Benchmark the performance of Asynchronous HP-aided MPC protocol (only online phase) in honest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/dm_async_asterisk_mpc`: Benchmark the performance of Asynchronous HP-aided MPC protocol (both offline and online phases) in dishonest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/dm_async_asterisk_offline`: Benchmark the performance of Asynchronous HP-aided MPC protocol (only offline phase) in dishonest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/dm_async_asterisk_online`: Benchmark the performance of Asynchronous HP-aided MPC protocol (only online phase) in dishonest-majority setiing by evaluating a circuit with a given depth and number of multiplication gates at each depth.
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
# repository root.
./benchmarks/hm_async_asterisk_mpc -p $PID --localhost -g 100 -d 10 -n 5

# The `async_asterisk_mpc` script in the repository root can be used to run the
# programs for all parties from the same terminal.
# For example, the previous benchmark can be run using the script as shown below.
# Here the last argument stands for either honest-majority setting (0) or
# dishonest majority setting (1) and the second last argument denotes the number 
# of delayed parties
./../async_asterisk_mpc.sh 100 10 5 0 0

# All other benchmark programs have similar options and behaviour. The '-h'
# option can be used for detailed usage information.

# Benchmark Asynchronous HP-aided MPC MPC for dishonest-majority setting.
./../async_asterisk_mpc.sh 100 10 5 0 1

# Benchmark Asynchronous HP-aided MPC offline phase for honest-majority setting.
./../async_asterisk_offline.sh 100 10 5 0 0

# Benchmark Asynchronous HP-aided MPC offline phase for dishonest-majority setting.
./../async_asterisk_offline.sh 100 10 5 0 1

# Benchmark Asynchronous HP-aided MPC online phase for honest-majority setting.
./../async_asterisk_online.sh 100 10 5 0 0

# Benchmark Asynchronous HP-aided MPC online phase for dishonest-majority setting.
./../async_asterisk_online.sh 100 10 5 0 1
```





```
#!/bin/bash
# set -eu

# Usage: ./run_sync_mpc.sh <g> <d> <players>
# g: number of multiplication gates at each level
# d: multiplication depth of the circuit
# players: total number of parties
# Example: ./run_sync_mpc.sh 10 10 5

# --- SINGULARITY AND PATH CONFIGURATION ---
# Define paths relative to the current directory (where the SIF file is)
SIF_FILE="./my-image.sif"
# Assuming project root is in a subdirectory named 'HP-aided-MPC-with-GOD'
ASTERISK_ROOT="./HP-aided-MPC-with-GOD" 

# Absolute path to the compiled binary on the host filesystem
BINARY_NAME=""
run_app="/home/jjb140/${BINARY_NAME}" 

# Check for essential files
if [ ! -f "$SIF_FILE" ]; then
    echo "FATAL: Singularity image not found at $SIF_FILE"
    exit 1
fi
if [ ! -f "$run_app" ]; then
    echo "FATAL: Compiled binary not found at $run_app. Did you compile with the SIF?"
    exit 1
fi
# --- END CONFIGURATION ---

# Arguments (using descriptive variable names)
gates=$1
depth=$2
players=$3

threads=64

echo "Running Synchronous Asterisk GOD MPC (Dishonest Majority)"
echo "*****************************************************************"
pkill -f "$BINARY_NAME"
dir=~/benchmark_data/${BINARY_NAME} # Logs still go to your home directory

# rm -rf $dir/*.log $dir/g*.json
mkdir -p "$dir"

num_repeat=1

for repeat in $(seq 1 $num_repeat)
do

for players in $players
do
    echo "Starting synchronous benchmark for N=$players parties..."

    # Launch parties 1 to N
    for party in $(seq 1 $players)
    do
        log="$dir/g_${gates}_d_${depth}_p${party}.log"
        json="$dir/g_${gates}_d_${depth}_p${party}.json"

        # Base command for the binary (using absolute path on the host)
        BASE_CMD="$run_app -p $party --localhost -g $gates -d $depth -n $players"
        
        # --- WRAP EXECUTION WITH SINGULARITY ---
        # The key is to use the -B bind mount option to ensure the container can find the binary path
        
        if test $party = 1
        then
            # Party 1: Output to log and save JSON
            singularity exec -B $PWD:$PWD $SIF_FILE $BASE_CMD -o "$json" 2>&1 >> "$log" &
        else
            # Other Parties: Suppress output
            singularity exec -B $PWD:$PWD $SIF_FILE $BASE_CMD 2>&1 > /dev/null &
        fi
        
        codes[$party]=$!
    done

    # --- WRAP EXECUTION FOR PARTY 0 ---
    # Party 0: Output to log and screen, set threads
    FULL_CMD="$run_app -p 0 --localhost -g $gates -d $depth -n $players -o $dir/g_${gates}_d_${depth}_p0.json -t $threads"
    singularity exec -B $PWD:$PWD $SIF_FILE $FULL_CMD 2>&1 | tee -a "$dir/g_${gates}_d_${depth}_p0.log" &
    codes[0]=$!

    for party in $(seq 0 $players)
    do
        wait ${codes[$party]} || exit 1 # FIX: Changed 'return 1' to 'exit 1'
    done
    echo "Benchmark for N=$players completed."

    pkill -f "$BINARY_NAME"
done

done
```