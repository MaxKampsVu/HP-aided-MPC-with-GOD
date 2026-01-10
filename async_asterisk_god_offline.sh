#!/bin/bash
# set -x

# Usage: ./../async_asterisk_mpc.sh <g> <d> <players> <delay_party_count>
# g: number of multiplication gates at each level
# d: multiplication depth of the circuit
# players: total number of parties
# delay_party_count: number of parties that will be delayed
# Example: ./../async_asterisk_mpc.sh 10 10 5 2

delay_party_count=$4
latency=200 # latency in milliseconds

threads=64

# Utility to check if TC rules are applied
verify_tc() {
    local expected=$1
    local actual_config=$(tc qdisc show dev lo)
    
    # Special case for LAN: check for 0.5ms OR 499us OR 500us
    if [[ "$expected" == "0.5ms" ]]; then
        if [[ "$actual_config" == *"0.5ms"* || "$actual_config" == *"499us"* || "$actual_config" == *"500us"* ]]; then
             echo "[SUCCESS] Verified LAN network delay (approximated as 0.5ms/500us)"
             return 0
        fi
    fi

    # Standard check for other values
    if [[ "$actual_config" == *"$expected"* ]]; then
        echo "[SUCCESS] Verified network delay: $expected"
    else
        echo "[ERROR] TC Configuration failed! Could not find $expected in config."
        echo "Current config: $actual_config"
        exit 1
    fi
}

function tc_off() {
    echo "Cleaning up network traffic control..."
    # Silence the error if no qdisc exists to delete
    sudo tc qdisc del dev lo root 2>/dev/null || true
}

function tc_lan() {
    tc_off
    echo "Applying LAN configuration (1Gbit, 0.5ms delay)..."
    sudo tc qdisc add dev lo root handle 1:0 htb default 10 && \
    sudo tc class add dev lo parent 1:0 classid 1:10 htb rate 1Gbit burst 15k && \
    sudo tc qdisc add dev lo parent 1:10 handle 10:0 netem delay 0.5ms 0.03ms 5% distribution normal
    
    verify_tc "0.5ms"
}

function tc_man() {
    tc_off
    echo "Applying MAN configuration (500Mbit, 10ms delay)..."
    sudo tc qdisc add dev lo root handle 1:0 htb default 10 && \
    sudo tc class add dev lo parent 1:0 classid 1:10 htb rate 500Mbit burst 15k && \
    sudo tc qdisc add dev lo parent 1:10 handle 10:0 netem delay 10.0ms
    
    verify_tc "10.0ms"
}

function tc_wan() {
    tc_off
    echo "Applying WAN configuration (100Mbit, 50ms delay)..."
    sudo tc qdisc add dev lo root handle 1:0 htb default 10 && \
    sudo tc class add dev lo parent 1:0 classid 1:10 htb rate 100Mbit burst 15k && \
    sudo tc qdisc add dev lo parent 1:10 handle 10:0 netem delay 50.0ms 3ms 25% distribution normal
    
    verify_tc "50.0ms"
}

echo "Running asynchronous Asterisk GOD offline phase (Dishonest Majority)"
echo "*****************************************************************"
pkill -f "dm_async_asterisk_god_offline"
run_app=./benchmarks/dm_async_asterisk_god_offline
dir=~/benchmark_data/dm_async_asterisk_god_offline

# rm -rf $dir/*.log $dir/g*.json
# mkdir -p $dir

num_repeat=1

for repeat in $(seq 1 $num_repeat)
do

# for players in {5,10,20,30}
for players in $3
do
    tc_lan

    for party in $(seq 1 $players)
    do
        log=$dir/g_$1_d_$2_$party.log
        json=$dir/g_$1_d_$2_$party.json

        if test $party -gt $(($players - $delay_party_count))
        then
            if test $party = $players
            then
                $run_app -p $party --localhost -g $1 -d $2 -n $players -l $latency 2>&1 >> $log &
            else
                $run_app -p $party --localhost -g $1 -d $2 -n $players -l $latency 2>&1 > /dev/null &
            fi
        else
            if test $party = $players
            then
                $run_app -p $party --localhost -g $1 -d $2 -n $players 2>&1 >> $log &
            else
                $run_app -p $party --localhost -g $1 -d $2 -n $players 2>&1 > /dev/null &
            fi
        fi
        
        codes[$party]=$!
    done

    $run_app -p 0 --localhost -g $1 -d $2 -n $players -t $threads 2>&1 | tee -a $dir/g_$1_d_$2_0.log &
    codes[0]=$!

    for party in $(seq 0 $players)
    do
        wait ${codes[$party]} || return 1
    done
    tc_off

    pkill -f $run_app
done

done
