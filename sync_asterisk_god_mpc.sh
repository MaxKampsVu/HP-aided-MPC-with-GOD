#!/bin/bash
# set -x

# Usage: ./../async_asterisk_mpc.sh <g> <d> <players> <delay_party_count>
# g: number of multiplication gates at each level
# d: multiplication depth of the circuit
# players: total number of parties
# Example: ./../async_asterisk_mpc.sh 10 10 5 2

threads=64

echo "Running Asynchronous Asterisk GOD MPC (Dishonest Majority)"
echo "*****************************************************************"
pkill -f "dm_sync_asterisk_god_mpc"
run_app=./benchmarks/dm_sync_asterisk_god_mpc
dir=~/benchmark_data/dm_sync_asterisk_god_mpc

# rm -rf $dir/*.log $dir/g*.json
mkdir -p $dir

num_repeat=1

for repeat in $(seq 1 $num_repeat)
do

# for players in {5,10,20,30}
for players in $3
do
    for party in $(seq 1 $players)
    do
        log=$dir/g_$1_d_$2_$party.log
        json=$dir/g_$1_d_$2_$party.json

            log=$dir/g_$1_d_$2_$party.log
        json=$dir/g_$1_d_$2_$party.json
        if test $party = 1
        then
            $run_app -p $party --localhost -g $1 -d $2 -n $players -o $json 2>&1 >> $log &
        else
            $run_app -p $party --localhost -g $1 -d $2 -n $players 2>&1 > /dev/null &
        fi
        
        codes[$party]=$!
    done

    $run_app -p 0 --localhost -g $1 -d $2 -n $players -o $dir/g_$1_d_$2_0.json 2>&1 | tee -a $dir/g_$1_d_$2_0.log &
    codes[0]=$!

    for party in $(seq 0 $players)
    do
        wait ${codes[$party]} || return 1
    done

    pkill -f $run_app
done

done
