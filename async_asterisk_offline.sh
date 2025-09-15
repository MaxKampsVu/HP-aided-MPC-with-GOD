#!/bin/bash
# set -x

# Usage: ./../async_asterisk_offline.sh <g> <d> <players> <delay_party_count> <run_opt>
# g: number of multiplication gates at each level
# d: multiplication depth of the circuit
# players: total number of parties
# delay_party_count: number of parties that will be delayed
# run_opt: 0 for Honest Majority, 1 for Dishonest Majority
# Example: ./../async_asterisk_offline.sh 10 10 5 2 0

delay_party_count=$4
run_opt=$5
latency=200 # latency in milliseconds

threads=64

if test $run_opt = 0
then
	echo "Running Asynchronous Asterisk offline phase (Honest Majority)"
	echo "*****************************************************************"
    pkill -f "hm_async_asterisk_offline"
	run_app=./benchmarks/hm_async_asterisk_offline
	dir=~/benchmark_data/hm_async_asterisk_offline
else
	echo "Running Asynchronous Asterisk offline phase (Dishonest Majority)"
	echo "*****************************************************************"
    pkill -f "dm_async_asterisk_offline"
	run_app=./benchmarks/dm_async_asterisk_offline
	dir=~/benchmark_data/dm_async_asterisk_offline
fi

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

    pkill -f $run_app
done

done
