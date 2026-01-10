#!/bin/bash
# set -x

delay_party_count=$4
latency=200 
threads=64

echo "Running Asynchronous Asterisk GOD online phase (Dishonest Majority)"
echo "*****************************************************************"

pkill -f "dm_sync_asterisk_god_offline"
run_app=./benchmarks/dm_sync_asterisk_god_online
dir=~/benchmark_data/dm_sync_asterisk_god_online

mkdir -p $dir

num_repeat=1

for repeat in $(seq 1 $num_repeat)
do
for players in $3
do
    # --- Launch Parties 1 to N ---
    for party in $(seq 1 $players)
    do
        log=$dir/g_$1_d_$2_$party.log
        if test $party -gt $(($players - $delay_party_count))
        then
            $run_app -p $party --localhost -g $1 -d $2 -n $players -l $latency > "$log" 2>&1 &
        else
            $run_app -p $party --localhost -g $1 -d $2 -n $players > "$log" 2>&1 &
        fi
        codes[$party]=$!
    done

    # --- Launch Party 0 (HP) ---
    $run_app -p 0 --localhost -g $1 -d $2 -n $players -t $threads 2>&1 | tee "$dir/g_$1_d_$2_0.log" &
    codes[0]=$!

    # --- Wait for completion ---
    for party in $(seq 0 $players)
    do
        wait ${codes[$party]}
    done

    # --- SUMMARY SECTION ---
    echo -e "\n--- Final Results for All Parties ---"
    
    fast_total_time=0
    fast_count=0
    all_total_time=0
    all_count=0

    # Determine the ID threshold for delayed parties
    # Delayed parties are those with ID > (players - delay_party_count)
    threshold=$(awk "BEGIN {print $players - $delay_party_count}")

    for party in $(seq 0 $players)
    do
        log_file=$dir/g_$1_d_$2_$party.log
        if [ -f "$log_file" ]; then
            echo "Party $party:"
            grep -E "pid:|time:|sent:" "$log_file" | sed 's/^/  /'
            
            # Extract time for averaging (excluding Party 0)
            if [ $party -ne 0 ]; then
                p_time=$(grep "time:" "$log_file" | awk '{print $2}')
                
                if [ ! -z "$p_time" ]; then
                    # Track all non-zero parties
                    all_total_time=$(awk "BEGIN {print $all_total_time + $p_time}")
                    all_count=$((all_count + 1))

                    # Track only non-delayed parties
                    if [ "$party" -le "$threshold" ]; then
                        fast_total_time=$(awk "BEGIN {print $fast_total_time + $p_time}")
                        fast_count=$((fast_count + 1))
                    fi
                fi
            fi
        fi
    done

    # --- PRINT AVERAGE ---
    echo -e "\n**************************************"
    if [ $fast_count -gt 0 ]; then
        # Average of non-delayed parties
        avg_time=$(awk "BEGIN {printf \"%.4f\", $fast_total_time / $fast_count}")
        echo "Average Time (Non-delayed Parties): $avg_time ms"
        echo "Counted $fast_count fast parties."
    elif [ $all_count -gt 0 ]; then
        # Fallback: All parties were delayed
        avg_time=$(awk "BEGIN {printf \"%.4f\", $all_total_time / $all_count}")
        echo "Average Time (All parties delayed): $avg_time ms"
        echo "Counted $all_count total parties."
    else
        echo "[!] No party timing data found."
    fi
    echo "**************************************"

    pkill -f $run_app
done
done