#!/bin/bash

for (( c=1000; c<=10000; c += 1000 ))
do
    # Generate pattern
    pattern=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
    # Generate compressed file
    seq  -f "a" -s '' $c > file
    compress -f file

    time_lz=$(echo $( TIMEFORMAT="%3U + %3S"; { time echo $pattern > ./build/lz_test; } 2>&1) "*1000" | bc -l)
    time_dp=$(echo $( TIMEFORMAT="%3U + %3S"; { time echo $pattern > ./build/dp_test; } 2>&1) "*1000" | bc -l)
    time_agrep=$(echo $( TIMEFORMAT="%3U + %3S"; { time echo $pattern > ./run_agrep.sh; } 2>&1) "*1000" | bc -l)

    echo $time_agrep $time_lz $time_dp
done