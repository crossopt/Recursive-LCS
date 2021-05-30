#!/bin/bash

for (( c=200000; c<=1000000000; c += 200000 ))
do
	# for cycle is too long, break when enough
	time_lz=0
	time_dp=0
	time_agrep=0
	lenf=0
	for (( i=0; i<=10; i++))
	do

	    # Generate pattern
	    pattern=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
	    # Generate compressed file TODO both types
	    printf 'a%.0s' $(seq 1 $c) > file
	    # Generation: fake texts
	    # echo $c $i | python3 zgen.py
	    # Generation: real texts
	    # st=$((($i * $c)))
	    # cat test | tail -c +$st | head -c $c > file
	    compress -f file
	    lenf=$(cat file.Z | wc -c)

	    
	    ts=$(date +%s%N) ; printf $pattern | ./build/dp_test; tt=$((($(date +%s%N) - $ts))) ; time_dp=$((time_dp + tt)); # echo dp $tt
	    ts=$(date +%s%N) ; echo $c $pattern | ./build/lz_test; tt=$((($(date +%s%N) - $ts))) ; time_lz=$((time_lz + tt)); # echo lz $tt
	    ts=$(date +%s%N) ; printf $pattern | ./run_agrep.sh; tt=$((($(date +%s%N) - $ts))) ; time_agrep=$((time_agrep + tt)); # echo agrep $tt
	done
	# / 20 / 1000000
	echo $c $lenf $time_agrep $time_lz $time_dp
done