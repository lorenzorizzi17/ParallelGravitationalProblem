#!/bin/bash

REP=1
for n in {40000,}; do
    for i in $(seq 1 $REP); do
        ./timer.out $n; 
    done
done
