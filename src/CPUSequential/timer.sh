#!/bin/bash

REP=15
for n in {4000,}; do
    for i in $(seq 1 $REP); do
        ./timer.out $n; 
    done
done
