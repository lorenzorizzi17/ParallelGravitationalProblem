#!/bin/bash

REP=10
for n in {10,100,200,500,1000}; do
    for i in $(seq 1 $REP); do
        ./timer.out $n; 
    done
done
