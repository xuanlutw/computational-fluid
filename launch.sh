#!/bin/bash

for x in {3..15}
do
    ./main 1.0 1.0 1.0 ${x} &
done
