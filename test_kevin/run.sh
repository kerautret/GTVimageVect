#!/bin/bash

ls input
sleep 5

i=0

for f in input/*; do
    echo -e "\033[1;32m ../build/bin/tv-triangulation-color -i $f -b 16 -D 16 -C output/$i.svg \033[0m";
    sleep 1;
    ../build/bin/tv-triangulation-color -i $f -b 16 -D 16 -C "output/${i}.svg";
    i=$((i+1))
done