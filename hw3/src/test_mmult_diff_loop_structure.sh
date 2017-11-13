#!/bin/bash
../bin/loop 100 100 > ../data/n100_r100_ordAll.csv
echo "100 done"
../bin/mmult 250 50 > ../data/n250_r50_ordAll.csv
echo "250 done"
../bin/mmult 1000 15 > ../data/n1000_r15_ordAll.csv
echo "1000 done"
