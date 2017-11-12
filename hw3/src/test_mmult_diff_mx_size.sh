#!/bin/bash
../bin/mmult 100 1000 > ../data/n100_r1000_ord1.csv
echo "100 done"
../bin/mmult 200 100 > ../data/n200_r100_ord1.csv
echo "200 done"
../bin/mmult 500 20 > ../data/n500_r20_ord1.csv
echo "500 done"
../bin/mmult 1000 20 > ../data/n1000_r20_ord1.csv
echo "1000 done"
../bin/mmult 2000 10 > ../data/n2000_r10_ord1.csv
echo "2000 done"
