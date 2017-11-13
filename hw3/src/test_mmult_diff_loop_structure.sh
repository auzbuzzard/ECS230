#!/bin/bash

# # loop is compiled with no optimization options
# ../bin/loop 100 100 > ../data/n100_r100_ordAll.csv
# echo "100 done"
# ../bin/loop 250 50 > ../data/n250_r50_ordAll.csv
# echo "250 done"

# # loop2 is compiled with the -O2 optimization option
# ../bin/loop2 100 1000 > ../data/n100_r1000_ordAll_O2opt.csv
# echo "100 done"
# ../bin/loop2 250 500 > ../data/n250_r500_ordAll_O2opt.csv
# echo "250 done"
# ../bin/loop2 500 200 > ../data/n500_r200_ordAll_O2opt.csv
# echo "500 done"
# ../bin/loop2 2000 100 > ../data/n2000_r100_ordAll_O2opt.csv
# echo "2000 done"

# loop3 is compiled with the -O3 optimization option
../bin/loop3 100 1000 > ../data/n100_r1000_ordAll_O3opt.csv
echo "100 done"
../bin/loop3 250 500 > ../data/n250_r500_ordAll_O3opt.csv
echo "250 done"
../bin/loop3 500 200 > ../data/n500_r200_ordAll_O3opt.csv
echo "500 done"
../bin/loop3 2000 100 > ../data/n2000_r100_ordAll_O3opt.csv
echo "2000 done"
