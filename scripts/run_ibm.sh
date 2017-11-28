#!/bin/bash

# run this script from the executable ./spicy directory
echo "Running IBM spice experiments"

make clean; make


echo "IMB1 Test"
make clean_outputs > /dev/null 2>&1
./spicy benchmarks/ibm/ibm1/ibmpg1.spice > /dev/null
./scripts/lc_sort.sh nodes_op_point_all.txt
./scripts/lc_sort.sh benchmarks/ibm/ibm1/ibmpg1.solution
python ./scripts/compare_sorted.py ibmpg1_lc_sorted.txt nodes_op_point_all_lc_sorted.txt

echo "IMB2 Test"
make clean_outputs > /dev/null 2>&1
./spicy benchmarks/ibm/ibm2/ibmpg2.spice > /dev/null
./scripts/lc_sort.sh nodes_op_point_all.txt
./scripts/lc_sort.sh benchmarks/ibm/ibm2/ibmpg2.solution
python ./scripts/compare_sorted.py ibmpg2_lc_sorted.txt nodes_op_point_all_lc_sorted.txt

#echo "IMB3 Test"
#make clean_outputs > /dev/null 2>&1
#./spicy benchmarks/ibm/ibm3/ibmpg3.spice > /dev/null
#./scripts/lc_sort.sh nodes_op_point_all.txt
#./scripts/lc_sort.sh benchmarks/ibm/ibm3/ibmpg3.solution
#python ./scripts/compare_sorted.py ibmpg3_lc_sorted.txt nodes_op_point_all_lc_sorted.txt

#echo "IMB4 Test"
#make clean_outputs > /dev/null 2>&1
#./spicy benchmarks/ibm/ibm4/ibmpg4.spice > /dev/null
#./scripts/lc_sort.sh nodes_op_point_all.txt
#./scripts/lc_sort.sh benchmarks/ibm/ibm4/ibmpg4.solution
#python ./scripts/compare_sorted.py ibmpg4_lc_sorted.txt nodes_op_point_all_lc_sorted.txt

#echo "IMB5 Test"
#make clean_outputs > /dev/null 2>&1
#./spicy benchmarks/ibm/ibm5/ibmpg5.spice > /dev/null
#./scripts/lc_sort.sh nodes_op_point_all.txt
#./scripts/lc_sort.sh benchmarks/ibm/ibm5/ibmpg5.solution
#python ./scripts/compare_sorted.py ibmpg5_lc_sorted.txt nodes_op_point_all_lc_sorted.txt

#echo "IMB6 Test"
#make clean_outputs > /dev/null 2>&1
#./spicy benchmarks/ibm/ibm6/ibmpg6.spice > /dev/null
#./scripts/lc_sort.sh nodes_op_point_all.txt
#./scripts/lc_sort.sh benchmarks/ibm/ibm6/ibmpg6.solution
#python ./scripts/compare_sorted.py ibmpg6_lc_sorted.txt nodes_op_point_all_lc_sorted.txt
