#!/bin/bash

# Creating result files if not already existing
if [ ! -f resultsNumerical.txt ]; then
	cat > resultsNumerical.txt <<-EOF
	"Beta" "Phi" "S" "Cost" "SegregationLoss" "Alpha" "N0" "S0" "TimeEst" "TimeMaxIntegrand" "ResolutionIntegrand" "pENumerical"
	EOF
fi

if [ ! -f resultsSimulation.txt ]; then
	cat > resultsSimulation.txt <<-EOF
	"Beta" "Phi" "S" "Cost" "SegregationLoss" "Alpha" "N0" "S0" "TimeEst" "NSim" "pESimulation"
	EOF
fi



cd numerical/src
make
cd ../../simulation/src
make
cd ../../



nb=$(python inputfileCreation.py 2>&1)
echo $nb
for i in $(LC_ALL=C seq 1 1 $nb)
do
	sbatch submit.sh $i
	# If results arrive too fast, there can be conflict in writing in result files.
	# In this case, use the sleep function:
	# sleep 0.2
done
