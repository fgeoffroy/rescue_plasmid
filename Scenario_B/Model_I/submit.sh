#!/bin/bash

#SBATCH -J scenario_B_model_I # job name
#SBATCH -t 0-100:00  # time (D-HH:MM)

# for numerical results:
numerical/bin/rescue_maintext inputfiles/input_$1.dat

# for stochastic simulations:
# simulation/bin/rescue_maintext inputfiles/input_$1.dat

# Both can be used at the same time, meaning a single node will compute the
# numerical result and then the stochastic simulation result.
