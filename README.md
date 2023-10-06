# Code for the paper "Limits to evolutionary rescue by conjugative plasmids"

The C++ scripts require functions from the [GNU Scientific Library](https://www.gnu.org/software/gsl/)


## Generate a single data point

#### Numerical results

To generate a numerical result, go to the desired folder, for instance `Scenario_A/Model_I/`. An input file `input.dat` with the desired parameters is required in the `inputfiles/` folder. A result file `resultsNumerical.txt` is also required. Then, in a terminal, and in the desired folder, type:

```bash
$ cd numerical/src
$ make
$ cd ../../
$ numerical/bin/rescue_maintext inputfiles/input.dat
```

#### Stochastic simulation results

To generate a stochastic simulation result, go to the desired folder, for instance `Scenario_A/Model_I/`. An input file `input.dat` with the desired parameters is required in the `inputfiles/` folder. A result file `resultsSimulation.txt` is also required. Then, in a terminal, and in the desired folder, type:

```bash
$ cd simulation/src
$ make
$ cd ../../
$ simulation/bin/rescue_maintext inputfiles/input.dat
```

The source code was tested on Ubuntu 20.04.

## Generate many data points on a cluster

We also provide additional files to generate multiple data points for different parameter values on a cluster managed with [SLURM](https://slurm.schedmd.com/). Note that the `run.sh` and `submit.sh` could be adapted to run on cluster with different management systems.
The different parameter values can be modified in the `inputfileCreation.py` file with generate a input file for each parameter set. In the `submit.sh`, one can decide to generate numerical and/or stochastic simulation results.
To generate results, go to the desired folder, for instance `Scenario_A/Model_I/`. Then, on the cluster terminal, type:

```bash
$ sh run.sh
```
