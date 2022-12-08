import os
import sys
import numpy as np


# Parameters
betaV = np.logspace(np.log10(1e-9), np.log10(1e-4), 100)[:-1]  # for numerical results
# betaV = np.logspace(np.log10(1e-9), np.log10(1e-7), 16)[:-1]   # for stochastic simulations
phi0V = [0, 1]
gammaV = [1]
sMV = [0.1]
cV = [0.09]
sigmaV = [0]
s0V = [-0.1]
alphaV = [1e-6]
N0V = [1e6]
u = 1e-6



# Deleting all previous inputfiles
dir = 'inputfiles/'
for f in os.listdir(dir):
    os.remove(os.path.join(dir, f))



# Writing new inputfiles
count = 1
for beta in betaV:
    for phi0 in phi0V:
        for gamma in gammaV:
            for sM in sMV:
                for c in cV:
                    for sigma in sigmaV:
                        for s0 in s0V:
                            for alpha in alphaV:
                                for N0 in N0V:
                                    filepath = dir + "input_" + str(count) + ".dat"
                                    with open(filepath, 'w') as f:
                                        f.write("%e #beta\n" % beta)
                                        f.write("%f #phi0\n" % phi0)
                                        f.write("%f #gamma\n" % gamma)
                                        f.write("%f #sM\n" % sM)
                                        f.write("%f #c\n" % c)
                                        f.write("%f #sigma\n" % sigma)
                                        f.write("%e #alpha\n" % alpha)
                                        f.write("%e #u\n" % u)
                                        f.write("%d #N0\n" % N0)
                                        f.write("%f #s0\n" % s0)
                                    count += 1



# Returning the number of files on cluster
print(count - 1)
