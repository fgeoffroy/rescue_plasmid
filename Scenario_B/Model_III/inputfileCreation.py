import os
import sys
import numpy as np


# Parameters
gammaV = [1]
sPV = [0.1]
sMV = [0.1]
cV = [0.001, 0.01, 0.05, 0.09]
sigmaV = [0]
alphaV = [1e-6]
s0V = [-0.1]
u = 1e-6
betaMin = 1e-9
betaMax = 1e-4
betaNbSteps = 20    # for numerical results (careful, also change in the loop below)
# betaNbSteps = 11    # for stochastic simulations (careful, also change in the loop below)










# Deleting all previous inputfiles
dir = 'inputfiles/'
for f in os.listdir(dir):
    os.remove(os.path.join(dir, f))



# Writing new inputfiles
count = 1
for gamma in gammaV:
    for alpha in alphaV:
        for sP in sPV:
            for s0 in s0V:
                for sM in sMV:
                    for c in cV:
                        for sigma in sigmaV:

                            ### For numerical results ###
                            betaFree = alpha * (c + sigma * (1 + sP)) / sP
                            betaV1 = np.logspace(np.log10(betaMin), np.log10(betaFree), betaNbSteps)
                            betaV1[0] = betaMin
                            betaV1[betaNbSteps-1] = betaFree
                            if sigma == 0:
                                betaWild = c * alpha / (sP - c)
                                betaV2 = np.logspace(np.log10(betaFree), np.log10(betaWild), betaNbSteps)
                                betaV3 = np.logspace(np.log10(betaWild), np.log10(betaMax), betaNbSteps)
                                betaV2[0] = betaFree
                                betaV3[0] = betaWild
                                betaV2[betaNbSteps-1] = betaWild
                                betaV3[betaNbSteps-1] = betaMax
                                betaV = np.concatenate((betaV1[:-1], betaV2[:-1], betaV3))
                            else:
                                betaV2 = np.logspace(np.log10(betaFree), np.log10(betaMax), betaNbSteps)
                                betaV2[0] = betaFree
                                betaV2[betaNbSteps-1] = betaMax
                                betaV = np.concatenate((betaV1[:-1], betaV2))

                            ### For stochastic simulations ###
                            # betaV = np.logspace(np.log10(betaMin), np.log10(betaMax), betaNbSteps)
                            # betaV[0] = betaMin
                            # betaV[betaNbSteps-1] = betaMax

                            ### For specific values of beta (Figures with sigma in SI) ###
                            # betaV = [1e-7, 1e-6, 1e-5, 1e-4]

                            for beta in betaV:
                                filepath = dir + "input_" + str(count) + ".dat"
                                with open(filepath, 'w') as f:
                                    f.write("%e #beta\n" % beta)
                                    f.write("%f #gamma\n" % gamma)
                                    f.write("%f #sP\n" % sP)
                                    f.write("%f #sM\n" % sM)
                                    f.write("%f #c\n" % c)
                                    f.write("%f #sigma\n" % sigma)
                                    f.write("%e #alpha\n" % alpha)
                                    f.write("%e #u\n" % u)
                                    f.write("%f #s0\n" % s0)
                                count += 1



# Returning the number of files on cluster
print(count - 1)
