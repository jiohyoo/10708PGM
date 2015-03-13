syntheticExp.m: toy dataset demo - 1st order optimization method
syntheticExp2.m: same as syntheticExp.m except optimization method - uses PNOPT as suggested in the Jason's paper

syntheticExp2_PGM.m: same toy dataset main script - we need to fill in this!
Three functions we need to fill in:
1. flhoodv5_PGM.m: negative log likelihood computation (w/o regularizers)
2. fgradlhoodv5_PGM.m: gradient computation for the parameters (w/o regularizers)
3. tfocsProxGroupv6_PGM.m: regularizer terms go here
