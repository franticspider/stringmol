

#Testing and Refactoring

Stringmol is research code, but the time has come to improve its organisation.
We need to standardise parameter loading on the following run types:

0: 1 on 1
1: ALife XII
2: Con pop
4: Comass GA

##Parameters

These are the values that are used in the config files (missing parameters are standardised already)

###RANDSEED

0: not set
1: initmyrand(-1) sets rseed
2: longinitmyrand() based on seedin sets rseed 
4: longinitmyrand() based on seedin sets rseed 

###NSTEPS / MAXNSTEPS

0: arg_load() sets A.nsteps
1: readordef_param_int() sets maxnsteps
2: A[0]->nsteps is set, but not used; MAXCONSTEPS is a hard coded const, value 10,000,000
4: setupSMol() sets R.maxnsteps using read_param_int()


###NTRIALS

0: not set
1: readordef_param_int() sets rlim 
2: not needed (?)
4:

###GAQNN

0: not needed
1: not needed
2: not needed
4: read_param_int

###GRANULAR

_This is used for the `granular_SM' experiment in the ALife Journal paper

0: not set
1: readordef_param_int() sets A.granular_1 
2: not set
4: not set


###NCON

0: not needed
1: not needed
2: NCON initialised to 4, then read_param_int("NCONTAINERS") 
4: not needed (but similar to popsize)


