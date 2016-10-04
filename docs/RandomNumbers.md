

#Introduction

Random numbers have been a problem in stringmol. One would ideally like a system where you can repeat any experiment 
regardless of whether you've specified a rand seed or not. This is achieved in Stringmol by having a variable called ```RANDSEED``` in the config file, which can be set to any number. If the option is missing, the system should detect it, and use /dev/random as the seed generator. 

The RNG used in Stringmol is the Mersenne Twister. It's seeding this that has been the problem. 


#Current Method 

##Users

*This was developed for spatial stringmol*

###Step 1:

Specify a random number seed in the config file like so:
```
RANDSEED 234
```
would use the number 234 as the seed. 

You can explicitly specify that a new seed needs to be generated like so: 
```
RANDSEED 0
```

The simulation will save config files at various points. The random number seed will be recorded in that file. 

##Current Method - Developers

###Initialising


###Storing





#Earlier Attempts: 