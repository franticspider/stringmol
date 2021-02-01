

# Intro

Stringmol has a tangled development history. It begain in the Plazzmid project, but has mostly been developed by me (Simon Hickinbotham) in my own time ever since (apart from 2 years further work during EvoEvo). When I started, I knew little about version control, tdd etc. 

# History

- Particle metabolome
- molecular microprograms
- Diversity from a monoculture
- reaction logic
- Degeneracy enriches
- conservation of matter
- maximize adjacent possible
- extinction in replicase models



# Stringmol flavours

- **Vanilla:** Used on most/all stringmol papers published before 2016
- **RUTSAC:** Developed in collaboration with Adam Nellis for the stringmol web app
- **Spatial:** Developd in collaboration with Paulien Hogeweg for R-P evolution - *Actively under development*
- **FPGA:** Developed in collaboartion with Matthew Rowlings/Martin Trezfer for running on their Centurion hardware - *Currently external to this repository*
- **RStringmol:** Analysis tools for spatial stringmol *Currently external to this repository*


# Tasks for release 0.2.3.6

## Functionality

- Turn off bind caching - have an option for this in the config file


## Compiler

convert all .c files to .cpp - makes compiling easier and improves compatability with e.g. Rstringmol

## Make targets

The makefile contains a number of different executable targets that need to be resolved, particularly `smspatial`. First, let's document what we have, and then we can think about refactoring

### Make targets in version 2.3.5

*file paths are relative to the `src` directory

- **`../release/stringmol`** the release version of stringmol
- **`../debug/stringmol`** the debug version of stringmol. 

## Tests

Testing began with version 0.2.3.4 but isn't really in a good place at the moment. 
Consider using https://github.com/catchorg/Catch2

## Warnings

DONE, 200607: Remove all compiler warnings

## Valgrind

Check that the spatial version doesn't leak

## Documentation

Commence summary documentation




