stringmol
=========

Automata Chemistry


url
===

https://github.com/franticspider/stringmol.git


Hello. You have found brief instructions to install, compile and run stringmol version 0.2.2. Please read carefully, and mail any comments or questions to sjh436@gmail.com.

If you are new to stringmol, please also see the technical report "smtr0.2.pdf" included in the same directory as this document

BUILDING STRINGMOL 
==================

THE QUICK WAY
-------------

Type the following at the command prompt:

    sh install
 
 
FURTHER INFORMATION
-------------------

Stringmol only builds in linux. Maybe one day we'll make it work on other platforms, but to be honest we are more interested in deploying it on the web and running experiments. email sjh436@gmail.com if this policy is causing you problems.

To compile stringmol on your system:

    cd src
    make all

should build an executable called stringmol in that directory. Do NOT run stringmol in there - it'll just make a mess. 
Instead, read the next section and do what it recommends. 

RUNNING STRINGMOL
=================

It is NOT a good idea to run stringmol from the directory it was compiled into, since a lot of files will be generated and it'll be hard to work out what was generated and what is already there. The best thing to do is to create a folder for runs, and execute stringmol by entering the relative path to it. 

Let's assume that we are starting from the directory containing this README. We recommend you try the following:

    mkdir test1
    cd test1
    ../src/stringmol

If all goes well, you'll get the following output:

	USAGE for TestSM

	The first argument after TestSM should be the trial type,
        followed by the arguments needed to run it

	TRIAL TYPES LIST (NUMBERS IN BRACKETS ARE EXPERIMENTAL!)
	NAME          TRIAL No:   ARGUMENTS

	1 on 1            0       .conf (.mtx)
	ALife XII         1       .conf (.mtx)
	Con Pop           2       .conf (.mtx)
	Origlife         (3)
	Comass GA        (4)
	Joinsplists      (5)
	Energetic ALXII  (6)
	swdist           (7)
	Comass ALXII     (8)
	speigmonst       (9)
	Check setup       10      .conf (.mtx)
	Finished!

This shows the current ways that stringmol can be run. More details below.

RUNNING STRINGMOL
=================

Firstly, it is important to know that not all of the trial types above actually work. If you know a little C/C++ programming, you can look through all the options if you want, or mail sjh436@gmail.com for more info. In this section, we are going to go through the "tested" options - those not in brackets in the list above.

The general method for running these trial is to type the following

PATH/stringmol TTYPE ARGS

where 'PATH' is the relative or absolute path to the strignmol-exe folder (e.g. '/home/simon/stringmol'), 'TTYPE' is the type of trial you want to run (e.g. '1'), and ARGS are the other arguments that are needed to run that particular trial. Arguments in brackets are optional.

For the currently working trials, the arguments are all the same. The first argument is the path to a config file that contains the details of the trial, and the second optional argument is to an 'mtx' file - the substitution matrices that are used in calculating alignments and creating mutations. A couple of examples of these files are contained in the 'config' subdirectory. We will one day make a guide to these files available. In the meantime, email sjh436@gmail.com with any questions. 

TTYPE 10: Check setup
--------------------

This option is there purely to check if your config file has been read
as you intended. It is a new feature in stringmol0.2.2. 

To test this out, navigate to the folder containing this README, and type the following:

$ mkdir test10; cd test10
$ ../src/stringmol 10 ../config/test_replicase.conf

The output should look like:

	Hello World! Let's run the stringmol tests!
	Argc = 3

	BEFORE loading the config, params are:
	Non-stringPM variables:
	NTRIALS     not set - the default value would be used if needed
	NSTEPS      not set - the default value would be used if needed
	CELLRAD     2500.000000 (vcellrad = 0.000000)
	AGRAD       10.000000
	ENERGY      0
	NSTEPS      -0.000012
	BLOSUM      0 size table loaded
	MUTATE      indelrate = 0.000000; subrate = 0.000000
	DECAY       0.000000
	MAXLEN      2000, (maxl0 = 2001)
	ESTEP       20
	..c'est ca!

	NTRIALS not specified;
	Setting NTRIALS to 1
	NTRIALS = 1
	Setting NSTEPS to 1200000000
	NSTEPS = 1200000000
	CELLRAD = 2500.000000
	AGRAD = 10.000000
	ENERGY = 0.000000
	NSTEPS = 1200000000.000000
	ERROR 60 on loading table
	No mutation rate found in config file
	Using ALifeXII values instead


	AFTER loading the config, params are:
	Non-stringPM variables:
	NTRIALS     1
	NSTEPS      1200000000
	CELLRAD     2500.000000 (vcellrad = 2500.000000)
	AGRAD       10.000000
	ENERGY      0
	NSTEPS      1200000000.000000
	BLOSUM      0 size table loaded
	MUTATE      indelrate = 0.000100; subrate = 0.000100
	DECAY       0.000237
	MAXLEN      2000, (maxl0 = 2001)
	ESTEP       20
	..c'est ca!

	Finished!

If you open the file config/test_replicase.conf, you'll see that the values in the config file correspond to the entries following the text "AFTER loading the config..". This is handy for checking any config files that you build. 

A new feature with version 0.2.2 is the ability to specify the substitution matrix at the command line. This was done to make stringmol available on youShare more easily - the .mtx file can be specified directly, rather than being fetched out of the .conf file. This may be easier for some manual configurations too. Note that since this modification was designed for youShare, runs of stringmol will not be interactive where the extra argument is used to specify the .mtx file - this is only important for TType 0, where the user presses the space bar to step through a reaction - if three arguments are present, the 1 on 1 simulation just runds until limits are reached (see below)

TTYPE 0: 1 on 1
---------------

This option allows you to see how one molecule reacts with another. To test it, starting from the stringmol directory where this README is located, type the following:

$ mkdir test;cd test
$ ../src/stringmol 0 ../config/test_1on1.conf ../config/ALXII.mtx > out.txt

This will create a file called 'out.txt' that will show each step of a reaction between the two molecules listed in the file test_1on1.conf, located in the config directory. 

Should you want to step through the reaction, you'll need to enter the correct file path for the .mtx file in the .conf file that you use - see the line beginning with the keyword "SUBMAT" in 'test_1on1.conf'. Once this is done, the command is:

$ mkdir test;cd test
$ ../src/stringmol 0 ../config/test_1on1.conf

you can then "step through" a reaction between molecules by repeatedly pressing the space bar. 

We'd like to improve this facility. Please send feedback to sjh436@gmail.com

TTYPE 1: Alife XXII
-------------------

This is the grand-daddy stringmol program - running a "container" of molecular reactions together. We have included a config file with the same configuration as for the paper "ALxii_dif.pdf" that is included in the tarball. 

$ mkdir test;cd test
$ ../src/stringmol 1 ../../config/ALXII.conf ../../config/ALXII.mtx

TTYPE 2: Con Pop
----------------

This is an unpublished trial type, in which a set of containers of molecules are run. If the population of molecules in a container goes extinct, half the contents of another randomly-selected container are moved into the empty one. 

You can test this trial type using the same configuration for TTYPE 1 using:

$ mkdir test;cd test
$ ../src/stringmol 2 ../../config/ALXII.conf ../../config/ALXII.mtx

Currently only four containers are run, but we have done trials where there are 16 (C programmers can set this by changing the value of NCON on line 1565 of the file TestSM.cpp, and recompiling). 

This trial runs fairly reliably, but there may need to be changes in the output files to make sense of what is going on. At the moment, the output to screen after the config files have loaded looks like:

	0	0	1	2	3
	Count:	150	150	150	150
	Printing species list
	1000	0	1	2	3
	Count:	165	172	172	164
	2000	0	1	2	3
	Count:	198	222	214	181
	3000	0	1	2	3

A list of species that the simulations have created is printed every 10000 time steps - the line "Printing species list" indicates this. Apart from this, there are alternating lines. The first line (e.g. "2000	0	1	2	3" indicates the timestep, and the number of each container.
The second line (e.g. "Count:	198	222	214	181") indicates the total number of molecules in each container.







