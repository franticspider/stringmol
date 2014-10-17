
# the name you want to give the executable
TARGET = conpop

# the code files you want to compile
# problem with this approach - you don't know if the objs have come from .c .cpp .cc files...s
#CFILES = RobotUpdater.o IrHandler.o WanderModule.o AvoidModule.o EpuckSensors.o Horst.o

CFILES = *.c
CPPFILES = *.cpp
CCFILES = *.cc

# the directory where all your code is stored
MAIN_DIR =.

# the folder where header files are stored
HEADERS = $(MAIN_DIR)

#############################################################
################## DON'T CHANGE CODE BELOW ##################
#############################################################


# which compiler to use 
C4C = gcc
CC = g++


# sjh utility functions - no dependencies:
 sjhutils: 		microbial_ga.o mt19937-2.o mathutil.o memoryutil.o params.o randutil.o hsort.o

hsort.o: 	hsort.c hsort.h 
	$(C4C) -O3 -Wall -c hsort.c

mathutil.o: 	mathutil.c mathutil.h 
	$(C4C) -O3 -Wall -c mathutil.c

memoryutil.o: 	memoryutil.c memoryutil.h 
	$(C4C) -O3 -Wall -c memoryutil.c

params.o:	params.c params.h 
	$(C4C) -O3 -Wall -c params.c

randutil.o: 	randutil.c randutil.h mt19937-2.h
	$(C4C) -O3 -Wall -c randutil.c

mt19937-2.o: 	mt19937-2.c mt19937-2.h
	$(C4C) -O3 -Wall -c mt19937-2.c

microbial_ga.o: microbial_ga.c microbial_ga.h
	$(C4C) -O3 -Wall -c microbial_ga.c

#########################
### Stringmol library ###
#########################
smlib: 			stringmanip.o	alignment.o	instructions.o\
				SMspp.o	rules.o agents_base.o stringPM.o\
				signalSM.o

#stringmol cpp files:


signalSM.o: 	signalSM.cpp signalSM.h
	$(CC) -O3 -Wall -c signalSM.cpp

stringmanip.o: 	stringmanip.cpp stringmanip.h
	$(CC) -O3 -Wall -c stringmanip.cpp

alignment.o: 	alignment.cpp alignment.h randutil.h
	$(CC) -O3 -Wall -c alignment.cpp

instructions.o: instructions.cpp instructions.h memoryutil.h randutil.h alignment.h stringmanip.h
	$(CC) -O3 -Wall -c instructions.cpp

SMspp.o:	SMspp.cpp SMspp.h
	$(CC) -O3 -Wall -c SMspp.cpp

rules.o:	rules.cpp rules.h memoryutil.h
	$(CC) -O3 -Wall -c rules.cpp

agents_base.o: agents_base.cpp agents_base.h memoryutil.h	randutil.h	params.h	rules.h
	$(CC) -O3 -Wall -c agents_base.cpp

stringPM.o:	stringPM.cpp stringPM.h \
		memoryutil.h	randutil.h	\
            	alignment.h		stringmanip.h	instructions.h	\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -O3 -Wall -c stringPM.cpp



#  g++ -o prog link-options file1.o file2.o ... other-libraries

#########################
### stringmol executable ###
#########################


setupSM.o:	setupSM.cpp \
		randutil.h	alignment.h \
		rules.h agents_base.h SMspp.h stringPM.h signalSM.h
	$(CC) -O3 -Wall -c setupSM.cpp 	
			


TestSM.o:  	TestSM.cpp \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h
	$(CC) -O3 -Wall -c TestSM.cpp 	


stringmol: TestSM.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o
	$(CC) -O3 -Wall -o stringmol -lm TestSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o\
			stringmanip.o	alignment.o	instructions.o\
			stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o


all: sjhutils smlib setupSM.o stringmol



#all: smlib sjhutils

# removes all the object files. Not automatically done by a call to make,
# need to type "make clean" (without the quotes) into the terminal.	
# or "make all clean"
clean:
	rm -f *.o
	rm stringmol

