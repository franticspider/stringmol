
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
#C4C = gcc -O0 -Wall
#CC = g++  -O0 -Wall
C4C = gcc -O0 -Wall -Wunused
CC = g++ -O0 -Wall -Wunused

#DBGC4C = $(C4C) -DDEBUG -g
#DBGCC = $(CC) -DDEBUG -g

debug: CXX += -DDEBUG -g
debug: CC += -DDEBUG -g



#############################################################
################## TOP LEVEL               ##################
#############################################################

# Notes on the all/debug thing. There's a lot of duplication here, but this helps to keep things clear. 
# 	We have to run a `cleano' before all/debug because otherwise there'll be .o files that may have been built with a different config

#Make everything from scratch:
all: clean release debug

release: cleano stringmol web_api

debug: cleano stringmol_dbg web_api_dbg

web_api: bind step

web_api_dbg: bind_dbg step_dbg

# removes all the object files. Not automatically done by a call to make,
# need to type "make clean" (without the quotes) into the terminal.	
# or "make all clean"
clean:
	rm -f *.o
	rm -f stringmol
	rm -f stringmol_dbg
	rm -f webapi_*.cgi
cleano:
	rm -f *.o

	
#############################################################
################## LIBS and EXES           ##################
#############################################################


# sjh utility functions - no dependencies:
sjhutils: 	microbial_ga.o mt19937-2.o mathutil.o \
			memoryutil.o params.o randutil.o hsort.o

# string manipulation library
smlib: 		stringmanip.o	alignment.o	instructions.o \
			SMspp.o	rules.o agents_base.o stringPM.o \
			signalSM.o


stringmol: 	stringmol.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o stringmol -lm stringmol.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
			

stringmol_dbg: 	stringmol.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o stringmol_dbg -lm stringmol.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o			
			
bind: 	webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o webapi_bind.cgi -lm webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o			
		

bind_dbg: 	webapi_bind_dbg.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o webapi_bind_dbg.cgi -lm webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o			
			
step: 	webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o webapi_step.cgi -lm webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o			
		

step_dbg: 	webapi_step_dbg.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o	mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o
	$(CC) -o webapi_step_dbg.cgi -lm webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o signalSM.o setupSM.o \
			hsort.o mathutil.o memoryutil.o params.o randutil.o mt19937-2.o microbial_ga.o \
			stringmanip.o	alignment.o	instructions.o tests.o			
		


	
#############################################################
################## C OBJECTS	           ##################
#############################################################


hsort.o: 	hsort.c hsort.h 
	$(C4C) -c hsort.c

mathutil.o: 	mathutil.c mathutil.h 
	$(C4C) -c mathutil.c

memoryutil.o: 	memoryutil.c memoryutil.h 
	$(C4C) -c memoryutil.c

params.o:	params.c params.h 
	$(C4C) -c params.c

randutil.o: 	randutil.c randutil.h mt19937-2.h
	$(C4C) -c randutil.c

mt19937-2.o: 	mt19937-2.c mt19937-2.h
	$(C4C) -c mt19937-2.c

microbial_ga.o: microbial_ga.c microbial_ga.h
	$(C4C) -c microbial_ga.c

#########################
### Stringmol library ###
#########################

#stringmol cpp files:

tests.o:	tests.cpp tests.h
	$(CC) -c tests.cpp

signalSM.o: 	signalSM.cpp signalSM.h
	$(CC) -c signalSM.cpp

stringmanip.o: 	stringmanip.cpp stringmanip.h
	$(CC) -c stringmanip.cpp

alignment.o: 	alignment.cpp alignment.h randutil.h
	$(CC) -c alignment.cpp

instructions.o: instructions.cpp instructions.h memoryutil.h randutil.h alignment.h stringmanip.h
	$(CC) -c instructions.cpp

SMspp.o:	SMspp.cpp SMspp.h
	$(CC) -c SMspp.cpp

rules.o:	rules.cpp rules.h memoryutil.h
	$(CC) -c rules.cpp

agents_base.o: agents_base.cpp agents_base.h memoryutil.h	randutil.h	params.h	rules.h
	$(CC) -c agents_base.cpp

stringPM.o:	stringPM.cpp stringPM.h \
		memoryutil.h	randutil.h	\
        alignment.h		stringmanip.h	instructions.h	\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c stringPM.cpp


webapi_util.o:	webapi_util.cpp webapi_util.h stringPM.h \
		randutil.h	params.h hsort.h\
        alignment.h	signalSM.h	setupSM.h\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c webapi_util.cpp


#  g++ -o prog link-options file1.o file2.o ... other-libraries

#########################
### stringmol executable ###
#########################


setupSM.o:	setupSM.cpp \
		randutil.h	alignment.h \
		rules.h agents_base.h SMspp.h stringPM.h signalSM.h
	$(CC) -c setupSM.cpp 	
			


stringmol.o:  	stringmol.cpp \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h \
			tests.h
	$(CC) -c stringmol.cpp 	

#########################
### web_api executables ###
#########################

webapi_bind.o:  webapi_bind.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_bind.cpp 

webapi_step.o:  webapi_step.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_step.cpp 		

#################################
### web_api debug executables ###
#################################


#WE NEED A "*_dbg.o" version for anything we want to use '#ifdef DEBUG' in...
webapi_bind_dbg.o:  webapi_bind.cpp \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h \
			tests.h
	$(CC) -c webapi_bind.cpp 		


#WE NEED A "*_dbg.o" version for anything we want to use '#ifdef DEBUG' in...
webapi_step_dbg.o:  webapi_step.cpp \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h signalSM.h \
			tests.h
	$(CC) -c webapi_step.cpp 



	
