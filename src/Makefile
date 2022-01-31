#############################################################
# This makefile will: 
# make ../release/stringmol
# make ../debug/stringmol
#
# C4C is the compiler for C files
# CC is the compiler for C++ files
#
# https://makefiletutorial.com/ is really helpful, *but*
# it seems current best practice is to use the "-M" command
# with gcc to automatically detect which header file to 
# include.. 
# see also 
# https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#
# Conclusion - we'll have to write this out longhand and 
# maintain it carefully - but at least we'll have a makefile
# that is legible to a range of skillsets
#############################################################


#############################################################
################## VARIABLES               ##################
#############################################################


# which compiler to use 
C4C = gcc -O3 -Wall -Wunused
CC = g++ -O3 -Wall -Wunused -std=c++11

#TODO: Handle profiling better...
#debug: C4C = gcc -O0 -Wall -Wunused -DDEBUG -g 
#debug: CC = g++ -O0 -Wall -Wunused -std=c++11 -DDEBUG -g


#If you want profiling, comment out the above and replace with: 
debug: C4C = gcc -O0 -Wall -Wunused -DDEBUG -g -pg
debug: CC = g++ -O0 -Wall -Wunused -std=c++11 -DDEBUG -g -pg


#build directories:
DEBUG_DIR =  ../debug
RELEASE_DIR =  ../release

#TODO: There's a lot of repetition of $(OUTPUT_DIR) atm...!!!
#OUTPUT_DIR = $(RELEASE_DIR)
#debug: OUTPUT_DIR = $(DEBUG_DIR)



#############################################################
################## TOP LEVEL TARGETS       ##################
#############################################################

# Notes on the all/debug thing. There's a lot of duplication here, but this 
# helps to keep things clear. We have to run a `cleano' before all/debug because
# otherwise there'll be .o files that may have been built with a different
# config

##Make everything from scratch:
#all: $(OUTPUT_DIR)/stringmol

##TODO: Figure out how to set the debug flag to call SDL-devel (or SDL-dev)
##debug: stringmol web_api
#debug: $(OUTPUT_DIR)/stringmol

all: debug release

release: $(RELEASE_DIR)/stringmol
debug: $(DEBUG_DIR)/stringmol

web_api: bind step

# removes all the object files. Not automatically done by a call to make,
# need to type "make clean" (without the quotes) into the terminal.	
# or "make all clean"
clean:
	rm -f $(RELEASE_DIR)/*.o
	rm -f stringmol
	rm -f stringmol_dbg
	rm -f webapi_*.cgi

#just remove the object files:
cleano:
	rm -f *.o
	rm -f ../debug/*.o
	rm -f ../release/*.o

	

#############################################################
################## C OBJECTS	           ##################
#############################################################

RELEASE_COBJECTS = $(RELEASE_DIR)/hsort.o \
		$(RELEASE_DIR)/mathutil.o \
		$(RELEASE_DIR)/memoryutil.o \
		$(RELEASE_DIR)/microbial_ga.o \
		$(RELEASE_DIR)/params.o \
		$(RELEASE_DIR)/randutil.o \
		$(RELEASE_DIR)/mt19937-2.o \
		
DEBUG_COBJECTS = $(DEBUG_DIR)/hsort.o \
		$(DEBUG_DIR)/mathutil.o \
		$(DEBUG_DIR)/memoryutil.o \
		$(DEBUG_DIR)/microbial_ga.o \
		$(DEBUG_DIR)/params.o \
		$(DEBUG_DIR)/randutil.o \
		$(DEBUG_DIR)/mt19937-2.o \


$(RELEASE_DIR)/hsort.o: 		hsort.c hsort.h 
	$(C4C) -c hsort.c -o $(RELEASE_DIR)/hsort.o
	
$(DEBUG_DIR)/hsort.o: 		hsort.c hsort.h 
	$(C4C) -c hsort.c -o $(DEBUG_DIR)/hsort.o


$(RELEASE_DIR)/mathutil.o: 	mathutil.c mathutil.h 
	$(C4C) -c mathutil.c -o $(RELEASE_DIR)/mathutil.o

$(DEBUG_DIR)/mathutil.o: 	mathutil.c mathutil.h 
	$(C4C) -c mathutil.c -o $(DEBUG_DIR)/mathutil.o


$(RELEASE_DIR)/memoryutil.o: 	memoryutil.c memoryutil.h 
	$(C4C) -c memoryutil.c -o $(RELEASE_DIR)/memoryutil.o

$(DEBUG_DIR)/memoryutil.o: 	memoryutil.c memoryutil.h 
	$(C4C) -c memoryutil.c -o $(DEBUG_DIR)/memoryutil.o
	

$(RELEASE_DIR)/microbial_ga.o: microbial_ga.c microbial_ga.h
	$(C4C) -c microbial_ga.c -o $(RELEASE_DIR)/microbial_ga.o
	
$(DEBUG_DIR)/microbial_ga.o: microbial_ga.c microbial_ga.h
	$(C4C) -c microbial_ga.c -o $(DEBUG_DIR)/microbial_ga.o


$(RELEASE_DIR)/params.o:		params.c params.h 
	$(C4C) -c params.c -o $(RELEASE_DIR)/params.o
	
$(DEBUG_DIR)/params.o:		params.c params.h 
	$(C4C) -c params.c -o $(DEBUG_DIR)/params.o


$(RELEASE_DIR)/mt19937-2.o: 	mt19937-2.c mt19937-2.h
	$(C4C) -c mt19937-2.c -o $(RELEASE_DIR)/mt19937-2.o

$(DEBUG_DIR)/mt19937-2.o: 	mt19937-2.c mt19937-2.h
	$(C4C) -c mt19937-2.c -o $(DEBUG_DIR)/mt19937-2.o
	

$(RELEASE_DIR)/randutil.o: 	randutil.c randutil.h mt19937-2.h
	$(C4C) -c randutil.c -o $(RELEASE_DIR)/randutil.o

$(DEBUG_DIR)/randutil.o: 	randutil.c randutil.h mt19937-2.h
	$(C4C) -c randutil.c -o $(DEBUG_DIR)/randutil.o



#############################################################
################## C++ OBJECTS	           ##################
#############################################################


RELEASE_CPPOBJECTS = $(RELEASE_DIR)/lodepng.o \
		$(RELEASE_DIR)/agents_base.o \
		$(RELEASE_DIR)/alignment.o \
		$(RELEASE_DIR)/instructions.o \
		$(RELEASE_DIR)/rules.o \
		$(RELEASE_DIR)/setupSM.o \
		$(RELEASE_DIR)/SMspp.o \
		$(RELEASE_DIR)/stringmanip.o \
		$(RELEASE_DIR)/stringPM.o \
		$(RELEASE_DIR)/tests.o \
		
DEBUG_CPPOBJECTS = $(DEBUG_DIR)/lodepng.o \
		$(DEBUG_DIR)/agents_base.o \
		$(DEBUG_DIR)/alignment.o \
		$(DEBUG_DIR)/instructions.o \
		$(DEBUG_DIR)/rules.o \
		$(DEBUG_DIR)/setupSM.o \
		$(DEBUG_DIR)/SMspp.o \
		$(DEBUG_DIR)/stringmanip.o \
		$(DEBUG_DIR)/stringPM.o \
		$(DEBUG_DIR)/tests.o \
		

###################
### PNG library ###
###################

$(RELEASE_DIR)/lodepng.o: lodepng.cpp lodepng.h
	$(CC) -c lodepng.cpp -o $(RELEASE_DIR)/lodepng.o

$(DEBUG_DIR)/lodepng.o: lodepng.cpp lodepng.h
	$(CC) -c lodepng.cpp -o $(DEBUG_DIR)/lodepng.o


#########################
### Stringmol library ###
#########################

#stringmol cpp files:

$(RELEASE_DIR)/agents_base.o: agents_base.cpp agents_base.h \
			memoryutil.h	randutil.h	params.h	rules.h
	$(CC) -c agents_base.cpp -o $(RELEASE_DIR)/agents_base.o

$(DEBUG_DIR)/agents_base.o: agents_base.cpp agents_base.h \
			memoryutil.h	randutil.h	params.h	rules.h
	$(CC) -c agents_base.cpp -o $(DEBUG_DIR)/agents_base.o



$(RELEASE_DIR)/alignment.o: 	alignment.cpp alignment.h randutil.h
	$(CC) -c alignment.cpp -o $(RELEASE_DIR)/alignment.o

$(DEBUG_DIR)/alignment.o: 	alignment.cpp alignment.h randutil.h
	$(CC) -c alignment.cpp -o $(DEBUG_DIR)/alignment.o
	
	

$(RELEASE_DIR)/instructions.o: instructions.cpp instructions.h \
			memoryutil.h randutil.h alignment.h stringmanip.h
	$(CC) -c instructions.cpp -o $(RELEASE_DIR)/instructions.o	

$(DEBUG_DIR)/instructions.o: instructions.cpp instructions.h \
			memoryutil.h randutil.h alignment.h stringmanip.h
	$(CC) -c instructions.cpp -o $(DEBUG_DIR)/instructions.o
	
	

$(RELEASE_DIR)/rules.o:	rules.cpp rules.h memoryutil.h
	$(CC) -c rules.cpp -o $(RELEASE_DIR)/rules.o

$(DEBUG_DIR)/rules.o:	rules.cpp rules.h memoryutil.h
	$(CC) -c rules.cpp -o $(DEBUG_DIR)/rules.o
	
	

$(RELEASE_DIR)/setupSM.o:	setupSM.cpp setupSM.h \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			lodepng.h mt19937-2.h
	$(CC) -c setupSM.cpp -o $(RELEASE_DIR)/setupSM.o	

$(DEBUG_DIR)/setupSM.o:	setupSM.cpp setupSM.h \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			lodepng.h mt19937-2.h
	$(CC) -c setupSM.cpp -o $(DEBUG_DIR)/setupSM.o
	
	

$(RELEASE_DIR)/SMspp.o:	SMspp.cpp SMspp.h
	$(CC) -c SMspp.cpp -o $(RELEASE_DIR)/SMspp.o

$(DEBUG_DIR)/SMspp.o:	SMspp.cpp SMspp.h
	$(CC) -c SMspp.cpp -o $(DEBUG_DIR)/SMspp.o
	
	

$(RELEASE_DIR)/stringmanip.o: 	stringmanip.cpp stringmanip.h
	$(CC) -c stringmanip.cpp -o $(RELEASE_DIR)/stringmanip.o

$(DEBUG_DIR)/stringmanip.o: 	stringmanip.cpp stringmanip.h
	$(CC) -c stringmanip.cpp -o $(DEBUG_DIR)/stringmanip.o



$(RELEASE_DIR)/stringPM.o:	stringPM.cpp stringPM.h \
		memoryutil.h	randutil.h	\
        alignment.h		stringmanip.h	instructions.h	\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c stringPM.cpp -o $(RELEASE_DIR)/stringPM.o

$(DEBUG_DIR)/stringPM.o:	stringPM.cpp stringPM.h \
		memoryutil.h	randutil.h	\
        alignment.h		stringmanip.h	instructions.h	\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c stringPM.cpp -o $(DEBUG_DIR)/stringPM.o
	
	

$(RELEASE_DIR)/tests.o:	tests.cpp tests.h
	$(CC) -c tests.cpp -o $(RELEASE_DIR)/tests.o

$(DEBUG_DIR)/tests.o:	tests.cpp tests.h
	$(CC) -c tests.cpp -o $(DEBUG_DIR)/tests.o


	
	
	
	
#############################################################
################## LIBS and EXECUTABLES    ##################
#############################################################


###########################################################
# stringmol
$(RELEASE_DIR)/stringmol: 	$(RELEASE_DIR)/stringmol.o \
			$(RELEASE_COBJECTS) $(RELEASE_CPPOBJECTS)
	$(CC) -o $(RELEASE_DIR)/stringmol -lm $(RELEASE_DIR)/stringmol.o \
			$(RELEASE_COBJECTS) $(RELEASE_CPPOBJECTS)\
			
			
$(DEBUG_DIR)/stringmol: 	$(DEBUG_DIR)/stringmol.o \
			$(DEBUG_COBJECTS) $(DEBUG_CPPOBJECTS)
	$(CC) -o $(DEBUG_DIR)/stringmol -lm $(DEBUG_DIR)/stringmol.o \
			$(DEBUG_COBJECTS) $(DEBUG_CPPOBJECTS)\
			
				

###########################################################
# stringmol web api			
			
bind: 	webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(RELEASE_DIR)/$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o lodepng.o
	$(CC) -o webapi_bind.cgi -lm webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(RELEASE_DIR)/$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o lodepng.o
			
			
step: 	webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(RELEASE_DIR)/$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o  lodepng.o
	$(CC) -o webapi_step.cgi -lm webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(RELEASE_DIR)/$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o	 lodepng.o		
		


#########################
### stringmol objects ###
#########################


$(RELEASE_DIR)/stringmol.o:  	stringmol.cpp \
			$(RELEASE_COBJECTS) $(RELEASE_CPPOBJECTS)
	$(CC) -c stringmol.cpp -o  	$(RELEASE_DIR)/stringmol.o
	
$(DEBUG_DIR)/stringmol.o:  	stringmol.cpp \
			$(DEBUG_COBJECTS) $(DEBUG_CPPOBJECTS)
	$(CC) -c stringmol.cpp -o  	$(DEBUG_DIR)/stringmol.o
			

#########################
### web_api objects ###
#########################


webapi_util.o:	webapi_util.cpp webapi_util.h stringPM.h \
		randutil.h	params.h hsort.h\
        alignment.h	setupSM.h\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c webapi_util.cpp
	
	
	
webapi_bind.o:  webapi_bind.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_bind.cpp 

webapi_step.o:  webapi_step.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_step.cpp 		

