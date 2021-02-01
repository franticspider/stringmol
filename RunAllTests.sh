#!/bin/bash


#g++ main.cpp -o travis_gcc_cpp98
#./travis_gcc_cpp98
#cppcheck --quiet --error-exitcode=1 main.cpp



echo "Running Tests.  Please Wait."
#g++ -Wall Shapes-Catch-Testing-Example/Source/Shapes-Catch-Testing-Example.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o main
#g++ -Wall Shapes-Catch-Testing-Example/Test/*.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o test
cd tests
echo "compiling..."
g++ -std=gnu++11 -Wall -o test  *.cpp ../release/mt19937-2.o ../release/randutil.o 
cd ..
echo ""
echo "  now testing.."
./tests/test
