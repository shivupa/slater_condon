#!/bin/bash

#icpc -std=c++11 -O3 -xAVX -ipo test.cpp ../nsubst1.x.cpp ../read_input.cpp -o testAVX
#icpc -std=c++11 -O3 -xSSE2 -ipo test.cpp ../nsubst1.x.cpp ../read_input.cpp -o testSSE2
g++ -O3 -mavx test.cpp ../read_input.cpp ../nsubst1.x.cpp -o testGNUAVX
g++ -O3 -msse2 test.cpp ../read_input.cpp ../nsubst1.x.cpp -o testGNUSSE2
g++ -O3 -msse3 test.cpp ../read_input.cpp ../nsubst1.x.cpp -o testGNUSSE3
g++ -O3 -msse4.1 test.cpp ../read_input.cpp ../nsubst1.x.cpp -o testGNUSSE4_1
g++ -O3 -msse4.2 test.cpp ../read_input.cpp ../nsubst1.x.cpp -o testGNUSSE4_2
