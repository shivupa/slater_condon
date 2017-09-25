#!/bin/bash

icpc -std=c++11 -O3 -xAVX -ipo test.cpp ../nsubst1.x.cpp ../read_input.cpp -o testAVX
icpc -std=c++11 -O3 -xSSE2 -ipo test.cpp ../nsubst1.x.cpp ../read_input.cpp -o testSSE2
