all:
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main.cpp -o Calculate.exe -O3 -lgsl -lgslcblas