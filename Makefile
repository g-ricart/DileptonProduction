all:
	g++ -std=c++11 -g -W -Wall -Wextra -pedantic -Wno-stringop-truncation src/Main.cpp -o Calculate.exe -O3 -lgsl -lgslcblas
