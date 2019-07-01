CXX:=g++
CXXFLAGS:= -O0 -g -fopenmp -ggdb
LD_INC_FLAGS:= -I./tinyFA -I./tinyFA/pliib/ -I./gfakluge/src -I./sparsepp/sparsepp -I./tinyVCF -I./tinyVCF/digestpp
LD_LIB_FLAGS:=

svaha2: main.cpp tinyFA/tinyfa.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)
