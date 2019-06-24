CXX:=g++
CXXFLAGS:= -O1 -g -fopenmp -mtune=native -march=native
LD_INC_FLAGS:= -I./tinyFA -I./tinyFA/pliib/ -I./gfakluge/src -I./sparsepp/sparsepp -I./tinyVCF
LD_LIB_FLAGS:=

svaha2: main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)
