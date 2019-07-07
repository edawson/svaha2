CXX:=g++
CXXFLAGS:= -O0 -g -fopenmp -ggdb -std=c++14
LD_INC_FLAGS:= -I./tinyFA -I./tinyFA/pliib/ -I./gfakluge/src -I./sparsepp/sparsepp -I./tinyVCF -I./tinyVCF/Hash-master/src
LD_LIB_FLAGS:=

EXEC:=svaha2

$(EXEC): main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp  Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

.PHONY: clean

clean:
	$(RM) $(EXEC)
