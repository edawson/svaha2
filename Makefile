CXX:=g++
CXXFLAGS:= -O0 -g -fopenmp -ggdb -std=c++14
LD_INC_FLAGS:=  -I./gfakluge/src -I./sparsepp/sparsepp -I./tinyVCF -I./tinyVCF/Hash-master/src -I./tinyFA -I./tinyFA/pliib/
LD_LIB_FLAGS:=

EXEC:=svaha2

$(EXEC): main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

fast: main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) -O3 -mtune=native -march=native -funroll-loops -std=c++14 -fopenmp -o $(EXEC) $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)


.PHONY: clean fast test

test: $(EXEC)
	cd tests/ && prove test.sh

clean:
	$(RM) $(EXEC)
