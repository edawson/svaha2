CXX:=g++
CXXFLAGS:= -O0 -g -fopenmp -ggdb -std=c++14
LD_INC_FLAGS:=  -I./gfakluge/src -I./sparsepp/sparsepp -I./tinyVCF -I./tinyVCF/Hash-master/src -I./tinyFA -I./tinyFA/pliib/
LD_LIB_FLAGS:=

EXEC:=svaha2

fast: main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) -O3 -mtune=native -march=native -funroll-loops -std=c++14 -fopenmp -o $(EXEC) $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

debug: main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) $(CXXFLAGS) -DDEBUG=1 -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)


cloud: main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) -O2 -mtune=haswell -march=haswell -o svaha2 $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

$(EXEC): main.cpp tinyFA/tinyFA.hpp gfakluge/src/gfakluge.hpp tinyVCF/tinyVCF.hpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

.PHONY: clean fast test debug

test: $(EXEC)
	cd tests/ && prove -v test.sh

clean:
	$(RM) $(EXEC)
