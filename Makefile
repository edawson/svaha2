CXX:=clang++
CXXFLAGS:= -O3 -march=native -ffast-math -std=c++17
LD_INC_FLAGS:=  -I./include
LD_LIB_FLAGS:=

EXEC:=svaha

# Ensure all headers are included in the dependency list
HEADERS:=include/svaha.hpp

fast: src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $(EXEC) $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

debug: src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -g -DDEBUG=1 -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

$(EXEC): src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

.PHONY: clean fast test debug

test: tests/test_svaha.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o test_svaha $< $(LD_INC_FLAGS)
	./test_svaha

clean:
	$(RM) $(EXEC) debug test_svaha
