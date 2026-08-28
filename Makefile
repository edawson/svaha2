CXX:=clang++
CXXFLAGS:= -O3 -march=native -ffast-math -std=c++17
LD_INC_FLAGS:=  -I./include
LD_LIB_FLAGS:=

EXEC:=svaha
VERSION:=0.2.0
DOCKER_IMAGE:=erictdawson/svaha2:v$(VERSION)
DOCKER_LATEST:=erictdawson/svaha2:latest

# Ensure all headers are included in the dependency list
HEADERS:=include/svaha.hpp

fast: src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $(EXEC) $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

debug: src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -g -DDEBUG=1 -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

$(EXEC): src/main.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

.PHONY: clean fast test debug docker-build docker-push

docker-build:
	docker build -t $(DOCKER_IMAGE) -t $(DOCKER_LATEST) .

docker-push: docker-build
	docker push $(DOCKER_IMAGE)
	docker push $(DOCKER_LATEST)

test: tests/test_svaha.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o test_svaha $< $(LD_INC_FLAGS)
	./test_svaha

clean:
	$(RM) $(EXEC) debug test_svaha
