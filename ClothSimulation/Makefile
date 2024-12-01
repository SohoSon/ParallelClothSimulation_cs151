#Makefile for cloth simulation
ASSIGN    = cloth_simulation
BREWPATH  = $(shell brew --prefix)
CXX       = /opt/homebrew/opt/llvm/bin/clang++
CXXFLAGS  = -std=c++11 -D_CRT_SECURE_NO_WARNINGS -DGL_SILENCE_DEPRECATION -Wno-macro-redefined -O2 -I./Headers -I$(BREWPATH)/include -fopenmp -arch arm64
LDFLAGS   = -L$(BREWPATH)/lib -lglfw -framework OpenGL -framework Cocoa -framework IOKit -arch arm64

all: $(ASSIGN)

$(ASSIGN): main.o glad.o
	$(CXX) $^ $(LDFLAGS) -o $@

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

glad.o: glad.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(ASSIGN) *.o *~ *.dSYM

.PHONY: clean all