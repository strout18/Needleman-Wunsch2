CXX=g++
CXXFLAGS=-Wall -Wextra -pedantic -std=c++17 -O0 -g
LDFLAGS=$(CXXFLAGS)
OBJ=$(SRC:.cc=.o)

all:	test_nw

test_nw:	test_nw.cc nw.cc
	$(CXX) $(LDFLAGS) -o $@ $^
