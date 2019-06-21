CXX=g++
CXXFLAGS=-Wall -Wextra -pedantic -std=c++17 -O0 -g
LDFLAGS=$(CXXFLAGS)
OBJ=$(SRC:.cc=.o)

all:	test_nw	test_alt_nw

test_nw:	test_nw.cc nw.cc nw_sim_mat.cc
	$(CXX) $(LDFLAGS) -o $@ $^

test_alt_nw:	test_alt_nw.cc alt_nw.cc
	$(CXX) $(LDFLAGS) -o $@ $^
