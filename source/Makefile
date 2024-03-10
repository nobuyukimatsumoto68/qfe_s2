CXX=g++
CXXFLAGS=-std=c++14 -O3 -Wall
CXXFLAGS+=-I/share/pkg.7/eigen/3.4.0/install/include/eigen3/

GRP_DIR=$(shell pwd)/../grp
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""

SRC=$(wildcard *.cc)
BIN_DIR=$(shell pwd)/../bin
BIN=$(basename $(SRC))
INC_DIR=$(shell pwd)/../include
INC=$(wildcard $(INC_DIR)/*.h)

.PHONY: clean

all: $(BIN)

%: %.cc $(INC)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) $< -o $(BIN_DIR)/$@

clean:
	$(RM) -r $(addprefix $(BIN_DIR)/, $(BIN))
