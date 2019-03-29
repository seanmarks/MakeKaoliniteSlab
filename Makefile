
# Output binary's name
PROJECT=MakeKaoliniteSlab

# Compilers
CC=gcc
CXX=g++


### !!! Shouldn't have to change anything below this line !!! ###


#############
### Flags ###
#############

SHELL := /bin/bash

# Compiler flags
# - Misc: -march=native -Wno-comment -Wno-sign-compare -DPLUMED_MODE
CXXFLAGS += -g -std=c++11 -DCPLUSPLUS -O3 -Wall

# Linking
LDFLAGS += -g
LIBS    += -lm


#############################
### Directories and Files ###
#############################

START_DIR := $(PWD)

SRC_DIR         := src
BUILD_DIR       := build
BUILD_BIN_DIR   := $(BUILD_DIR)/bin
INSTALL_BIN_DIR := $(INSTALL_DIR)/bin

# Get source and target objects
SRC_FILES = $(shell find $(SRC_DIR) -name '*.cpp')
SRC_DIRS  = $(shell find $(SRC_DIR) -type d | sed 's/$(SRC_DIR)/./g' )
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


#############
#### Make ###
#############

all : buildrepo $(PROJECT)

# Link project
$(PROJECT) : $(OBJECTS)
	$(CXX) -o $(BUILD_BIN_DIR)/$@ $(LDFLAGS) $(LIBS) $(OBJECTS)

# Compile
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $<

.PHONY : clean	
clean :
	rm -rf $(BUILD_DIR)/*

.PHONY : clean_install
clean_install :
	rm -rf $(INSTALL_DIR)/*

# Make directories 
buildrepo :
	@{ \
	for dir in $(BUILD_DIR) $(BUILD_BIN_DIR) $(INSTALL_DIR)/bin ; do \
		if [[ ! -d $$dir ]]; then \
			echo "Making directory $$dir ..." ;\
			mkdir -p $$dir ;\
		fi ;\
	done ;\
	}

# Run regression test(s)
.PHONY: test
test:
	@{ \
	cd test/grid_2_2_2 ;\
	./run.sh ;\
	./compare.sh ;\
	}

