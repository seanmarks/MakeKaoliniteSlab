# NOTES:
# -g	Extra symbolic debugging info for use with the gdb debugger
# -o	Specify output file name
# -c	Compile into object file
# -Wall	Print all warning messages
#
# $@	Name of the target
# $<	First dependency of the target
# $^	All dependencies of the target
#
# Linking occurs when the object files (.o) are put together into an executable (.exe)

# Output program
PROJECT=MakeKaoliniteSlab

# Compiler
CC= g++

# Compiler flags 
CXXFLAGS= -Wall -g -std=c++11 -DCPLUSPLUS -O3

# Linking
LFLAGS= -g


#######


# Make all
.PHONY: all
all : $(PROJECT)

# Linking
$(PROJECT) : 
	$(CC) $(LFLAGS) -o bin/$@ $(CXXFLAGS) src/main.cpp

clean:
	rm -rf bin/*

.PHONY: test
test:
	@{ \
	cd test/grid_2_2_2 ;\
	./run.sh ;\
	./compare.sh ;\
	}
