## GNU Makefile


## Compilers
# CXX : Serial compiler.
export CXX 						= g++

## Build flags
# CPPFLAGS : Compile flags for all components. -Wall (all warnings), and -O3
# (all optimisations) by default.
# LFLAGS : Linker flags for all components.
# LIBS : Libraries to link to.
export CPPFLAGS				= -Wall -O3
export LFLAGS					=
export LIBS						=

## Paths
SRC									= BasicCatheter/BasicCatheter

all: Catheter

.PHONY: Catheter
Catheter:
	$(MAKE) --directory=$(SRC) -e

.PHONY: clean
clean:
	$(MAKE) --directory=$(SRC) clean



