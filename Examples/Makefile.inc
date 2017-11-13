# Parameters for Makefile

# !!!! MAKE SURE TO USE C++11 OR HIGHER !!!!

CC = g++
ARCH = ar -rUu 

# Path to StoichPack directory
STOICHPACK = ../../StoichPack
# Path to eigen directory
EIGEN = ../eigen/Eigen

# all projects will include this
DEFAULT_INCLUDE = -I$(STOICHPACK) -I$(EIGEN) -I..

# default compiler flags
# g++: -std=c++0x assures C++11
DEFAULT_STD_FLAGS = -ansi -pedantic -Wall -Werror -std=c++0x

# debug flags
DEFAULT_DEGUG_FLAGS = -g -O0 -fsanitize=address -fno-omit-frame-pointer

# flags for max. performance
DEFAULT_OPTIMIZE_FLAGS = -O3 -DNDEBUG

# debug or performance?
# USE THIS OPTION FOR DEBUGING
#DEFAULT_OPTIMIZE_DEBUG_FLAGS = $(DEFAULT_DEGUG_FLAGS)
# USE THIS OPTION FOR PERFORMANCE
DEFAULT_OPTIMIZE_DEBUG_FLAGS = $(DEFAULT_OPTIMIZE_FLAGS)

# use these flags everytime
DEFAULT_CC_FLAGS = $(DEFAULT_OPTIMIZE_DEBUG_FLAGS) $(DEFAULT_INCLUDE) $(DEFAULT_STD_FLAGS)

# link this in every link step
# actually this is not necessary if g++ is used
DEFAULT_LINK = -lstdc++

