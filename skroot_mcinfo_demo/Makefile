# vim: filetype=Makefile #
# first things first: pull in all the SKOFL stuff
include    $(SKOFL_ROOT)/config.gmk

# programs used during compilation
##################################
CC          = g++
CXX         = g++
SHELL      := /bin/bash

## pull out source files. Match .cc, .cpp, .hh and .h
#####################################################
SRC_DIR     = src
INC_DIR     = include
SRC         = $(wildcard $(SRC_DIR)/*.cc)
SRC         += $(wildcard $(SRC_DIR)/*.cpp)
OBJ1         = $(SRC:%.cc=%.o)
OBJ         += $(OBJ1:%.cpp=%.o)
HDR         = $(wildcard $(INC_DIR)/*.hh)
HDR         += $(wildcard $(INC_DIR)/*.h)
HDR_NAMES   = $(patsubst $(INC_DIR)/%,%,$(HDR))

# includes, libraries and compilation flags
###########################################
# user include paths
LOCAL_INC   = -I$(INC_DIR)
# user library paths
LDFLAGS = -L/HOME/LinuxSystemFiles/SK/stllibs
# user libraries
LOCAL_LIBS = -lRootStl
# C++ compiler flags - XXX config.gmk sets this already, so APPEND ONLY XXX
CXXFLAGS    += -std=c++11 -fdiagnostics-color=always -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable
# -MMD -MP

# targets
#########
EXECUTABLE=show_neutrons
TESTER=test_show_neutrons
LIBRARY=$(SRC_DIR)/libShowNeutrons.so

all: $(EXECUTABLE) test_show_neutrons

# generic stuff
###############
OUTPUT_OPTION = -o $@
.PHONY: clean all

# THE RULES
###########
$(OBJ): $(SRC)
	@echo "my making $@"
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXECUTABLE) : $(OBJ)
	@echo "making $@"
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(OUTPUT_OPTION) $(LDLIBS) -Wl,--as-needed
#	$(RM) $(OBJ)

# a simple reader to check output files have been generated
# since this doesn't use any of the SKROOT classes we don't need to specify SKOFL/ATMPD via LD_RUN_PATH
$(TESTER): $(SRC_DIR)/test_show_neutrons.cc
	@echo "making $@"
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(OUTPUT_OPTION) $(LDLIBS) -Wl,--as-needed

# surprisingly this also works
$(LIBRARY):  $(OBJ)
	$(CXX) $(CXXFLAGS) -shared $(LDFLAGS) $^ $(OUTPUT_OPTION)

clean:
	@$(RM) $(OBJ) \
	       $(EXECUTABLE) \
	       $(LIBRARY) \
	       $(TESTER)

#-include $(SRC:%.cc=%.d)

## REFERENCE NOTES:
# we do not need to write our own rule for %.o : %.cpp etc, if we use make's implicit rules. 
# These will compile .o files (with g++) for the '.cpp' or '.cc' extension, but NOT for the '.cxx' extension
# We can pass flags to the implicit rules by using the standard flag sets: 
#   CPPFLAGS for pre-processor flags (e.g. 'include' paths)
#   CXXFLAGS for other c++ compiler flags (including -shared)
#   LDFLAGS for library paths
#   LDLIBS for libraries (-lMyLib)
#   OUTPUT_OPTION for defining target naming rules. (Perhaps this could somehow be used to put the objects in a build dir, but then the implicit rules do not get triggered for $(BUILD_DIR)/%.o files...)
# compiling and linking against ROOT dictionary generated with rootcint is necessary for classes that depend on ROOT, otherwise you'll be left with undefined references to Streamer, Class, IsA... etc
# rootcint calls '$CC', which unless not otherwised defined will resolve to 'gcc' NOT 'g++' and will not recognise the -cstd=c++1y flag (and some CXX flags) and will fall over.
# rootcint also requires the include directory to be specified WITHOUT A SPACE AFTER THE -I, which would normally be ok. To allow the same CPPFLAGS to be used for both, ensure CPPFLAGS uses no space!


#######################################################
# simple template with automatic dependency inclusion #
# SRC = $(wildcard *.cpp)                             #
# all: main                                           #
#                                                     #
# main: $(SRC:%.cpp=%.o)                              #
#	g++ -std=c++1y -MD -MP -o $@ $^               #
#                                                     #
# -include $(SRC:%.cpp=%.d)                           #
#######################################################
