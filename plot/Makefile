########################################################################
########################################################################
#                                                                      #
# This is an automatic 'Makefile'. One generally does not have to edit #
# this file. Just put the program sources (source files with the       #
# program body) into the 'src' directory, put the library sources      #
# (source files to be linked to the programs, or the compiled          #
# libraries or shared libraries themselves) into the 'lib'             #
# directory, put the header sources (source files to be included in)   #
# into the 'inc' directory. You can directly refer to include files in #
# 'inc' directory, as they are added to the include search path.       #
# Then, by saying 'make' (or 'make bin/targetfile', where              #
# 'targetfile' is substituted by an appropriate name) source is        #
# compiled.                                                            #
#                                                                      #
# File naming conventions:                                             #
# - header files: '*.h' for C, '*.h', '*.H' for C++,                   #
#   '*.inc' for FORTRAN,                                               #
# - source files: '*.c' for C, '*.cc', '*.cpp', '*.cxx', '*.c++',      #
#   '*.C' for C++, '*.f', '*.F' for FORTRAN.                           #
# - library files: 'lib*.a' for static libraries, 'lib*.so' for shared #
#   libraries.                                                         #
#                                                                      #
# The compilers are specified by the 'C', 'CPP' and 'F' variables for  #
# C, C++ and FORTRAN. The linker is specified by 'L' variable.         #
#                                                                      #
# If you want to pass some global switches while compilation and       #
# linking, edit the variable 'ARGS'. If you want to change             #
# optimization settings, edit the variable 'OPT'.                      #
#                                                                      #
# If you want to pass arguments to the linker, then edit 'LARGS'       #
# variable.                                                            #
# If you want to pass arguments to the C, C++ or FORTRAN compilers,    #
# respectively, then edit 'CARGS', 'CPPARGS', 'FARGS'.                 #
#                                                                      #
# To compile in debug mode, say 'make debug' instead of 'make'.        #
#                                                                      #
########################################################################
########################################################################


####################### User definitions section: ######################


### C compiler:
ifeq "$(strip $(USEDC))" ""
 USEDC		= gcc
endif
### C++ compiler:
ifeq "$(strip $(USEDCPP))" ""
 USEDCPP	= g++
endif
### FORTRAN compiler:
ifeq "$(strip $(USEDF))" ""
 USEDF		= g77
endif
### Linker:
USEDL		= $(USEDCPP)

### Global arguments:
ARGS		= 
### Global optimization arguments:
OPT		= -O3
#OPT		= -O0
#OPT		= -O2
#OPT		= -O3 -fno-fast-math -ffloat-store

### Arguments to be passed when C compiling. E.g. -Idir switches:
CARGS		= -Wall
### Arguments to be passed when C++ compiling. E.g. -Idir switches:
CPPARGS		= -Wall $(ROOTCFLAGS)
### Arguments to be passed when FORTRAN compiling. E.g. -Idir switches:
FARGS		= -Wall
### Arguments to be passed while linking. E.g. -Llibdir and -llib 
### switches:
LARGS		= -Wall $(ROOTGLIBS) $(ROOTLIBS)

### Debug flags if debug mode is switched on:
ifeq "$(strip $(DEBUG))" "TRUE"
 CARGS += -g
 CPPARGS += -g
 FARGS += -g
 LARGS += -g
endif

### Static flags if static mode is switched on:
ifeq "$(strip $(STATIC))" "TRUE"
 LARGS += -static
endif

### Other stuff:
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2
ROOTGLIBS    := $(shell root-config --glibs) -lMinuit2

#################### End of user definitions section. ##################


### Extract source file names etc. from directories:
CHEADERS	= $(filter %.h,$(wildcard inc/*))
CPPHEADERS	= $(filter %.h %.H,$(wildcard inc/*))
FHEADERS	= $(filter %.inc,$(wildcard inc/*))
CSOURCES	= $(filter %.c,$(wildcard src/*))
CPPSOURCES	= $(filter %.cc %.cpp %.cxx %.c++ %.C,$(wildcard src/*))
FSOURCES	= $(filter %.f %.F,$(wildcard src/*))
CLIBSOURCES	= $(filter %.c,$(wildcard lib/*))
CPPLIBSOURCES	= $(filter %.cc %.cpp %.cxx %.c++ %.C,$(wildcard lib/*))
FLIBSOURCES	= $(filter %.f %.F,$(wildcard lib/*))
OBJECTS		= $(patsubst %.c,%.o,$(CSOURCES)) $(patsubst %.cc,%.o,$(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.c++,%.o,$(patsubst %.C,%.o,$(CPPSOURCES)))))) $(patsubst %.f,%.o,$(patsubst %.F,%.o,$(FSOURCES)))
LIBOBJECTS	= $(patsubst %.c,%.o,$(CLIBSOURCES)) $(patsubst %.cc,%.o,$(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.c++,%.o,$(patsubst %.C,%.o,$(CPPLIBSOURCES)))))) $(patsubst %.f,%.o,$(patsubst %.F,%.o,$(FLIBSOURCES)))
LIBS		= $(filter lib/lib%.a,$(wildcard lib/*))
BINARIES	= $(patsubst src/%.o,bin/%,$(OBJECTS))
CINCFLAGS	= -Iinc
CPPINCFLAGS	= -Iinc
FINCFLAGS	= -Iinc
LIBFLAGS	= -Llib $(patsubst lib/lib%.a,-l%,$(LIBS))

### SVN revision:
REVISIONN       = $(shell svn info 2>/dev/null | grep "^Revision:" | cut -f 2 -d ':' | tr -d ' \t')
REVISIONM       = $(shell svn status 2>/dev/null | grep -v "^?" | head -n 1 | cut -c 1 | tr '[:alnum:]' 'm')
ifneq "$(strip $(REVISIONN))" ""
 CINCFLAGS += -DREVISIONSTRING=\"$(REVISIONN)$(REVISIONM)\"
 CPPINCFLAGS += -DREVISIONSTRING=\"$(REVISIONN)$(REVISIONM)\"
 FINCFLAGS += -DREVISIONSTRING=\"$(REVISIONN)$(REVISIONM)\"
endif

### These entries do not correspond to file names:
.PHONY		: all debug clean Clean

### These are the secondary files:
.SECONDARY	: $(OBJECTS) $(LIBOBJECTS)

### Make everything. This is the default:
all		: $(BINARIES)

### Make everything, with debug mode:
debug		:
	$(MAKE) DEBUG=TRUE

### Make everything, dynamic linking:
dynamic		:
	rm -f $(BINARIES)
	$(MAKE) STATIC=FALSE

### Make everything, static linking:
static		:
	rm -f $(BINARIES)
	$(MAKE) STATIC=TRUE

### Link objects to create binaries:
bin/%		: src/%.o $(LIBOBJECTS) $(LIBS) Makefile
	$(USEDL) $(ARGS) $(OPT) -o $@ $< $(LIBOBJECTS) $(LIBFLAGS) $(LARGS)

### Compile library sources:
lib/%.o		: lib/%.c $(CHEADERS) Makefile
	$(USEDC) $(ARGS) $(OPT) -c -o $@ $< $(CINCFLAGS) $(CARGS)
lib/%.o		: lib/%.cc $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
lib/%.o		: lib/%.cpp $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
lib/%.o		: lib/%.cxx $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
lib/%.o		: lib/%.c++ $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
lib/%.o		: lib/%.C $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
lib/%.o		: lib/%.f $(FHEADERS) Makefile
	$(USEDF) $(ARGS) $(OPT) -c -o $@ $< $(FINCFLAGS) $(FARGS)
lib/%.o		: lib/%.F $(FHEADERS) Makefile
	$(USEDF) $(ARGS) $(OPT) -c -o $@ $< $(FINCFLAGS) $(FARGS)

### Compile program sources:
src/%.o		: src/%.c $(CHEADERS) Makefile
	$(USEDC) $(ARGS) $(OPT) -c -o $@ $< $(CINCFLAGS) $(CARGS)
src/%.o		: src/%.cc $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
src/%.o		: src/%.cpp $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
src/%.o		: src/%.cxx $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
src/%.o		: src/%.c++ $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
src/%.o		: src/%.C $(CPPHEADERS) Makefile
	$(USEDCPP) $(ARGS) $(OPT) -c -o $@ $< $(CPPINCFLAGS) $(CPPARGS)
src/%.o		: src/%.f $(FHEADERS) Makefile
	$(USEDF) $(ARGS) $(OPT) -c -o $@ $< $(FINCFLAGS) $(FARGS)
src/%.o		: src/%.F $(FHEADERS) Makefile
	$(USEDF) $(ARGS) $(OPT) -c -o $@ $< $(FINCFLAGS) $(FARGS)

### Clean up sources:
clean		:
	rm -f $(OBJECTS) $(BINARIES)
Clean		:
	rm -f $(OBJECTS) $(BINARIES) $(LIBOBJECTS)

