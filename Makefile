
CC = g++
# Set this next flag to enable fractional powers for Units
# This needs a compiler with C++20 support AND MAY INCREASE compilation time
# If you don't need fractional powers but want other c++20 support, set to 0
# but set DEFAULT_STD to c++20
USE_FRACTIONAL_POWERS = 1
DEFAULT_STD = c++17

#Set this for compilers which dont support templating friendship, such as GCC pre 13
#NO_SUPPORT_TEMPLATE_FRIENDSHIP = 1

INCLUDE = -I ./include/
SRCDIR = src
OBJDIR = obj

CFLAGS = -g -c $(INCLUDE) -O3
CFLAGS += -DNO_NARROWING_CONVERSIONS
ifeq ($(strip $(USE_FRACTIONAL_POWERS)),1)
  CFLAGS += -std=c++20 -DUSE_FRACTIONAL_POWERS
  $(info *****Enabled fractional powers*****)
else
  CFLAGS += -std=$(DEFAULT_STD)
  $(info *****Disabled fractional powers*****)
endif

DEBUG = -W -Wall -O0 -pedantic -D_GLIBCXX_DEBUG
#DEBUG+= -Wno-sign-compare
DEBUG+= -Wno-unused-parameter

EXENAME = Test
DEPS = -MM

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG) -DDEBUG

endif
ifeq ($(strip $(TIMER)),true)
  CFLAGS += -DRUN_TIMER_TEST

endif
#list of all header and cpp pairs. Add new files here.....
INCLS =

#make lists of source and object files, all headers plus main
SOURCE := $(INCLS:.h=.cpp)
SOURCE += test.cpp
OBJS := $(SOURCE:.cpp=.o)

#header files only (no .cpp)
INCLS += helper.hpp SimpleFrac.hpp UnitCheckedType.hpp StorageTypes.hpp

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
INCLS := $(addprefix include/, $(INCLS))

VPATH = include/

default : $(OBJS)
	$(CC) $(OBJS) -o $(EXENAME)

.PHONY: deps
deps: clean
	$(CC) $(CFLAGS) $(DEPS) $(SOURCE)


#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

obj/%.o:./src/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

#Dependencies
$(OBJDIR)/test.o : $(SRCDIR)/test.cpp $(INCLS)

.PHONY : tar clean veryclean

clean:
	$(RM) $(EXENAME) $(OBJS)

veryclean:
	$(RM) $(EXENAME) $(OBJS)
	$(RM) -r $(OBJDIR)

tar:
	tar -cvf $(EXENAME).tar $(SOURCE) $(INCLS) Makefile

