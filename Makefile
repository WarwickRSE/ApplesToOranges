
CC = g++
INCLUDE = -I ./include/
SRCDIR = src
OBJDIR = obj
CFLAGS = -g -c $(INCLUDE) -std=c++20 -O0
DEBUG = -W -Wall -O0 -pedantic -D_GLIBCXX_DEBUG
#DEBUG+= -Wno-sign-compare
DEBUG+= -Wno-unused-parameter

EXENAME = Test
DEPS = -MM

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
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

.PHONY : tar clean veryclean

clean:
	$(RM) $(EXENAME) $(OBJS)

veryclean:
	$(RM) $(EXENAME) $(OBJS)
	$(RM) -r $(OBJDIR)

tar:
	tar -cvf $(EXENAME).tar $(SOURCE) $(INCLS) Makefile

