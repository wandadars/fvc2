PROG := msh2fvc
SRCS := msh2fvc.f90
OBJS := $(SRCS:.f90=.o)
DEPS :=

FC ?= gfortran
FCFLAGS ?= -g
FCFLAGS += $(OPTFLAGS) $(DBGFLAGS)
LDFLAGS ?=
LDLIBS ?=

RM ?= rm -f
CP ?= cp

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.f90 $(DEPS)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	-$(RM) $(OBJS) $(PROG) *.mod *.smod

.PHONY: all clean
