PROG := fvc2_solver
SRCS := data.f90 utils.f90 thermo.f90 user.f90 \
	solver_ops.f90 convection.f90 io.f90 main.f90
OBJS := $(SRCS:.f90=.o)
DEPS :=

FC := gfortran
BASEFLAGS ?= -Wno-argument-mismatch -ffree-line-length-none
OPTFLAGS ?= -O2
DBGFLAGS ?= -g -O0 -fcheck=all -fbacktrace
LDFLAGS ?=
LDLIBS ?=

RM ?= rm -f

all: release

release: FCFLAGS = $(BASEFLAGS) $(OPTFLAGS)
release: $(PROG)

debug: FCFLAGS = $(BASEFLAGS) $(DBGFLAGS)
debug: clean $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.f90 $(DEPS)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	-$(RM) $(OBJS) $(PROG)  profiles* *.mod *.smod

clean_all:
	-$(RM) $(OBJS) $(PROG)  profiles* *.mod *.smod

.PHONY: all release debug clean clean_all
