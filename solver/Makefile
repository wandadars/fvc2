PROG := fvc2_solver
SRCS := main.f solve.f convection.f user.f io.f
OBJS := $(SRCS:.f=.o)
DEPS := MYDATA

FC ?= gfortran
FCFLAGS ?= -Wno-argument-mismatch -g -ffree-line-length-none -free
FCFLAGS += $(OPTFLAGS) $(DBGFLAGS)
LDFLAGS ?=
LDLIBS ?=

RM ?= rm -f

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.f $(DEPS)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	-$(RM) $(OBJS) $(PROG) soln* profiles*

clean_all:
	-$(RM) $(OBJS) $(PROG) soln* edgl elml elnorm nodel profiles*

.PHONY: all clean clean_all
