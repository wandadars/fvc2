PROG := sort_grid
SRCS := main.f grid.f
OBJS := $(SRCS:.f=.o)
DEPS := MYDATA

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

%.o: %.f $(DEPS)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	-$(RM) $(OBJS) $(PROG) soln* elml elnorm *.vtk nodel edgl centroids edges

move_all:
	-$(CP) edgl nodel elml elnorm ../solver/

.PHONY: all clean move_all
