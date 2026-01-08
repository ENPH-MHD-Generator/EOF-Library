SOLVER = MHDSolve
SRC    = MHDSolve.F90
LIBDIR = /usr/local/lib/elmersolver

all:
	elmerf90 -shared $(SRC) -o lib$(SOLVER).so
	cp lib$(SOLVER).so $(LIBDIR)

clean:
	rm -f lib$(SOLVER).so
