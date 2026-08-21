CF=gfortran

DESTLIBDIR=${FOAM_USER_LIBBIN}
DESTBINDIR=${FOAM_USER_APPBIN}

OBJ = aamods.o in.o out.o props.o react.o section.o solve.o species.o main.o 

all: shared program install clean

%.o: %.f90 ; $(CF) -c -fPIC -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none aamods.f90 in.f90 out.f90 props.f90 react.f90 section.f90 solve.f90 species.f90 main.f90 

shared: $(OBJ) ; $(CF) -shared -o libgenpol.so *.o

program:
# program: $(OBJ) ; $(CF) -o genpol *.o

install: ; cp libgenpol.so $(DESTLIBDIR); # cp genpol $(DESTBINDIR)

clean: ; rm -f *.o *.so *.mod *~

