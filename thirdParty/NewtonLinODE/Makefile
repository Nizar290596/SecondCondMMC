LIBS=-lm

#CF=gfortran
CF=f95

DESTLIBDIR=${FOAM_USER_LIBBIN}
DESTBINDIR=${FOAM_USER_APPBIN}

OBJ = dummyNewtonLinODE.o

#all: shared program install clean
all: shared install clean

#%.o: %.F90 ; $(CF) -c -fPIC -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none dummyNewtonLinODE.F90
%.o: %.F90 ; $(CF) -c -fPIC dummyNewtonLinODE.F90

shared: $(OBJ) ; $(CF) -shared -o libNewtonLinF90.so *.o

install: ; cp libNewtonLinF90.so $(DESTLIBDIR) 

clean: ; rm -f *.o *.so *.mod *~

