# CHOOCH Makefile
# Created on Oct 30 1999 by Gwyndaf Evans
#
# Edit the three first directory definitions to specify
# a) Edit the host type
# b)  the directory where GSL (Gnu Scientific Library) is kept
# c)  the directory where Cgraph (PS plotting library) is kept
# d)  and the directory where you would like you executables to go (BINDIR).
#
ARCH   = Linux
#ARCH   = OSF1
#ARCH   = SunOS
#
PGPLOTDIR  = /usr/local/pgplot
DISLIN = /usr/local/dislin
GSLDIR = /usr/local/lib
CGRAPHDIR = /usr/local/lib
INCLUDE   = /usr/local/pgplot
X11LIBDIR  = /usr/X11R6/lib
######################################
#
VERSION = 5.0.6
CGRAPH = -lcgraph
LIBS = -lgsl -lgslcblas 
PGLIBS =  -lcpgplot -lpgplot
DISLINLIBS= -ldislnc
EXE    = chooch-$(VERSION).$(ARCH)
EXEPG    = chooch-$(VERSION)-pg.$(ARCH)

PYTHON_INCLUDE_DIR ?= $(shell python -c 'from distutils.sysconfig import get_python_inc; print(get_python_inc())')
PYTHON_DEST_DIR ?= $(shell python -c 'import site; print(site.getsitepackages()[0])')

#
# How to compile and link
#
include Makefile.$(ARCH)
#
# Basic definitions
#
RM    = /bin/rm
MV    = /bin/mv
CP    = /bin/cp
#
#
#OBJECTS = main.o      fluread.o printbanner.o minmax.o  spline.o \
#          mucal.o     fdprime.o smooth.o      fits.o    normalize.o \
#          checks.o    usage.o   integrate.o   psplot.o  selwavel.o \
#          copyright.o toplot.o  license.c     pngplt.o  savwin.o
#
#
PYOBJECTS = printbanner.o minmax.o  spline.o \
        mucal.o     fdprime.o smooth.o      fits.o    normalize.o \
        checks.o    integrate.o   selwavel.o \
        toplot.o  savwin.o
OBJECTS = main.o      fluread.o printbanner.o minmax.o  spline.o \
	mucal.o     fdprime.o smooth.o      fits.o    normalize.o \
	checks.o    usage.o   integrate.o   selwavel.o \
	copyright.o toplot.o  license.c     savwin.o

chooch : clean ${OBJECTS} Makefile
	$(CC) -o ${EXE} ${CFLAGS} ${OBJECTS} $(LDFLAGS)

pychooch : clean ${OBJECTS} 
	$(CC) -shared -o PyChooch.so PyChooch.c ${CFLAGS} ${PYOBJECTS} $(LDFLAGS)

chooch-pg : 
	make chooch-with-pgplot "CFLAGS = -I$(INCLUDE) $(CFLAGS) -DPGPLOT"

chooch-with-pgplot : clean ${OBJECTS} Makefile
	$(FC) -v $(CFLAGS) -o ${EXEPG} ${OBJECTS} $(LDFLAGS)
#
all: chooch chooch-pg
#
install : pychooch
	$(MV) PyChooch.so $(PYTHON_DEST_DIR)
#
clean :
	${RM} -f *.o
#
# End
#
