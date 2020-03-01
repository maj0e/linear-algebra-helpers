#####################################################################
## MAKEFILE to compile the Linear Algebra Helpers Library   	   ##
#####################################################################

# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

LIBRARY = liblah

include ../Makefile.inc

RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f
I = -I../Include

AR_TARGET = $(LIBRARY).a
SO_TARGET = $(LIBRARY).so

all: install

LAH_OBJ = \
 lah_util.o \
 lah_noise.o \
 lah_matUpdate.o \
 lah_matTranspose.o \
 lah_matMul.o \
 lah_matAdd.o \
 lah_forwardSub.o \
 lah_chol.o \
 lah_cholUpdate.o \
 lah_SVD.o \
 lah_matSqrt.o \
 lah_eigenValue.o \
 lah_matConstruct.o \
 lah_QR.o \
 lah_LU.o \
 lah_solveLU.o \
 lah_properties.o
 


$(LAH_OBJ): ../Include/lah.h Makefile

%.o: ../Source/%.c ../Include/lah.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(LAH_OBJ)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@
	
$(SO_TARGET): $(LAH_OBJ)
	$(CC) $(CF) $(LAH_OBJ) -o $@ $(LDFLAGS) -shared

# install archive in this directory
install: $(AR_TARGET) $(SO_TARGET)

.PHONY: clean purge

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

