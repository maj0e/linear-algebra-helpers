include ../Makefile.inc

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/liblah.a $(LDLIBS)
TEST_INC =-I../Include/ -I./Include
BUILD_DIR = ./build


TESTS = test_matLoad test_MACROS test_matAdd test_matMul test_chol test_forwardSub test_cholUpdate test_matConstruct
LAPACK_ONLY_TESTS = test_LU test_solveLU test_SVD test_eigenValue test_matSqrt
OLD_TESTS := test_SimpleMatOperations test_Cholesky

ifdef HAVE_LAPACK
	ALL_TESTS = $(TESTS) $(LAPACK_ONLY_TESTS)
else
	ALL_TESTS = $(TESTS)
endif

TEST_PROGRAMS = $(addsuffix .test, $(ALL_TESTS))

all: $(ALL_TESTS) $(OLD_TESTS) check

$(ALL_TESTS): %: %.o
		$(CC) $(CF) $(TEST_INC) -o $(BUILD_DIR)/$@.test $(BUILD_DIR)/$< $(CS)
		
#Can be deleted later when all tests transitioned to minunit
$(OLD_TESTS): %: %.o
		$(CC) $(CF) $(TEST_INC) -o $(BUILD_DIR)/$@.test $(BUILD_DIR)/$< $(CS)

%.o : %.c
		$(CC) $(CF) $(TEST_INC) -o $(BUILD_DIR)/$@ -c $<
%.test:
	$(BUILD_DIR)/$@	
	
.PHONY: clean purge check %.test

check: $(TEST_PROGRAMS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r ./build \
    *.a *.dSYM *.obj *.dll
	- mkdir ./build
