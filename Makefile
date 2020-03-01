#------------------------------------------------------------------------------
# Main MakeFile for koopman project 
#------------------------------------------------------------------------------

C:
	( cd Lib ; $(MAKE) )
	( cd Test ; $(MAKE))

all: C 

library:
	( cd Lib ; $(MAKE) )

.PHONY: clean purge

clean:
	( cd Lib ; $(MAKE) clean )
	( cd Test ; $(MAKE) clean)

purge:
	( cd Lib ; $(MAKE) purge )
	( cd Test ; $(MAKE) purge )


