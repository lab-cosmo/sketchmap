.PHONY: all clean libs tools distclean

all: libs tools

tools: libs
	$(MAKE) -C tools all
	
libs:
	$(MAKE) -C libs all
	
clean: 
	$(MAKE) -C libs clean
	$(MAKE) -C tools clean
	
distclean: clean
	rm -f bin/*
