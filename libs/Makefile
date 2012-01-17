-include ../../make.in
.PHONY: all clean 
SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDS=$(SOURCES:.cpp=.d)

all: libtoolbox.a 

libtoolbox.a: $(OBJECTS)
	ar rcs libtoolbox.a $(OBJECTS)

-include $(DEPENDS)
%.d: %.cpp
	$(MAKEDEPEND) $< > $@.$$$$;                      \
	   sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
	   rm -f $@.$$$$

clean: 
	rm -f $(OBJECTS) $(DEPENDS)
