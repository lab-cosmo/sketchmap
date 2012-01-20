-include ../make.in

SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDS=$(SOURCES:.cpp=.d)
TOOLS=dimdist dimlandmark dimproj dimred 
.PHONY: all obj clean $(TOOLS)

all: $(TOOLS)

clean:
	rm -f $(OBJECTS) $(DEPENDS)

-include $(DEPENDS)
%.d: %.cpp
	$(MAKEDEPEND) $< > $@.$$$$;                      \
           sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
           rm -f $@.$$$$


# generates -- and targets -- executables in another directory
$(TOOLS): % : ../bin/%
$(TOOLS:%=../bin/%): ../bin/%: %.o libdimred.o
	$(CXX) $(CXXFLAGS) -o $@ $< libdimred.o $(LDFLAGS)
