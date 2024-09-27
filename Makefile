xCXX 	= g++ 
OPT	= -g -fPIE
CFLAGS 	= -I./ $(shell /usr/local/bin/root-config --cflags --auxcflags  )
LDFLAGS	= $(shell /usr/local/bin/root-config --ldflags )
LIBS    = $(shell /usr/local/bin/root-config --libs --glibs --auxlibs )  -lMinuit -lFoam

$(info $(LDFLAGS))

TARGETS =   cover2 t12limit discover discover2 sensit3 pval



.PHONY: all clean

all: $(TARGETS)

clean: 
	-rm $(TARGETS)
	-rm *.o
	-rm *~

pval: pval.o pf.o
	$(CXX)  $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 

sensit3: sensit3.o pf.o
	$(CXX)  $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 

sensit: sensit.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 
cover: 	cover.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 
cover2: 	cover2.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 
t12limit: 	t12limit.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 
discover: 	discover.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 
discover2: 	discover2.o pf.o
	$(CXX) $(OPT) $(CFLAGS) -o $@  $(LDFLAGS) $^  $(LIBS) 

%.o: %.C
	$(CXX)  $(OPT) -o $@ -c $< $(CFLAGS) 

%.o: %.cxx
	$(CXX) $(OPT) -o $@ -c $< $(CFLAGS) 

%.cxx: %.h

