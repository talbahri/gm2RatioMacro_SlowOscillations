Target  = wiggle.exe wiggle_fit.exe RatioMethod.exe RatioMethodFit.exe
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
all:$(Target)

wiggle.exe: wiggle.cc
	g++ -o $@ wiggle.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

wiggle_fit.exe: wiggle_fit.cc
	g++ -std=c++11 -o $@ wiggle_fit.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

RatioMethod.exe: RatioMethod.cc
	g++ -o $@ RatioMethod.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

RatioMethodFit.exe: RatioMethodFit.cc
	g++ -o $@ RatioMethodFit.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
