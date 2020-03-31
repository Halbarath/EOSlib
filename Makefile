objects = EOSlib.o ../ANEOSmaterial/ANEOSmaterial.o ../ANEOSmaterial/interpBilinear.o ../tillotson/tillotson.o ../tillotson/tillinitlookup.o ../tillotson/tillsplint.o ../tillotson/interpol/brent.o

defs = -DTILL_OUTPUT_ALL_WARNINGS -DTILL_VERBOSE

execs = testEOSmaterial calcTemperature calcColdcurveEnergy

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas

CFLAGS ?= -O3

FFLAGS ?= $(CFLAGS)
LIBS ?= -lm $(GSL_LIB)

default: 
	@echo "Please specify which tool you want to make."	

all: default

testEOSmaterial: testEOSmaterial.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcTemperature: calcTemperature.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcColdcurveEnergy: calcColdcurveEnergy.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects)

cleanall:
	rm -f $(execs) *.txt *.mat *.o *.pdf
