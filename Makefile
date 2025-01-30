objects = EOSlib.o igeos.o ../ANEOSmaterial/ANEOSmaterial.o ../ANEOSmaterial/interpBilinear.o ../tillotson/tillotson.o ../tillotson/tillinitlookup.o ../tillotson/tillsplint.o ../tillotson/interpol/brent.o ../reos3/reos3.o ../scvh/scvheos.o

defs = -DTILL_OUTPUT_ALL_WARNINGS -DTILL_VERBOSE -EOSLIB_VERBOSE

execs = testEOSmaterial calcTemperature calcColdcurveEnergy timingTest writeColdCurve calcPressure calcSoundSpeed

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas

# Enable the use of the external libraries
INCLUDES = -DHAVE_ANEOSMATERIAL_H -DHAVE_TILLOTSON_H -DHAVE_REOS3_H -DHAVE_SCVHEOS_H -I../ANEOSmaterial -I../tillotson -I../reos3 -I../scvh

CFLAGS ?= -O3 -Wall $(INCLUDES)

FFLAGS ?= $(CFLAGS)
LIBS ?= -lm $(GSL_LIB)

default: 
	@echo "Please specify which tool you want to make."	

all: default

testEOSmaterial: testEOSmaterial.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testigeos: testigeos.o $(objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

timingTest: timingTest.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcTemperature: calcTemperature.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	
calcPressureRhoU: calcPressureRhoU.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcPressureRhoT: calcPressureRhoT.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcEnergyRhoT: calcEnergyRhoT.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
    
calcSoundSpeed: calcSoundSpeed.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcColdcurveEnergy: calcColdcurveEnergy.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	
writeColdCurve: writeColdCurve.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

callEOSlibforMixing: callEOSlibforMixing.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

callEOSlibRhoofPT: callEOSlibRhoofPT.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects)

cleanall:
	rm -f $(execs) *.txt *.mat *.o *.pdf
