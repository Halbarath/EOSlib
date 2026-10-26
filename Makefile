objects = EOSlib.o igeos.o ../ANEOSmaterial/ANEOSmaterial.o ../tillotson/tillotson.o ../tillotson/tillinitlookup.o ../tillotson/tillsplint.o ../reos3/reos3.o ../scvh/scvheos.o

# Position-independent copies of the same objects for the shared library.
pic_objects = $(objects:.o=.pic.o)

# -DEOSLIB_VERBOSE makes EOSlib and every sublibrary print a diagnostic for
# each failed lookup (see "Failure contract" in EOSlib.h). Off by default:
# a failed lookup is counted in nLookupFail either way.
DEFS =

execs = testEOSmaterial calcTemperature calcColdcurveEnergy timingTest writeColdCurve calcPressure calcSoundSpeed

# Enable the use of the external libraries
INCLUDES = -DHAVE_ANEOSMATERIAL_H -DHAVE_TILLOTSON_H -DHAVE_REOS3_H -DHAVE_SCVHEOS_H -I../ANEOSmaterial -I../tillotson -I../reos3 -I../scvh

CFLAGS ?= -O3 -Wall $(DEFS) $(INCLUDES)

FFLAGS ?= $(CFLAGS)
LIBS ?= -lm

default: 
	@echo "Please specify which tool you want to make."	

all: default

# Shared library for the Python wrapper (eoslib.py):  make libEOSlib
libEOSlib: libEOSlib.so

libEOSlib.so: $(pic_objects)
	$(CC) $(LDFLAGS) -shared -o $@ $^ $(LIBS)

%.pic.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

testEOSmaterial: testEOSmaterial.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testigeos: testigeos.o $(objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

timingTest: timingTest.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	
calcPressureRhoU: calcPressureRhoU.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcTempRhoU: calcTempRhoU.o $(objects)
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

.PHONY: default all libEOSlib clean cleanall

clean:
	rm -f $(execs) $(objects) $(pic_objects) libEOSlib.so

cleanall:
	rm -f $(execs) *.txt *.mat *.o *.pdf libEOSlib.so $(pic_objects)