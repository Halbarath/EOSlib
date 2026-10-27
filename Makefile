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

# Header dependencies, so that no object compiled against an older header survives:
#  - DEPFLAGS: the compiler writes a .d file next to every object (foo.o -> foo.d, foo.pic.o ->
#    foo.pic.d) listing the headers it actually included; -MP adds empty rules for them, so a
#    deleted or renamed header does not stop make.
#  - headers: every header of EOSlib and the sublibraries (including the vendored interpBilinear
#    copies). All objects depend on all of them, which also covers an object that was compiled
#    elsewhere (e.g. by a sublibrary's own Makefile) and has no .d file yet. The libraries are
#    small, so rebuilding everything after a header change costs little.
#  - The Makefile itself is a prerequisite as well, so changed flags rebuild the objects.
DEPFLAGS = -MMD -MP
headers = $(wildcard *.h \
                     ../ANEOSmaterial/*.h ../ANEOSmaterial/interpBilinear/*.h \
                     ../tillotson/*.h \
                     ../reos3/*.h ../reos3/interpBilinear/*.h \
                     ../scvh/*.h ../scvh/interpBilinear/*.h)

default:
	@echo "Please specify which tool you want to make."

all: default

# Shared library for the Python wrapper (eoslib.py):  make libEOSlib
libEOSlib: libEOSlib.so

libEOSlib.so: $(pic_objects)
	$(CC) $(LDFLAGS) -shared -o $@ $^ $(LIBS)

# Objects are built next to their sources (no separate build directory).
%.o: %.c $(headers) Makefile
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

%.pic.o: %.c $(headers) Makefile
	$(CC) $(CFLAGS) $(DEPFLAGS) -fPIC -c $< -o $@

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

# The dependency files written by DEPFLAGS (missing ones are ignored).
-include $(objects:.o=.d) $(pic_objects:.o=.d) $(wildcard *.d)

.PHONY: default all libEOSlib clean cleanall

clean:
	rm -f $(execs) $(objects) $(pic_objects) libEOSlib.so $(objects:.o=.d) $(pic_objects:.o=.d) *.d

cleanall:
	rm -f $(execs) *.txt *.mat *.o *.d *.pdf libEOSlib.so $(pic_objects) $(objects) $(objects:.o=.d) $(pic_objects:.o=.d)