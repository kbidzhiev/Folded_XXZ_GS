#################################################################
#################################################################
#################################################################


include ./this_dir.mk
include ./options.mk


#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/core.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)



#Rules ------------------

%.o: %.cpp $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cpp $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

3siteHam_GS.exe: 3siteHam_GS.cpp tech_param.cpp Folded_XXZ.cpp XXZ.cpp observables_GS.cpp quantum_measurement.cpp $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 3siteHam_GS.cpp tech_param.cpp Folded_XXZ.cpp XXZ.cpp observables_GS.cpp quantum_measurement.cpp -o $@ $(LIBFLAGS)


all: 3siteHam_GS.exe
clean:
	rm *.exe *.o

