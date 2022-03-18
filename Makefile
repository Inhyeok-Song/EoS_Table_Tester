SOURCES=driver.cpp readtable.cpp linterp_some.cpp cubinterp_some.cpp \
	NuclearEoS.cpp findtoreps.cpp findtemp_NR_bisection.cpp \
	test_temp_self.cpp test_temp_self_2.cpp test_engy_self.cpp test_engy_temp_err.cpp
INCLUDES=Macro.h NuclearEoS.h
HDF5DIR=/cluster/software/hdf5-parallel/1.8.21/gcc--8.3.0/openmpi--3.1.4
HDF5INCS=-I$(HDF5DIR)/include
MPIINCS=-I/cluster/software/openmpi/3.1.4/gcc--8.3.0/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5
OBJECTS=$(SOURCES:.cpp=.o )

CXX=mpicxx
CFLAGS=-g -O3 -std=c++03
EXTRALIBS=-lm

driver: $(OBJECTS) $(INCLUDES)
	$(CXX) $(CFLAGS) -o driver $(OBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.cpp $(INCLUDES)
	$(CXX) $(CFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@

clean:
	rm -f *.o *.txt driver
