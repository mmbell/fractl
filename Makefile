CXX = /usr/bin/g++
CXX_FLAGS = -gdwarf-4 -std=c++11

TDRP_GEN = /usr/local/lrose/bin/tdrp_gen

GEOGRAPHIC_LIB_DIR = /home/bpmelli/.linuxbrew/Cellar/geographiclib/1.48/lib
ZIP_LIB_DIR        = /home/bpmelli/.linuxbrew/Cellar/libzip/1.3.0/lib
LROSE_LIB_DIR      = /usr/local/lrose-20170929.x86_64/lib

GEOGRAPHIC_LIB = $(GEOGRAPHIC_LIB_DIR)/libGeographic.so.17.1.1
ZIP_LIB        =  $(ZIP_LIB_DIR)/libzip.so.5.0.0

INCLUDES = -I/home/bpmelli/.linuxbrew/include \
	-I/usr/local/lrose-20170929.x86_64/include \
	-I/home/bpmelli/.linuxbrew/include/eigen3 
LIBS = -L$(LROSE_LIB_DIR) \
	-ltdrp -lkd -lNcxx -lnetcdf_c++ -lRadx -lhdf5_cpp -lhdf5 \
	$(GEOGRAPHIC_LIB) $(ZIP_LIB)

LINK_FLAGS= \
	-Wl,-rpath,$(GEOGRAPHIC_LIB_DIR) \
	-Wl,-rpath,$(ZIP_LIB_DIR) \
	-Wl,-rpath,$(LROSE_LIB_DIR)

OBJS = Args.o Fractl.o readRadx.o readPreGridded.o calcWinds.o \
       writeNetCdf.o Params.o

all: Fractl

Fractl: $(OBJS)
	$(CXX) $(CXX_FLAGS) -o Fractl $(OBJS) $(LIBS) $(LINK_FLAGS)

%.o : %.cc
	$(CXX) $(CXX_FLAGS) $(INCLUDES) -c $< -o $@

Args.cc: Params.hh
Params.o: Params.cc

Params.hh: paramdef.Fractl
	$(TDRP_GEN) -c++ -f paramdef.Fractl

clean:
	rm -f *.o Fractal Params.hh Params.cc
