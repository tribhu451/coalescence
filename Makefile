
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

GSLLIBS      := $(shell gsl-config --libs)

EXTRA_FLAGS  = -D TABLE 

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) $(GSLLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp events.cpp read_input_file.cpp particles.cpp observables.cpp\
             decay_channel.cpp decay_table.cpp  pdg_properties.cpp  read_pdg.cpp\
             reso_decays.cpp inparams.cpp
                
             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = observables
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

cleanoutput : 
	rm -rf outputs/*

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
