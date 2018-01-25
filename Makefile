##################################################
# makefile for the Coverage & Anisotropy Toolkit #
##################################################

SYSTEM = $(shell uname)
CC = gcc 
CXX = g++
OPTFLAGS = -O2
#DBGFLAGS = -g3
WFLAGS = -D__USE_FIXED_PROTOTYPES__ -Wall
OBJ = ./

.PHONY: clean
#------------------- defs ------------------------------
CFITSIO_INCDIR = $(CFITSIO_DIR)/include
CFITSIO_LIBDIR = $(CFITSIO_DIR)/lib

CHEALPIX_INCDIR = $(HEALPIX_DIR)/include
CHEALPIX_LIBDIR = $(HEALPIX_DIR)/lib

STC_INCDIR = $(STC_DIR)
STC_LIBDIR = $(STC_DIR)

INCDIR = -I$(CFITSIO_INCDIR) \
	-I$(CHEALPIX_INCDIR) \
	-I$(STC_INCDIR) \
	$(shell root-config --cflags)

LIBDIR = -L${CHEALPIX_LIBDIR} -lchealpix -lm \
	-L${CFITSIO_LIBDIR} -lcfitsio \
	-L${STC_LIBDIR} -lSTCoordinates \
	$(shell root-config --glibs)
#-------------------------------------------------------

#------- alias -----------------------------------------
libobjs = \
	healpixmap.o \
	harmotools.o \
	lobe.o \
	maptools.o \
	projmap.o \
	coverage.o \
	angdist.o \
	timemod.o \
	events.o \
	simuevents.o \
	Cl.o \
	lima.o \
	correlation.o \
	common.o \
	dipoleHiRes.o \
	fitdipole.o \
	rayleigh.o \
	userfcn.o \
	agntools.o \
	scanner.o \
	crosscorrelation.o

execs = \
	example_blindsearch.exe \
	example_Cl.exe \
	example_bllac.exe \
	example_correlation.exe \
	compute_lobe.exe \
	example_simuDipole.exe \
	example_simuSources.exe \
	example_simuAcceptance.exe \
	example_dipoleHiRes.exe \
	example_fitDipole.exe \
	example_rayleigh.exe \
	example_phiMod.exe \
	example_simuAccDip.exe \
	example_prescription.exe \
	example_deflexion.exe \
	example_spectrum.exe \
	example_spectrumSimu.exe \
	example_plot2MASSvsVCV.exe \
	example_percolation.exe \
	example_simuAnisotropic.exe \
	example_drawSources.exe \
	example_ironPrescription.exe \
	example_MC.exe \
	example_analyzeMC.exe \
	example_mapMC.exe \
	example_crossCorrelation.exe \
	example_LemoineWaxman.exe \
	example_testLiMa.exe \
	ScanWithoutCenA.exe \
	Magnification.exe \
	ProducePixFileForGMF.exe 

exeobjs = $(patsubst %.exe,%.o,$(execs))

HEADERS = $(patsubst %.o,%.h,$(libobjs))

thelib = libtkit.a
#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the library sources
$(OBJ)%.o:%.cc %.h
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<

# rules for the executable sources
$(OBJ)%.o:%.cc
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all : lib $(execs)
lib: $(thelib) userfcn.h
clean :
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(thelib) $(execs) $(libobjs) $(exeobjs)
	@echo "Done"
#-------------------------------------------------------

#-------- specific rules -------------------------------
$(thelib): $(libobjs)
	@echo "Making "$(thelib)
	@ar r $@ $^
	@ranlib $@
	@echo "Done."

example_simuAcceptance.exe : example_simuAcceptance.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_blindsearch.exe : example_blindsearch.o $(thelib)
	$(CXX) $(DBGFLAGS) -o $@ $^ $(LIBDIR)
example_simuSources.exe : example_simuSources.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_simuDipole.exe : example_simuDipole.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR) -lMinuit
example_Cl.exe : example_Cl.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_bllac.exe : example_bllac.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_correlation.exe : example_correlation.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
compute_lobe.exe : compute_lobe.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_dipoleHiRes.exe : example_dipoleHiRes.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_fitDipole.exe : example_fitDipole.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR) -lMinuit
example_rayleigh.exe : example_rayleigh.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_phiMod.exe : example_phiMod.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_simuAccDip.exe : example_simuAccDip.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR) -lMinuit
example_prescription.exe : example_prescription.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_deflexion.exe : example_deflexion.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_spectrum.exe : example_spectrum.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_spectrumSimu.exe : example_spectrumSimu.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_plot2MASSvsVCV.exe: example_plot2MASSvsVCV.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_percolation.exe: example_percolation.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_simuAnisotropic.exe: example_simuAnisotropic.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_drawSources.exe: example_drawSources.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_ironPrescription.exe: example_ironPrescription.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_MC.exe: example_MC.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_mapMC.exe: example_mapMC.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_analyzeMC.exe: example_analyzeMC.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_crossCorrelation.exe: example_crossCorrelation.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_LemoineWaxman.exe: example_LemoineWaxman.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_testLiMa.exe: example_testLiMa.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_GalacticSources.exe: example_GalacticSources.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
ScanWithoutCenA.exe: ScanWithoutCenA.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
Magnification.exe: Magnification.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)  
ProducePixFileForGMF.exe: ProducePixFileForGMF.o $(thelib) 
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------

