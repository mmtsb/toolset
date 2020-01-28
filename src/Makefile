SHELL = /bin/sh

CC  = gcc -Wno-unused-result -Wno-write-strings
CXX = g++ -Wno-unused-result -Wno-write-strings

# try the following on MacOS:
#CC  = gcc-mp-8 -Wno-unused-result -Wno-write-strings
#CXX = g++-mp-8 -Wno-unused-result -Wno-write-strings

FC  = gfortran
F90 = gfortran

FLAGS = -DLINUX 
OPTFLAGS = -O2
CFLAGS = $(FLAGS)
COPTFLAGS = $(OPTFLAGS)
CXXFLAGS = $(FLAGS)
CXXOPTFLAGS = $(OPTFLAGS)
FFLAGS = $(FLAGS) 
FOPTFLAGS = $(OPTFLAGS)

BINDIR  = ../bin
ODIR    = obj
DATADIR = ../data

UTIL = util
NR   = numrec

NONE :
	@echo ""
	@echo "please specify target:"
	@echo ""
	@echo "    all"
	@echo "    mfold mener rebuild rebprimo mchain gdt"
	@echo "    jclust kclust spear solvate "
	@echo "    reducescwrlib scwrlbinary scwrllib"
	@echo ""


NROBJS = $(ODIR)/brent.o $(ODIR)/f1dim.o $(ODIR)/frprmn.o\
         $(ODIR)/linmin.o $(ODIR)/mnbrak.o $(ODIR)/nrutil.o 

all : $(BINDIR) $(ODIR) $(DATADIR) mfold mener rebuild rebprimo mchain gdt kclust jclust spear reducescwrlib scwrlbinary scwrllib solvate;

$(BINDIR):
	mkdir -p $(BINDIR)

$(ODIR):
	mkdir -p $(ODIR)

$(DATADIR):
	mkdir -p $(DATADIR)

rebuild : $(BINDIR)/rebuild;
$(BINDIR)/rebuild : rebuild.c\
            $(UTIL)/vector.h $(UTIL)/aa.h $(UTIL)/pdb.h $(UTIL)/field.h \
            $(ODIR)/aa.o $(ODIR)/pdb.o $(ODIR)/field.o $(NROBJS)
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(NR) -I$(UTIL) \
	-o $@ rebuild.c $(NROBJS) $(ODIR)/aa.o $(ODIR)/pdb.o $(ODIR)/field.o -lm 

rebprimo : $(BINDIR)/rebprimo;
$(BINDIR)/rebprimo : rebprimo.c\
            $(UTIL)/vector.h $(UTIL)/aa.h $(UTIL)/pdb.h $(UTIL)/field.h \
            $(ODIR)/aa.o $(ODIR)/pdb.o $(ODIR)/field.o 
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ rebprimo.c $(ODIR)/aa.o $(ODIR)/pdb.o $(ODIR)/field.o -lm

mchain : $(BINDIR)/mchain;
$(BINDIR)/mchain : mchain.c\
                $(UTIL)/vector.h $(UTIL)/pdb.h $(ODIR)/pdb.o \
		$(UTIL)/field.h $(ODIR)/field.o
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ mchain.c $(ODIR)/pdb.o $(ODIR)/field.o

jclust : $(BINDIR)/jclust;
$(BINDIR)/jclust : jclust.c\
                $(UTIL)/vector.h $(UTIL)/pdb.h $(ODIR)/pdb.o\
                $(ODIR)/charmmlsq.o $(ODIR)/field.o $(UTIL)/field.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ jclust.c $(ODIR)/pdb.o $(ODIR)/charmmlsq.o $(ODIR)/field.o $(FORLINK) -lm -lgfortran

kclust : $(BINDIR)/kclust;
$(BINDIR)/kclust : kclust.c\
                $(UTIL)/vector.h $(UTIL)/pdb.h $(ODIR)/pdb.o\
                $(ODIR)/charmmlsq.o $(ODIR)/field.o $(UTIL)/field.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ kclust.c $(ODIR)/pdb.o $(ODIR)/charmmlsq.o $(ODIR)/field.o $(FORLINK) -lm -lgfortran


gdt : $(BINDIR)/gdt;
$(BINDIR)/gdt : gdt.c\
         $(ODIR)/charmmlsq.o $(UTIL)/vector.h $(UTIL)/pdb.h $(ODIR)/pdb.o\
         $(ODIR)/field.o $(UTIL)/field.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ gdt.c $(ODIR)/pdb.o $(ODIR)/charmmlsq.o $(ODIR)/field.o $(FORLINK) -lm -lgfortran


solvate : $(BINDIR)/solvate;
$(BINDIR)/solvate : solvate.c\
         $(UTIL)/vector.h $(UTIL)/pdb.h $(ODIR)/pdb.o\
         $(ODIR)/field.o $(UTIL)/field.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ solvate.c $(ODIR)/pdb.o $(ODIR)/field.o $(FORLINK) -lm

spear : $(BINDIR)/spear;
$(BINDIR)/spear : spear.c
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -o $@ spear.c -lm

reducescwrlib : $(BINDIR)/reducescwrlib;
$(BINDIR)/reducescwrlib : reducescwrlib.c
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) \
	-o $@ reducescwrlib.c

scwrlbinary : $(BINDIR)/scwrlbinary;
$(BINDIR)/scwrlbinary : scwrlbinary.c
	$(CC) $(CFLAGS) $(LDFLAGS) \
	-o $@ scwrlbinary.c

scwrllib : $(DATADIR)/scwrlbin34.red.lib $(DATADIR)/scwrlbin34.lib;
$(DATADIR)/scwrlbin34.lib: $(BINDIR)/scwrlbinary
	gunzip -c ../data/bbdep01.Jul.sortlib.gz | $(BINDIR)/scwrlbinary; \
	mv scwrlbin34.lib $(DATADIR); rm README.TO.CHECK.LIBRARY
$(DATADIR)/scwrlbin34.red.lib: $(BINDIR)/reducescwrlib $(DATADIR)/scwrlbin34.lib
	$(BINDIR)/reducescwrlib $(DATADIR)/scwrlbin34.lib $(DATADIR)/scwrlbin34.red.lib

mfold : $(BINDIR)/mfold;
$(BINDIR)/mfold : mfold.f $(ODIR)/pclient.o
	$(FC) $(FFLAGS) $(FOPTFLAGS) $(LDFLAGS) -o $@ $(ODIR)/pclient.o mfold.f

mener : $(BINDIR)/mener;
$(BINDIR)/mener : mener.c $(ODIR)/melib.o $(ODIR)/field.o
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) $(LDFLAGS) -I$(UTIL) \
	-o $@ $(ODIR)/melib.o mener.c $(FORLINK) $(ODIR)/field.o -lm -lgfortran

$(ODIR)/melib.o : $(UTIL)/melib.f
	$(FC) $(FFLAGS) $(FOPTFLAGS) -o $@ -c $(UTIL)/melib.f

$(ODIR)/pclient.o : $(UTIL)/pclient.c
	$(CC) $(CFLAGS) -o $@ -c $(UTIL)/pclient.c

$(ODIR)/brent.o : $(NR)/brent.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/brent.c 

$(ODIR)/f1dim.o : $(NR)/f1dim.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/f1dim.c

$(ODIR)/frprmn.o : $(NR)/frprmn.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/frprmn.c

$(ODIR)/linmin.o : $(NR)/linmin.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/linmin.c

$(ODIR)/mnbrak.o : $(NR)/mnbrak.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/mnbrak.c

$(ODIR)/nrutil.o : $(NR)/nrutil.c $(NR)/nrutil.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(NR) -o $@ -c $(NR)/nrutil.c

$(ODIR)/aa.o : $(UTIL)/aa.c $(UTIL)/aa.h $(UTIL)/vector.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(UTIL) \
	-o $@ -c $(UTIL)/aa.c

$(ODIR)/field.o : $(UTIL)/field.c $(UTIL)/field.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(UTIL) \
	-o $@ -c $(UTIL)/field.c

$(ODIR)/pdb.o : $(UTIL)/vector.h $(UTIL)/pdb.c $(UTIL)/pdb.h
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -I$(UTIL) \
	-o $@ -c $(UTIL)/pdb.c

$(ODIR)/charmmlsq.o : $(UTIL)/charmmlsq.f
	$(FC) $(FFLAGS) $(FOPTFLAGS) -o $@ -c $(UTIL)/charmmlsq.f

