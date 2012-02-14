SYSTEM     = x86_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio_Academic123/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio_Academic123/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0
CC  = gcc -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m32 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
COPT  = -m32 -fPIC 

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -m32 -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -m32 -lm -pthread
JAVA      = java  -Djava.library.path=$(CPLEXDIR)/bin/x86_sles10_4.1 -classpath $(CPLEXJARDIR):


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)

C_EX = lpex1 

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CX_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
#
# The examples
#

lptest1: lptest1.o
	$(CC) $(CFLAGS) lptest1.o -o lptest1 $(CLNFLAGS)

lptest1.o: ./lptest1.c
	$(CC) -c $(CFLAGS) ./lptest1.c -o lptest1.o

formInicial: formInicial.o
	$(CCC) $(CFLAGS) formInicial.o -o formInicial $(CLNFLAGS)

formInicial.o: ./formInicial.cpp
	$(CCC) -c $(CFLAGS) ./formInicial.cpp -o formInicial.o

