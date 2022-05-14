SRCS=	random.c main.c energy.c force.c init.c simonsparam.c mc.c constructnn.c analysis.c print.c switchfunction.c quaternion.c ffs.c storage.c

OBJS=	random.o main.o energy.o force.o init.o simonsparam.o mc.o constructnn.o analysis.o print.o switchfunction.o quaternion.o ffs.o storage.o

HFILES = glut.h

#

CF =   gcc
CF2 =   gcc-10
FFLAGS = -Wno-deprecated-declarations   
PARALLEL = -fopenmp
LIBS =   -lm -lGL   
OUT =	./run_test/runfile.run



opt:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -o $(OUT) $(SRCS:.c=.o) $(LIBS)

prof:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -pg -o $(OUT) $(SRCS) $(LIBS)

opt3:	$(SRCS)
	$(CF) -Wno-deprecated-declarations  -O3 -o $(OUT) $(SRCS) $(LIBS)

opt4:	$(SRCS)
	$(CF2) -Wno-deprecated-declarations  -O3 -o $(OUT) $(SRCS) $(LIBS) $(PARALLEL)

deb:	$(SRCS)
	$(CF) -g -O2 -o rods.deb $(SRCS) $(LIBS)

$(SRCS):	$(HFILES)
$(OBJS):	$(HFILES)	







