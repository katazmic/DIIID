#implying that GSLHOME=/usr/pppl/gsl/1.13
#     gcc  -I$(GSLHOME)/include -L/usr/pppl/gsl/1.13/lib/ -lgslcblas -lgsl -lm DIIID.o modules.o -o DIIID
ifdef FCNNG
	GSLH	=	/usr/lib64/
	GSLI	=	/usr/include/gsl/
else
	GSLH	=	$(GSLHOME)/lib/
	GSLI	=	$(GSLHOME)/include/
endif
all: DIIID clean		

clean:
	rm *.o

DIIID: DIIID.o modules.o
	gcc -pedantic -I$(GSLI) -L$(GSLH) -lgslcblas -lgsl -lm  DIIID.o modules.o -o DIIID

DIIID.o: DIIID.c modules.h
	gcc -c DIIID.c

modules.o: modules.c modules.h
	gcc -c modules.c 


