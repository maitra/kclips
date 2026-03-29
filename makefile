CC =	gcc

CFLAGS = -Wall -pedantic -O3 -std=c11
#CFLAGS = -Wall -pedantic -g -std=c11
LDFLAGS = -lm -lpthread

# Local variables
srcs = $(wildcard *.c rmath-stalone/*.c)
hds = $(wildcard *.h)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.d)

objsrun = $(objs:run_kclips.o=) #program to be excluded is run_kclips.c

run_kclips_estk:	$(objsrun)
	$(CC) -o run_kclips_estk $(objsrun) $(CFLAGS) $(LDFLAGS)

archive: 
	tar czvf kmeans_init.tar.gz $(srcs) $(hds) Makefile

clean:	
	-rm $(objs) *~ run_kclips_estk 2>/dev/null

clear:
	rm $(objs) 

dummy:
	@echo $(objs)
