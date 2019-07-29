# Makefile for AncesBin 
#CC= gcc
CC= gcc
CFLAGS= -O2 -march=x86-64 -mtune=generic
LFLAGS= -lm -pthread -lz

SOURCES=bin10X.c binHiC.c binLong.c seqbin_fastq.c seqbin_split.c seqbin_reads.c seqbin_bcodes.c seqbin_binreads.c seqbin_binlong.c seqbin_getlong.c seqbin_rename.c


OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,seqbin-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) 


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe seqbin-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp seqbin-bin/binHiC .
	@cp seqbin-bin/bin10X .
	@cp seqbin-bin/binLong .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
