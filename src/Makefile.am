# File generated by Zed2Soft Project Manager at Tue Jan 22 10:15:41 2019


bin_PROGRAMS = lfsrintruder deconv
lfsrintruder_CFLAGS = -I. -I../util @GLOBAL_CFLAGS@ @fftw3_CFLAGS@
lfsrintruder_LDFLAGS = @GLOBAL_LDFLAGS@

lfsrintruder_LDADD = ../util/libutil.la  @FFTW3_EXTRA_LIBS@ @fftw3_LIBS@ @GLOBAL_LDFLAGS@

lfsrintruder_SOURCES = correlator.c correlator.h lfsr.c lfsr.h lfsrdesc.c lfsrdesc.h main.c lfsrintruder.h


deconv_CFLAGS = -I. -I../util @GLOBAL_CFLAGS@
deconv_LDFLAGS = @GLOBAL_LDFLAGS@

deconv_LDADD = ../util/libutil.la  @GLOBAL_LDFLAGS@

deconv_SOURCES = lfsr.c lfsr.h lfsrdesc.c lfsrdesc.h deconv.c viterbi.c viterbi.h
