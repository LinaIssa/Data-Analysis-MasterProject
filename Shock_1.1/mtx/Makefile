# Lapack makefile
LIB=libmtx.a
LIB_TOOL = ar crs

PACKAGE_DIR = .
# compilo='intel'
compilo='gfortran'
ifeq ($(compilo), 'intel')
FC=ifort
FLAGS= -fixed -132
FREEFLAGS=-free -r8
else
FC =  gfortran 
FLAGS =  -ffixed-form -ffixed-line-length-132
FREEFLAGS=  -ffree-form  -fdefault-double-8 -fdefault-real-8 
endif
SRC_DIR = $(PACKAGE_DIR)/src
INCLUDES = -I$(PACKAGE_DIR)/include
#OTHER_INCLUDES = -I$(FRAME_DIR)/include

LIB_OBJS =  caxpy.o ccopy.o cdotc.o cdotu.o cgbmv.o cgemm.o cgemv.o cgerc.o \
 cgeru.o chbmv.o chemm.o chemv.o cher.o cher2.o cher2k.o cherk.o \
 chpmv.o chpr.o chpr2.o crotg.o cscal.o csrot.o csscal.o cswap.o \
 csymm.o csyr2k.o csyrk.o ctbmv.o ctbsv.o ctpmv.o ctpsv.o ctrmm.o \
 ctrmv.o ctrsm.o ctrsv.o dasum.o daxpy.o dbdsdc.o dbdsqr.o \
 dcabs1.o dcopy.o ddisna.o ddot.o dgbbrd.o dgbcon.o dgbequ.o \
 dgbmv.o dgbrfs.o dgbsv.o dgbsvx.o dgbtf2.o dgbtrf.o dgbtrs.o \
 dgebak.o dgebal.o dgebd2.o dgebrd.o dgecon.o dgeequ.o dgees.o \
 dgeesx.o dgeev.o dgeevx.o dgegs.o dgegv.o dgehd2.o dgehrd.o \
 dgelq2.o dgelqf.o dgels.o dgelsd.o dgelss.o dgelsx.o dgelsy.o \
 dgemm.o dgemv.o dgeql2.o dgeqlf.o dgeqp3.o dgeqpf.o dgeqr2.o \
 dgeqrf.o dger.o dgerfs.o dgerq2.o dgerqf.o dgesc2.o dgesdd.o \
 dgesv.o dgesvd.o dgesvx.o dgetc2.o dgetf2.o dgetrf.o dgetri.o  \
 dgetrs.o dggbak.o dggbal.o dgges.o dggesx.o dggev.o dggevx.o  \
 dggglm.o dgghrd.o dgglse.o dggqrf.o dggrqf.o dggsvd.o dggsvp.o  \
 dgtcon.o dgtrfs.o dgtsv.o dgtsvx.o dgttrf.o dgttrs.o dgtts2.o  \
 dhgeqz.o dhsein.o dhseqr.o disnan.o dlabad.o dlabrd.o dlacn2.o  \
 dlacon.o dlacpy.o dladiv.o dlae2.o dlaebz.o dlaed0.o dlaed1.o  \
 dlaed2.o dlaed3.o dlaed4.o dlaed5.o dlaed6.o dlaed7.o dlaed8.o  \
 dlaed9.o dlaeda.o dlaein.o dlaev2.o dlaexc.o dlag2.o dlag2s.o  \
 dlags2.o dlagtf.o dlagtm.o dlagts.o dlagv2.o dlahqr.o dlahr2.o  \
 dlahrd.o dlaic1.o dlaisnan.o dlaln2.o dlals0.o dlalsa.o  \
 dlalsd.o dlamch.o dlamrg.o dlaneg.o dlangb.o dlange.o dlangt.o  \
 dlanhs.o dlansb.o dlansp.o dlanst.o dlansy.o dlantb.o dlantp.o  \
 dlantr.o dlanv2.o dlapll.o dlapmt.o dlapy2.o dlapy3.o dlaqgb.o  \
 dlaqge.o dlaqp2.o dlaqps.o dlaqr0.o dlaqr1.o dlaqr2.o dlaqr3.o  \
 dlaqr4.o dlaqr5.o dlaqsb.o dlaqsp.o dlaqsy.o dlaqtr.o dlar1v.o  \
 dlar2v.o dlarf.o dlarfb.o dlarfg.o dlarft.o dlarfx.o dlargv.o  \
 dlarnv.o dlarra.o dlarrb.o dlarrc.o dlarrd.o dlarre.o  \
 dlarrf.o dlarrj.o dlarrk.o dlarrr.o dlarrv.o dlartg.o  \
 dlartv.o dlaruv.o dlarz.o dlarzb.o dlarzt.o dlas2.o dlascl.o  \
 dlasd0.o dlasd1.o dlasd2.o dlasd3.o dlasd4.o dlasd5.o  \
 dlasd6.o dlasd7.o dlasd8.o dlasda.o dlasdq.o dlasdt.o  \
 dlaset.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o  \
 dlasq6.o dlasr.o dlasrt.o dlassq.o dlasv2.o dlaswp.o  \
 dlasy2.o dlasyf.o dlatbs.o dlatdf.o dlatps.o dlatrd.o  \
 dlatrs.o dlatrz.o dlatzm.o dlauu2.o dlauum.o dlazq3.o  \
 dlazq4.o dnrm2.o dopgtr.o dopmtr.o dorg2l.o dorg2r.o  \
 dorgbr.o dorghr.o dorgl2.o dorglq.o dorgql.o dorgqr.o  \
 dorgr2.o dorgrq.o dorgtr.o dorm2l.o dorm2r.o dormbr.o  \
 dormhr.o dorml2.o dormlq.o dormql.o dormqr.o dormr2.o  \
 dormr3.o dormrq.o dormrz.o dormtr.o dpbcon.o dpbequ.o  \
 dpbrfs.o dpbstf.o dpbsv.o dpbsvx.o dpbtf2.o dpbtrf.o  \
 dpbtrs.o dpocon.o dpoequ.o dporfs.o dposv.o dposvx.o  \
 dpotf2.o dpotrf.o dpotri.o dpotrs.o dppcon.o dppequ.o  \
 dpprfs.o dppsv.o dppsvx.o dpptrf.o dpptri.o dpptrs.o  \
 dptcon.o dpteqr.o dptrfs.o dptsv.o dptsvx.o dpttrf.o  \
 dpttrs.o dptts2.o drot.o drotg.o drotm.o drotmg.o drscl.o  \
 dsbev.o dsbevd.o dsbevx.o dsbgst.o dsbgv.o dsbgvd.o  \
 dsbgvx.o dsbmv.o dsbtrd.o dscal.o dsdot.o dsgesv.o dspcon.o  \
 dspev.o dspevd.o dspevx.o dspgst.o dspgv.o dspgvd.o dspgvx.o  \
 dspmv.o dspr.o dspr2.o dsprfs.o dspsv.o dspsvx.o dsptrd.o  \
 dsptrf.o dsptri.o dsptrs.o dstebz.o dstedc.o dstegr.o  \
 dstein.o dstemr.o dsteqr.o dsterf.o dstev.o dstevd.o dstevr.o  \
 dstevx.o dswap.o dsycon.o dsyev.o dsyevd.o dsyevr.o dsyevx.o  \
 dsygs2.o dsygst.o dsygv.o dsygvd.o dsygvx.o dsymm.o dsymv.o  \
 dsyr.o dsyr2.o dsyr2k.o dsyrfs.o dsyrk.o dsysv.o dsysvx.o  \
 dsytd2.o dsytf2.o dsytrd.o dsytrf.o dsytri.o dsytrs.o dtbcon.o  \
 dtbmv.o dtbrfs.o dtbsv.o dtbtrs.o dtgevc.o dtgex2.o dtgexc.o  \
 dtgsen.o dtgsja.o dtgsna.o dtgsy2.o dtgsyl.o dtpcon.o dtpmv.o  \
 dtprfs.o dtpsv.o dtptri.o dtptrs.o dtrcon.o dtrevc.o dtrexc.o  \
 dtrmm.o dtrmv.o dtrrfs.o dtrsen.o dtrsm.o dtrsna.o dtrsv.o  \
 dtrsyl.o dtrti2.o dtrtri.o dtrtrs.o dtzrqf.o dtzrzf.o dzasum.o  \
 dznrm2.o dzsum1.o icamax.o idamax.o ieeeck.o ilaenv.o isamax.o  \
 izamax.o lsame.o sasum.o saxpy.o scabs1.o scasum.o scnrm2.o  \
 scopy.o sdot.o sdsdot.o sgbmv.o sgemm.o sgemv.o sger.o snrm2.o  \
 srot.o srotg.o srotm.o srotmg.o ssbmv.o sscal.o sspmv.o  \
 sspr.o sspr2.o sswap.o ssymm.o ssymv.o ssyr.o ssyr2.o ssyr2k.o  \
 ssyrk.o stbmv.o stbsv.o stpmv.o stpsv.o strmm.o strmv.o strsm.o  \
 strsv.o xerbla.o zaxpy.o zcopy.o zdotc.o zdotu.o zdrot.o zdscal.o  \
 zgbmv.o zgemm.o zgemv.o zgerc.o zgeru.o zhbmv.o zhemm.o zhemv.o  \
 zher.o zher2.o zher2k.o zherk.o zhpmv.o zhpr.o zhpr2.o zrotg.o zscal.o \
 zswap.o zsymm.o zsyr2k.o zsyrk.o ztbmv.o ztbsv.o ztpmv.o  \
 ztpsv.o ztrmm.o ztrmv.o ztrsm.o ztrsv.o square_solve.o

OTHER_LIBS =

$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)	

%.o:  $(SRC_DIR)/%.f
#	$(FC) $(FLAGS) $(INCLUDES) -c -O2 -w $<
	$(FC) $(FLAGS) -c -O2 -w $<

%.o:  $(SRC_DIR)/%.f90
#	$(FC) $(FREEFLAGS) $(INCLUDES) -c -O2 -w $<
	$(FC) $(FREEFLAGS) -c -O2 -w $<

clean:
	rm -f *.o *.mod $(LIB)
