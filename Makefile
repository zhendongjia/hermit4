.KEEP_STATE:
FFLAGS = -O3
#FFLAGS = -O3 -ffortran-bounds-check
FC = gfortran

SOURCE = \
hermit4.f bhinit.f bhint.f bhlist.f bhpert.f bhreg.f \
bhterm.f block.f bodies.f cmf.f cputim.f data.f \
energy.f fpoly1.f fpert.f iblock.f inext.f input.f \
intgrt.f mydump.f nbint.f output.f ran2.f remove.f \
resolv.f search.f start.f stepi.f stepk.f steps.f \
tstep.f xvpred.f zero.f

OBJECTS = $(SOURCE:.f=.o)

hermit4:  $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o hermit4
 
clean:
	rm -f $(OBJECTS) hermit4

print:
	@- \rm -f HERMIT4.TEXT
	@cat $(SOURCE) > HERMIT4.TEXT
