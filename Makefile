.KEEP_STATE:
#FFLAGS = -pg -g -ffpe-trap=invalid,zero,overflow -ffixed-line-length-none
FFLAGS = -ffixed-line-length-none -O3 -march=sandybridge -static


FC = gfortran

SOURCE = \
hermit4.f bhinit.f bhint.f bhlist.f bhpert.f bhreg.f \
bhterm.f block.f bodies.f cmf.f cputim.f data.f \
energy.f fpoly1.f fpert.f iblock.f inext.f input.f \
intgrt.f mydump.f nbint.f output.f ran2.f remove.f \
resolv.f search.f start.f stepi.f stepk.f steps.f \
tstep.f xvpred.f zero.f get_precession.f damping.f \
gas_potential.f gr.f rm_from_list.f \
gas_potential_add_inner_pot.f get_outer_gravity.f \
gas_potential_3d.f gas_potential_thommes08.f \
gas_potential_full.f \
gas_potential_red.f \
gas_potential_sg.f \
damping_a.f \
gas_potential_j.f gas_potential_p.f gas_potential_s.f \
load_data.f


OBJECTS = $(SOURCE:.f=.o)

hermit4:  $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o hermit4

clean:
	rm -f $(OBJECTS) hermit4

print:
	@- \rm -f HERMIT4.TEXT
	@cat $(SOURCE) > HERMIT4.TEXT
