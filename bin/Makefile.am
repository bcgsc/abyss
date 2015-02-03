if !HAVE_LIBMPI
nompi=ABYSS-P
endif

ABYSS-P: $(top_builddir)/bin/abyss-nompi
	cp $< $@

dist_bin_SCRIPTS = \
	abyss-bowtie \
	abyss-bowtie2 \
	abyss-bwa \
	abyss-bwasw \
	abyss-bwamem \
	abyss-fatoagp \
	abyss-kaligner \
	abyss-pe \
	abyss-samtoafg \
	abyss-tabtomd \
	abyss-bloom-dist.mk \
	abyss-dida
	$(nompi)
dist_noinst_SCRIPTS = \
	abyss-adjtodot.pl \
	abyss-cstont \
	abyss-fac.pl \
	abyss-joindist \
	abyss-nompi
