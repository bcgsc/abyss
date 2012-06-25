if HAVE_GHC
bin_PROGRAMS = abyss-samtobreak
endif

abyss_samtobreak_SOURCES = samtobreak.hs

CLEANFILES=*.hi

GHC_FLAGS=-O2

abyss-samtobreak: $(abyss_samtobreak_SOURCES)
	$(GHC) $(GHC_FLAGS) -hidir $(dir $@) -odir $(dir $@) -o $@ $^
