if HAVE_LIBMPI
Parallel=Parallel
endif

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
dist_pkgdoc_DATA = ChangeLog LICENSE README

SUBDIRS = Align Common DataLayer ExpPaired Graph ABYSS $(Parallel) \
	bin doc fasta2psq \
	AdjList DistanceEst KAligner MergePaths Overlap ParseAligns \
	SimpleGraph
