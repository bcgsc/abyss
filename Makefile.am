if HAVE_LIBMPI
Parallel=Parallel
endif

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
dist_pkgdoc_DATA = ChangeLog LICENSE README

SUBDIRS = Align Assembly Common DataLayer ExpPaired Graph \
	ABYSS $(Parallel) \
	bin doc fasta2psq kmerprint \
	AdjList DistanceEst KAligner MergePaths Overlap ParseAligns \
	SimpleGraph Consensus PathOverlap
