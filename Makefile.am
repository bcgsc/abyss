if HAVE_LIBMPI
Parallel=Parallel
endif

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
dist_pkgdoc_DATA = ChangeLog LICENSE README

SUBDIRS = Align Assembly Common DataLayer ExpPaired Graph \
	ABYSS $(Parallel) \
	bin doc kmerprint \
	AdjList DistanceEst KAligner MergePaths Overlap \
	ParseAligns PopBubbles SimpleGraph Consensus PathOverlap
