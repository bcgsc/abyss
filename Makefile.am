if HAVE_LIBMPI
Parallel=Parallel
endif

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
dist_pkgdoc_DATA = ChangeLog COPYRIGHT LICENSE README

EXTRA_DIST=doxygen.conf

SUBDIRS = \
	bin \
	doc \
	Assembly \
	Common \
	DataLayer \
	FMIndex \
	Graph \
	dialign \
	Align \
	ABYSS $(Parallel) \
	AdjList \
	Consensus \
	DAssembler \
	DistanceEst \
	KAligner \
	MergePaths \
	Overlap \
	ParseAligns \
	PathOverlap \
	PopBubbles \
	Scaffold \
	SimpleGraph \
	kmerprint
