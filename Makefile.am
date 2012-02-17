if HAVE_LIBMPI
Parallel=Parallel
endif

dist_doc_DATA = ChangeLog COPYRIGHT LICENSE \
	README README.css README.html

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
	Map \
	MergePaths \
	Overlap \
	ParseAligns \
	PathOverlap \
	PopBubbles \
	Scaffold \
	SimpleGraph \
	kmerprint \
	FilterGraph

README.html: %.html: %
	multimarkdown $< >$@
