if HAVE_LIBMPI
Parallel=Parallel
endif

if SAM_SEQ_QUAL
GapFiller=GapFiller
endif

if HAVE_GTEST
Unittest=Unittest
else
check:
	@echo "error: To use 'make check' command, please setup"; \
	echo "       gtest header&library as shown on the Wiki"; \
	false
endif

dist_doc_DATA = ChangeLog COPYRIGHT LICENSE \
	README.css README.html README.md

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
	Layout \
	Map \
	MergePaths \
	Misc \
	Overlap \
	ParseAligns \
	PathOverlap \
	PopBubbles \
	Scaffold \
	SimpleGraph \
	kmerprint \
	FilterGraph \
	$(Unittest) \
	$(GapFiller)

%.html: %.md
	-multimarkdown $< >$@
