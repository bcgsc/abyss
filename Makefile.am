if HAVE_LIBMPI
Parallel=Parallel
endif

if HAVE_GTEST
Unittest=Unittest
else
check:
	@echo "error: To use 'make check' command, please setup"; \
	echo "       gtest header&library as shown on the Wiki"; \
	false
endif

if HAVE_SQLITE
DataBase=DataBase
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
	$(DataBase) \
	ABYSS $(Parallel) \
	AdjList \
	Bloom \
	Konnector \
	Consensus \
	DAssembler \
	DistanceEst \
	KAligner \
	Layout \
	LogKmerCount \
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
	GapFiller \
	$(Unittest)

%.html: %.md
	-multimarkdown $< >$@
