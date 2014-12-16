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

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib CITATION.md \
	COPYRIGHT \
	LICENSE \
	README.css README.html README.md

EXTRA_DIST=doxygen.conf

SUBDIRS = \
	bin \
	doc \
	Common \
	Graph \
	DataLayer \
	$(DataBase) \
	FMIndex \
	Assembly \
	dialign \
	Align \
	ABYSS $(Parallel) \
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
	PairedDBG \
	ParseAligns \
	PathOverlap \
	PopBubbles \
	Scaffold \
	SimpleGraph \
	kmerprint \
	FilterGraph \
	GapFiller \
	Sealer \
	AdjList \
	$(Unittest)

%.html: %.md
	-multimarkdown $< >$@

clean-local:
	rm -f README.html
