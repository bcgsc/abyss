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
<<<<<<< HEAD
	Sealer \
=======
	GapFiller2 \
>>>>>>> 026d310249f124d1ded1f280ec36c622c52d72cf
	$(Unittest)

%.html: %.md
	-multimarkdown $< >$@
