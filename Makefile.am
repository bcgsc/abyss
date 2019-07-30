if HAVE_LIBMPI
Parallel=Parallel
endif

if HAVE_PTHREAD
GTest=lib/gtest-1.7.0
UnitTest=Unittest
endif

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib CITATION.md \
	COPYRIGHT \
	LICENSE \
	README.md

if HAVE_PANDOC
dist_doc_DATA += README.html
endif

EXTRA_DIST=autogen.sh doxygen.conf

SUBDIRS = \
	bin \
	doc \
	Common \
	Graph \
	DataLayer \
	DataBase \
	FMIndex \
	Assembly \
	dialign \
	Align \
	ABYSS $(Parallel) \
	Bloom \
	BloomDBG \
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
	lib \
	vendor/nthash \
	$(GTest) \
	$(UnitTest)

%.html: $(srcdir)/%.md
	-pandoc -s -o $@ $<

clean-local:
	rm -f README.html

clang-format:
	for i in Bloom/RollingBloomDBGVisitor.h Bloom/bloom.cc  BloomDBG/BloomIO.h \
	BloomDBG/Checkpoint.h  BloomDBG/HashAgnosticCascadingBloom.h BloomDBG/bloom-dbg.* \
	Unittest/BloomDBG/HashAgnosticCascadingBloomTest.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in Bloom/RollingBloomDBGVisitor.h Bloom/bloom.cc  BloomDBG/BloomIO.h \
	BloomDBG/Checkpoint.h  BloomDBG/HashAgnosticCascadingBloom.h BloomDBG/bloom-dbg.* \
	Unittest/BloomDBG/HashAgnosticCascadingBloomTest.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
