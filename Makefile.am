if HAVE_LIBMPI
Parallel=Parallel
endif

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
pkgdoc_DATA = ChangeLog LICENSE README

SUBDIRS = Align Common DataLayer ExpPaired Graph ABYSS $(Parallel) \
	bin doc fasta2psq \
	AdjList DistanceEst SimpleGraph MergePaths KAligner ParseAligns
