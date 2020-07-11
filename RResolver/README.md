RResolver
===================

RResolver, which stands for repeat or read resolver, whichever you prefer, improves genome assemblies at unitig stage by using a sliding window at read size level across repeats to determine which paths are correct.

Synopsis
===================

`abyss-rresolver-short " [OPTION]... <contigs> <graph> [<reads1> <reads2> ...]`

Options
===================

 * `b`: read Bloom filter size.
        Unit suffixes 'K' (kilobytes), 'M' (megabytes),
        or 'G' (gigabytes) may be used. [`required`]
 * `k`: assembly k-mer size [`required`]
 * `g`: write the contig adjacency graph to FILE [`required`]
 * `c`: write the contigs to FILE [`required`]
 * `j`: use N parallel threads [`1`]
 * `h`: write the algorithm histograms with the given prefix.
        Histograms are omitted if no prefix is given.
 * `t`: set path support threshold to N. [`5`]
 * `m`: set minimum number of sliding window moves to N. [`20`]
 * `r`: set maximum number of branching paths to N. [`75`]
 * `S`: write supported paths to FILE.
 * `U`: write unsupported paths to FILE.

Authors
===================
+ [**Vladimir Nikolic**](https://github.com/schutzekatze)
+ Supervised by [**Dr. Inanc Birol**](https://www.bcgsc.ca/people/inanc-birol).

Copyright 2020 Canada's Michael Smith Genome Sciences Centre