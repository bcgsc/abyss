/**
 * Place each contig on a one-dimesional coordinate system using
 * distance estimates.
 * Written by Shaun Jackman.
 */

/** Disable Boost uBLAS runtime sanity checks. */
#define BOOST_UBLAS_NDEBUG 1

#include "ContigProperties.h"
#include "Estimate.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "cholesky.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <string>

using namespace std;
using boost::tie;
namespace ublas = boost::numeric::ublas;

#define PROGRAM "abyss-drawgraph"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FASTA|OVERLAP DIST...\n"
"Place each contig on a one-dimensional coordinate system using\n"
"distance estimates and output a DOT graph with coordinates.\n"
"\n"
" Arguments:\n"
"\n"
"  FASTA    contigs in FASTA format\n"
"  OVERLAP  the contig overlap graph\n"
"  DIST     estimates of the distance between contigs\n"
"\n"
" Options:\n"
"\n"
"      --l2=REAL         set the L2 regularizer to REAL [1e-323]\n"
"  -x, --xscale=N        set the x scale to N nt/inch [100e3]\n"
"      --gv, --dot       output in GraphViz DOT format [default]\n"
"      --tsv             output in TSV format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** The L2 regularizer. */
	double l2 = 1e-323;

	/** The x scale. */
	double xscale = 100e3; // nt/inch

	/** Verbose output. */
	int verbose; // used by PopBubbles

	/** Output format */
	int format = DOT; // used by DistanceEst
}

static const char shortopts[] = "x:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_L2 };

static const struct option longopts[] = {
	{ "dot", no_argument, &opt::format, DOT },
	{ "gv", no_argument, &opt::format, DOT },
	{ "tsv", no_argument, &opt::format, TSV },
	{ "l2", required_argument, NULL, OPT_L2 },
	{ "xscale", required_argument, NULL, 'x' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A distance estimate graph. */
typedef DirectedGraph<Length, DistanceEst> DG;
typedef ContigGraph<DG> Graph;

/** A matrix. */
typedef ublas::matrix<double> Matrix;

/** A vector. */
typedef ublas::vector<double> Vector;

/** Read a graph from the specified file. */
static void readGraph(const string& path, Graph& g)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";
	ifstream fin(path.c_str());
	istream& in = path == "-" ? cin : fin;
	assert_good(in, path);
	read_graph(in, g, BetterDistanceEst());
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);
	g_contigNames.lock();
}

/** Solve Ax = b for x using Cholesky decomposition.
 * The matrix A is symmetric and positive-definite.
 * @param[in,out] a input matrix A, which is clobbered
 * @param[in,out] b input vector b and output vector x
 */
static void solve(Matrix& a, Vector& b)
{
	int ret = cholesky_decompose(a);
	if (ret > 0) {
		cerr << PROGRAM ": error: The graph matrix is singular. "
			"It may have multiple connected components. "
			"Try increasing the L2 regularizer, for example --l2=1e-300\n";
		exit(EXIT_FAILURE);
	}
	cholesky_solve(a, b, ublas::lower());
}

/** Run abyss-drawgraph. */
int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true;
			break;
		  case OPT_L2:
			arg >> opt::l2;
			break;
		  case 'x':
			arg >> opt::xscale;
			break;
		  case 'v':
			opt::verbose++;
			break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (argc - optind < 0) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::edge_iterator Eit;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::vertex_iterator Vit;

	Graph g;
	if (optind < argc) {
		for (; optind < argc; optind++)
			readGraph(argv[optind], g);
	} else
		readGraph("-", g);

	// Add any missing complementary edges.
	addComplementaryEdges(g);

	size_t n = num_vertices(g);
	Matrix a(n, n);
	Vector b = ublas::zero_vector<double>(n);

	// Set the origin of the layout. Pin the first contig and its reverse
	// complement to an extreme negative and positive position, so that the two
	// components should not overlap. The two vertices should be in separate
	// components. If the graph has multiple components, the matrix should be
	// regularized by adding lambda to all values on the diagonal. The
	// components will overlap be can be separated afterward.
	const double MAX_SCAFFOLD_SIZE = 10e6;
	a(0, 0) = 1;
	b[0] = -MAX_SCAFFOLD_SIZE;
	a(1, 1) = 1;
	b[1] = MAX_SCAFFOLD_SIZE;

	// Add the L2 regularizer to the matrix.
	if (opt::l2 > 0)
		for (unsigned i = 0; i < n; ++i)
			a(i, i) += opt::l2;

	// Build the information matrix.
	Eit eit, elast;
	for (tie(eit, elast) = edges(g); eit != elast; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		size_t ui = get(vertex_index, g, u);
		size_t vi = get(vertex_index, g, v);
		int d = get(edge_weight, g, e);
		double err = max(0.1, (double)g[e].stdDev);
		double weight = 1 / err;
		a(ui, ui) += weight;
		a(ui, vi) -= weight;
		a(vi, ui) -= weight;
		a(vi, vi) += weight;
		b[ui] -= d * weight;
		b[vi] += d * weight;
	}

	// Solve the equation Ax = b for x.
	solve(a, b);

	// Determine the origin and width of the layout.
	double min_pos = std::numeric_limits<double>::max();
	double max_pos = -std::numeric_limits<double>::max();
	double min_rc = std::numeric_limits<double>::max();
	for (unsigned i = 0; i < n; ++i) {
		if (b[i] == 0) {
			// A disconnected vertex.
		} else if (b[i] < 0) {
			min_pos = min(min_pos, b[i] - g[vertex(i, g)].length);
			max_pos = max(max_pos, b[i]);
		} else
			min_rc = min(min_rc, b[i] - g[vertex(i, g)].length);
	}
	assert(min_pos != std::numeric_limits<double>::max());
	assert(max_pos != -std::numeric_limits<double>::max());
	assert(min_rc != std::numeric_limits<double>::max());

	// Set the origin of the positive and negative components.
	for (unsigned i = 0; i < n; ++i) {
		if (b[i] == 0)
			// A disconnected vertex.
			b[i] = NAN;
		else if (b[i] < 0)
			b[i] = b[i] - min_pos;
		else
			b[i] = b[i] - min_rc + (max_pos - min_pos);
	}

	// Output the coordinates of each contig.
	if (opt::format == TSV) {
		std::cout << "Name\tSize\tLeftPos\tRightPos\n";
		Vit uit, ulast;
		for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit) {
			V u = *uit;
			size_t ui = get(vertex_index, g, u);
			std::cout << get(vertex_name, g, u) << '\t' << g[u].length;
			if (std::isnan(b[ui])) {
				// A disconnected vertex.
				std::cout << "\tNA\tNA\n";
			} else {
				ssize_t x1 = (ssize_t)round(b[ui]);
				ssize_t x0 = x1 - g[u].length;
				std::cout << '\t' << x0 << '\t' << x1 << '\n';
			}
		}
		exit(EXIT_SUCCESS);
	}
	assert(opt::format == DOT);

	// Sort the contigs by their right coordinate.
	std::vector< std::pair<double, V> > sorted;
	sorted.reserve(n);
	Vit uit, ulast;
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit) {
		V u = *uit;
		size_t ui = get(vertex_index, g, u);
		double x1 = std::isnan(b[ui]) ? 0 : b[ui];
		sorted.push_back(std::make_pair(x1, u));
	}
	sort(sorted.begin(), sorted.end());

	// Write the graph.
	cout << "digraph g {\n"
		"node [shape=\"box\" height=0.3 fixedsize=1]\n";
	assert_good(cout, "stdout");

	// Write the vertices.
	double pos0 = -numeric_limits<double>::max();
	size_t yi = 0;
	for (size_t i = 0; i < n; ++i) {
		V u = sorted[i].second;
		double pos = sorted[i].first;
		size_t l = g[u].length;

		if (pos - l < pos0)
			yi++;
		pos0 = pos;

		double dpi = 72;
		const double W = dpi / opt::xscale; // pt/nt
		const double H = 32; // pt
		double w = l * W;
		double x = pos * W - w / 2;
		double y = yi * H;
		cout << '"' << get(vertex_name, g, u) << "\""
			" [pos=\"" << x << ",-" << y << "\""
			" width=" << w / dpi << "]\n";
	}

	// Write the edges.
	for (tie(eit, elast) = edges(g); eit != elast; ++eit)
		cout << get(edge_name, g, *eit) << '\n';

	cout << "}\n";
	assert_good(cout, "stdout");

	return 0;
}
