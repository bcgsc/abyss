#include "config.h"
#include "Uncompress.h"
#include "UnorderedMap.h"
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "FastaReader.h"
#include "Rotation.h"
#include "RotatedRead.h"

using namespace std;

#define PROGRAM "DAssembler"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Rod Docking.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

namespace opt {
    static unsigned max_overlap = 10;
    static unsigned max_mismatch = 2;
    static unsigned min_coverage = 2;
    static unsigned read_length = 50;
    static int verbose = 0;
}

static const char shortopts[] = "o:m:c:r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "max_overlap",            required_argument, NULL, 'o' },
    { "max_mismatch",           required_argument, NULL, 'm' },
    { "min_coverage",           required_argument, NULL, 'c' },
    { "read_length",            required_argument, NULL, 'r' },
    { "verbose",                no_argument,       NULL, 'v' },
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [READS]\n"
"Assemble a single contig from reads in a single orientation.\n"
"\n"
" Arguments:\n"
"\n"
"  READS  fasta-formatted reads file: the first read is used as the seed.\n"
"\n"
" Options:\n"
"\n"
"  -o, --max_overlap=INT            maximum tier overlap for consensus calling"
" [10]\n"
"  -m, --max_mismatch=INT           maximum mismatches allowed for consensus"
" calling [2]\n"
"  -c, --min_coverage=INT           minimum coverage to call a consensus base"
" [2]\n"
"  -r, --read_length=INT            read length\n"
"  -v, --verbose                    display verbose output\n"
"      --help                       display this help and exit\n"
"\n";

/* A small struct for holding overlap information:
    just a sequence and offset from the focal read */
struct Overlap{
    string seq;
    unsigned offset;
    Overlap(const string& seq, unsigned offset)
        : seq(seq), offset(offset) { }
};

/* Additional sort of Overlap objects used for pretty-printing alignments*/
static bool offset_sort(const Overlap& x, const Overlap& y)
{
    return x.offset < y.offset;
}

/* Calculate the tier overlap between two rotated reads */
static int tier_overlap(const string& seq1, const string& seq2,
    bool allow_mismatch = false);

/* From a rotated read, return the original sequence */
static string original_read_from_rotated_read(
        const string& rotated_read)
{
    size_t dollar_pos = rotated_read.find('$');
    string orig_seq = rotated_read.substr(
        dollar_pos+1,rotated_read.size()-dollar_pos) +
        rotated_read.substr(0,dollar_pos);

    return orig_seq;
}

/* Struct for holding base counts at a given position */
struct BaseCount {
    unsigned x[4];
    BaseCount() { fill(x, x + 4, 0); }

    /** Return the number of reads at this position. */
    unsigned sum() const { return accumulate(x, x+4, 0); }

    friend ostream& operator <<(ostream& out, const BaseCount& base)
    {
        out << base.x[0];
        for (int x = 1; x < 4; x++)
            out << '\t' << base.x[x];
        return out;
    }
};

/* Call a consensus base */
static char call_consensus_base(BaseCount counts, char orig_base)
{
    unsigned coverage = accumulate(counts.x, counts.x+4, 0);

    /*If we're below minimum coverage but have an already-called base
          (i.e., for the first few bases)*/
    if (coverage < opt::min_coverage)
        return orig_base;

    // Call a consensus base
    unsigned *maxIt = max_element(counts.x, counts.x+4);

    // Very verbose debugging output
    // char base_to_return = float(*maxIt) <= float(coverage)*0.60 ? orig_base : codeToBase(maxIt - counts.x);
    // char new_base_to_return = *maxIt < opt::min_coverage ? orig_base : codeToBase(maxIt - counts.x);
    // float call_threshold = float(coverage)*0.60;
    // float max_base_percent = float(*maxIt) / float(coverage);
    // cerr << "Coverage: " << coverage << "\tOriginal: " << orig_base <<
    //  "\tReturning: " << base_to_return << "\tBases: " << counts <<
    //  "\tCall threshold: " << call_threshold <<
    //  "\tMax base percent: " << max_base_percent << 
    //  "\tNew base to return: " << new_base_to_return << endl;

    // Original version with hard-coded coverage frequency of 60%
    // return float(*maxIt) <= float(coverage)*0.60 ? orig_base
    //     : codeToBase(maxIt - counts.x);
    
    // Return the most-frequent base, so long as that base has
    //  coverage >= min_coverage
    return *maxIt < opt::min_coverage ? orig_base 
        : codeToBase(maxIt - counts.x);
    
}

/* Return the frequency of the most-common base*/
static float most_common_base_frequency(BaseCount counts)
{
    // Calculate coverage as before
    unsigned coverage = accumulate(counts.x, counts.x+4, 0);
    
    // Call a consensus base
    unsigned *maxIt = max_element(counts.x, counts.x+4);

    // Return the frequency of the most-common base
    return float(*maxIt) / float(coverage);
    
}

typedef vector<Rotation> Rotations;

// Find all overlaps with the given focal read and call a consensus sequence
static string find_complex_overlap(const RotatedRead& f,
        const Rotations& r,
        vector<RotatedRead>& rl)
{
    // A vector for tracking all the overlaps, seeded with the initial read
    vector<Overlap> o;
    o.push_back(Overlap(f.seq, 0));

    // The pre-pended string to use to seed the search
    const string& seq1 = '$' + f.seq;

    //Find it in the sorted list - NOTE: if the flank read doesn't correspond
    // to a real read, this iterator will not be used
    Rotations::const_iterator rt = lower_bound(r.begin(), r.end(), seq1);

    /*
    Continue down the sorted list, checking for other matches
     - for real reads (seq1 == rt->seq), continue from the position
    of seq1 in the list
     - otherwise, just start at the beginning
     */
    for(Rotations::const_iterator st = (seq1 == rt->seq) ? rt+1 : r.begin();
            st != r.end(); ++st)
    {
        // Check for an overlap between the two sequences,
        //  allowing for mismatches
        const string& seq2 = st->seq;
        unsigned new_overlap = tier_overlap(seq1, seq2, true);

        // Continue if there's no match
        if (new_overlap == 0
                || new_overlap > opt::max_overlap)
            continue;

        // Add a new overlap object for each appropriate overlap found
        o.push_back(Overlap(
            original_read_from_rotated_read(seq2), new_overlap));
    }

    // Counters for calculating coverage
    // Vector size should be something like "read_length + maximum tier"
    // THIS WILL BREAK WITH LONGER READS
    vector<BaseCount> counts(300);

    // Pretty-print the alignment for verbose only
    if(opt::verbose){
        cerr << endl;
        sort(o.begin(), o.end(), offset_sort);
    }

    // Go through each overlap to count bases offset by the appropriate amount
    for(vector<Overlap>::const_iterator ot = o.begin();
            ot != o.end(); ++ot){

        // Pretty-print each found read aligned with the focal read
        if (opt::verbose)
            cerr << string(ot->offset, ' ') << ot->seq << " t:" << ot->offset;

        // Retrieve the original RotatedRead object to get the
        //  count for each read
        vector<RotatedRead>::const_iterator rt = lower_bound(
            rl.begin(), rl.end(), ot->seq);
        if(opt::verbose){cerr << " x" << rt->count <<
            " used: " << rt->used << endl;}

        // Continue if we've marked this read as used already
        if(rt->used == true){
            continue;
        }

        // Increment the coverage lists appropriately
        for (size_t i = 0; i < opt::read_length; ++i){
            if ((ot->seq[i] == 'X') || (ot->seq[i] == 'N')){
                continue;
            }
            counts[i+ot->offset].x[baseToCode(ot->seq[i])]
                += rt->count;
        }
    }

    // Call consensus bases until we run out of coverage
    ostringstream new_contig;
    char new_base = '*';
    float current_consensus_freq = 1.0;
    float next_consensus_freq = 1.0;
    for (unsigned i = 0; new_base != 'X'; i++) {
        
        // Retrieve the original base, or 'X' if we're past the end
        //  of the original flank
        char orig_base = i < opt::read_length ? f.seq[i] : 'X';

        // Call a new consensus base if possible
        new_base = call_consensus_base(counts[i], orig_base);

        // Check the frequency of the most-common base
        current_consensus_freq = most_common_base_frequency(counts[i]);
        next_consensus_freq = most_common_base_frequency(counts[i+1]);
        //cerr << "Current: " << current_consensus_freq << " Next: " << next_consensus_freq << endl;
        
        // Bail out if we encounter two SNPs in a row
        // Set the current base to 'X' and trim the last one
        if ((current_consensus_freq <= 0.8) && (next_consensus_freq <= 0.8))
            new_base = 'X';
        
        // If we've found a new base, add it to the growing consensus
        if (new_base != 'X') 
            new_contig << new_base;
        
    }
    if (opt::verbose)
        cerr << new_contig.str() << " (consensus) " << endl;

    // Mark reads that shouldn't be used again
    unsigned growth = new_contig.str().size() - opt::read_length;
    for(vector<Overlap>::const_iterator ot = o.begin();
            ot != o.end(); ++ot){
        // Reads are used if they don't extend to the end of the consensus
        if (ot->offset <= (growth-1)){
            //Find the correct RotatedRead object and mark it as used
            vector<RotatedRead>::iterator rt = lower_bound(
                rl.begin(), rl.end(), ot->seq);
            if(rt->seq == ot->seq)
                rt->used = true;
        }
    }

    /*The sequence returned here contains the original focal
     read plus any extension
    The main routine is responsible for trimming back the growing contig */
    return new_contig.str();
}

// Calculate the tier overlap between two reads
static int tier_overlap(const string& seq1, const string& seq2,
        bool allow_mismatch)
{
    assert(seq1 != seq2);

    //Find the position of the '$' character in both reads
    unsigned first_dollar_pos = seq1.find('$');
    unsigned second_dollar_pos = seq2.find('$');
    unsigned earliest_dollar_pos = first_dollar_pos <= second_dollar_pos ?
        first_dollar_pos : second_dollar_pos;
    unsigned latest_dollar_pos = first_dollar_pos > second_dollar_pos ?
        first_dollar_pos : second_dollar_pos;

    //If the two strings are equal outside the dollar signs,
    //  return the tier - this is a no-mismatch overlap
    if( (seq1.substr(0, earliest_dollar_pos) ==
         seq2.substr(0, earliest_dollar_pos)) &&
        (seq1.substr(latest_dollar_pos+1,
        (opt::read_length+1)-latest_dollar_pos+1) ==
        seq2.substr(latest_dollar_pos+1,
        (opt::read_length+1)-latest_dollar_pos+1)) ){
           return latest_dollar_pos - earliest_dollar_pos;
    }

    //Otherwise, if mismatches are allowed, calculate that overlap
    if (allow_mismatch){
        unsigned num_mismatch = 0;

        for(unsigned i = 0; i < (opt::read_length+1); ++i)
        {
            if ((i >= earliest_dollar_pos) && (i <= latest_dollar_pos)){
                continue;
            }else if (seq1[i] != seq2[i]){
                num_mismatch++;
            }
        }

        //NOTE: this is also checking that the second read is
        //  DOWNSTREAM from the first
        if ((num_mismatch <= opt::max_mismatch) &&
            (second_dollar_pos > first_dollar_pos)){
            return latest_dollar_pos - earliest_dollar_pos;
        }
    }

    //Otherwise...
    // NOTE: we're currently defining both "no overlap"
    //  and "no offset but mismatches" as "0"
    //  ==> This could probably be changed
    return 0;
}

// From a sorted list of RotatedRead objects, generate a sorted
//  vector of Rotation objects
// This function generates the main vector we traverse when looking
//  for simple extensions
static Rotations generate_rotations_list(
    const vector<RotatedRead>& read_list)
{
    // Each rotation has a sequence (string) and overlap (int)
    // The overlap refers to the overlap between the current sequence
    //  and the next sequence in the list
    Rotations s;

    // For each distinct read, add all the rotations to the overlap list
    for (vector<RotatedRead>::const_iterator it = read_list.begin();
            it != read_list.end(); ++it)
        for(vector<string>::size_type i = 0; i != (*it).rotations.size(); ++i)
            // Initialize object with sequence
            s.push_back(Rotation((*it).rotations[i]));

    // Sort the list of rotated reads
    sort(s.begin(), s.end());

    // Find the 0-mismatch tier overlap between each pair of reads
    //  in the sorted list
    for (Rotations::iterator rt = s.begin(); rt != s.end()-1; ++rt)
        rt->overlap = tier_overlap(rt->seq, rt[1].seq);

    // Define last entry in the list as 0
    s[s.size()-1].overlap = 0;

    return s;
}

// Main control flow function
int main(int argc, char** argv)
{
    // Parse command-line options
    bool die = false;
    for (int c; (c = getopt_long(argc, argv,
                    shortopts, longopts, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?': die = true; break;
            case 'o': arg >> opt::max_overlap; break;
            case 'm': arg >> opt::max_mismatch; break;
            case 'c': arg >> opt::min_coverage; break;
            case 'r': arg >> opt::read_length; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                cout << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                cout << VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) {
        cerr << PROGRAM ": missing arguments\n";
        die = true;
    } else if (argc - optind > 1) {
        cerr << PROGRAM ": too many arguments\n";
        die = true;
    } else if (opt::max_overlap > (opt::read_length-1)){
        cerr << PROGRAM ": max_overlap cannot be larger than (read_length-1)\n";
        die = true;
    }

    if (die) {
        cerr << "Try `" << PROGRAM
            << " --help' for more information.\n";
        exit(1);
    }

    const char* fasta_file = argv[optind++];

    if(opt::verbose){
        cerr << PROGRAM <<
            "\n  max_overlap:         " << opt::max_overlap <<
            "\n  max_mismatch:        " << opt::max_mismatch <<
            "\n  min_coverage:        " << opt::min_coverage <<
            "\n  read_length:         " << opt::read_length << endl;
    }

    // Use the ABySS FastaReader class to read in a fasta file of reads
    // Assume the first read in the file is the seed for the assembly
    if(opt::verbose){cerr << "Reading `" << fasta_file << "'...  ";}

    typedef unordered_map<string, unsigned> ReadMap;
    ReadMap read_map;

    bool first_read = true;
    string contig;

    FastaReader in(fasta_file,
         FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {

        string read_seq = rec.seq;

        if (first_read){
            contig = read_seq;
            first_read = false;
        }
        // Count the reads as we collect them here...
        read_map[read_seq]++;
    }

    vector<RotatedRead> read_list;
    // ... Then put them into the vector of RotatedRead objects
    for (ReadMap::iterator i = read_map.begin();
            i != read_map.end(); ++i)
        read_list.push_back(RotatedRead(i->first, i->second));
    read_map.clear();
    sort(read_list.begin(), read_list.end());

    if(opt::verbose){cerr << "finished reading fasta file with " <<
        read_list.size() << " distinct reads.\n\n" << endl;}

    // Generate the sorted lists
    Rotations rotation_list = generate_rotations_list(read_list);

    // Main assembly loop
    if(opt::verbose){cerr << contig << " (seed)" << endl;}
    bool time_to_die = false;
    int hard_cap = 0;

    while(! time_to_die){

        // Temporary hard-cap to prevent runaway execution
        hard_cap++;
        if (hard_cap >= 500){
            time_to_die = true;
            //cerr << "Hard cap hit - I give up!" << endl;
            //cerr << contig.size();
            // cout << ">contig (" << contig.size() << "bp)" << endl
            //      << contig << endl;
            //cout << contig.size();
            cout << contig << endl;
            exit(1);
        }

        // Another break if the contig grows too long
        if (contig.size() >= 1500){
            //cerr << contig.size();
            cout << contig << endl;
            exit(1);
        }

        // Extract the flanking sequence and fetch rotations of that read
        string flank = contig.substr(contig.size()-opt::read_length);

        // Retrieve the flanking read
        vector<RotatedRead>::iterator low = lower_bound(
            read_list.begin(), read_list.end(), flank);

        //TODO - This search sometimes causes a segfault at higher -o values
        // Figure out why!
        RotatedRead flank_read = (*low);

        // If the flank sequence doesn't correspond to a real read:
        //  - generate a temporary RotatedRead object for the complex search
        if (flank != flank_read.seq){
            if (opt::verbose) cerr <<
                "Flank doesn't correspond to a real read" << endl;
            flank_read = RotatedRead(flank, 1);
        }

        bool found_complex_overlap = false;
        string read_to_add;

        string extension_seq = find_complex_overlap(
            flank_read, rotation_list, read_list);

        if (! (extension_seq == flank_read.seq)){
            found_complex_overlap = true;
            // The new contig = old contig - flank read + extension sequence
            // (the extension sequence contains the flank read)
            contig = contig.substr(0, contig.size()-opt::read_length) +
             extension_seq;

            // This is very verbose - prints out a fasta sequence for each
            //  step of the assembly process
            if(opt::verbose){
                cout << ">p" << opt::max_overlap << "_" <<
            contig.size() << "bp_complex\n" << contig << endl;
            }
        }

        // If the search fails, stop extension
        if (! found_complex_overlap){
            time_to_die = true;
        }
    }

    // Output the final contig
    cout << contig << endl;    
    
    return 0;
}
