/*    sortmapped: a program for sorting mapped read files
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "MappedReadsFile.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

typedef vector<GenomicRegion>::const_iterator GenomicRegionPointer;

struct region_pointer_less {
  bool operator()(const GenomicRegionPointer a, 
		  const GenomicRegionPointer b) const {
    return (*a) < (*b);
  }
};

struct region_pointer_name_less {
  bool operator()(const GenomicRegionPointer a, 
		  const GenomicRegionPointer b) const {
    return a->get_name() < b->get_name();
  }
};

void
sort_regions(const bool SORT_ON_NAME,
	     const vector<GenomicRegion> &regions, 
	     const vector<string> &sequences, const vector<string> &scores, 
	     const string outfile) {
  vector<GenomicRegionPointer> sorter;
  for (vector<GenomicRegion>::const_iterator i = regions.begin(); 
       i != regions.end(); ++i) sorter.push_back(i);
  
  if (SORT_ON_NAME)
    sort(sorter.begin(), sorter.end(), region_pointer_name_less());
  else 
    sort(sorter.begin(), sorter.end(), region_pointer_less());
  
  ostream* out = (outfile.empty()) ? 
    &std::cout : new std::ofstream(outfile.c_str());
  for (vector<GenomicRegionPointer>::const_iterator i(sorter.begin());
       i != sorter.end(); ++i) {
    const size_t j = std::distance(regions.begin(), *i);
    *out << *(*i) << '\t' << sequences[j] << '\t' << scores[j] << '\n';
  }
  if (out != &std::cout) delete out;
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool SORT_ON_NAME = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("sortmapped", "A program for sorting mapped read files",
			   "<mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("name", 'N', "Sort by the names of reads", 
		      false , SORT_ON_NAME);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/
    
    vector<GenomicRegion> regions;
    vector<string> sequences;
    vector<string> scores;
    LoadMappedReadsFile(input_file_name, regions, sequences, scores);
    
    if (!check_sorted(regions) || SORT_ON_NAME)
      sort_regions(SORT_ON_NAME, regions, sequences, scores, outfile);
    else {
      ostream* out = (outfile.empty()) ? 
	&std::cout : new std::ofstream(outfile.c_str());
      for (size_t i = 0; i < regions.size(); ++i)
	*out << regions[i] << '\t' << sequences[i] << '\t' << scores[i] << '\n';
      if (out != &std::cout) delete out;
    }
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
