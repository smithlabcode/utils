/*    sigoverlap: 
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

using std::ofstream;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;

using std::unordered_map;


template <class T> 
size_t
genomic_region_intersection(const std::vector<T>& regions_a, 
			    const std::vector<T>& regions_b) {
  typename std::vector<T>::const_iterator a(regions_a.begin());
  typename std::vector<T>::const_iterator a_lim(regions_a.end());
  typename std::vector<T>::const_iterator b(regions_b.begin());
  typename std::vector<T>::const_iterator b_lim(regions_b.end());
  size_t c = 0;
  while (a != a_lim && b != b_lim) {
    if (a->overlaps(*b))
      ++c;
    if (a == b) {++a; ++b;}
    else if (*a < *b) ++a;
    else ++b; //  if (*b < *a)
  }
  return c;
}


int
main(int argc, const char **argv) {

  /* FILES */
  string target_regions_file;
  string outfile;
  
  bool VERBOSE = false;

  double genome_size = 0.0;

  /****************** GET COMMAND LINE ARGUMENTS ***************************/
  OptionParser opt_parse(argv[0], "",
			 "<bed-format-file>");
  opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		    false , outfile);
  opt_parse.add_opt("verbose", 'v', "print more run info", 
		    false , VERBOSE);
  opt_parse.add_opt("target", 't', "target regions file", 
		    true, target_regions_file);
  opt_parse.add_opt("size", 's', "genome size", 
		    true, genome_size);
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
  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string input_file_name = leftover_args.back();
  /**********************************************************************/
  
  try {
    
    vector<GenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);
    if (!check_sorted(regions)) {
      cerr << "ERROR: regions in \"" << input_file_name
	   << "\" not sorted" << endl;
      return EXIT_FAILURE;
    }
    
    vector<GenomicRegion> other_regions;
    ReadBEDFile(target_regions_file, other_regions);
    if (!check_sorted(other_regions)) {
      cerr << "ERROR: regions in \"" << target_regions_file
	   << "\" not sorted" << endl;
      return EXIT_FAILURE;
    }
    double other_regions_size = 0.0;
    for (size_t i = 0; i < other_regions.size(); ++i)
      other_regions_size += other_regions[i].get_width();
    
    
    double expected = 0.0;
    for (size_t i = 0; i < regions.size(); ++i) {
      const double p = (other_regions_size + regions[i].get_width()*(other_regions.size() - 1))/genome_size;
      genome_size -= regions[i].get_width();
      other_regions_size -= p*regions[i].get_width();
      expected += p;
    }

    typedef vector<GenomicRegion>::const_iterator region_itr;
    size_t overlap = 0;
    for (size_t i = 0; i < regions.size(); ++i) {
      region_itr closest(find_closest(other_regions, regions[i]));
      if (closest->overlaps(regions[i]))
	++overlap;
    }
    // size_t c = genomic_region_intersection(other_regions, regions);
    
    cout << other_regions.size() << "\t"
	 << regions.size() << "\t" 
	 << overlap << "\t" 
	 << expected << "\t"
	 << overlap/expected << endl;
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
