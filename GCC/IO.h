/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <vector>
// using namespace std;

class IO
{
private:
public:
	IO(void);
	~IO(void);

	static int file2code(std::string datafile_in, std::string codefile_out);
	static int file2vector(std::string datafile_in, std::vector<std::string>& vector_out);
	static std::vector<std::string> string2vector(
	    std::string line); // split string line to pieces, and store in a vector

	static std::vector<int> formula2index(std::string formula);
	static std::vector<int> formula2index_total(std::string formula);

	static std::vector<int> formula2index_define(std::string formula, std::vector<std::string> Chemical_Element);

	static bool file_compare(std::string, std::string);
	static std::vector<std::string> vector_reduce(std::vector<std::string>);
	static void entrance(void);
};
