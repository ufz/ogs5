/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "misc.h"
#include "msh_mesh.h"

void DisplayOption()
{
	std::string opt
	    = "Options:\n"
	      "    --version:  display version number.\n"
	      "    --help:     display help info.\n"
	      "    --output-directory: set output directory.\n";
	std::cout << opt << std::endl;
}
void DisplayMessage()
{
	std::string info
	    = "Pre-process tool for OGS using mHM recharge data.\n "
	      "Command: mHM2OGS [file name] [option]\n";
	std::cout << info << std::endl;
}

void DisplayMessageConsole()
{
	DisplayMessage();

	std::cout << "Run in console.\nInput file name or with option:" << std::endl;
}
void DisplayVersion()
{
	std::string ver = "Version: 1.0.";
	std::cout << ver << std::endl;
}

int main(int argc, char* argv[])
{
	std::vector<std::string> arg_strings;
	if (argc > 1)
	{
		for (int i = 1; i < argc; i++)
		{
			arg_strings.push_back(std::string(argv[i]));
		}
	}
	else
	{
		DisplayMessageConsole();

		std::string s_buff;
		std::stringstream ss;
		getline(std::cin, s_buff);
		ss.str(s_buff);
		while (!ss.eof())
		{
			ss >> s_buff;
			arg_strings.push_back(s_buff);
		}
	}

	std::string file_name;
	std::string o_path;
	for (std::size_t i = 0; i < arg_strings.size(); i++)
	{
		const std::string anArg = arg_strings[i];
		if (anArg == "--help" || anArg == "-h")
		{
			DisplayMessage();
			DisplayOption();
			exit(0);
		}
		if (anArg == "--version")
		{
			DisplayVersion();
			exit(0);
		}
		if (anArg == "--output-directory")
		{
			if (i + 1 >= arg_strings.size())
			{
				std::cerr << "Error: Parameter " << anArg << " needs a path for output files" << std::endl;
				std::exit(EXIT_FAILURE);
			}
			std::string path = arg_strings[++i];

			if (!path.empty())
				o_path = path;
			continue;
		}
		else
		{
			file_name = arg_strings[i];
		}
	}

	if (argc > 1)
	{
		DisplayMessage();
	}

	std::basic_string<char>::size_type indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = file_name.find_last_of('\\');
	indexChLinux = file_name.find_last_of('/');
	//
	std::string file_path;
	if (indexChWin != std::string::npos)
		file_path = file_name.substr(0, indexChWin) + "\\";
	else if (indexChLinux != std::string::npos)
		file_path = file_name.substr(0, indexChLinux) + "/";
	if (o_path.empty())
		o_path = file_path;

	MeshLib::CFEMesh* a_mesh = NULL;

	std::string fname = file_name + ".msh";
	std::ifstream is_mesh(fname.c_str(), std::ios::in);

	std::string s_buff;
	std::getline(is_mesh, s_buff);

	if (s_buff.find("#FEM_MSH") != std::string::npos)
	{
		a_mesh = new MeshLib::CFEMesh(NULL, &file_name);
		a_mesh->Read(&is_mesh);
	}
	else
	{
		std::cout << "Cannot open mesh file " << fname << std::endl;
		return EXIT_FAILURE;
	}

	a_mesh->mHM2NeumannBC();
	if (a_mesh)
	{
		delete a_mesh;
		a_mesh = NULL;
	}
	return EXIT_SUCCESS;
}
