/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "mHMPreprocessor.h"

#include "Base/FileTools.h"

void DisplayOption()
{
    std::string opt =
        "Options:\n"
        "    --version:  display version number.\n"
        "    --help:     display help info.\n"
        "    --output-directory: set output directory.\n";
    std::cout << opt << std::endl;
}
void DisplayMessage()
{
    std::string info =
        "Pre-process tool for OGS5 using mHM recharge data.\n "
        "Command: mHM2OGS [file name] [option]\n";
    std::cout << info << std::endl;
}

void DisplayMessageConsole()
{
    DisplayMessage();

    std::cout << "Run in console.\nInput file name or with option:"
              << std::endl;
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
                std::cerr << "Error: Parameter " << anArg
                          << " needs a path for output files" << std::endl;
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

    //
    const std::string file_path = pathDirname(file_name);
    if (o_path.empty())
        o_path = file_path;

    std::string fname = file_name + ".msh";
    std::ifstream is_mesh(fname.c_str(), std::ios::in);

    std::string s_buff;
    std::getline(is_mesh, s_buff);

    MeshLib::mHMPreprocessor* mHM_preprocessor = NULL;
    if (s_buff.find("#FEM_MSH") != std::string::npos)
    {
        mHM_preprocessor = new MeshLib::mHMPreprocessor(&file_name);
        mHM_preprocessor->Read(&is_mesh);
    }
    else
    {
        std::cout << "Cannot open mesh file " << fname << std::endl;
        return EXIT_FAILURE;
    }

    mHM_preprocessor->transform_mHMData(o_path);

    delete mHM_preprocessor;

    std::cout << "Terminate normally ^O^ ." << std::endl;
    return EXIT_SUCCESS;
}
