/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// expre_new_Operator.cpp
// compile with: /EHsc
#include <cmath>

#include <fstream>

#include "misc.h"
#include "msh_mesh.h"

using namespace std;
using namespace MeshLib;

void DisplayMessage()
{
   string _infor =
      "\nData pre-process tool for the finite element analysis of subsurface flow problems.\n "
      "By Wenqing Wang\n\n"
      "Command: fem_toolkits -opt [Option] [File name with extension]\n"
      "Option:\n"
      "  1:  Generation OGS Neumman BC from the recharge data in GIS raster file,\n"
      "  2:  Top surface integration for 3D mesh,\n"
      "  3:  Convert GIS raster cells into FE mesh\n";
   cout<<_infor<<endl;
}

int main(int argc, char* argv[])
{
   string s_buff;
   string opt;
   stringstream ss;
   string file_name;
   if(argc>1)
   {
      for(int i=1; i<argc; i++)
      {
         s_buff = argv[i];
         if(s_buff.find("--version")!=string::npos)
         {
            DisplayMessage();
            exit(0);
         }

         if( s_buff.find("-opt") != string::npos)
         {
            opt = argv[i+1];
         }

         if(   !(s_buff.find("-")!=string::npos || s_buff.find("--")!=string::npos) )
         {
            file_name = s_buff;
         }

      }
   }
   else //terminal
   {
      DisplayMessage();
      cout<<"\nInput options and file name (non extension):\n ";

      getline(cin, s_buff);
      ss.str(s_buff);
      while(!ss.eof())
      {
         ss>>s_buff;
         if(s_buff.find("-opt")!=string::npos)
         {
            ss >> opt;
         }
         if( s_buff[0] != '-')
         {
            file_name = s_buff;
         }
      }
      ss.clear();

    }

    int option;
    ss.str(opt);
    ss >> option;
    ss.clear();

	basic_string <char>::size_type indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = file_name.find_last_of('\\');
	indexChLinux = file_name.find_last_of('/');
	//
    string file_path;
	if(indexChWin != string::npos)
		file_path = file_name.substr(0,indexChWin) + "\\";
	else if(indexChLinux != string::npos)
		file_path = file_name.substr(0,indexChLinux) + "/";

	CFEMesh* a_mesh = NULL;

    if(option != 3)
	{
       	string fname = file_name + ".msh";
		ifstream is_mesh(fname.c_str(), ios::in);

		std::getline(is_mesh, s_buff);

        if(s_buff.find("#FEM_MSH") != std::string::npos)
		{
           a_mesh = new CFEMesh(NULL, &file_name);
           a_mesh->Read(&is_mesh);
		}
		else
		{
            std::cout<<"Cannot open mesh file "<< fname << endl;
            return EXIT_FAILURE;
		}
	}

	switch(option)
	{
       case 1:
          a_mesh->mHM2NeumannBC();
          break;
       case 2:
          a_mesh->TopSurfaceIntegration();
          break;
       case 3:
          a_mesh = new CFEMesh();
          a_mesh->ConvertShapeCells(file_name + ".asc");
          break;
       default:
          break;
	}

	if(a_mesh)
	{
       delete a_mesh;
	   a_mesh = NULL;
	}
	return EXIT_SUCCESS;
}
