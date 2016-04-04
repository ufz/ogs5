/**
 * \file StationIO.cpp
 * 23/03/2010 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "StationIO.h"
#include "StringTools.h"
#include <fstream>

int StationIO::readStationFile(const std::string& path,
                               std::string& name,
                               std::vector<GEOLIB::Point*>* stations,
                               GEOLIB::Station::StationType type)
{
	int returnValue = 1;
	std::string line;

	std::ifstream in(path.c_str());

	if (!in.is_open())
	{
		std::cout << "StationIO::readFromFile() - Could not open file...\n";
		return 0;
	}

	/* try to find a name for the list in the first line of the file */
	getline(in, line);
	if ((line.substr(0, 1)).compare("!") == 0)
		name = line.substr(1, line.length() - 1);
	else
		in.seekg(0); // sets stream to the beginning of the file

	/* read all stations */
	while (getline(in, line))
	{
		GEOLIB::Station* newStation(NULL);
		if (type == GEOLIB::Station::STATION)
			newStation = GEOLIB::Station::createStation(line);
		if (type == GEOLIB::Station::BOREHOLE)
			newStation = GEOLIB::StationBorehole::createStation(line);

		if (newStation)
			// newStation->setList(name);
			stations->push_back(newStation);
		else
			returnValue = -1;
	}

	in.close();

	return returnValue;
}

void StationIO::writeStratigraphyTable(const std::vector<GEOLIB::Point*>* boreholes, const std::string& filename)
{
	size_t maxIterations = 1;
	std::ofstream out(filename.c_str(), std::ios::out);

	for (size_t j = 0; j < boreholes->size(); j++)
		out << static_cast<GEOLIB::StationBorehole*>((*boreholes)[j])->getName() << "\t";
	out << "\n";

	for (size_t i = 0; i < maxIterations; i++)
	{
		for (size_t j = 0; j < boreholes->size(); j++)
		{
			std::vector<std::string> soilNames = static_cast<GEOLIB::StationBorehole*>((*boreholes)[j])->getSoilNames();
			if (!soilNames.empty())
			{
				if (i == 0 && maxIterations < soilNames.size())
					maxIterations = soilNames.size();
				if (i < soilNames.size())
					out << soilNames[i];
			}
			out << "\t";
		}
		out << "\n";
	}

	out.close();
}
