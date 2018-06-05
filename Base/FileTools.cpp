/**
 * \file FileTools.h
 * 26/4/2010 LB Initial implementation
 *
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FileTools.h"
#include "StringTools.h"

#include <iostream>
#include <fstream>
#include <cstdio>

#include <sys/stat.h>

// for getCwd
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(std::string const& strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(), &stFileInfo);

	if (intStat == 0)
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	else
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;

	return blnReturn;
}

bool HasCRInLineEnding(std::string const& strFilename)
{
	std::ifstream is(strFilename.c_str(), std::ios::in | std::ios::binary);
	if (!is)
	{
		std::cout << "*** error: could not open " << strFilename.data() << std::endl;
		return false;
	}

	bool foundCR = false;
	while (is.good())
	{
		char c;
		is.read(&c, sizeof(c));
		if (c == '\r')
		{
			foundCR = true;
			break;
		}
		else if (c == EOF || c == '\n')
		{
			break;
		}
	}

	is.close();

	return foundCR;
}

char getDirSep()
{
#ifdef WIN32
	return '\\';
#else
	return '/';
#endif
}

std::string pathJoin(const std::string& path1, const std::string& path2)
{
	if (path1.empty())
		return path2;
	if (path2.empty())
		return path1;

	const char dirSep = getDirSep();

	const std::string s = rtrim(path1, dirSep) + dirSep + ltrim(path2, dirSep);

	return s;
}

std::string pathBasename(const std::string& path)
{
	if (path.empty())
		return "";

	const char dirSep = getDirSep();
	const std::string p = rtrim(path, dirSep);

	const size_t idx = p.find_last_of(dirSep);
	if (idx == std::string::npos)
	{
		return path; // no dirSep in path
	}
	else
	{
		return p.substr(idx + 1);
	}
}

std::string pathDirname(const std::string& path)
{
	if (path.empty())
		return ".";

	const char dirSep = getDirSep();
	const std::string p = rtrim(path, dirSep);

	const size_t idx = p.find_last_of(dirSep);
	if (idx == std::string::npos)
	{
		return "."; // no dirSep in path
	}
	else if (idx == 0)
	{
		return std::string(1, dirSep); // only one dirSep at the beginning of path
	}
	else
	{
		return p.substr(0, idx);
	}
}

std::string getCwd()
{
	char cwd[FILENAME_MAX];

#ifdef WIN32
	_getcwd(cwd, FILENAME_MAX);
#else
	getcwd(cwd, FILENAME_MAX);
#endif
	return cwd;
}
