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

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <string>

/**
 * Returns true if given file exists. From
 * http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(std::string const& strFilename);

/// Returns true if given file includes CR
bool HasCRInLineEnding(std::string const& strFilename);

/**
 * @brief computes the basename of the given path, i.e. the component after the
 * last diretory separator (/ or \).
 */
std::string pathBasename(const std::string& path);

/**
 * @brief computes the dirname of the given path, i.e. the component before the
 * last diretory separator (/ or \).
 */
std::string pathDirname(const std::string& path);

/**
 * @brief joins two paths using the correct directory separator.
 *
 * trailing and preceding (back)slashes at the join point are ignored.
 *
 * @returns a string: path1/path2. if any of the paths is empty, only the other
 * one is returned
 */
std::string pathJoin(const std::string& path1, const std::string& path2);

/// returns the current process working directory
std::string getCwd();

char getDirSep();

#endif  // FILETOOLS_H
