/*
 * CRSIO.h
 *
 *  Created on: Feb 7, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CRSIO_H_
#define CRSIO_H_

namespace FileIO
{
template <class T>
void writeCompressedStorageFmt(std::ostream& os, unsigned n, unsigned* iA, unsigned* jA, T* A)
{
	os.write((char*)&n, sizeof(unsigned));
	os.write((char*)iA, (n + 1) * sizeof(unsigned));
	os.write((char*)jA, iA[n] * sizeof(unsigned));
	os.write((char*)A, iA[n] * sizeof(T));
}

template <class T>
void readCompressedStorageFmt(std::istream& is, unsigned& n, unsigned*& iA, unsigned*& jA, T*& A)
{
	is.read((char*)&n, sizeof(unsigned));
	if (iA != NULL)
	{
		delete[] iA;
		delete[] jA;
		delete[] A;
	}
	iA = new unsigned[n + 1];
	assert(iA != NULL);
	is.read((char*)iA, (n + 1) * sizeof(unsigned));

	jA = new unsigned[iA[n]];
	assert(jA != NULL);
	is.read((char*)jA, iA[n] * sizeof(unsigned));

	A = new T[iA[n]];
	assert(A != NULL);
	is.read((char*)A, iA[n] * sizeof(T));

#ifndef NDEBUG
	// do simple checks
	if (iA[0] != 0)
		std::cerr << "\n"
		          << "CRS matrix: array iA doesn't start with 0"
		          << "\n";

	unsigned i = 0;
	while (i < iA[n] && jA[i] < n)
		++i;
	if (i < iA[n])
		std::cerr << "\n"
		          << "CRS matrix: the " << i << "th entry of jA has the value " << jA[i] << ", which is out of bounds."
		          << "\n";
#endif
}

template <class T>
void readCompressedStorageFmt(std::istream& is, long& n, long*& iA, long*& jA, T*& A, T*& rhs)
{
	is.read((char*)&n, sizeof(long));
	if (iA != NULL)
	{
		delete[] iA;
		delete[] jA;
		delete[] A;
	}

	iA = new long[n + 1];
	assert(iA != NULL);
	is.read((char*)iA, (n + 1) * sizeof(long));

	jA = new long[iA[n]];
	assert(jA != NULL);
	is.read((char*)jA, iA[n] * sizeof(long));

	A = new T[iA[n]];
	assert(A != NULL);
	is.read((char*)A, iA[n] * sizeof(T));

	rhs = new T[n];
	assert(rhs != NULL);
	is.read((char*)rhs, n * sizeof(T));

#ifndef NDEBUG
	// do simple checks
	if (iA[0] != 0)
		std::cerr << "\n"
		          << "CRS matrix: array iA doesn't start with 0"
		          << "\n";

	long i = 0;
	while (i < iA[n] && jA[i] < n)
		++i;
	if (i < iA[n])
		std::cerr << "\n"
		          << "CRS matrix: the " << i << "th entry of jA has the value " << jA[i] << ", which is out of bounds."
		          << "\n";
#endif
}
} // end namespace FileIO

#endif /* CRSIO_H_ */
