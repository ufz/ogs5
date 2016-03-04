/*
 * SparseMatrixDOK.h
 *
 *  Created on: Oct 6, 2011
 *      Author: NW
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef SPARSEMATRIXDOK_H_
#define SPARSEMATRIXDOK_H_

#include <vector>
#include <map>
#include <set>


namespace Math_Group {

/** Sparse matrix (type: Dictionary of keys) */
class SparseMatrixDOK {
public:
	//#define USE_HASHMAP
#ifdef USE_HASHMAP
	typedef std::vector<hash_map<size_t , double> > mat_t;
	typedef stdext::hash_map<size_t, double> col_t;
	typedef std::vector<std::set<long>> mat_id_t;
	typedef std::set<long> col_id_t;
	typedef col_id_t::const_iterator col_id_itr;
#else
	typedef std::vector<std::map<size_t, double> > mat_t;
	typedef std::map<size_t, double> col_t;
#endif
	typedef mat_t::iterator row_iter;
	typedef col_t::iterator col_iter;
	typedef std::vector<std::vector<size_t> > mat_col_t;
protected:
	size_t nrows, nrows0;
	size_t ncols, ncols0;
	size_t size;
	size_t non_zero_entry_size;
	bool Sym;
	bool is_constructed;
private:
	double dummy_zero;
	mat_t mat_row;
	mat_col_t mat_col;
#ifdef USE_HASHMAP
	mat_id_t set_col_id;
#endif

public:
	SparseMatrixDOK(size_t dim);
	SparseMatrixDOK(size_t rows, size_t cols);
	SparseMatrixDOK();
	explicit SparseMatrixDOK(const SparseMatrixDOK& m);

	const mat_t& GetRawData() const
	{
		return mat_row;
	}

	mat_t& GetRawData()
	{
		return mat_row;
	}
	;
	size_t SizeOfNonZeroEntries();
	void CalculateNonZeroEntries();

	//void resize(const int dim);

	virtual ~SparseMatrixDOK();

	// Operators
	SparseMatrixDOK& operator =(double a);
	void operator *=(double a);
	void operator +=(double a);
	void operator =(const SparseMatrixDOK& m);
	void operator +=(const SparseMatrixDOK& m);
	void operator -=(const SparseMatrixDOK& m);
	void multiVec(double *vec_s, double *vec_r);
	void LimitSize(size_t nRows, size_t nCols = 1);
	void Diagonize(size_t idiag, const double b_given, double *b);

	// Access to members
	double& operator()(size_t i, size_t j); //const;
	double operator() (size_t i, size_t j) const; //const;
	double& operator()(size_t i); //const;

	size_t Rows() const
	{
		return nrows;
	}
	size_t Cols() const
	{
		return ncols;
	}
	size_t Size() const
	{
		return size;
	}
	size_t Dim() const
	{
		return nrows;
	}

	void Write(std::ostream &os = std::cout, int format = 0);
	bool IsSymmetry();

#if defined(LIS) || defined(MKL)  // These two pointers are in need for Compressed Row Storage
	int* ptr;
	int* col_idx;
	int* entry_index;
	int GetCRSValue(double* value);
	void ConstructCRSstructure();
#endif
#ifdef USE_HASHMAP
	void ConstructSortedColumnID();
#endif
};

} // end namespace Math_Group

#endif /* SPARSEMATRIXDOK_H_ */
