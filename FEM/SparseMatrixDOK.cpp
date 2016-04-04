/*
 * SparseMatrixDOK.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: NW
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <iostream>
#include <cmath>
#include <iomanip>
#include "SparseMatrixDOK.h"

namespace Math_Group
{
// SparseMatrixDOK
//:Matrix(0)
SparseMatrixDOK::SparseMatrixDOK(size_t _nrows, size_t _ncols)
{
	nrows = _nrows;
	ncols = _ncols;
	size = nrows * ncols;
	// data = new double[dim];
	nrows0 = nrows;
	ncols0 = ncols;
	//    for(int i=0; i<size; i++) data[i] = 0.0;
	mat_row.resize(_nrows);
	mat_col.resize(_ncols);
	dummy_zero = 0.0;
	non_zero_entry_size = 0;
	is_constructed = false;
#if defined(LIS) || defined(MKL)
	ptr = NULL;
	col_idx = NULL;
	entry_index = NULL;
#endif
}

SparseMatrixDOK::SparseMatrixDOK(size_t dim) //:Matrix(0)
{
	nrows = ncols = dim;
	size = dim;
	//    data = new double[dim];
	nrows0 = ncols0 = dim;
	//    for(int i=0; i<size; i++) data[i] = 0.0;
	mat_row.resize(dim);
	mat_col.resize(dim);
	dummy_zero = 0.0;
	non_zero_entry_size = 0;
	is_constructed = false;
#if defined(LIS) || defined(MKL)
	ptr = NULL;
	col_idx = NULL;
	entry_index = NULL;
#endif
}

SparseMatrixDOK::SparseMatrixDOK() //:Matrix(0)
{
	//   Sym = true;
	nrows = 0;
	ncols = 0;
	nrows0 = 0;
	ncols0 = 0;
	size = 0;
	//   data = 0;
	dummy_zero = 0.0;
	non_zero_entry_size = 0;
	is_constructed = false;
#if defined(LIS) || defined(MKL)
	ptr = NULL;
	col_idx = NULL;
	entry_index = NULL;
#endif
}

SparseMatrixDOK::~SparseMatrixDOK()
{
#if defined(LIS) || defined(MKL)
	if (ptr)
		delete[] ptr;
	ptr = NULL;
	if (col_idx)
		delete[] col_idx;
	col_idx = NULL;
	if (entry_index)
		delete[] entry_index;
	entry_index = NULL;
#endif
}

//:Matrix(0)
SparseMatrixDOK::SparseMatrixDOK(const SparseMatrixDOK& m)
{
	Sym = m.Sym;
	nrows = m.nrows;
	ncols = m.ncols;
	nrows0 = m.nrows0;
	ncols0 = m.ncols0;
	size = m.size;
	//   data = new double[size];
	// for(int i=0; i<size; i++) data[i] = 0.0;
	dummy_zero = 0.0;
	non_zero_entry_size = m.non_zero_entry_size;
#if defined(LIS) || defined(MKL)
	ptr = m.ptr;
	col_idx = m.col_idx;
	entry_index = m.entry_index;
#endif
}

SparseMatrixDOK& SparseMatrixDOK::operator=(double a)
{
	//    for(int i=0; i<size; i++) data[i] = a;

	row_iter ii;
	col_iter jj;

	const row_iter row_end = this->mat_row.end();

	for (ii = this->mat_row.begin(); ii != row_end; ii++)
	{
		const col_iter col_end = (*ii).end();
		for (jj = (*ii).begin(); jj != col_end; jj++)
		{
			//        for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
			(*jj).second = a;
		}
	}

	return *this;
}
void SparseMatrixDOK::operator*=(double a)
{
	row_iter ii;
	col_iter jj;

	for (ii = this->mat_row.begin(); ii != this->mat_row.end(); ii++)
	{
		for (jj = (*ii).begin(); jj != (*ii).end(); jj++)
		{
			//        for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
			(*jj).second *= a;
		}
	}
}
void SparseMatrixDOK::operator+=(double a)
{
	row_iter ii;
	col_iter jj;
	const row_iter row_end = this->mat_row.end();

	for (ii = this->mat_row.begin(); ii != row_end; ii++)
	{
		const col_iter col_end = (*ii).end();
		for (jj = (*ii).begin(); jj != col_end; jj++)
		{
			//        for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
			(*jj).second += a;
		}
	}
}

//
void SparseMatrixDOK::operator=(const SparseMatrixDOK& m)
{
#ifdef gDEBUG
	if (nrows != m.Rows() || ncols != m.Cols())
	{
		cout << "\n The sizes of the two matrices are not matched"
		     << "\n";
		abort();
	}
#endif
	const mat_t& tmp_mat = m.mat_row;
	col_t* col;
	mat_t::const_iterator ii;
	col_t::const_iterator jj;

	// for(ii=tmp_mat.begin(); ii!=tmp_mat.end(); ii++){
	//    for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
	//      this->mat[(*ii).first][(*jj).first] = (*jj).second;
	//    }
	//}

	for (size_t i = 0; i < this->nrows; i++)
	{
		col = const_cast<col_t*>(&tmp_mat[i]);
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			this->mat_row[i][(*jj).first] = (*jj).second;
		}
	}
}

//
void SparseMatrixDOK::operator+=(const SparseMatrixDOK& m)
{
#ifdef gDEBUG
	if (nrows != m.Rows())
	{
		cout << "\n The sizes of the two matrices are not matched"
		     << "\n";
		abort();
	}
#endif
	const mat_t& tmp_mat = m.mat_row;
	// col_t *col;
	// mat_t::const_iterator ii;
	// col_t::const_iterator jj;

	// for(ii=m.mat.begin(); ii!=m.mat.end(); ii++){
	//    for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
	//      this->mat[(*ii).first][(*jj).first] += (*jj).second;
	//    }
	//}

	const long n_rows = (long)this->nrows;
#ifdef _OPENMP
#pragma omp parallel for // default(none) //shared(tmp_mat)
#endif
	for (long i = 0; i < n_rows; i++)
	{
		col_t& this_col = this->mat_row[i];
		const col_t* col = const_cast<col_t*>(&tmp_mat[i]);
		const col_t::const_iterator& col_end = col->end();

		for (col_t::const_iterator jj = col->begin(); jj != col_end; jj++)
		{
			this_col[(*jj).first] += (*jj).second;
			// this->mat_row[i][(*jj).first] += (*jj).second;
		}
	}
}

//
void SparseMatrixDOK::operator-=(const SparseMatrixDOK& m)
{
#ifdef gDEBUG
	if (nrows != m.Rows()) // Assertion, will be removed
	{
		cout << "\n The sizes of the two matrices are not matched"
		     << "\n";
		abort();
	}
#endif
	const mat_t& tmp_mat = m.mat_row;
	col_t* col;
	mat_t::const_iterator ii;
	col_t::const_iterator jj;

	// for(ii=m.mat.begin(); ii!=m.mat.end(); ii++){
	//    for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
	//      this->mat[(*ii).first][(*jj).first] -= (*jj).second;
	//    }
	//}

	for (size_t i = 0; i < this->nrows; i++)
	{
		col = const_cast<col_t*>(&tmp_mat[i]);
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			this->mat_row[i][(*jj).first] -= (*jj).second;
		}
	}
}
//
// const
double& SparseMatrixDOK::operator()(size_t i, size_t j)
{
#ifdef gDEBUG
	if (i >= nrows || j >= nrows)
	{
		cout << "\n Index exceeds the size of the matrix"
		     << "\n";
		abort();
	}
#endif

	return this->mat_row[i][j];
}

double SparseMatrixDOK::operator()(size_t i, size_t j) const
{
#ifdef gDEBUG
	if (i >= nrows || j >= nrows)
	{
		cout << "\n Index exceeds the size of the matrix"
		     << "\n";
		abort();
	}
#endif

	col_t::const_iterator ii = mat_row[i].find(j);
	if (ii != mat_row[i].end())
		return ii->second;
	else
		return .0;
}

double& SparseMatrixDOK::operator()(size_t i) // const
{
#ifdef gDEBUG
	if (i >= size)
	{
		cout << "\n Index exceeds the size of the matrix"
		     << "\n";
		abort();
	}
#endif

	return mat_row[i][i];
}

void SparseMatrixDOK::LimitSize(size_t nRows, size_t nCols)
{
#ifdef gDEBUG
	if (nRows > nrows0 || nCols > ncols0)
	{
		cout << "\n Given size exceeds the original size of the matrix"
		     << "\n";
		abort();
	}
#endif
	nrows = nRows;
	ncols = nCols;
	size = nRows * nCols;
}

size_t SparseMatrixDOK::SizeOfNonZeroEntries()
{
	return non_zero_entry_size;
}

void SparseMatrixDOK::CalculateNonZeroEntries()
{
	row_iter ii;
	col_iter jj;

	long cnt = 0;

	for (ii = this->mat_row.begin(); ii != this->mat_row.end(); ii++)
	{
		cnt += (*ii).size();
	}

	non_zero_entry_size = cnt;
}

/*****************************************************************/ /**
  Set
  A(ii,ii) = x_i,
  A(ii, j) = 0., j!=ii
  A(i, ii) = 0., i!=ii
  b_i -= A(i,k)b_k  // b_k is given
  Programm:
  10/2007 WW
  ********************************************************************/
void SparseMatrixDOK::Diagonize(size_t idiag, const double b_given, double* b)
{
	row_iter ii;
	col_iter jj;

	if (!this->is_constructed)
	{
		std::cout << "-> Constructing colmun info in SparseMatrixDOK::Diagonize()"
		          << "\n";
		const size_t n_rows = this->mat_row.size();
		for (size_t i = 0; i < n_rows; i++)
		{
			col_t& i_mat = this->mat_row[i];
			col_t::const_iterator itr_end = i_mat.end();
			//#ifdef USE_HASHMAP
			//      const col_id_itr colid_end = set_col_id[cnt_rows].end();
			//      for(col_id_itr kk=set_col_id[cnt_rows].begin(); kk!=colid_end; kk++){
			//        jj = ii->find(*kk);
			//#else
			for (jj = i_mat.begin(); jj != itr_end; jj++)
			{
				//#endif
				this->mat_col[(*jj).first].push_back(i);
			}
		}
		is_constructed = true;
	}

	double vdiag = .0;

	col_t& i_mat = this->mat_row[idiag];
	for (jj = i_mat.begin(); jj != i_mat.end(); jj++)
	{
		if ((*jj).first == idiag)
		{
			vdiag = (*jj).second;
		}
		else
		{
			(*jj).second = 0.0;
		}
	}

	mat_col_t::iterator itr_col_diag = this->mat_col.begin() + idiag;
	const size_t n_itr_col_diag = itr_col_diag->size();
	for (size_t i = 0; i < n_itr_col_diag; i++)
	{
		size_t row_id = (*itr_col_diag)[i];
		if (row_id == idiag)
			continue;
		jj = this->mat_row[row_id].find(idiag);
		b[row_id] -= (*jj).second * b_given;
		(*jj).second = 0.0;
	}

	// for (size_t i=0; i<this->mat_row.size(); i++) {
	//  if (i==idiag) continue;
	//  col_t &i_mat = this->mat_row[i];
	//  for(jj=i_mat.begin(); jj!=i_mat.end(); jj++){
	//    if ((*jj).first==idiag) {
	//      b[i] -= (*jj).second*b_given;
	//      (*jj).second = 0.0;
	//    } else if ((*jj).first > idiag) {
	//      break;
	//    }
	//  }
	//}

	b[idiag] = vdiag * b_given;
}

void SparseMatrixDOK::multiVec(double* vec_s, double* vec_r)
{
	(void)vec_s; // unused
	(void)vec_r; // unused
	std::cout << "***ERROR: SparseMatrixDOK::multiVec() is not implemented yet."
	          << "\n";
}

void SparseMatrixDOK::Write(std::ostream& os, int format)
{
	row_iter ii;
	col_iter jj;

#ifdef USE_HASHMAP
	if (this->set_col_id.size() == 0)
		this->ConstructSortedColumnID();
#endif

	//
	if (format == 0)
	{
		os << "*** Non-zero entries of matrix:  "
		   << "\n";
		// os.width(25);
		// os.precision(10);

		os.setf(std::ios_base::scientific, std::ios_base::floatfield);
		os.precision(20);

		for (size_t i = 0; i < this->mat_row.size(); i++)
		{
#ifdef USE_HASHMAP
			for (col_id_itr kk = set_col_id[i].begin(); kk != set_col_id[i].end(); kk++)
			{
				jj = this->mat_row[i].find(*kk);
#else
			for (jj = this->mat_row[i].begin(); jj != this->mat_row[i].end(); jj++)
			{
#endif
				os << std::setw(10) << i + 1 << " " << std::setw(10) << (*jj).first + 1 << " " << std::setw(15)
				   << (*jj).second << "\n";
			}
		}

		os.unsetf(std::ios_base::scientific);
		//
	} //
	else if (format == 1)
	{
		os << Dim() << "\n";
		os.setf(std::ios::scientific);
		// os.width(25);
		os.precision(10);

		for (size_t i = 0; i < this->mat_row.size(); i++)
		{
			for (size_t j = 0; j < this->mat_row.size(); j++)
			{
				jj = this->mat_row[i].find(j);
				double v = 0.0;
				if (jj != this->mat_row[i].end())
				{
					v = jj->second;
				}
				os << v << " ";
			}

			os << "\n";
		} //
	}
}

bool SparseMatrixDOK::IsSymmetry()
{
	row_iter ii;
	col_iter jj, jj2;
#define ZERO_TOLERANCE 1.E-6
	//
	for (size_t i = 0; i < this->mat_row.size(); i++)
	{
		for (jj = this->mat_row[i].begin(); jj != this->mat_row[i].end(); jj++)
		{
			if (jj->first < i)
				continue;
			jj2 = this->mat_row[jj->first].find(i);
			if (jj2 != this->mat_row[jj->first].end())
			{
				double diff = jj->second - jj2->second;
				// if (jj->second != jj2->second) {
				if (fabs(diff) > ZERO_TOLERANCE)
				{
					std::cout << "->unsymmetry: " << i << " - " << jj->first << "\n";
					return false;
				}
				//} else if (jj->second!=0.0) {
			}
			else
			{
				std::cout << "->unsymmetry: " << i << " - " << jj->first << "\n";
				return false;
			}
		}
	}

	return true;
}

#ifdef USE_HASHMAP
void SparseMatrixDOK::ConstructSortedColumnID()
{
	if (this->set_col_id.size() > 0)
	{
		cout << "->SparseMatrixDOK::ConstructSortedColumnID() - Already sorted"
		     << "\n";
		return;
	}
	cout << "->SparseMatrixDOK::ConstructSortedColumnID()"
	     << "\n";

	// Construct list of sorted col id
	const size_t n_row = this->mat_row.size();
	this->set_col_id.resize(n_row);
	for (size_t i = 0; i < n_row; i++)
	{
		col_t::const_iterator col_end = this->mat_row[i].end();
		col_id_t& colid = this->set_col_id[i];
		for (col_iter jj = this->mat_row[i].begin(); jj != col_end; jj++)
		{
			colid.insert((*jj).first);
		}
	}
}
#endif

#if defined(LIS) || defined(MKL)
void SparseMatrixDOK::ConstructCRSstructure()
{
//
#ifdef USE_HASHMAP
	if (this->set_col_id.size() == 0)
		this->ConstructSortedColumnID();
#endif

	// ptr:         an integer array with a length of n + 1, which stores the starting
	//              points of the rows of the arrays value and index.
	// col_idx:     an integer array with a length of nnz, which stores the column
	//              numbers of the nonzero elements stored in the array value.
	// entry_index: ?
	const long total_matrix_length = nrows;
	const long total_matrix_entry_size = non_zero_entry_size;

	this->ptr = new int[total_matrix_length + 1];
	this->col_idx = new int[total_matrix_entry_size];
	this->entry_index = new int[total_matrix_entry_size];

	long counter_ptr = 0, counter_col_idx = 0;
	long J, K;

	row_iter ii;
	col_iter jj;
	long cnt_row = 0;

	for (ii = this->mat_row.begin(); ii != this->mat_row.end(); ii++)
	{
		ptr[cnt_row++] = counter_ptr; // starting point of the row
#ifdef USE_HASHMAP
		col_id_t& colid = this->set_col_id[cnt_row - 1];
		for (col_id_itr kk = colid.begin(); kk != colid.end(); kk++)
		{
			jj = ii->find(*kk);
#else
		for (jj = (*ii).begin(); jj != (*ii).end(); jj++)
		{
#endif
			J = (*jj).first; // column in global matrix
			K = counter_ptr; // index in entry

			this->col_idx[counter_col_idx] = J;
			this->entry_index[counter_col_idx] = K;

			++counter_ptr;
			++counter_col_idx;
		}
	}
	ptr[total_matrix_length] = counter_ptr;
}

int SparseMatrixDOK::GetCRSValue(double* v)
{
	int success = 1;

	row_iter ii;
	col_iter jj;
	long cnt = 0;
	long cnt_rows = 0;

	const row_iter row_end = this->mat_row.end();

	for (ii = this->mat_row.begin(); ii != row_end; ii++)
	{
#ifdef USE_HASHMAP
		const col_id_itr colid_end = set_col_id[cnt_rows].end();
		for (col_id_itr kk = set_col_id[cnt_rows].begin(); kk != colid_end; kk++)
		{
			jj = ii->find(*kk);
#else
		const col_iter col_end = (*ii).end();
		for (jj = (*ii).begin(); jj != col_end; jj++)
		{
#endif
			v[cnt++] = (*jj).second;
		}
		cnt_rows++;
	}

	return success;
}
#endif

} // end namespace Math_Group
