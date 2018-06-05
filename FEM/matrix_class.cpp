/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*========================================================================
   GeoSys - class Matrix (Definition)
          class vec
   Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
   Function:   See the definition below
   programming:
   22/08/2004  WW
   ==========================================================================*/

/// Matrix
#include <cfloat>
#include <cmath>
#include <iomanip>
//
#include "mathlib.h"
#include "matrix_class.h"

#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
#ifdef NEW_EQS
#include "msh_mesh.h"
#include "par_ddc.h"
#endif
#endif

namespace Math_Group
{
MatrixBase::MatrixBase(size_t rows, size_t cols, size_t size)
    : nrows(rows), nrows0(rows), ncols(cols), ncols0(cols), size(size), data(size > 0 ? new double[size] : NULL)
{
	for (size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

MatrixBase::MatrixBase(const MatrixBase& m)
    : nrows(m.nrows), nrows0(m.nrows0), ncols(m.ncols), ncols0(m.ncols0), size(m.size),
      data(size > 0 ? new double[size] : NULL)
{
	for (size_t i = 0; i < size; i++)
		data[i] = m.data[i];
}

MatrixBase::~MatrixBase()
{
	delete[] data;
	data = NULL;
}

// 06.2010. WW
void MatrixBase::ReleaseMemory()
{
	delete[] data;
	data = NULL;
}

// m_results = this*m. m_results must be initialized
void MatrixBase::multi(const MatrixBase& m, MatrixBase& m_result, double fac)
{
#ifdef gDEBUG
	if (ncols != m.Rows() && nrows != m_result.Rows() && m.Cols() != m_result.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	for (size_t i = 0; i < m_result.Rows(); i++)
		for (size_t j = 0; j < m_result.Cols(); j++)
		{
			// m_result(i,j) = 0.0;
			for (size_t k = 0; k < ncols; k++)
				//            m_result(i,j) += fac*data[i*ncols+k]*m(k,j);
				m_result(i, j) += fac * (*this)(i, k) * m(k, j);
		}
}

//
// m_results = this*m1*m2. m_results must be  initialized
void MatrixBase::multi(const MatrixBase& m1, const MatrixBase& m2, MatrixBase& m_result)
{
#ifdef gDEBUG
	if (ncols != m1.Rows() && m1.Cols() != m2.Rows() && m2.Cols() != m_result.Cols() && nrows != m_result.Rows())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	for (size_t i = 0; i < m_result.Rows(); i++)
		for (size_t j = 0; j < m_result.Cols(); j++)
		{
			// m_result(i,j) = 0.0;
			for (size_t k = 0; k < ncols; k++)
				for (size_t l = 0; l < m2.Rows(); l++)
					//                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
					m_result(i, j) += (*this)(i, k) * m1(k, l) * m2(l, j);
		}
}
// vec_result = This*vec. vec_result must be  initialized
void MatrixBase::multi(const double* vec, double* vec_result, double fac)
{
	for (int i = 0; (size_t)i < nrows; i++)
		for (int j = 0; (size_t)j < ncols; j++)
			vec_result[i] += fac * (*this)(i, j) * vec[j];
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2004 WW Implementation
   02/2005 WW Change name
**************************************************************************/
void MatrixBase::Write(std::ostream& os)
{
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(12);

	for (size_t i = 0; i < nrows; i++)
	{
		os << "| ";
		for (size_t j = 0; j < ncols; j++)
			os << (*this)(i, j) << " ";
		os << "| "
		   << "\n";
	}
	os << "\n";
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
   03/2010 TF write whole matrix in one chunk
**************************************************************************/
void MatrixBase::Write_BIN(std::fstream& os)
{
	os.write((char*)data, size * sizeof(double));
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
**************************************************************************/
void MatrixBase::Read_BIN(std::fstream& is)
{
	is.read((char*)data, size * sizeof(double));
}

// Constructors
Matrix::Matrix(size_t rows, size_t cols) : MatrixBase(rows, cols, rows * cols)
{
}

Matrix::Matrix() : MatrixBase(0, 0, 0)
{
}

Matrix::Matrix(const Matrix& m) : MatrixBase(m)
{
}

void Matrix::resize(size_t rows, size_t cols)
{
	if (size > 0)
	{
		delete[] data;
		data = NULL;
	}

	nrows = rows;
	ncols = cols;
	nrows0 = rows;
	ncols0 = ncols;
	size = nrows * ncols;
	data = new double[size];
	for (size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

Matrix::~Matrix()
{
}

//
void Matrix::GetTranspose(Matrix& m)
{
#ifdef gDEBUG
	if (ncols != m.Rows() && nrows != m.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif

	double* m_data = m.getEntryArray();
	const size_t mrows = m.Rows();
	const size_t mcols = m.Cols();
	for (size_t i = 0; i < mrows; i++)
	{
		double* row_m_data = &m_data[i * mcols];
		for (size_t j = 0; j < mcols; j++)
		{
			row_m_data[j] = data[j * ncols + i];
		}
	}
}

// m_results = this*m. m_results must be initialized
void Matrix::multi(const Matrix& m, Matrix& m_result, double fac)
{
#ifdef gDEBUG
	if (ncols != m.Rows() && nrows != m_result.Rows() && m.Cols() != m_result.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	const double* m_data = m.getEntryArray();
	const size_t mcols = m.Cols();
	double* r_data = m_result.getEntryArray();
	const size_t r_rows = m_result.Rows();
	const size_t r_cols = m_result.Cols();

	for (size_t i = 0; i < r_rows; i++)
	{
		const double* row_data = &data[i * ncols];
		double* r_row_data = &r_data[i * r_cols];
		for (size_t j = 0; j < r_cols; j++)
		{
			// r_row_data[j] = 0.0;
			double val = 0.;
			for (size_t k = 0; k < ncols; k++)
			{
				val += row_data[k] * m_data[k * mcols + j];
			}
			r_row_data[j] += val * fac;
		}
	}
}

//
// m_results = this*m1*m2. m_results must be  initialized
void Matrix::multi(const Matrix& m1, const Matrix& m2, Matrix& m_result)
{
#ifdef gDEBUG
	if (ncols != m1.Rows() && m1.Cols() != m2.Rows() && m2.Cols() != m_result.Cols() && nrows != m_result.Rows())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	for (size_t i = 0; i < m_result.Rows(); i++)
		for (size_t j = 0; j < m_result.Cols(); j++)
		{
			// m_result(i,j) = 0.0;
			for (size_t k = 0; k < ncols; k++)
				for (size_t l = 0; l < m2.Rows(); l++)
					//                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
					m_result(i, j) += (*this)(i, k) * m1(k, l) * m2(l, j);
		}
}
// vec_result = This*vec. vec_result must be  initialized
void Matrix::multi(const double* vec, double* vec_result, double fac)
{
	for (size_t i = 0; i < nrows; i++)
	{
		double val = 0.;
		const double* row_data = &data[i * ncols];
		for (size_t j = 0; j < ncols; j++)
		{
			val += row_data[j] * vec[j];
		}
		vec_result[i] += fac * val;
	}
}

double& Matrix::operator()(size_t i, size_t j) const
{
#ifdef gDEBUG
	if (i >= nrows || j >= ncols)
	{
		std::cout << "\n Index exceeds the size of the matrix"
		          << "\n";
		abort();
	}
#endif
	return data[i * ncols + j];
}
void Matrix::LimitSize(size_t nRows, size_t nCols)
{
#ifdef gDEBUG
	if (nRows > nrows0 || nCols > ncols0)
	{
		std::cout << "\n Given size exceeds the original size of the matrix"
		          << "\n";
		abort();
	}
#endif
	nrows = nRows;
	ncols = nCols;
	size = nrows * ncols;
}

//-----------------------------------------------------
// Symmetrical matrix
SymMatrix::SymMatrix(size_t dim) : MatrixBase(dim, dim, (size_t)dim * (dim + 1) / 2)
{
}

SymMatrix::SymMatrix() : MatrixBase(0, 0, 0)
{
}
SymMatrix::SymMatrix(const SymMatrix& m) : MatrixBase(m)
{
}

void SymMatrix::resize(size_t dim)
{
	if (size > 0)
	{
		delete[] data;
		data = NULL;
	}

	nrows = ncols = dim;
	size = (int)nrows * (nrows + 1) / 2;
	data = new double[size];
	nrows0 = ncols0 = dim;
	for (size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

//
void SymMatrix::LimitSize(size_t dim)
{
#ifdef gDEBUG
	if (dim > nrows0)
	{
		std::cout << "\n Given size exceeds the original size of the matrix"
		          << "\n";
		abort();
	}
#endif
	nrows = ncols = dim;
	size = nrows * (nrows + 1) / 2;
}

// m_results = this*m. m_results must be initialized
void SymMatrix::multi(const SymMatrix& m, Matrix& m_result, double fac)
{
#ifdef gDEBUG
	if (ncols != m.Rows() && nrows != m_result.Rows() && m.Cols() != m_result.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	double* r_data = m_result.getEntryArray();
	const size_t r_rows = m_result.Rows();
	const size_t r_cols = m_result.Cols();

	for (size_t i = 0; i < r_rows; i++)
	{
		const double* row_data = &data[(i * (i + 1) / 2)];
		double* row_data_r = &r_data[i * r_cols];

		for (size_t j = 0; j < r_cols; j++)
		{
			// row_data_r[j] = 0.0;
			double val = 0.;
			for (size_t k = 0; k <= i; k++)
			{
				val += row_data[k] * m(k, j);
			}
			for (size_t k = i + 1; k < ncols; k++)
			{
				val += data[(k * (k + 1) / 2) + i] * m(k, j);
			}
			row_data_r[j] += val * fac;
		}
	}
}

//
// m_results = this*m1*m2. m_results must be  initialized
void SymMatrix::multi(const SymMatrix& m1, const Matrix& m2, Matrix& m_result)
{
#ifdef gDEBUG
	if (ncols != m1.Rows() && m1.Cols() != m2.Rows() && m2.Cols() != m_result.Cols() && nrows != m_result.Rows())
	{
		std::cout << "\n The sizes of the two matrices are not matched"
		          << "\n";
		abort();
	}
#endif
	double* r_data = m_result.getEntryArray();
	const size_t r_rows = m_result.Rows();
	const size_t r_cols = m_result.Cols();

	for (size_t i = 0; i < r_rows; i++)
	{
		double* row_data_r = &r_data[i * r_cols];
		for (size_t j = 0; j < r_cols; j++)
		{
			// m_result(i,j) = 0.0;
			double val = 0.;
			for (size_t k = 0; k < ncols; k++)
			{
				const double entry_of_this = data[getArrayIndex(i, k)];
				for (size_t l = 0; l < m2.Rows(); l++)
					val += entry_of_this * m1(k, l) * m2(l, j);
			}
			row_data_r[j] += val;
		}
	}
}
// vec_result = This*vec. vec_result must be  initialized
void SymMatrix::multi(const double* vec, double* vec_result, double fac)
{
	for (size_t i = 0; i < nrows; i++)
	{
		double val = 0.;

		const double* row_data = &data[static_cast<size_t>(i * (i + 1) / 2)];
		for (size_t j = 0; j <= i; j++)
		{
			val += row_data[j] * vec[j];
		}

		for (size_t j = i + 1; j < ncols; j++)
		{
			val += data[static_cast<size_t>(j * (j + 1) / 2) + i] * vec[j];
		}

		vec_result[i] += fac * val;
	}
}

//-----------------------------------------------------
// Diagonal matrix
DiagonalMatrix::DiagonalMatrix(size_t dim) : MatrixBase(dim, dim, dim)
{
	dummy_zero = 0.0;
}

DiagonalMatrix::DiagonalMatrix() : MatrixBase(0, 0, 0)
{
	dummy_zero = 0.0;
}
DiagonalMatrix::DiagonalMatrix(const DiagonalMatrix& m) : MatrixBase(m)
{
	dummy_zero = 0.0;
}

void DiagonalMatrix::resize(size_t dim)
{
	if (size > 0)
	{
		delete[] data;
		data = NULL;
	}

	nrows = ncols = dim;
	size = dim;
	data = new double[size];
	nrows0 = ncols0 = dim;
	for (size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

void DiagonalMatrix::LimitSize(size_t dim)
{
#ifdef gDEBUG
	if (dim > nrows0)
	{
		cout << "\n Given size exceeds the original size of the matrix"
		     << "\n";
		abort();
	}
#endif
	nrows = ncols = dim;
	size = dim;
}

/*========================================================================
   MathLib-Method:
   Task:       Carry out vector operation
   Function:   See the declaration below
   programming:
   05/2005  WW
   ==========================================================================*/
// 1.
template <class T>
vec<T>::vec(int argSize) : _size(argSize), _entry(new T[argSize])
{
#ifdef gDEBUG
	if (!_entry)
	{
		std::cout << "\n*** failed in memory allocating for vec ";
		abort();
	}
#endif
}

template <class T>
vec<T>::vec(const vec<T>& v) : _size(v.Size()), _entry(new T[v.Size()])
{
#ifdef gDEBUG
	if (!_entry)
	{
		std::cout << "\n*** failed in memory allocatiing for vec ";
		abort();
	}
#endif
	for (size_t i = 0; i < _size; i++)
		_entry[i] = v._entry[i];
}

template <class T>
vec<T>::~vec()
{
	delete[] _entry;
	_entry = 0;
}

template <class T>
void vec<T>::resize(const int argSize)
{
	if (_size > 0)
	{
		delete[] _entry;
		_entry = NULL;
	}
	_size = argSize;
	_entry = new T[argSize];
#ifdef gDEBUG
	if (!_entry)
	{
		std::cout << "\n*** failed in memory allocatiing for vec ";
		abort();
	}
#endif
}

template <class T>
void vec<T>::operator=(const vec<T>& v)
{
#ifdef gDEBUG
	if (_size != v.Size())
	{
		std::cout << "\n*** Sizes do not match in vec ";
		abort();
	}
#endif
	for (size_t i = 0; i < _size; i++)
		_entry[i] = v[i];
}

template <class T>
void vec<T>::Write(std::ostream& os) const
{
	for (size_t i = 0; i < _size; i++)
		os << _entry[i] << "  ";
	os << "\n";
}

// 2.
vec<void*>::vec(const int argSize) : _size(argSize)
{
	_entry = new void*[argSize];
#ifdef gDEBUG
	if (!_entry)
	{
		cout << "\n*** failed in memory allocating for vec ";
		abort();
	}
#endif
}

vec<void*>::vec(const vec<void*>& v)
{
	_size = v.Size();
	resize(_size);

#ifdef gDEBUG
	if (!_entry)
	{
		std::cout << "\n*** failed in memory allocating for vec ";
		abort();
	}
#endif
	for (size_t i = 0; i < _size; i++)
		_entry[i] = v._entry[i];
}

vec<void*>::~vec()
{
	delete[] _entry;
	_entry = 0;
}
void vec<void*>::resize(const int argSize)
{
	if (_size > 0)
	{
		delete[] _entry;
		_entry = NULL;
	}
	_size = argSize;
	_entry = new void*[argSize];
#ifdef gDEBUG
	if (!_entry)
	{
		std::cout << "\n*** failed in memory allocating for vec ";
		abort();
	}
#endif
}

void vec<void*>::Write(std::ostream& os) const
{
	for (size_t i = 0; i < _size; i++)
		os << _entry[i] << "  ";
	os << "\n";
}

void vec<void*>::operator=(const vec<void*>& v)
{
#ifdef gDEBUG
	if (_size != v.Size())
	{
		std::cout << "\n*** Sizes do not match in vec ";
		abort();
	}
#endif
	for (size_t i = 0; i < _size; i++)
		_entry[i] = v._entry[i];
}

// 3.
template <class T>
void vec<T*>::operator=(const vec<T*>& v)
{
#ifdef gDEBUG
	if (_size != v.Size())
	{
		std::cout << "\n*** Sizes do not match in vec ";
		abort();
	}
#endif
	for (int i = 0; i < _size; i++)
		_entry[i] = v._entry[i];
}

////////////////////////////////////////////////////////////
#ifdef NEW_EQS
/*\!
 ********************************************************************
   Create sparse matrix table
   01/2006 WW
   08/2007 WW
   10/2007 WW
   02/2008 PCH  Compressed Row Storage with LIS option
   03/2010 WW: CRS storage for matrix algbraic
 ********************************************************************
 */
SparseTable::SparseTable(CFEMesh* a_mesh, bool quadratic, bool symm, StorageType stype)
    : symmetry(symm), storage_type(stype)
{
	long i = 0, j = 0, ii = 0, jj = 0;
	long lbuff0 = 0, lbuff1 = 0;
	long** larraybuffer;
	larraybuffer = NULL;
	//
	// In sparse table, = number of nodes
	rows = a_mesh->GetNodesNumber(quadratic);
	size_entry_column = 0;
	diag_entry = new long[rows];

	if (storage_type == JDS)
	{
		row_index_mapping_n2o = new long[rows];
		row_index_mapping_o2n = new long[rows];
	}
	else if (storage_type == CRS)
	{
		row_index_mapping_n2o = NULL;
		row_index_mapping_o2n = NULL;
	}

	if (symmetry)
	{
		larraybuffer = new long*[rows];
		for (i = 0; i < rows; i++)
		{
			if (storage_type == JDS)
				row_index_mapping_n2o[i] = i;
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			lbuff1 = (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
			larraybuffer[i] = new long[lbuff1 + 1];
			//
			larraybuffer[i][0] = lbuff1;
			for (j = 0; j < lbuff1; j++)
				larraybuffer[i][j + 1] = a_mesh->nod_vector[i]->getConnectedNodes()[j];
			a_mesh->nod_vector[i]->getConnectedNodes().clear();
			for (j = 0; j < lbuff1; j++)
			{
				jj = larraybuffer[i][j + 1];
				if (i <= jj)
					a_mesh->nod_vector[i]->getConnectedNodes().push_back(jj);
			}
		}
	}

	/// CRS storage
	if (storage_type == CRS)
	{
		/// num_column_entries saves vector ptr of CRS
		num_column_entries = new long[rows + 1];

		std::vector<long> A_index;
		long col_index;

		for (i = 0; i < rows; i++)
		{
			num_column_entries[i] = (long)A_index.size();

			for (j = 0; j < (long)a_mesh->nod_vector[i]->getConnectedNodes().size(); j++)
			{
				col_index = a_mesh->nod_vector[i]->getConnectedNodes()[j];

				/// If linear element is used
				if ((!quadratic) && (col_index >= rows))
					continue;

				if (i == col_index)
					diag_entry[i] = (long)A_index.size();
				A_index.push_back(col_index);
			}
		}

		size_entry_column = (long)A_index.size();
		num_column_entries[rows] = size_entry_column;

		entry_column = new long[size_entry_column];
		for (i = 0; i < size_entry_column; i++)
			entry_column[i] = A_index[i];
	}
	else if (storage_type == JDS)
	{
		//
		//--- Sort, from that has maximum connect nodes to that has minimum connect nodes
		//
		for (i = 0; i < rows; i++)
		{
			row_index_mapping_n2o[i] = i;
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			diag_entry[i] = (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
			if (!quadratic)
			{
				lbuff0 = 0;
				for (j = 0; j < diag_entry[i]; j++)
					if (a_mesh->nod_vector[i]->getConnectedNodes()[j] < static_cast<size_t>(rows))
						lbuff0++;
				diag_entry[i] = lbuff0;
			}
			size_entry_column += diag_entry[i];
		}

		//
		for (i = 0; i < rows; i++)
		{
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			lbuff0 = diag_entry[i]; // Nodes to this row
			lbuff1 = row_index_mapping_n2o[i];
			j = i;
			while ((j > 0) && (diag_entry[j - 1] < lbuff0))
			{
				diag_entry[j] = diag_entry[j - 1];
				row_index_mapping_n2o[j] = row_index_mapping_n2o[j - 1];
				j = j - 1;
			}
			diag_entry[j] = lbuff0;
			row_index_mapping_n2o[j] = lbuff1;
		}
		// Old index to new one
		for (i = 0; i < rows; i++)
			row_index_mapping_o2n[row_index_mapping_n2o[i]] = i;
		// Maximum number of columns in the sparse table
		max_columns = diag_entry[0];
		//--- End of sorting
		//
		//--- Create sparse table
		//
		num_column_entries = new long[max_columns];
		entry_column = new long[size_entry_column];
		// 1. Count entries in each column in sparse table
		for (i = 0; i < max_columns; i++)
			num_column_entries[i] = 0;
		for (i = 0; i < rows; i++)
			// 'diag_entry' still is used as a temporary array
			// it stores that numbers of nodes connect to this nodes
			for (j = 0; j < diag_entry[i]; j++)
				num_column_entries[j]++;

		// 2. Fill the sparse table, i.e. store all its entries to
		//    entry_column
		lbuff0 = 0;

		for (i = 0; i < max_columns; i++)
			for (j = 0; j < num_column_entries[i]; j++)
			{
				// ii is the real row index of this entry in matrix
				ii = row_index_mapping_n2o[j];
				// jj is the real column index of this entry in matrix
				jj = a_mesh->nod_vector[ii]->getConnectedNodes()[i];
				entry_column[lbuff0] = jj;

				// Till to this stage, 'diag_entry' is really used to store indices of the diagonal entries.
				// Hereby, 'index' refers to the index in entry_column array.
				if (ii == jj)
					diag_entry[ii] = lbuff0;
				//
				lbuff0++;
			}
	}

	// For the case of symmetry matrix
	if (symmetry)
	{
		for (i = 0; i < rows; i++)
		{
			lbuff0 = larraybuffer[i][0];
			a_mesh->nod_vector[i]->getConnectedNodes().resize(lbuff0);
			//
			for (j = 0; j < lbuff0; j++)
				a_mesh->nod_vector[i]->getConnectedNodes()[j] = larraybuffer[i][j + 1];
		}
		for (i = 0; i < rows; i++)
		{
			delete[] larraybuffer[i];
			larraybuffer[i] = 0;
		}
		delete[] larraybuffer;
		larraybuffer = 0;
	}
}
/*\!
 ********************************************************************
   Create sparse matrix table for each domain
   12/2007 WW
 ********************************************************************
 */
SparseTable::SparseTable(CPARDomain& m_dom, bool quadratic, bool symm) : symmetry(symm)
{
	long i = 0, j = 0, ii = 0, jj = 0;
	long lbuff0 = 0, lbuff1 = 0;
	storage_type = JDS;
	//
	rows = m_dom.GetDomainNodes(quadratic);
	size_entry_column = 0;
	//
	row_index_mapping_n2o = new long[rows];
	row_index_mapping_o2n = new long[rows];
	diag_entry = new long[rows];

	if (symmetry)
	{
		std::vector<long> conc;

		for (i = 0; i < rows; i++)
		{
			row_index_mapping_n2o[i] = i;
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			lbuff1 = m_dom.num_nodes2_node[i];
			//
			for (j = 0; j < lbuff1; j++)
			{
				jj = m_dom.node_conneted_nodes[i][j];
				if (i <= jj)
					conc.push_back(jj);
				m_dom.node_conneted_nodes[i][j] = 0;
			}
			// Number of nodes connected to this node.
			m_dom.num_nodes2_node[i] = (long)conc.size();
			// New
			for (j = 0; j < m_dom.num_nodes2_node[i]; j++)
				m_dom.node_conneted_nodes[i][j] = conc[j];
		}
	}
	//
	//--- Sort, from that has maximum connect nodes to that has minimum connect nodes
	//
	for (i = 0; i < rows; i++)
	{
		row_index_mapping_n2o[i] = i;
		// 'diag_entry' used as a temporary array
		// to store the number of nodes connected to this node
		diag_entry[i] = m_dom.num_nodes2_node[i];
		if (!quadratic)
		{
			lbuff0 = 0;
			for (j = 0; j < diag_entry[i]; j++)
				if (m_dom.node_conneted_nodes[i][j] < rows)
					lbuff0++;
			diag_entry[i] = lbuff0;
		}
		size_entry_column += diag_entry[i];
	}
	//
	for (i = 0; i < rows; i++)
	{
		// 'diag_entry' used as a temporary array
		// to store the number of nodes connected to this node
		lbuff0 = diag_entry[i]; // Nodes to this row
		lbuff1 = row_index_mapping_n2o[i];
		j = i;
		while ((j > 0) && (diag_entry[j - 1] < lbuff0))
		{
			diag_entry[j] = diag_entry[j - 1];
			row_index_mapping_n2o[j] = row_index_mapping_n2o[j - 1];
			j = j - 1;
		}
		diag_entry[j] = lbuff0;
		row_index_mapping_n2o[j] = lbuff1;
	}
	// Old index to new one
	for (i = 0; i < rows; i++)
		row_index_mapping_o2n[row_index_mapping_n2o[i]] = i;
	// Maximum number of columns in the sparse table
	max_columns = diag_entry[0];
	//--- End of sorting
	//
	//--- Create sparse table
	//
	num_column_entries = new long[max_columns];
	entry_column = new long[size_entry_column];
	// 1. Count entries in each column in sparse table
	for (i = 0; i < max_columns; i++)
		num_column_entries[i] = 0;
	for (i = 0; i < rows; i++)
		// 'diag_entry' still is used as a temporary array
		// it stores that numbers of nodes connect to this nodes
		for (j = 0; j < diag_entry[i]; j++)
			num_column_entries[j]++;
	// 2. Fill the sparse table, i.e. store all its entries to
	//    entry_column
	lbuff0 = 0;
	for (i = 0; i < max_columns; i++)
		for (j = 0; j < num_column_entries[i]; j++)
		{
			// ii is the real row index of this entry in matrix
			ii = row_index_mapping_n2o[j];
			// jj is the real column index of this entry in matrix
			jj = m_dom.node_conneted_nodes[ii][i];
			entry_column[lbuff0] = jj;
			// Till to this stage, 'diag_entry' is really used to store indices of the diagonal entries.
			// Hereby, 'index' refers to the index in entry_column array.
			if (ii == jj)
				diag_entry[ii] = lbuff0;
			//
			lbuff0++;
		}
}
/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
void SparseTable::Write(std::ostream& os)
{
	long i, k, counter = 0;

	os.width(10);
	os << "Symmetry: " << symmetry << "\n";
	os << "\n*** Row index  "
	   << "\n";

	if (storage_type == CRS)
	{
		os << "\n*** Sparse entry  "
		   << "\n";
		for (i = 0; i < rows; i++)
		{
			for (k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
				os << entry_column[k] + 1 << " ";
			os << "\n";
		}
	}
	else if (storage_type == JDS)
	{
		for (i = 0; i < rows; i++)
			os << row_index_mapping_n2o[i] + 1 << "\n";
		//
		os << "\n*** Sparse entry  "
		   << "\n";
		for (k = 0; k < max_columns; k++)
		{
			os << "--Column: " << k + 1 << "\n";
			for (i = 0; i < num_column_entries[k]; i++)
			{
				os << entry_column[counter] + 1 << "\n";
				counter++;
			}
			os << "\n";
		}
	}
}

/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
SparseTable::~SparseTable()
{
	if (entry_column)
		delete[] entry_column;
	if (num_column_entries)
		delete[] num_column_entries;
	if (row_index_mapping_n2o)
		delete[] row_index_mapping_n2o;
	if (row_index_mapping_o2n)
		delete[] row_index_mapping_o2n;
	if (diag_entry)
		delete[] diag_entry;
	entry_column = NULL;
	num_column_entries = NULL;
	row_index_mapping_n2o = NULL;
	row_index_mapping_o2n = NULL;
	diag_entry = NULL;
}
/*\!
 ********************************************************************
   Constructor of sparse matrix
   Arguments:
      sparse_table: Sparse graph
      dof:  Degree of freedom given by PDE
   08/2007 WW
   10/2007 WW
   02/2008 PCH Compressed Row Storage
 ********************************************************************/
CSparseMatrix::CSparseMatrix(const SparseTable& sparse_table, const int dof) : DOF(dof)
{
	symmetry = sparse_table.symmetry;
	size_entry_column = sparse_table.size_entry_column;
	max_columns = sparse_table.max_columns;
	rows = sparse_table.rows;
	storage_type = sparse_table.storage_type; // WW
	// Topology mapping from data array to matrix
	// Only refer address
	entry_column = sparse_table.entry_column;
	num_column_entries = sparse_table.num_column_entries;
	row_index_mapping_n2o = sparse_table.row_index_mapping_n2o;
	row_index_mapping_o2n = sparse_table.row_index_mapping_o2n;
	diag_entry = sparse_table.diag_entry;
	// Values of all sparse entries
	entry = new double[dof * dof * size_entry_column + 1];
	entry[dof * dof * size_entry_column] = 0.;
	zero_e = 0.;
//
#if defined(LIS) || defined(MKL) // PCH
	int counter, counter_ptr = 0, counter_col_idx = 0;
	int i, k, ii, jj, J, K;
	int row_in_sparse_table;

	ptr = new int[rows * dof + 1];
	col_idx = new int[dof * dof * size_entry_column];
	entry_index = new int[dof * dof * size_entry_column];

	for (ii = 0; ii < DOF; ii++)
		for (i = 0; i < rows; i++)
		{
			// Store ptr arrary for CRS
			ptr[i + rows * ii] = counter_ptr;
			row_in_sparse_table = row_index_mapping_o2n[i];
			for (jj = 0; jj < DOF; jj++)
			{
				counter = row_in_sparse_table;
				for (k = 0; k < max_columns; k++)
				{
					if (row_in_sparse_table < num_column_entries[k])
					{
						// I = ii * rows + i; // row in global matrix
						// column in global matrix
						J = jj * rows + entry_column[counter];
						K = (ii * DOF + jj) * size_entry_column + counter;

						// Store column index for CRS
						col_idx[counter_col_idx] = J;
						entry_index[counter_col_idx] = K;

						++counter_ptr;
						++counter_col_idx;
						counter += num_column_entries[k];
					}
					else
						break;
				}
			}
		}
	ptr[i + rows * (ii - 1)] = counter_ptr;
#endif
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
CSparseMatrix::~CSparseMatrix()
{
	delete[] entry;
	entry = NULL;

#if defined(LIS) || defined(MKL) // PCH
	delete[] ptr;
	ptr = NULL;
	delete[] col_idx;
	col_idx = NULL;
	delete[] entry_index;
	entry_index = NULL;
#endif
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
double& CSparseMatrix::operator()(const long i, const long j) const
{
#ifdef gDEBUG
	if (i >= rows * DOF || j >= rows * DOF)
	{
		std::cout << "\n Index exceeds the dimension of the matrix"
		          << "\n";
		abort();
	}
#endif
	long ii, jj, ir, jr, k;
	ii = i;
	jj = j;
	if (symmetry)
		if (i > j)
		{
			ii = j;
			jj = i;
		}
	ir = ii % rows;
	jr = jj % rows;
	ii /= rows;
	jj /= rows;
	//
	k = -1;

	if (storage_type == JDS)
	{
		long row_in_parse_table, counter;
		row_in_parse_table = row_index_mapping_o2n[ir];
		counter = row_in_parse_table;
		for (k = 0; k < max_columns; k++)
		{
			if (row_in_parse_table >= num_column_entries[k])
				return zero_e;
			if (entry_column[counter] == jr)
				break; // Found the entry
			counter += num_column_entries[k];
		}
		if (counter >= size_entry_column)
			return zero_e;
		//  Zero entry;
		k = (ii * DOF + jj) * size_entry_column + counter;
	}
	else if (storage_type == CRS)
	{
		/// Left boundary of this row: num_column_entries[ir]
		/// Right boundary of this row: num_column_entries[ir+1]
		/// Search target is jr
		k = binarySearch(entry_column, jr, num_column_entries[ir], num_column_entries[ir + 1]);
		if (k == -1)
			return zero_e;

		k = (ii * DOF + jj) * size_entry_column + k;
	}

	return entry[k]; //
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] = a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator*=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] *= a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator+=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] += a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] = m.entry[i];
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator+=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] += m.entry[i];
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator-=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] -= m.entry[i];
}
/*\!
 ********************************************************************
   Output sparse matrix
   08/2007 WW
   10/2007 WW
   03/2011 WW  CRS
 ********************************************************************/
void CSparseMatrix::Write(std::ostream& os)
{
	//
	long i, k, ii, jj, row_in_parse_table, counter;
	os << "*** Non-zero entries of matrix:  "
	   << "\n";
	os.width(14);
	os.precision(8);
	//
	if (storage_type == CRS)
		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
				for (jj = 0; jj < DOF; jj++)
					for (k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
						// TEST
						// if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN) //DBL_EPSILON)
						os << std::setw(10) << ii * rows + i << " " << std::setw(10) << jj * rows + entry_column[k]
						   << " " << std::setw(15) << entry[(ii * DOF + jj) * size_entry_column + k] << "\n";

	else if (storage_type == JDS)
	{
		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
			{
				row_in_parse_table = row_index_mapping_o2n[i];
				for (jj = 0; jj < DOF; jj++)
				{
					counter = row_in_parse_table;
					for (k = 0; k < max_columns; k++)
					{
						if (row_in_parse_table < num_column_entries[k])
						{
							// TEST
							// if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN) //DBL_EPSILON)
							os << std::setw(10) << ii * rows + i << " " << std::setw(10)
							   << jj * rows + entry_column[counter] << " " << std::setw(15)
							   << entry[(ii * DOF + jj) * size_entry_column + counter] << "\n";
							counter += num_column_entries[k];
						}
						else
							break;
					}
				}
			}
	}
}
//--------------------------------------------------------------
/*!
   \brief Write matrix to a binary file

   03.2011. WW
 */
void CSparseMatrix::Write_BIN(std::ostream& os)
{
	if (storage_type == JDS)
		return;
	//
	if (DOF == 1)
	{
		os.write((char*)&rows, sizeof(long));
		os.write((char*)num_column_entries, (rows + 1) * sizeof(long));
		os.write((char*)entry_column, num_column_entries[rows] * sizeof(long));
		os.write((char*)entry, num_column_entries[rows] * sizeof(double));
	}
	else
	{
		long i, k, ii, jj, size;
		long* ptr;
		long* A_index;
		double* A_value;

		ptr = new long[DOF * rows + 1];
		size = DOF * DOF * num_column_entries[rows];
		A_index = new long[size];
		A_value = new double[size];

		long counter = 0;

		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
			{
				ptr[ii * rows + i] = counter;
				for (jj = 0; jj < DOF; jj++)
					for (k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
					{
						A_index[counter] = jj * rows + entry_column[k];
						A_value[counter] = entry[(ii * DOF + jj) * size_entry_column + k];
						counter++;
					}
			}
		ptr[DOF * rows] = counter;

		ii = DOF * rows;
		os.write((char*)&ii, sizeof(long));
		os.write((char*)ptr, (ii + 1) * sizeof(long));
		os.write((char*)A_index, size * sizeof(long));
		os.write((char*)A_value, size * sizeof(double));

		delete[] ptr;
		delete[] A_index;
		delete[] A_value;
	}
}

/********************************************************************
   Perform A*x
   Arguments:
     vec_sr: M*vec_s-->vec_r
   01/2006 WW
   08/2007 WW
   10/2007 WW
   03/2011 WW      CRS storage
********************************************************************/
void CSparseMatrix::multiVec(double* vec_s, double* vec_r)
{
	long i, j, k, ii, jj, kk, ll, idof, jdof, counter;
	for (i = 0; i < rows * DOF; i++)
		vec_r[i] = 0.0;
	//
	counter = 0;
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1]; j++)
				{
					jj = entry_column[j];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							k = (idof * DOF + jdof) * size_entry_column + j;
							vec_r[kk] += entry[k] * vec_s[ll];
							if (symmetry & (kk != ll))
								vec_r[ll] += entry[k] * vec_s[kk];
						}
					}
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							j = (idof * DOF + jdof) * size_entry_column + counter;
							vec_r[kk] += entry[j] * vec_s[ll];
							if (symmetry & (kk != ll))
								vec_r[ll] += entry[j] * vec_s[kk];
						}
					}
					counter++;
				}
		}
	}
	else // DOF = 1
	{
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1]; j++)
				{
					jj = entry_column[j];
					vec_r[ii] += entry[j] * vec_s[jj];
					if (symmetry & (ii != jj))
						vec_r[jj] += entry[j] * vec_s[ii];
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					vec_r[ii] += entry[counter] * vec_s[jj];
					if (symmetry & (ii != jj))
						vec_r[jj] += entry[counter] * vec_s[ii];
					counter++;
				}
		}
	}
}

/*\!
 ********************************************************************
   Perform A^T*x
   Arguments:
      vec_sr: M^T*vec_s-->vec_r
   10/2010 WW
   03/2011 WW      CRS storage
 ********************************************************************/
void CSparseMatrix::Trans_MultiVec(double* vec_s, double* vec_r)
{
	long i, j, k, ii, jj, kk, ll, idof, jdof, counter;
	for (i = 0; i < rows * DOF; i++)
		vec_r[i] = 0.0;
	//
	counter = 0;
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1]; j++)
				{
					jj = entry_column[j];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							k = (idof * DOF + jdof) * size_entry_column + j;
							vec_r[ll] += entry[k] * vec_s[kk];
							if (symmetry & (kk != ll))
								vec_r[kk] += entry[k] * vec_s[ll];
						}
					}
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							j = (idof * DOF + jdof) * size_entry_column + counter;
							vec_r[ll] += entry[j] * vec_s[kk];
							if (symmetry & (kk != ll))
								vec_r[kk] += entry[j] * vec_s[ll];
						}
					}
					counter++;
				}
		}
	}
	else // DOF = 1
	{
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1]; j++)
				{
					jj = entry_column[j];
					vec_r[jj] += entry[j] * vec_s[ii];
					if (symmetry & (ii != jj))
						vec_r[ii] += entry[j] * vec_s[jj];
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					vec_r[jj] += entry[counter] * vec_s[ii];
					if (symmetry & (ii != jj))
						vec_r[ii] += entry[counter] * vec_s[jj];
					counter++;
				}
		}
	}
}
/*\!
 ********************************************************************
   Set
        A(ii,ii) = x_i,
        A(ii, j) = 0., j!=ii
        A(i, ii) = 0., i!=ii
        b_i -= A(i,k)b_k  // b_k is given
   Programm:
   10/2007 WW
   03/2011 WW  CRS storage
 ********************************************************************/
void CSparseMatrix::Diagonize(const long idiag, const double b_given, double* b)
{
	//
	double vdiag = 0.;
	long j, k, ii, jj, j0;
	long id = idiag % rows;

	ii = idiag / rows;

	if (storage_type == CRS)
	{
		const long row_end = num_column_entries[id + 1];
		/// Diagonal entry and the row where the diagonal entry exists
		j = diag_entry[id];
		vdiag = entry[(ii * DOF + ii) * size_entry_column + j];
		/// Row where the diagonal entry exists
		for (jj = 0; jj < DOF; jj++)
		{
			const long ij = (ii * DOF + jj) * size_entry_column;
			for (k = num_column_entries[id]; k < row_end; k++)
			{
				j0 = entry_column[k];
				if (id == j0 && jj == ii) // Diagonal entry
					continue;
				entry[ij + k] = 0.;
			}
		}
#ifdef colDEBUG
		/// Clean column id
		for (i = 0; i < rows; i++)
		{
			j = binarySearch(entry_column, id, num_column_entries[i], num_column_entries[i + 1]);
			if (j == -1)
				continue;
			j0 = entry_column[j];

			for (jj = 0; jj < DOF; jj++)
			{
				if (i == j0 && ii == jj)
					continue;
				k = (jj * DOF + ii) * size_entry_column + j;
				b[jj * rows + i] -= entry[k] * b_given;
				entry[k] = 0.;
				// Room for symmetry case
			}
		}
#endif
	}
	else if (storage_type == JDS)
	{
		const long kk = ii * DOF;
		long row_in_parse_table, counter;

		// Row is zero
		row_in_parse_table = row_index_mapping_o2n[id];
		counter = row_in_parse_table;
		for (k = 0; k < max_columns; k++)
		{
			if (row_in_parse_table < num_column_entries[k])
			{
				j0 = entry_column[counter];
				for (jj = 0; jj < DOF; jj++)
				{
					if (id == j0 && jj == ii)
					{
						vdiag = entry[(kk + jj) * size_entry_column + counter];
					}
					else
					{
						entry[(kk + jj) * size_entry_column + counter] = 0.;
					}
				}
				counter += num_column_entries[k];
			}
			else
				break;
		}
#ifdef colDEBUG
		//
		counter = 0;
		for (k = 0; k < max_columns; k++)
			for (i = 0; i < num_column_entries[k]; i++)
			{
				i0 = row_index_mapping_n2o[i];
				/*
				   if(i0 == id)
				   {
				   counter++;
				   continue;
				   }
				 */
				j0 = entry_column[counter];
				if (j0 == id)
					for (jj = 0; jj < DOF; jj++)
					{
						if (i0 == j0 && ii == jj)
							continue;
						j = (jj * DOF + ii) * size_entry_column + counter;
						b[jj * rows + i0] -= entry[j] * b_given;
						entry[j] = 0.;
						// Room for symmetry case
					}
				//
				counter++;
			}
#endif
	}
	b[idiag] = vdiag * b_given;
}

/*\!
 ********************************************************************
   M^{-1}*A

          a_ij  i=j
   M = {
          0     i!=j
   Programm:
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::Precond_Jacobi(double* vec_s, double* vec_r)
{
	long i, idof;
	double diag = 0.;
	//
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		for (i = 0; i < rows; i++)
			for (idof = 0; idof < DOF; idof++)
			{
				diag = entry[(idof * DOF + idof) * size_entry_column + diag_entry[i]];
				if (fabs(diag) < DBL_MIN)
					//        if(fabs(diag)<DBL_EPSILON)
					diag = 1.0;
				//  std::cout<<"Diagonal entry is zero. Abort simulation!!  " <<"\n";
				vec_r[idof * rows + i] = vec_s[idof * rows + i] / diag;
			}
		//
	}
	else // DOF = 1

		for (i = 0; i < rows; i++)
		{
			diag = entry[diag_entry[i]];
			// if(fabs(diag)<DBL_EPSILON)
			if (fabs(diag) < DBL_MIN)
				diag = 1.0;
			//   std::cout<<"Diagonal entry is zero. Abort simulation!!  " <<"\n";
			//
			vec_r[i] = vec_s[i] / diag;
		}
}
#if defined(USE_MPI)
/*\!
 ********************************************************************
   Get diagonal entries of the matrix
   Programm:
   12/2007 WW
 ********************************************************************/
void CSparseMatrix::DiagonalEntries(double* diag_e)
{
	long i, idof;
	//
	if (DOF > 1)
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		for (i = 0; i < rows; i++)
			for (idof = 0; idof < DOF; idof++)
				diag_e[idof * rows + i] = entry[(idof * DOF + idof) * size_entry_column + diag_entry[i]];
	//
	else // DOF = 1

		for (i = 0; i < rows; i++)
			diag_e[i] = entry[diag_entry[i]];
}
#endif // USE_MPI

#if defined(LIS) || defined(MKL)
/********************************************************************
   Get sparse matrix values in compressed row storage
   Programm:
   02/2008 PCH
********************************************************************/
int CSparseMatrix::GetCRSValue(double* value)
{
	int success = 1;
	int i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < size_entry_column * DOF * DOF; ++i)
		value[i] = entry[entry_index[i]];

	return success;
}
#endif // LIS
#endif // NEW_EQS
///////////////////////////////////////////////////////////
} // Namespace

using Math_Group::vec;
using Math_Group::SymMatrix;

template class vec<int>;
template class vec<long>;
template class vec<double>;

// End of class Matrix
//==========================================================================
