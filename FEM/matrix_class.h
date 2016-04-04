/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*========================================================================
   GeoSys - class Matrix, Sparse matrix (Declaration)
   Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
   Function:   See the declaration below
   Design and programm WW
   03/2010 some improvements TF
   03/2014 Rewritten A*B and A*x functions WW
   ==========================================================================*/
#ifndef matrix_class_INC
#define matrix_class_INC

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#ifdef NEW_EQS
namespace MeshLib
{
class CFEMesh;
}
// 08.2007 WW
class CPARDomain;
#endif
//#define OverLoadNEW_DELETE

namespace Math_Group
{
/// Base class for matrix objects
class MatrixBase
{
public:
	/**
	 * Constructor
	 * @param rows  the number of rows
	 * @param cols  the number of columns
	 * @param size  data size
	 */
	MatrixBase(size_t rows, size_t cols, size_t size);

	MatrixBase(const MatrixBase& m);

	virtual ~MatrixBase();

	void ReleaseMemory();
	size_t Rows() const { return nrows; }
	size_t Cols() const { return ncols; }
	size_t Size() const { return size; }
	// Access to members
	virtual double& operator()(const size_t i, const size_t j = 0) const = 0;
	double* getEntryArray() { return data; }
	const double* getEntryArray() const { return data; }
	double& operator[](const size_t i) const { return data[i]; }
	// Operators
	inline void operator=(double a)
	{
		for (size_t i = 0; i < size; i++)
			data[i] = a;
	}

	inline void operator*=(double a)
	{
		for (size_t i = 0; i < size; i++)
			data[i] *= a;
	}

	inline void operator/=(double a)
	{
		for (size_t i = 0; i < size; i++)
			data[i] /= a;
	}

	inline void operator+=(double a)
	{
		for (size_t i = 0; i < size; i++)
			data[i] += a;
	}

	inline void operator-=(double a)
	{
		for (size_t i = 0; i < size; i++)
			data[i] -= a;
	}

	inline void operator=(const MatrixBase& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < nrows; i++)
			for (size_t j = 0; j < ncols; j++)
				data[i * ncols + j] = m(i, j);
	}

	inline void operator+=(const MatrixBase& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < nrows; i++)
			for (size_t j = 0; j < ncols; j++)
				data[i * ncols + j] += m(i, j);
	}

	inline void operator-=(const MatrixBase& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols()) // Assertion, will be removed
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < nrows; i++)
			for (size_t j = 0; j < ncols; j++)
				data[i * ncols + j] -= m(i, j);
	}

	// vec_result = This*vec. vec_result must be initialized.
	virtual void multi(const double* vec, double* vec_result, double fac = 1.0);
	// m_result = this*m. m_result must be initialized.
	virtual void multi(const MatrixBase& m, MatrixBase& m_result, double fac = 1.0);
	// m_result = this*m1*m2. m_result must be initialized.To be removed
	virtual void multi(const MatrixBase& m1, const MatrixBase& m2, MatrixBase& m_result);

	// Print
	void Write(std::ostream& os = std::cout);
	void Write_BIN(std::fstream& os);
	void Read_BIN(std::fstream& is);

protected:
	size_t nrows, nrows0;
	size_t ncols, ncols0;
	size_t size;
	double* data;
};

#ifdef _MSC_VER
#pragma warning(disable : 4522)
#endif

/// Dense matrix
class Matrix : public MatrixBase
{
public:
	using MatrixBase::operator=;
	using MatrixBase::operator+=;
	using MatrixBase::operator-=;
	using MatrixBase::operator*=;
	using MatrixBase::operator/=;

	explicit Matrix(size_t rows, size_t cols = 1);
	Matrix();
	Matrix(const Matrix& m);
	//
	void resize(size_t rows, size_t cols = 1);
	//
	virtual ~Matrix();

	// Operators
	inline void operator=(const Matrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		const double* m_data = m.getEntryArray();
		for (size_t i = 0; i < size; i++)
			data[i] = m_data[i];
	}

	inline void operator+=(const Matrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		const double* m_data = m.getEntryArray();
		for (size_t i = 0; i < size; i++)
			data[i] += m_data[i];
	}

	inline void operator-=(const Matrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols()) // Assertion, will be removed
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		const double* m_data = m.getEntryArray();
		for (size_t i = 0; i < size; i++)
			data[i] -= m_data[i];
	}

	void GetTranspose(Matrix& m);

	// vec_result = This*vec. vec_result must be initialized.
	virtual void multi(const double* vec, double* vec_result, double fac = 1.0);
	// m_result = this*m. m_result must be initialized.
	virtual void multi(const Matrix& m, Matrix& m_result, double fac = 1.0);
	// m_result = this*m1*m2. m_result must be initialized.To be removed
	virtual void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

	// Access to members
	virtual double& operator()(const size_t i, const size_t j = 0) const;
	void LimitSize(size_t nRows, size_t nCols = 1);
};

// Symmetrical matrix. 12-01-2005. WW
class SymMatrix : public MatrixBase
{
public:
	using MatrixBase::operator=;
	using MatrixBase::operator+=;
	using MatrixBase::operator-=;
	using MatrixBase::operator*=;
	using MatrixBase::operator/=;

	explicit SymMatrix(size_t dim);
	SymMatrix();
	SymMatrix(const SymMatrix& m);

	void resize(size_t dim);
	virtual ~SymMatrix() {}
	// Operators
	inline void operator=(const SymMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] = m.data[i];
	}

	inline void operator+=(const SymMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows())
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] += m.data[i];
	}

	inline void operator-=(const SymMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows()) // Assertion, will be removed
		{
			std::cout << "\n The sizes of the two matrices are not matched"
			          << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] -= m.data[i];
	}

	// Access the element
	virtual inline double& operator()(const size_t i, const size_t j = 1) const { return data[getArrayIndex(i, j)]; }
	// Access to members
	void LimitSize(size_t dim);

	// vec_result = This*vec. vec_result must be initialized
	virtual void multi(const double* vec, double* vec_result, double fac = 1.0);
	// m_result = this*m. m_result must be initialized. m_result must be a full stored matrix
	virtual void multi(const SymMatrix& m, Matrix& m_result, double fac = 1.0);
	// m_result = this*m1*m2. m_result must be initialized.  m_result must be a full stored matrix
	virtual void multi(const SymMatrix& m1, const Matrix& m2, Matrix& m_result);

	inline size_t getArrayIndex(const size_t i, const size_t j) const
	{
		if (i >= j)
			return static_cast<size_t>(i * (i + 1) / 2) + j;
		else
			return static_cast<size_t>(j * (j + 1) / 2) + i;
	}
};

class DiagonalMatrix : public MatrixBase
{
private:
	mutable double dummy_zero;

public:
	using MatrixBase::operator=;
	using MatrixBase::operator+=;
	using MatrixBase::operator-=;
	using MatrixBase::operator*=;
	using MatrixBase::operator/=;

	explicit DiagonalMatrix(size_t dim);
	DiagonalMatrix();
	DiagonalMatrix(const DiagonalMatrix& m);

	void resize(size_t dim);

	virtual ~DiagonalMatrix() {}
	// Forwarded operators
	inline void operator=(const DiagonalMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows() || ncols != m.Cols())
		{
			cout << "\n The sizes of the two matrices are not matched"
			     << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] = m.data[i];
	}

	inline void operator+=(const DiagonalMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows())
		{
			cout << "\n The sizes of the two matrices are not matched"
			     << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] += m.data[i];
	}

	inline void operator-=(const DiagonalMatrix& m)
	{
#ifdef gDEBUG
		if (nrows != m.Rows()) // Assertion, will be removed
		{
			cout << "\n The sizes of the two matrices are not matched"
			     << "\n";
			abort();
		}
#endif
		for (size_t i = 0; i < size; i++)
			data[i] -= m.data[i];
	}

	void LimitSize(size_t dim);

	// Access to members
	inline double& operator()(const size_t i, const size_t j) const
	{
#ifdef gDEBUG
		if (i >= nrows || j >= nrows)
		{
			cout << "\n Index exceeds the size of the matrix"
			     << "\n";
			abort();
		}
#endif
		if (i == j)
			return data[i]; // temporary
		else
			return dummy_zero;
	}

	inline double& operator()(const size_t i) const
	{
#ifdef gDEBUG
		if (i >= size)
		{
			cout << "\n Index exceeds the size of the matrix"
			     << "\n";
			abort();
		}
#endif

		return data[i];
	}
};

#ifdef _MSC_VER
#pragma warning(default : 4522)
#endif

typedef Matrix Vec;

/*========================================================================
   GeoSys - class my_vector (Declaration)
   Task:       Carry out vector operation
   Function:   See the declaration below
   programming:
   05/2005  WW
   ==========================================================================*/
template <class T>
class vec
{
public:
	vec(int argSize);
	vec() : _size(0), _entry(NULL) {}
	explicit vec(const vec<T>& v);

	virtual ~vec();
	// Operator
	virtual void operator=(T v)
	{
		for (size_t i = 0; i < _size; i++)
			_entry[i] = v;
	}
	virtual void operator=(const vec<T>&);
	virtual void resize(int newh);
	virtual T& operator[](size_t i) { return (T&)_entry[i]; }
	virtual const T& operator[](size_t i) const { return (const T&)_entry[i]; }
	virtual size_t Size() const { return _size; }
	T* Entry() { return _entry; }
	T* Entry() const { return _entry; }
	virtual void Write(std::ostream& os = std::cout) const;

protected:
	size_t _size;
	T* _entry;
};

template <>
class vec<void*>
{
public:
	vec(int argSize);
	vec() : _entry(NULL), _size(0) {}
	explicit vec(const vec<void*>& v);

	virtual ~vec();
	// Operator
	void operator=(void* v)
	{
		for (size_t i = 0; i < _size; i++)
			_entry[i] = v;
	}
	void operator=(const vec<void*>& v);
	void*& operator[](size_t i) { return _entry[i]; }
	const void*& operator[](size_t i) const { return (const void*&)_entry[i]; }
	// Access to memebers
	void** Entry() { return _entry; }
	const void** Entry() const { return (const void**)_entry; }
	virtual void resize(int newh);
	virtual size_t Size() const { return _size; }
	virtual void Write(std::ostream& os = std::cout) const;

protected:
	void** _entry;
	size_t _size;
};

template <class T>
class vec<T*> : public vec<void*>
{
public:
	vec(int Size) : vec<void*>(Size) {}
	vec() : vec<void*>() {}
	explicit vec(const vec<T*>& v) : vec<void*>(v) {}
	~vec() {}
	// Operator
	void operator=(T* v)
	{
		for (size_t i = 0; i < _size; i++)
			_entry[i] = v;
	}
	void operator=(const vec<T*>& v);
	T*& operator[](size_t i) { return (T*&)_entry[i]; }
	const T*& operator[](size_t i) const { return (const T*&)_entry[i]; }
	T** Entry() { return _entry; }
	T** Entry() const { return (const T**)_entry; }
};

#ifdef NEW_EQS
//
/// Sparse matrix storage type //04.2011. WW
enum StorageType
{
	CRS,
	JDS
};
class SparseTable
{
public:
	SparseTable(MeshLib::CFEMesh* a_mesh, bool quadratic, bool symm = false, StorageType stype = JDS);
	SparseTable(CPARDomain& m_dom, bool quadratic, bool symm = false);
	~SparseTable();
	void Write(std::ostream& os = std::cout);

private:
	bool symmetry;
	// Topology mapping from data array to matrix
	long* entry_column;
	long* num_column_entries; // number of entries of each columns in sparse table
	long* row_index_mapping_n2o; // Row index of sparse table to row index of matrix
	long* row_index_mapping_o2n; // Inverse of last
	long* diag_entry; // Global index to the index of  entry_column
	long size_entry_column;
	long max_columns;
	long rows;
	StorageType storage_type; // 04.2011. WW
	friend class CSparseMatrix;
};
// 08.2007 WW
// Jagged Diagonal Storage
class CSparseMatrix
{
public:
	CSparseMatrix(const SparseTable& sparse_table, const int dof);
	~CSparseMatrix();
	// Preconditioner
	void Precond_Jacobi(double* vec_s, double* vec_r);
	// TEMP
	void Precond_ILU(double* /*vec_s*/, double* /*vec_r*/) {}
	// Operator
	void operator=(const double a);
	void operator*=(const double a);
	void operator+=(const double a);
	void operator=(const CSparseMatrix& m);
	void operator+=(const CSparseMatrix& m);
	void operator-=(const CSparseMatrix& m);
	// Vector pass through augment and bring results back.
	void multiVec(double* vec_s, double* vec_r);
	void Trans_MultiVec(double* vec_s, double* vec_r);
	void Diagonize(const long idiag, const double b_given, double* b);
	//
	// Access to members
	double& operator()(const long i, const long j = 0) const;
	//
	StorageType GetStorageType() const { return storage_type; } // 05.2011. WW
	long Dim() const { return DOF * rows; }
	int Dof() const { return DOF; }
	void SetDOF(const int dof_n) //_new. 02/2010. WW
	{
		DOF = dof_n;
	}
	long Size() const { return rows; }
#if defined(LIS) || defined(MKL) // These two pointers are in need for Compressed Row Storage
	int nnz() const // PCH
	{
		return DOF * DOF * size_entry_column;
	}
	int* ptr;
	int* col_idx;
	int* entry_index;
	int GetCRSValue(double* value);
#endif
	// Print
	void Write(std::ostream& os = std::cout);
	void Write_BIN(std::ostream& os);
// Domain decomposition
#if defined(USE_MPI)
	void DiagonalEntries(double* diag_e);
#endif
private:
	// Data
	double* entry;
	mutable double zero_e;
	/// 0. 03.2011. WW
	StorageType storage_type;
	//
	bool symmetry;
	// Topology mapping from data array to matrix. All are only pointers to the
	// correpinding members in SparseTable, and no memory are allocated for them
	long* entry_column;
	long* num_column_entries; // number of entries of each columns in sparse table
	long* row_index_mapping_n2o; // Row index of sparse table to row index of matrix
	long* row_index_mapping_o2n; // Inverse of last
	long* diag_entry;
	long size_entry_column;
	long max_columns;
	long rows;
	//
	int DOF;
};
// Since the pointer to member funtions gives lower performance
#endif

//
// Cross production x^y. WW 12.01.2005
// const Vec& operator ^ (Vec& x,  Vec& y);

// End of class Matrix
}

//==========================================================================
#endif
