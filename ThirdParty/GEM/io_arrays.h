//-------------------------------------------------------------------
// $Id: io_arrays.h 761 2012-11-30 09:12:20Z dmitrieva $
/// \file io_arrays.h
/// Various service functions for writing/reading arrays in files
//
// Copyright (C) 2006-2012 S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include  <fstream>

#include "verror.h"

struct outField /// Internal descriptions of fields
 {
   gstring name; ///< name of field in structure
   long int alws;    ///< 1 - must be read, 0 - default values can be used
   long int readed;  ///< 0; set to 1 after reading the field from input file
   long int indexation;  ///< 1 - static object; 0 - undefined; <0 type of indexation, >1 number of elements in array
   gstring comment;

};

class TRWArrays  /// Basic class for red/write fields of structure
 {
 protected:
    fstream& ff;
    long int numFlds; ///< Size of array flds
    outField* flds;   ///< Array of permissible fields

 public:

    /// Constructor
	 TRWArrays( short aNumFlds, outField* aFlds, fstream& fin ):
    	ff( fin ), numFlds(aNumFlds), flds(aFlds)
    {}

    /// Find field by name
    virtual  long int findFld( const char *Name );

     /// Set the data object can be skipped from the file
     /// and default value(s) can be used
     /// \param ii index in array flds
    void  setNoAlws( long int ii )
    {  flds[ii].alws = 0; }

    /// Set the data object can be skipped from the file
    /// and default value(s) can be used
    /// \param Name of field in array flds
    void  setNoAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setNoAlws(ii);
    }

    /// Set the data object  must be always present in the file
    /// \param ii index in array flds
    void  setAlws( long int ii )
    {  flds[ii].alws = 1; }

    /// Set the data object Name must be always present in the file
    /// \param Name of field in array flds
    void  setAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setAlws(ii);
    }

    /// Test the data object  must be always present in the file
    /// \param ii index in array flds
    bool  getAlws( long int ii )
    {  return (flds[ii].alws == 1); }

    /// Test the data object Name must be always present in the file
    /// \param Name of field in array flds
    bool  getAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
           return getAlws(ii);
         else
        return false;	 
    }

};

/// Print fields of structure outField
class TPrintArrays: public  TRWArrays
{

public:

    /*inline*/ void writeValue(float val);
    /*inline*/ void writeValue(double val);
    /*inline*/ void writeValue(long val);

    /// Writes long field to a text file.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeField(long f_num, long value, bool with_comments, bool brief_mode  );

    /// Writes short field to a text file.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeField(long f_num, short value, bool with_comments, bool brief_mode  );

    /// Writes char field to a text file.
    /// <flds[f_num].name> 'value'
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeField(long f_num, char value, bool with_comments, bool brief_mode  );

    /// Writes double field to a text file.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeField(long f_num, double value, bool with_comments, bool brief_mode  );

    /// Writes array to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num,  double* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false );

    /// Writes array to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num, long* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false );

    /// Writes array to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num, short* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false );

    /// Writes char array to a text file.
    /// <flds[f_num].name> "arr[0]" ... "arr[size-1]"
    /// \param l_size - Setup number of characters in one element
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArrayF( long f_num, char* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false );


    /// Constructor
    TPrintArrays( short aNumFlds, outField* aFlds, fstream& fout ):
        TRWArrays( aNumFlds, aFlds, fout)
    {}

    /// Writes char array to a text file.
    void writeArray( const char *name, char*   arr, long int size, long int arr_size );
    /// Writes float array to a text file.
    void writeArray( const char *name, float*  arr, long int size, long int l_size=-1L );
    /// Writes double array to a text file.
    void writeArray( const char *name, double* arr, long int size, long int l_size=-1L );
    /// Writes long array to a text file.
    void writeArray( const char *name, long* arr, long int size, long int l_size=-1L  );

    /// Writes char array to a text file.
    void writeArray( const char *name, char*   arr, int size, int arr_size );
    /// Writes float array to a text file.
    void writeArray( const char *name, float*  arr, int size, int l_size=-1 );
    /// Writes double array to a text file.
    void writeArray( const char *name, double* arr, int size, int l_size=-1 );
    /// Writes short array to a text file.
    void writeArray( const char *name, short* arr, int size, int l_size=-1  );

    /// Writes selected elements from float array to a text file.
    void writeArray( const char *name, float*  arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );
    /// Writes selected elements from double array to a text file.
    void writeArray( const char *name, double* arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );
    /// Writes selected elements from long array to a text file.
    void writeArray( const char *name, long* arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );

    /// Writes selected elements from float array to a text file.
    void writeArray( const char *name, float*  arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );
    /// Writes selected elements from double array to a text file.
    void writeArray( const char *name, double* arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );
    /// Writes selected elements from short array to a text file.
    void writeArray( const char *name, short* arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );

};


 class TReadArrays : public  TRWArrays /// Read fields of structure
 {
    gstring curArray;

 protected:
    /// Reads value from a text file.
    inline void readValue(float& val);
    /// Reads value from a text file.
    inline void readValue(double& val);
    inline void setCurrentArray( const char* name, long int size );
 
 public:

    /// Constructor
    TReadArrays( short aNumFlds, outField* aFlds, fstream& fin ):
        TRWArrays( aNumFlds, aFlds, fin ), curArray("")
    {}

    void  skipSpace();
    void reset();  ///< Reset to 0 all flags (readed)

    long int findFld( const char *Name ); ///< Find field by name
    long int findNext();  ///< Read next name from file and find in fields list
    void  readNext( const char* label); ///< Read next name from file

    gstring testRead();   ///< Test for reading all fields must be always present in the file

    /// Reads array from a text file.
    void readArray( const char *name, short* arr, long int size );
    /// Reads array from a text file.
    void readArray( const char *name, int* arr, long int size );
    /// Reads array from a text file.
    void readArray( const char *name, long int* arr, long int size );
    /// Reads array from a text file.
    void readArray( const char *name, float* arr, long int size );
    /// Reads array from a text file.
    void readArray( const char *name, double* arr, long int size );
    /// Reads array from a text file.
    void readArray( const char *name, char* arr, long int size, long int el_size );

};

//=============================================================================
