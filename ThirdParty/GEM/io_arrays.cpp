//-------------------------------------------------------------------
// $Id: io_arrays.cpp 776 2012-12-16 18:46:33Z dmitrieva $
//
/// \file io_arrays.cpp
/// Implementation of service functions for writing/reading arrays in files
//
// Copyright (c) 2006-2012 S.Dmytriyeva
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

#include <iomanip>
#include  <iostream>

#include "io_arrays.h"
#include "v_user.h"

#ifdef IPMGEMPLUGIN

istream& f_getline(istream& is, gstring& str, char delim);

//    This constants should be 'defined' to satisfy all compilers
#define SHORT_EMPTY 	   -32768
#define LONG_EMPTY             -2147483648L
#define FLOAT_EMPTY	          1.17549435e-38F
#define DOUBLE_EMPTY         2.2250738585072014e-308
#define CHAR_EMPTY   	     '`'

inline bool IsFloatEmpty( const float v )
{
    return ( v>0. && v <= FLOAT_EMPTY);
}
inline bool IsDoubleEmpty( const double v )
{
    return ( v>0. && v <= DOUBLE_EMPTY);
}

#else

#include "v_vals.h"

#endif

long int TRWArrays::findFld( const char *Name )
  {
  	long int ii;

   for( ii=0; ii < numFlds; ii++ )
    if( !( strcmp( flds[ii].name.c_str(), Name ) ))
      return ii;
   return -1;
  }

//---------------------------------------------------------//
// print Arrays ( fields of structure )

/// Write long value to file
/*inline*/ void TPrintArrays::writeValue(long val)
    {

       ff << val << " ";
    }

/// Write float value to file
/*inline*/ void TPrintArrays::writeValue(float val)
    {
      if( IsFloatEmpty( val ))
        ff << CHAR_EMPTY << " ";
      else
    //    ff << setprecision(10) << scientific << arr[ii] << " ";
       ff << setprecision(7) << val << " ";
    }

/// Write double value to file
/*inline*/ void TPrintArrays::writeValue(double val)
    {
      if( IsDoubleEmpty( val ))
        ff << CHAR_EMPTY << " ";
      else
//    ff << setprecision(18) << scientific << arr[ii] << " ";
    ff << setprecision(15) << val << " ";
    }

 void TPrintArrays::writeField(long f_num, long value, bool with_comments, bool brief_mode  )
    {
      if(!brief_mode || getAlws( f_num ))
      {  if( with_comments && flds[f_num].comment.length()>1)
          ff << endl << flds[f_num].comment.c_str();
        ff << endl << "<" << flds[f_num].name.c_str() << ">  ";
        ff << /*left << setw(17)  <<*/  value;
      }
    }

 void TPrintArrays::writeField(long f_num, short value, bool with_comments, bool brief_mode  )
    {
      if(!brief_mode || getAlws( f_num ))
      {  if( with_comments && flds[f_num].comment.length()>1)
          ff <<  endl << flds[f_num].comment.c_str();
        ff << endl << "<" << flds[f_num].name.c_str() << ">  ";
        ff << /*left << setw(17)  <<*/  value;
      }
    }

 void TPrintArrays::writeField(long f_num, char value, bool with_comments, bool brief_mode  )
    {
      if(!brief_mode || getAlws( f_num ))
      {  if( with_comments && flds[f_num].comment.length()>1)
          ff <<  endl <<  flds[f_num].comment.c_str();
        ff << endl << "<" << flds[f_num].name.c_str() << ">  ";
        ff << "\'" << value << "\'";
      }
    }

 void TPrintArrays::writeField(long f_num, double value, bool with_comments, bool brief_mode  )
  {
     if(!brief_mode || getAlws(f_num ))
     { if( with_comments && flds[f_num].comment.length()>1 )
            ff <<  endl <<  flds[f_num].comment.c_str();
         ff << endl << "<" << flds[f_num].name.c_str() << ">  ";
         ff << /*left << setw(7) <<*/  value;
     }
  }

 void TPrintArrays::writeArray( long f_num,  double* arr,
                 long int size, long int l_size, bool with_comments, bool brief_mode )
 {

   if(!brief_mode || getAlws(f_num ))
   { if( with_comments )
          ff <<  endl << flds[f_num].comment.c_str();
     writeArray( flds[f_num].name.c_str(),  arr,size, l_size);
   }
 }

 void TPrintArrays::writeArray( long f_num, long* arr,
               long int size, long int l_size, bool with_comments, bool brief_mode  )
 {
     if(!brief_mode || getAlws(f_num ))
     { if( with_comments )
            ff <<  endl << flds[f_num].comment.c_str();
        writeArray( flds[f_num].name.c_str(),  arr,size, l_size);
     }
 }

 void TPrintArrays::writeArray( long f_num, short* arr,
               long int size, long int l_size, bool with_comments, bool brief_mode  )
 {
     if(!brief_mode || getAlws(f_num ))
     { if( with_comments )
            ff <<  endl << flds[f_num].comment.c_str();
        writeArray( flds[f_num].name.c_str(),  arr,size, l_size);
     }
 }

 void TPrintArrays::writeArrayF( long f_num, char* arr,
               long int size, long int l_size, bool with_comments, bool brief_mode  )
 {
     if(!brief_mode || getAlws(f_num ))
     { if( with_comments )
            ff <<  endl << flds[f_num].comment.c_str();
        writeArray( flds[f_num].name.c_str(),  arr,size, l_size);
     }
 }

/// If the first parameter is given as NULL then the char array
/// will be printed as a comment
void TPrintArrays::writeArray( const char *name, char* arr,
		long int size, long int arr_siz )
{
 bool isComment = false;

 if( name )
     ff << endl << "<" << name << ">" << endl;
 else
 { ff << endl << "#  ";
   isComment = true;
 }
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 40 )
    { jj=0;  ff << endl;
      if(isComment)
    	  ff << "#  ";
    }
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
#ifdef IPMGEMPLUGIN // 24/08/2010
    strip(str);
#else
    str.strip();
#endif
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void TPrintArrays::writeArray( const char *name, char* arr,
		int size, int arr_siz )
{
 bool isComment = false;

 if( name )
     ff << endl << "<" << name << ">" << endl;
 else
 { ff << endl << "#  ";
   isComment = true;
 }
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 40 )
    { jj=0;  ff << endl;
      if(isComment)
    	  ff << "#  ";
    }
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
#ifdef IPMGEMPLUGIN // 24/08/2010
    strip(str);
#else
    str.strip();
#endif
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		long int size, long int l_size )
{
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		int size, int l_size )
{
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		long int size, long int l_size )
{
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		int size, int l_size )
{
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name, long* arr, long int size, long int l_size  )
{
 long int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);//ff << arr[ii] << " ";
 }
}


void TPrintArrays::writeArray( const char *name, short* arr, int size, int l_size  )
{
 int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << arr[ii] << " ";
 }
}

//-------------------------------------------------------------------------

void TPrintArrays::writeArray( const char *name,  float* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
 if(!arr) return;
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
    for(long int cc=0; cc<nColumns; cc++ )
    {
    	if(jj == sz)
	    { jj=0;  ff << endl;}
         writeValue(arr[selArr[ii]*nColumns+cc]);
        jj++;
   }
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
 if(!arr) return;
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
   for( int cc=0; cc<nColumns; cc++ )
    {
    	if(jj == sz)
	    { jj=0;  ff << endl;}
       writeValue(arr[selArr[ii]*nColumns+cc]);
        jj++;
    }
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
 if(!arr) return;
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
    for(long int cc=0; cc<nColumns; cc++ )
        {
              if(jj == sz)
	        { jj=0;  ff << endl;}
              writeValue(arr[selArr[ii]*nColumns+cc]);
              jj++;
        }
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
 if(!arr) return;
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
     for( int cc=0; cc<nColumns; cc++ )
        {
           if(jj == sz)
           { jj=0;  ff << endl;}
             writeValue(arr[selArr[ii]*nColumns+cc]);
             jj++;
           }
 }
}


void TPrintArrays::writeArray( const char *name, long* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
  if(!arr) return;
  long int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
	for(long int cc=0; cc<nColumns; cc++ )
    {
		if(jj == sz)
        { jj=0;  ff << endl;}
    	ff << arr[selArr[ii]*nColumns+cc] << " ";
    	jj++;
    }
 }
}

void TPrintArrays::writeArray( const char *name, short* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
  if(!arr) return;
  int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
	for( int cc=0; cc<nColumns; cc++ )
    {
		if(jj == sz)
        { jj=0;  ff << endl;}
    	ff << arr[selArr[ii]*nColumns+cc] << " ";
    	jj++;
    }
 }
}

//-------------------------------------------------------------------------
//------------------------------------------------------------------

 inline void TReadArrays::readValue(float& val)
 {
   char input;
   skipSpace();

   ff.get( input );
   if( input == CHAR_EMPTY )
      val = FLOAT_EMPTY;
   else
      {
        ff.putback(input);
        ff >> val;
      }
  }

 inline void TReadArrays::setCurrentArray( const char* name, long int size )
 {
   char buf[200];
   sprintf( buf, "After successfully read <%s> %ld data items", name, size);
   curArray = buf;
 }

 inline void TReadArrays::readValue(double& val)
 {
   char input;
   skipSpace();

   ff.get( input );
   if( input == CHAR_EMPTY )
      val = DOUBLE_EMPTY;
   else
      {
        ff.putback(input);
        ff >> val;
      }
 }

// skip  ' ',  '\n', '\t' and comments (from '#' to end of line)
void  TReadArrays::skipSpace()
{
  char input;
  ff.get( input );
  while( input == '#' || input == ' ' ||
        input == '\n' || input == '\t')
 {
   if( input == '#' )
    do{
         ff.get( input );
      }while( input != '\n' && input != '\0' && !ff.eof());
   if( input == '\0' || ff.eof() )
     return;
   ff.get( input );
  }
 ff.putback(input);
}

void TReadArrays::reset()
{
 for(long int ii=0; ii < numFlds; ii++ )
    flds[ii].readed = 0;
}

long int TReadArrays::findFld( const char *Name )
{
	long int ii;
 gstring str = Name;
 size_t len = str.find('>');
 str = str.substr(0, len );

 for( ii=0; ii < numFlds; ii++ )
  if( !( strcmp( flds[ii].name.c_str(), str.c_str() ) ))
    return ii;
 return -1;
}

long int TReadArrays::findNext()
{
 char buf[200];
 skipSpace();

 if( ff.eof() )
   return -3;

 ff >> buf;

 if( !( memcmp( "END_DIM", buf+1, 7 )) )
  return -2;

 long int ii = findFld( buf+1 );
 if(  ii < 0 )
 {
    gstring msg = buf;
          msg += " - Invalid label of data.\n";
          msg += curArray;
     Error( "Formatted read error 01", msg );
 }

 flds[ii].readed = 1;
 return ii;
}


void TReadArrays::readNext( const char* label)
{
 char buf[200];
 gstring msg;
 skipSpace();

 if( ff.eof() )
 {
    msg = label;
    msg += " - No data where expected.\n";
    msg += curArray;
    Error( "Formatted read error 02", msg );
 }

 ff >> buf;
 gstring str = buf+1;
 size_t len = str.find('>');
 str = str.substr(0, len );

 if( !( strcmp( label, str.c_str() ) ))
     return;

 msg = buf;
 msg += " - Invalid label of data.\n";
 msg += curArray;
 Error( "Formatted read error 03", msg );

}

void TReadArrays::readArray( const char* name, short* arr, long int size )
{
 setCurrentArray( name, size);
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}

void TReadArrays::readArray( const char* name, int* arr, long int size )
{
 setCurrentArray( name, size);
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}

void TReadArrays::readArray( const char* name, long int* arr, long int size )
{
 setCurrentArray( name, size);
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}
void TReadArrays::readArray( const char* name, float* arr, long int size )
{
 setCurrentArray( name, size);
 for( long int ii=0; ii<size; ii++  )
    readValue(arr[ii]);
}

void TReadArrays::readArray( const char* name, double* arr, long int size )
{
 setCurrentArray( name, size);
 //ff << setprecision(15);
 for( long int ii=0; ii<size; ii++  )
    readValue(arr[ii]);
}

void TReadArrays::readArray( const char* name, char* arr, long int size, long int el_size )
{
 char ch;
 char buf[200];

 setCurrentArray( name, size);
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff.get(ch);
//   while( ff.good() && ch != '\'' )
//       ff.get(ch);
   ff.getline( buf, el_size+1, '\'');
   copyValues( arr +(ii*el_size), buf, el_size );
 }

}

gstring TReadArrays::testRead()
{
 gstring ret = "";
 for(long int ii=0; ii < numFlds; ii++ )
  if( flds[ii].alws==1 && flds[ii].readed != 1 )
  {  if( !ret.empty() )
       ret += ", ";
     ret += flds[ii].name;
  }
 return ret;
}

//=============================================================================
// io_arrays.cpp
