#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* sprintf */
#include "phrqtype.h"
#include "IPhreeqc.h"

#include "fwrap.h"

char *
f2cstring(char* fstring, int len)
{
    char *cstr, *str;
    int  i;

    str = fstring;
    for (i = len - 1; i >= 0 && !isgraph((int)str[i]); i--);
    cstr = (char *) malloc((size_t) (i + 2));
    if (!cstr) return 0;

    cstr[i + 1] = '\0';
    memcpy(cstr,str,i+1);
    return cstr;
}

void
padfstring(char *dest, const char *src, unsigned int len)
{
    unsigned int sofar;

    for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
        *dest++ = *src++;

    while (sofar++ < len)
        *dest++ = ' ';
}

IPQ_RESULT
AccumulateLineF(int *id, char *line, unsigned int line_length)
{
	IPQ_RESULT n;
	char* cline;

	cline = f2cstring(line, line_length);
	if (!cline)
	{
		::AddError(*id, "AccumulateLine: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AccumulateLine(*id, cline);
	free(cline);
	return n;
}

int
AddErrorF(int *id, char *error_msg, unsigned int len)
{
	int n;
	char* cmsg;

	cmsg = f2cstring(error_msg, len);
	if (!cmsg)
	{
		::AddError(*id, "AddError: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AddError(*id, cmsg);
	free(cmsg);
	return n;
}

int
AddWarningF(int *id, char *warn_msg, unsigned int len)
{
	int n;
	char* cmsg;

	cmsg = f2cstring(warn_msg, len);
	if (!cmsg)
	{
		::AddError(*id, "AddWarning: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AddWarning(*id, cmsg);
	free(cmsg);
	return n;
}

IPQ_RESULT
ClearAccumulatedLinesF(int *id)
{
	return ::ClearAccumulatedLines(*id);
}

int
CreateIPhreeqcF(void)
{
	return ::CreateIPhreeqc();
}

int
DestroyIPhreeqcF(int *id)
{
	return ::DestroyIPhreeqc(*id);
}

int
GetComponentCountF(int *id)
{
	return ::GetComponentCount(*id);
}

void
GetComponentF(int *id, int *n, char* comp, unsigned int line_length)
{
	padfstring(comp, ::GetComponent(*id, (*n) - 1), line_length);
}

int
GetCurrentSelectedOutputUserNumberF(int *id)
{
	return ::GetCurrentSelectedOutputUserNumber(*id);
}

void
GetDumpFileNameF(int *id, char* fname, unsigned int fname_length)
{
	padfstring(fname, ::GetDumpFileName(*id), fname_length);
}

int
GetDumpFileOnF(int *id)
{
	return ::GetDumpFileOn(*id);
}

/*
GetDumpStringF
*/

int
GetDumpStringLineCountF(int *id)
{
	return ::GetDumpStringLineCount(*id);
}

void
GetDumpStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetDumpStringLine(*id, (*n) - 1), line_length);
}

int
GetDumpStringOnF(int *id)
{
	return ::GetDumpStringOn(*id);
}

void
GetErrorFileNameF(int *id, char* fname, unsigned int fname_length)
{
	padfstring(fname, ::GetErrorFileName(*id), fname_length);
}

int
GetErrorFileOnF(int *id)
{
	return ::GetErrorFileOn(*id);
}

/*
GetErrorStringF
*/

int
GetErrorStringLineCountF(int *id)
{
	return ::GetErrorStringLineCount(*id);
}

void
GetErrorStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetErrorStringLine(*id, (*n) - 1), line_length);
}

int
GetErrorStringOnF(int *id)
{
	return ::GetErrorStringOn(*id);
}

void 
GetLogFileNameF(int *id, char* fname, unsigned int fname_length)
{
	padfstring(fname, ::GetLogFileName(*id), fname_length);
}

int
GetLogFileOnF(int *id)
{
	return ::GetLogFileOn(*id);
}

int
GetLogStringLineCountF(int *id)
{
	return ::GetLogStringLineCount(*id);
}

void
GetLogStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetLogStringLine(*id, (*n) - 1), line_length);
}

int
GetLogStringOnF(int *id)
{
	return ::GetLogStringOn(*id);
}

int
GetNthSelectedOutputUserNumberF(int *id, int* n)
{
	return ::GetNthSelectedOutputUserNumber(*id, (*n) - 1);
}

void
GetOutputFileNameF(int *id, char* fname, unsigned int fname_length)
{
	padfstring(fname, ::GetOutputFileName(*id), fname_length);
}

int
GetOutputStringLineCountF(int *id)
{
	return ::GetOutputStringLineCount(*id);
}

void
GetOutputStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetOutputStringLine(*id, (*n) - 1), line_length);
}

int
GetOutputStringOnF(int *id)
{
	return ::GetOutputStringOn(*id);
}

int
GetOutputFileOnF(int *id)
{
	return ::GetOutputFileOn(*id);
}

int
GetSelectedOutputColumnCountF(int *id)
{
	return ::GetSelectedOutputColumnCount(*id);
}

int
GetSelectedOutputCountF(int *id)
{
	return ::GetSelectedOutputCount(*id);
}

void
GetSelectedOutputFileNameF(int *id, char* fname, unsigned int fname_length)
{
	padfstring(fname, ::GetSelectedOutputFileName(*id), fname_length);
}

int
GetSelectedOutputFileOnF(int *id)
{
	return ::GetSelectedOutputFileOn(*id);
}

/*
GetSelectedOutputStringF
*/

int
GetSelectedOutputStringLineCountF(int *id)
{
	return ::GetSelectedOutputStringLineCount(*id);
}

void
GetSelectedOutputStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetSelectedOutputStringLine(*id, (*n) - 1), line_length);
}

int
GetSelectedOutputStringOnF(int *id)
{
	return ::GetSelectedOutputStringOn(*id);
}

int
GetSelectedOutputRowCountF(int *id)
{
	int rows = ::GetSelectedOutputRowCount(*id);
	if (rows > 0)
	{
		rows -= 1;
	}
	return rows;
}

IPQ_RESULT
GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	IPQ_RESULT result;
	VAR v;
	VarInit(&v);
	char buffer[100];

	int adjcol = *col - 1;
	result = ::GetSelectedOutputValue(*id, *row, adjcol, &v);

	switch (v.type)
	{
	case TT_EMPTY:
		*vtype = v.type;
		break;
	case TT_ERROR:
		*vtype = v.type;
		break;
	case TT_LONG:
		*vtype = TT_DOUBLE;
		*dvalue = (double)v.lVal;
		::sprintf(buffer, "%ld", v.lVal);
		padfstring(svalue, buffer, svalue_length);
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		::sprintf(buffer, "%23.15e", v.dVal);
		padfstring(svalue, buffer, svalue_length);
		break;
	case TT_STRING:
		*vtype = v.type;
		padfstring(svalue, v.sVal, svalue_length);
		break;
	default:
		assert(0);
	}
	::VarClear(&v);
	return result;
}

void
GetVersionStringF(char* version, unsigned int version_length)
{
	padfstring(version, ::GetVersionString(), version_length);
}

/*
GetWarningStringF
*/

int
GetWarningStringLineCountF(int *id)
{
	return ::GetWarningStringLineCount(*id);
}

void
GetWarningStringLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetWarningStringLine(*id, (*n) - 1), line_length);
}

int
LoadDatabaseF(int *id, char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabase(*id, cfilename);
	free(cfilename);
	return n;
}

int
LoadDatabaseStringF(int *id, char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "LoadDatabaseString: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseString(*id, cinput);
	free(cinput);
	return n;
}

void
OutputAccumulatedLinesF(int *id)
{
	::OutputAccumulatedLines(*id);
}

void
OutputErrorStringF(int *id)
{
	::OutputErrorString(*id);
}

void
OutputWarningStringF(int *id)
{
	::OutputWarningString(*id);
}

int
RunAccumulatedF(int *id)
{
	return ::RunAccumulated(*id);
}

int
RunFileF(int *id, char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunFile(*id, cfilename);
	free(cfilename);
	return n;
}

int
RunStringF(int *id, char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "RunString: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunString(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetBasicFortranCallbackF(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return ::SetBasicFortranCallback(*id, fcn);
}

IPQ_RESULT
SetCurrentSelectedOutputUserNumberF(int *id, int *n)
{
	return ::SetCurrentSelectedOutputUserNumber(*id, *n);
}

IPQ_RESULT
SetDumpFileNameF(int *id, char* fname, unsigned int fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetDumpFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetDumpFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetDumpFileOnF(int *id, int* dump_on)
{
	return ::SetDumpFileOn(*id, *dump_on);
}

IPQ_RESULT
SetDumpStringOnF(int *id, int* dump_string_on)
{
	return ::SetDumpStringOn(*id, *dump_string_on);
}

IPQ_RESULT
SetErrorFileNameF(int *id, char* fname, unsigned int fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetErrorFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetErrorFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetErrorFileOnF(int *id, int* error_file_on)
{
	return ::SetErrorFileOn(*id, *error_file_on);
}

IPQ_RESULT
SetErrorStringOnF(int *id, int* error_string_on)
{
	return ::SetErrorStringOn(*id, *error_string_on);
}

IPQ_RESULT
SetLogFileNameF(int *id, char* fname, unsigned int fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetLogFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetLogFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetLogFileOnF(int *id, int* log_on)
{
	return ::SetLogFileOn(*id, *log_on);
}

IPQ_RESULT
SetLogStringOnF(int *id, int* log_string_on)
{
	return ::SetLogStringOn(*id, *log_string_on);
}

IPQ_RESULT
SetOutputFileNameF(int *id, char* fname, unsigned int fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetOutputFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetOutputFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetOutputFileOnF(int *id, int* output_on)
{
	return ::SetOutputFileOn(*id, *output_on);
}

IPQ_RESULT
SetOutputStringOnF(int *id, int* output_string_on)
{
	return ::SetOutputStringOn(*id, *output_string_on);
}

IPQ_RESULT
SetSelectedOutputFileNameF(int *id, char* fname, unsigned int fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetSelectedOutputFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetSelectedOutputFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetSelectedOutputFileOnF(int *id, int* sel_on)
{
	return ::SetSelectedOutputFileOn(*id, *sel_on);
}

IPQ_RESULT
SetSelectedOutputStringOnF(int *id, int* selected_output_string_on)
{
	return ::SetSelectedOutputStringOn(*id, *selected_output_string_on);
}

#if defined(_WIN32) && !defined(_M_AMD64)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
IPQ_DLL_EXPORT int  __stdcall ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  __stdcall ADDERROR(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall ADDWARNING(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall CLEARACCUMULATEDLINES(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  __stdcall CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  __stdcall DESTROYIPHREEQC(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void __stdcall GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETCURRENTSELECTEDOUTPUTUSERNUMBER(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void __stdcall GETDUMPFILENAME(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPFILEON(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void __stdcall GETDUMPSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPSTRINGLINECOUNT(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPSTRINGON(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall GETERRORFILENAME(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORFILEON(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void __stdcall GETERRORSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORSTRINGLINECOUNT(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORSTRINGON(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall GETLOGFILENAME(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGFILEON(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void __stdcall GETLOGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGSTRINGLINECOUNT(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGSTRINGON(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETNTHSELECTEDOUTPUTUSERNUMBER(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void __stdcall GETOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTFILEON(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void __stdcall GETOUTPUTSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTSTRINGLINECOUNT(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTSTRINGON(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTCOUNT(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void __stdcall GETSELECTEDOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTFILEON(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void __stdcall GETSELECTEDOUTPUTSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTSTRINGLINECOUNT(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTSTRINGON(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void __stdcall GETVERSIONSTRING(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void __stdcall GETWARNINGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETWARNINGSTRINGLINECOUNT(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall LOADDATABASESTRING(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTACCUMULATEDLINES(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTERRORSTRING(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTWARNINGSTRING(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  __stdcall RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  __stdcall RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  __stdcall SETBASICFORTRANCALLBACK(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  __stdcall SETCURRENTSELECTEDOUTPUTUSERNUMBER(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPFILENAME(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPFILEON(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORFILENAME(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORFILEON(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORSTRINGON(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGFILENAME(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGFILEON(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGSTRINGON(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTFILEON(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTSTRINGON(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTFILEON(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTSTRINGON(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}
#if defined(__cplusplus)
}
#endif

#endif // _WIN32

