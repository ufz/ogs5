#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32) && !defined(_M_AMD64)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
// Intel Fortran compiler 9.1 /iface:stdref /names:uppercase /assume:underscore
//
IPQ_DLL_EXPORT int  __stdcall ACCUMULATELINE_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  __stdcall ADDERROR_(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall ADDWARNING_(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall CLEARACCUMULATEDLINES_(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  __stdcall CREATEIPHREEQC_(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  __stdcall DESTROYIPHREEQC_(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void __stdcall GETCOMPONENT_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETCOMPONENTCOUNT_(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETCURRENTSELECTEDOUTPUTUSERNUMBER_(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void __stdcall GETDUMPFILENAME_(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPFILEON_(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void __stdcall GETDUMPSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPSTRINGLINECOUNT_(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETDUMPSTRINGON_(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall GETERRORFILENAME_(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORFILEON_(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void __stdcall GETERRORSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORSTRINGLINECOUNT_(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETERRORSTRINGON_(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall GETLOGFILENAME_(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGFILEON_(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void __stdcall GETLOGSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGSTRINGLINECOUNT_(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETLOGSTRINGON_(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETNTHSELECTEDOUTPUTUSERNUMBER_(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void __stdcall GETOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTFILEON_(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void __stdcall GETOUTPUTSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTSTRINGLINECOUNT_(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETOUTPUTSTRINGON_(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTCOLUMNCOUNT_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTCOUNT_(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void __stdcall GETSELECTEDOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTFILEON_(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTROWCOUNT_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void __stdcall GETSELECTEDOUTPUTSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTSTRINGLINECOUNT_(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTSTRINGON_(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall GETSELECTEDOUTPUTVALUE_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void __stdcall GETVERSIONSTRING_(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void __stdcall GETWARNINGSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall GETWARNINGSTRINGLINECOUNT_(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall LOADDATABASE_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall LOADDATABASESTRING_(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTACCUMULATEDLINES_(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTERRORSTRING_(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void __stdcall OUTPUTWARNINGSTRING_(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  __stdcall RUNACCUMULATED_(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  __stdcall RUNFILE_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall RUNSTRING_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  __stdcall SETBASICFORTRANCALLBACK_(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  __stdcall SETCURRENTSELECTEDOUTPUTUSERNUMBER_(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPFILEON_(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  __stdcall SETDUMPSTRINGON_(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORFILEON_(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  __stdcall SETERRORSTRINGON_(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGFILEON_(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall SETLOGSTRINGON_(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTFILEON_(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall SETOUTPUTSTRINGON_(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTFILEON_(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  __stdcall SETSELECTEDOUTPUTSTRINGON_(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}
#if defined(__cplusplus)
}
#endif

#endif // _WIN32

