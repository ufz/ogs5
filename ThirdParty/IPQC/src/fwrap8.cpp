#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32) && !defined(_M_AMD64)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:stdref /names:lowercase
//
IPQ_DLL_EXPORT int  __stdcall accumulateline(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  __stdcall adderror(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall addwarning(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  __stdcall clearaccumulatedlines(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  __stdcall createiphreeqc(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  __stdcall destroyiphreeqc(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void __stdcall getcomponent(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getcomponentcount(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getcurrentselectedoutputusernumber(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void __stdcall getdumpfilename(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall getdumpfileon(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void __stdcall getdumpstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getdumpstringlinecount(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getdumpstringon(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall geterrorfilename(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall geterrorfileon(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void __stdcall geterrorstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall geterrorstringlinecount(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall geterrorstringon(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void __stdcall getlogfilename(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall getlogfileon(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void __stdcall getlogstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getlogstringlinecount(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getlogstringon(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall getnthselectedoutputusernumber(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void __stdcall getoutputfilename(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall getoutputfileon(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void __stdcall getoutputstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getoutputstringlinecount(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getoutputstringon(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputcolumncount(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputcount(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void __stdcall getselectedoutputfilename(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputfileon(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputrowcount(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void __stdcall getselectedoutputstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputstringlinecount(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputstringon(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  __stdcall getselectedoutputvalue(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void __stdcall getversionstring(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void __stdcall getwarningstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  __stdcall getwarningstringlinecount(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  __stdcall loaddatabase(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall loaddatabasestring(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void __stdcall outputaccumulatedlines(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void __stdcall outputerrorstring(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void __stdcall outputwarningstring(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  __stdcall runaccumulated(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  __stdcall runfile(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall runstring(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  __stdcall setbasicfortrancallback(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  __stdcall setcurrentselectedoutputusernumber(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  __stdcall setdumpfilename(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall setdumpfileon(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  __stdcall setdumpstringon(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  __stdcall seterrorfilename(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall seterrorfileon(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  __stdcall seterrorstringon(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  __stdcall setlogfilename(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall setlogfileon(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall setlogstringon(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  __stdcall setoutputfilename(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall setoutputfileon(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall setoutputstringon(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  __stdcall setselectedoutputfilename(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  __stdcall setselectedoutputfileon(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  __stdcall setselectedoutputstringon(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}
#if defined(__cplusplus)
}
#endif

#endif // _WIN32

