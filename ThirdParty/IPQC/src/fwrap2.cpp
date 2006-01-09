#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)


#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
IPQ_DLL_EXPORT int  ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  ADDERROR(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  ADDWARNING(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  CLEARACCUMULATEDLINES(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  DESTROYIPHREEQC(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  GETCURRENTSELECTEDOUTPUTUSERNUMBER(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void GETDUMPFILENAME(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETDUMPFILEON(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void GETDUMPSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGLINECOUNT(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGON(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void GETERRORFILENAME(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETERRORFILEON(int *id)
{
	return GetErrorFileOnF(id);
}
IPQ_DLL_EXPORT void GETERRORSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETERRORSTRINGLINECOUNT(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETERRORSTRINGON(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void GETLOGFILENAME(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETLOGFILEON(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void GETLOGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETLOGSTRINGON(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  GETLOGSTRINGLINECOUNT(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETNTHSELECTEDOUTPUTUSERNUMBER(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void GETOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETOUTPUTFILEON(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void GETOUTPUTSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETOUTPUTSTRINGLINECOUNT(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETOUTPUTSTRINGON(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTCOUNT(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void GETSELECTEDOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTFILEON(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void GETSELECTEDOUTPUTSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTSTRINGLINECOUNT(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTSTRINGON(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void GETVERSIONSTRING(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
IPQ_DLL_EXPORT void GETWARNINGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETWARNINGSTRINGLINECOUNT(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  LOADDATABASESTRING(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void OUTPUTACCUMULATEDLINES(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void OUTPUTERRORSTRING(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void OUTPUTWARNINGSTRING(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  SETBASICFORTRANCALLBACK(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  SETCURRENTSELECTEDOUTPUTUSERNUMBER(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  SETDUMPFILENAME(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETDUMPFILEON(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  SETERRORFILENAME(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETERRORFILEON(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  SETERRORSTRINGON(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  SETLOGFILENAME(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETLOGFILEON(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  SETLOGSTRINGON(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  SETOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETOUTPUTFILEON(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  SETOUTPUTSTRINGON(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTFILENAME(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTFILEON(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTSTRINGON(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}

#if defined(__cplusplus)
}
#endif

#endif