#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:default /names:lowercase
//
IPQ_DLL_EXPORT int  accumulateline(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  adderror(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  addwarning(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  clearaccumulatedlines(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  createiphreeqc(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  destroyiphreeqc(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void getcomponent(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getcomponentcount(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  getcurrentselectedoutputusernumber(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void getdumpfilename(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getdumpfileon(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void getdumpstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getdumpstringlinecount(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getdumpstringon(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void geterrorfilename(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  geterrorfileon(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void geterrorstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  geterrorstringlinecount(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  geterrorstringon(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void getlogfilename(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getlogfileon(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void getlogstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getlogstringlinecount(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getlogstringon(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  getnthselectedoutputusernumber(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void getoutputfilename(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getoutputfileon(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void getoutputstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getoutputstringlinecount(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getoutputstringon(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputcolumncount(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputcount(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void getselectedoutputfilename(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getselectedoutputfileon(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputrowcount(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void getselectedoutputstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getselectedoutputstringlinecount(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputstringon(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputvalue(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void getversionstring(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void getwarningstringline(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getwarningstringlinecount(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  loaddatabase(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  loaddatabasestring(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void outputaccumulatedlines(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void outputerrorstring(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void outputwarningstring(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  runaccumulated(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  runfile(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  runstring(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  setbasicfortrancallback(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  setcurrentselectedoutputusernumber(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  setdumpfilename(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setdumpfileon(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  setdumpstringon(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  seterrorfilename(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  seterrorfileon(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  seterrorstringon(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  setlogfilename(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setlogfileon(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  setlogstringon(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  setoutputfilename(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setoutputfileon(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  setoutputstringon(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  setselectedoutputfilename(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setselectedoutputfileon(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  setselectedoutputstringon(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}

#if defined(__cplusplus)
}
#endif

#endif
