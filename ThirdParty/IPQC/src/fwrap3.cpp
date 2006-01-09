#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cref /assume:underscore
//
IPQ_DLL_EXPORT int  accumulateline_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  adderror_(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  addwarning_(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  clearaccumulatedlines_(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  createiphreeqc_(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  destroyiphreeqc_(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void getcomponent_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getcomponentcount_(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  getcurrentselectedoutputusernumber_(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void getdumpfilename_(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getdumpfileon_(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void getdumpstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getdumpstringlinecount_(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getdumpstringon_(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void geterrorfilename_(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  geterrorfileon_(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void geterrorstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  geterrorstringlinecount_(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  geterrorstringon_(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void getlogfilename_(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getlogfileon_(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void getlogstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getlogstringlinecount_(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getlogstringon_(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  getnthselectedoutputusernumber_(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void getoutputfilename_(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getoutputfileon_(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void getoutputstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getoutputstringlinecount_(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getoutputstringon_(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputcolumncount_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputcount_(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void getselectedoutputfilename_(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  getselectedoutputfileon_(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputrowcount_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void getselectedoutputstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getselectedoutputstringlinecount_(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputstringon_(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  getselectedoutputvalue_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void getversionstring_(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void getwarningstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  getwarningstringlinecount_(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  loaddatabase_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  loaddatabasestring_(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void outputaccumulatedlines_(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void outputerrorstring_(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void outputwarningstring_(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  runaccumulated_(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  runfile_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  runstring_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  setbasicfortrancallback_(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  setcurrentselectedoutputusernumber_(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  setdumpfilename_(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setdumpfileon_(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  setdumpstringon_(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  seterrorfilename_(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  seterrorfileon_(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  seterrorstringon_(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  setlogfilename_(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setlogfileon_(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  setlogstringon_(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  setoutputfilename_(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setoutputfileon_(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  setoutputstringon_(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  setselectedoutputfilename_(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  setselectedoutputfileon_(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  setselectedoutputstringon_(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}

#if defined(__cplusplus)
}
#endif

#endif
