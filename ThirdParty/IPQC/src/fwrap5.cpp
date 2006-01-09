#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif


//
// Intel Fortran compiler 9.1 /iface:default /names:default /assume:underscore
//
IPQ_DLL_EXPORT int  ACCUMULATELINE_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  ADDERROR_(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  ADDWARNING_(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  CLEARACCUMULATEDLINES_(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  CREATEIPHREEQC_(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  DESTROYIPHREEQC_(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void GETCOMPONENT_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETCOMPONENTCOUNT_(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  GETCURRENTSELECTEDOUTPUTUSERNUMBER_(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void GETDUMPFILENAME_(int *id, char *filename, unsigned int len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETDUMPFILEON_(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void GETDUMPSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGLINECOUNT_(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGON_(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void GETERRORFILENAME_(int *id, char *filename, unsigned int len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETERRORFILEON_(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void GETERRORSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETERRORSTRINGLINECOUNT_(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETERRORSTRINGON_(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void GETLOGFILENAM_E(int *id, char *filename, unsigned int len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETLOGFILEON_(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void GETLOGSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETLOGSTRINGLINECOUNT_(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETLOGSTRINGON_(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  GETNTHSELECTEDOUTPUTUSERNUMBER_(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void GETOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETOUTPUTFILEON_(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void GETOUTPUTSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETOUTPUTSTRINGLINECOUNT_(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETOUTPUTSTRINGON_(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTCOLUMNCOUNT_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTCOUNT_(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void GETSELECTEDOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTFILEO_N(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTROWCOUNT_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void GETSELECTEDOUTPUTSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTSTRINGLINECOUNT_(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTSTRINGON_(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTVALUE_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void GETVERSIONSTRING_(char* version, unsigned int version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void GETWARNINGSTRINGLINE_(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETWARNINGSTRINGLINECOUNT_(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  LOADDATABASE_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  LOADDATABASESTRING_(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void OUTPUTACCUMULATEDLINES_(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void OUTPUTERRORSTRING_(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void OUTPUTWARNINGSTRING_(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  RUNACCUMULATED_(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  RUNFILE_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  RUNSTRING_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  SETBASICFORTRANCALLBACK_(int *id, double (*fcn)(double *x1, double *x2, char *str, int l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  SETCURRENTSELECTEDOUTPUTUSERNUMBER_(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  SETDUMPFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETDUMPFILEON_(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  SETDUMPSTRINGON_(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  SETERRORFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETERRORFILEON_(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  SETERRORSTRINGON_(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  SETLOGFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETLOGFILEON_(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  SETLOGSTRINGON_(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  SETOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETOUTPUTFILEON_(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  SETOUTPUTSTRINGON_(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTFILENAME_(int *id, char *filename, unsigned int len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTFILEON_(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  SETSELECTEDOUTPUTSTRINGON_(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}

#if defined(__cplusplus)
}
#endif

#endif
