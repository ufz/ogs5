#ifndef LOGGING_INCLUDE
#define LOGGING_INCLUDE

#ifndef NERROR
#define PRINT_ERROR(x)
#else
#define PRINT_ERROR(x) \
  cout << "Error! : " << x ;  // Note: No endl, for flexibility
#endif

#ifndef NDEBUG
#define PRINT_DEBUG(x)
#else
#define PRINT_DEBUG(x) \
  cout << "Debug : " << #x << ":\t" << x << endl;
#endif

#ifndef NWARNING
#define PRINT_WARNING(x)
#else
#define PRINT_WARNING(x) \
  cout << "Warning! : " << x ;  // Note: No endl, for flexibility
#endif

#ifndef NINFO
#define PRINT_INFO(x)
#else
#define PRINT_INFO(x) \
  cout << "Info :"x ;  // Note: No endl, for flexibility
#endif

#endif // LOGGING_INCLUDE
