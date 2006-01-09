#define EXTERNAL
#include "global.h"
#include "output.h"
#include "phrqproto.h"
#include "input.h"
#include <stdlib.h>
#include "phqalloc.h"

// MDL: P_run_simulation
int P_run_simulations (PFN_READ_CALLBACK pfn, void *cookie);
int P_return_err(int errors, char * tok);
