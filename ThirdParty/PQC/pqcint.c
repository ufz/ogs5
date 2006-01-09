#include "pqcint.h"
#include "phreeqc.h"

static char const svnid[] =
  "$Id: pqcint.c 13 2009-04-22 13:08:13Z delucia $";

// MDL:  entry function for the library, modified from the original main() function 
int
Phreeqcmain ( int argc, char* argv[], int nsim, int punch_dim, 
	      int inbuffer_dim, char** inbuffer, double* out)
{
  /* --------------- arguments -------------- */
  // MDL: argv must be an array of 5:
  // [0] ->  name of calling program {default: "libphreeqc"}
  // [1] ->  name of input file, irrelevant in this implementation {"FIXME"}
  // [2] ->  name of optional output file {"file_out.out"}
  // [3] ->  name of database file {"phreeqc.dat"}
  // [4] ->  flag to activate output print on file (if "T") or not (if everything else). 
  //         Defaults to "F", that is no file wile be printed.

  // nsim -> number of simulations in the current run
  // punch_dim -> number of punched entities
  // inbuffer_dim -> number of lines in the input
  // inbuffer -> char** buffer replacing the input file
  // out -> pointer to the array wich will contain the output (punch) of all simulations
  //        out must be allocated outside and has length nsim*punch_dim.

  int errors, i;
  int Perrors;
  void *db_cookie = NULL;
  void *input_cookie = NULL;
  char* dum = "T";

  // MDL: not defined on windows
  //#ifndef DOS
  // M: SVN_REV is defined in the Makefile
  // printf(":: libphreeqc rev %s :: \n", SVN_REV); 
  //  printf(":: SVN id: %s :: \n",svnid);
  //#endif

  // MDL: set "global" variables and local copies
  Pnsim = nsim;
  Ppunch_dim = punch_dim;
  Pbuff_dim = inbuffer_dim;
  Pout_dim = Pnsim * Ppunch_dim;
  PCountLine = 0;
  phast = FALSE;

  // MDL: prepare the vectors for output 
  Pout = (double *) malloc(Pout_dim * sizeof(double));
  if (Pout == NULL)
      malloc_error();

  Ppunch = (double *) malloc (Ppunch_dim * sizeof(double));
  if (Ppunch == NULL)
      malloc_error();

  // MDL: print/punch default
  Pfileprint = FALSE;
  Pdo_punch = FALSE;

  // MDL: Pbuff is the (global visible) pointer to array of arrays containing the input
  Pbuff = inbuffer;
  printf("::                  libPHREEQC\n");
  printf(":: called with %i simulations, output will be %i values\n", 
	 Pnsim, Pout_dim);

  if (strncmp(argv[4],dum, 1)==0) {
    Pfileprint = TRUE; 
    printf(":: Ouput written to file %s\n",argv[2]);
  } 

  
#ifdef MDL_DEBUG
  printf(": deb:: process_file_names con 1: %s, 2: %s, 3: %s, 4: %s\n",
	 argv[0],argv[1],argv[2],argv[3]);
#endif


/*
 *   Add callbacks for error_msg and warning_msg
 */
  if (add_output_callback (phreeqc_handler, NULL) != OK)
  {
    return(P_return_err(-1, "malloc"));
  }

/*
 *   Open input/output files
 */ 

#ifdef MDL_DEBUG
  printf(": deb :: Open io files\n");
#endif

  errors = process_file_names (argc, argv, &db_cookie, &input_cookie, FALSE);
  if (errors != 0)
    {
      clean_up ();
      return(P_return_err(errors, "process_file_names"));
    }
 
#ifdef MDL_DEBUG
  printf(": deb :: io files ok\n");
#endif

  /*
   *   Initialize arrays
   */
  errors = do_initialize ();
  if (errors != 0)
    {
      clean_up ();
      return(P_return_err(errors, "initialize"));
    }
  
  /*
   *   Load database into memory
   */ 
#ifdef MDL_DEBUG
  printf(": deb :: reading database from file ... \n");
#endif

  errors = read_database (getc_callback, db_cookie);
  if (errors != 0)
    {
      clean_up ();
      return(P_return_err(errors, "database"));
    }  
#ifdef MDL_DEBUG
  printf(": deb :: database loaded.\n");
#endif

/*
 *   Read input data for simulation
 */

  if (Pnsim > 1) 
    {
#ifdef MDL_DEBUG
      printf(": deb :: P_run_simulations (multiple) ...\n");
#endif
      Perrors = P_run_simulations (getc_callback, input_cookie);
      
#ifdef MDL_DEBUG
      printf(": deb :: P_run_simulations (multiple) OK, Perrors %i\n", Perrors);
#endif

/*
 *   Display successful status
 */
      errors = do_status ();
      if (errors != 0)
	{
	  clean_up ();
	  return(P_return_err(errors, "do_status"));
	}

#ifdef MDL_DEBUG
      printf(": deb :: do_status OK\n");
#endif


    }
  else
    { 
#ifdef MDL_DEBUG
      printf(": deb :: normal run_simulations ...\n");
#endif
      errors = run_simulations (getc_callback, input_cookie);
      if (errors != 0)
	{
	  clean_up ();
	  return(P_return_err(errors, "run_simulation"));
	}
      

#ifdef MDL_DEBUG
      printf(": deb :: normal run_simulations OK\n");
#endif

/*
 *   Display successful status
 */
      errors = do_status ();
      if (errors != 0)
	{
	  clean_up ();
	  return(P_return_err(errors, "do_status"));
	}
      
#ifdef MDL_DEBUG
      printf(": deb :: do_status OK\n");
#endif
    }


  clean_up ();
    
#ifdef MDL_DEBUG
  printf(": deb :: clean_up OK\n");
#endif

  close_input_files ();
#ifdef MDL_DEBUG
  printf(": deb :: close_input OK\n");
#endif


  close_output_files ();
#ifdef MDL_DEBUG
  printf(": deb :: close_output OK\n");
#endif

#ifdef MDL_DEBUG
  printf(": deb :: copy and set pointer for output\n");
#endif

  for (i=0;i<Pout_dim;i++)
    {
      out[i]=*(Pout+i);
      //      printf("punch_tot[%i]: %g\n",i,out[i]);
    }

  // MDL: deallocations
  free(Pout);
  free(Ppunch);

  printf(":: Seems ok, bye!\n\n");
  return(OK);
}


//
// MDL: P_return_err
//

int
P_return_err(int errors, char * tok)
{
  printf(":: BUM! phreeqc exits on error %s with number: %i\n",tok,errors);
  return(errors);
}
