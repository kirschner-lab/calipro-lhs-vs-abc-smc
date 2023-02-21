/* Using SUNDIALS CVODES.

   To learn how to use the SUNDIALS CVODES API, read the bundled
   examples:

   source /path/to/your/spack/git/directory/
   spack cd -i sundials
   cd examples/cvodes/serial/

   The function usage is online:
   https://sundials.readthedocs.io/en/latest/cvodes/Usage/SIM.html
   https://sundials.readthedocs.io/en/latest/nvectors/NVector_package_links.html#nvector-functions-used-by-cvode
   https://sundials.readthedocs.io/en/latest/sunnonlinsol/SUNNonlinSol_API_link.html
   https://sundials.readthedocs.io/en/latest/sunnonlinsol/SUNNonlinSol_package_links.html#cvodes-sunnonlinearsolver-interface

   ... as well as the constants:
   https://sundials.readthedocs.io/en/latest/cvodes/Constants_link.html
*/

#include <stdio.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_version.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#define CHECK_RETURN(expected, function_call)		\
  for (int actual = function_call;			\
       expected != actual;)				\
    {							\
      fprintf(stderr,					\
	      __FILE__ ":%d:%s "			\
	     "expected %d but got %d from %s\n",	\
	     __LINE__, __func__,			\
	     expected, actual,				\
	     #function_call);				\
      return(1);					\
    }							\

#define CHECK_POINTER(pointer)				\
  if (pointer == NULL)					\
    {							\
      fprintf(stderr,					\
	      __FILE__ ":%d:%s "			\
	      "failed to allocate memory for %s\n",	\
	      __LINE__, __func__,			\
	      #pointer);				\
      return(1);					\
    }							\

/* SUNDIALS simulation context. */
SUNContext sunctx;

/* ODE equations wrapper passed to CVodes. */
static int dydt(realtype t, N_Vector y, N_Vector dy, void *userdata) {
  /* Suppress compiler warning of unused function argument. */
  (void) t;
  /* Cast parameters to realtype. */
  realtype *parameters = (realtype *) userdata;
  /* Reference individual parameters from the pointer. */
  realtype a = parameters[0];
  realtype b = parameters[1];
  realtype c = parameters[2];
  realtype d = parameters[3];
  /* Reference the current solution values. */
  realtype u = NV_Ith_S(y, 0);
  realtype v = NV_Ith_S(y, 1);
  /* Equation RHS. */
  NV_Ith_S(dy, 0) = a*u - b*u*v;
  NV_Ith_S(dy, 1) = -c*v + d*b*u*v;
  /* Return success. */
  return 0;
}

typedef struct {
  realtype a, b, c, d;
} userdata_t;

int main()
{
  /* Logger object. */
  SUNLogger logger;

  /* Solver settings and objects. */
  sunindextype num_equations = 2;
  N_Vector y = NULL;		/* Solution vector. */
  void *cvode_mem = NULL;
  SUNLinearSolver ls = NULL;
  SUNNonlinearSolver nls = NULL;
  realtype reltol = RCONST(1.0e-6);
  realtype abstol = RCONST(1.0e-10);
  realtype endtol = 1.0e-15;
 
  /* Parameters. */
  userdata_t parameters = {
    1.0,			/* a */
    0.1,			/* b */
    1.5,			/* c */
    0.75			/* d */
  };
  /* Initial values. */
  realtype y0_guess[] = {
    RCONST(10.0),
    RCONST(5.0)
  };
  realtype t0 = RCONST(0.0);
  realtype tn = RCONST(15.0);
  realtype dt = RCONST((tn - t0) / 100);

  /* -----------------------------------------------------------------
     Initialization.
     -----------------------------------------------------------------

     Note that some return codes are stored in the first element of
     the initialized empty objects!
  */

  /* Logging. */
  CHECK_RETURN(CV_SUCCESS, SUNLogger_Create(NULL, -1, &logger));
  CHECK_RETURN(CV_SUCCESS, SUNLogger_SetErrorFilename(logger, "stderr"));
  CHECK_RETURN(CV_SUCCESS, SUNLogger_SetWarningFilename(logger, "stderr"));
  CHECK_RETURN(CV_SUCCESS, SUNLogger_SetInfoFilename(logger, "stdout"));
  CHECK_RETURN(CV_SUCCESS, SUNLogger_SetDebugFilename(logger, "stdout"));

  /* Version. */
  char version[25];
  CHECK_RETURN(CV_SUCCESS, SUNDIALSGetVersion(version, 25));
  printf("SUNDIALS version: %s\n", version);
  printf("SUNDIALS logging level: %d\n", SUNDIALS_LOGGING_LEVEL);

  /* SUNDIALS simulation context that all SUNDIALS objects require. */
  CHECK_RETURN(CV_SUCCESS, SUNContext_Create(NULL, &sunctx));
  CHECK_RETURN(CV_SUCCESS, SUNContext_SetLogger(sunctx, logger));

  /* Solution vector. */
  y = N_VMake_Serial(num_equations, y0_guess, sunctx);
  CHECK_POINTER(y);

  /* Integrator algorithm, memory, and tolerance settings. */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  CHECK_POINTER(cvode_mem);
  CHECK_RETURN(CV_SUCCESS, CVodeSetUserData(cvode_mem, &parameters));
  CHECK_RETURN(CV_SUCCESS, CVodeInit(cvode_mem, dydt, t0, y));
  CHECK_RETURN(CV_SUCCESS, CVodeSStolerances(cvode_mem, reltol, abstol));

  /* Nonlinear fixed point solver. */
  nls = SUNNonlinSol_FixedPoint(y, 0, sunctx);
  CHECK_POINTER(nls);
  CHECK_RETURN(CV_SUCCESS, CVodeSetNonlinearSolver(cvode_mem, nls));
  CHECK_RETURN(CV_SUCCESS, SUNNonlinSolInitialize(nls));

  /* -----------------------------------------------------------------
     Problem.
     -----------------------------------------------------------------
  */
  printf("\nInitial conditions guessed:\n");
  for (int i = 0; i < num_equations; ++i) {
    printf("%f\n", y0_guess[i]);
  }
  printf("\n");

  /* -----------------------------------------------------------------
     Solution.
     -----------------------------------------------------------------

     Solver loop.
  */
  realtype t = t0;
  realtype ti = t0 + dt;
  for (int ret_cvode = 0;
       tn - t > endtol && ret_cvode >= 0;
       ) {
    /* Integrate. */
    ret_cvode = CVode(cvode_mem, ti, y, &t, CV_NORMAL);
    if (ret_cvode >= 0) {	/* Success. */
      ti += dt;
      ti = (ti > tn) ? tn : ti;
    } else {			/* Failure. */
      fprintf(stderr, "Solver failure: stopping integration.\n");
    }
  }

  /* Save solver output. */
  printf("Initial conditions solved:\n");
  realtype *y0 = N_VGetArrayPointer(y);
  for (int i = 0; i < num_equations; ++i) {
    printf("%f\n", y0[i]);
  }
  printf("\nSolver statistics:\n");
  CHECK_RETURN(CV_SUCCESS, CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE));

  /* -----------------------------------------------------------------
     Clean up memory.
     -----------------------------------------------------------------
  */
  SUNNonlinSolFree(nls);
  SUNLinSolFree(ls);
  CVodeFree(&cvode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);
  SUNLogger_Flush(logger, SUN_LOGLEVEL_ALL);
  SUNLogger_Destroy(&logger);

  return 0;
}
