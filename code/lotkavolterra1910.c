/* Using SUNDIALS CVODES.

   To learn how to use the SUNDIALS CVODES API, read the bundled
   examples:
   
   source /path/to/your/spack/git/directory/
   spack cd -i sundials
   cd examples/cvodes/serial/

   Also the function reference is online:

   https://sundials.readthedocs.io/en/latest/cvodes/Usage/SIM.html
*/

#include <stdio.h>

#include <cvodes/cvodes.h>           /*  */
#include <sundials/sundials_types.h> /* SUNOutputFormat */

#define CHECK(expected, function_call)			\
  int actual = function_call;				\
  if (expected != actual)				\
    {							\
      printf(__FILE__ ":%d:%s "				\
	     "expected %d but got %d from %s\n",	\
	     __LINE__, __func__,			\
	     expected, actual,				\
	     #function_call);				\
      return(1);					\
    }							\

/* SUNDIALS simulation context */
static SUNContext sunctx;

/* ODE equations wrapper passed to CVodes. */
static int dydt(realtype t, N_Vector dy, void *data);

int main()
{
  /* Parameters. */

  /* Variables. */

  /* SUNDIALS simulation context that all SUNDIALS objects require. */
  CHECK(CV_SUCCESS, SUNContext_Create(NULL, &sunctx));

  /* Save output. */
  /* SUN_OUTPUTFORMAT_CSV; */

  /* Clean up memory. */
  SUNContext_Free(&sunctx);

  return 0;
}
