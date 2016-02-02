#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define SELDON_WITH_MKL

#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE
#define VERDANDI_WITH_PETSC
#define VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
#define VERDANDI_OBSERVATION_ERROR_SPARSE

#define VERDANDI_ROUKF_PARALLEL_INNOVATION

//#define VERDANDI_ROUKF_DEBUG_OUTPUT

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS
//#define VERDANDI_WITH_MPI
//#define VERDANDI_WITH_OMP

//#define VERDANDI_LOGGING_LEVEL -7
#include "Python.h"

#include <fstream>
#include <iomanip>

//#if defined(VERDANDI_WITH_MPI)
//#include <mpi.h>
//#endif
//
//#if defined(VERDANDI_WITH_OMP)
//#include <omp.h>
//#endif

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

#if defined(SELDON_WITH_MKL)
#include "mkl_service.h"
#endif

#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"
#include "seldon/vector/PetscVector.cxx"
#include "seldon/matrix/PetscMatrix.cxx"

#include "SimvascularVerdandiModel.cxx"
#include "SimvascularObservationManager.cxx"
#include "method/ReducedOrderUnscentedKalmanFilter.cxx"

#include "ROUKFModified.cxx"
#include "autoGeneratedVersionNumber.hxx"
#include "dateTools.hxx"
#include "CRIMSONPython.hxx"
#include "common_c.h"

using namespace Verdandi;

static char help[] = "ROUKF driver.\n\n";

int main(int argc, char** argv)
{
    kalmanFilterActive.kalmanFilterOn = true;
    char buildNumber[100];
    char buildTime[100];
    getBuildNumber(buildNumber);
    getBuildTime(buildTime);
    std::cout << "This is Simvascular version " << buildNumber << ", built at " << buildTime << "." << std::endl;

    // Expiry date check (uncomment enableExpiryDate() call below to enable):
    expiryDate expiry = expiryDate();
    expiry.setExpiryDayOfMonth(14);
    expiry.setExpiryMonthOfYear(11);
    expiry.setExpiryYear(2014);
    // UNCOMMENT TO DO A BUILD WITH AN EXPIRY DATE!
    // expiry.enableExpiryDate();
    expiry.checkWhetherExpiryDatePassed();

	TRY;

	if (argc < 2) {
		string mesg = "Usage:\n";
		mesg +=
				string("  ") + argv[0]
						+ " Nproc [assimilation method configuration file]";
		cout << mesg << endl;
		return 1;
	}
	else if (argc >= 2) {
           // Debugger snare:
           const char* debuggerFlag = "1";
           for(int ii=1; ii<argc; ii++)
           {
             // Look for a single "1" on the command line, indicating that we should
             // wait for the debugger...
             if(!strcmp(argv[ii], debuggerFlag))
             {
                 static volatile int debuggerPresent =0;
                 std::cout << "Debug flag spotted on the command line. Pausing to await debugger connection..." << std::endl;
                 while (!debuggerPresent ); // assign debuggerPresent=1
             }
           }
	}

#if defined(SELDON_WITH_MKL)
	mkl_cbwr_set(MKL_CBWR_AUTO);
#endif

	PetscInitialize(&argc, &argv, (char *)0, help);
    
    // Initialise Python integration
    initialisePython();

    std::string current_directory = getenv("PWD");
    std::string input_filename = argv[1];

    cout << input_filename.substr(input_filename.find_last_of("\\/")+1,input_filename.length()-1) << endl;

    ROUKFModified<double,SimvascularVerdandiModel, SimvascularObservationManager> driver;

    driver.Initialize(argv[1], true);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();

        driver.Forward();
        driver.Analyze();

        driver.FinalizeStep();
    }

    driver.Finalize();

    END;

    // Terminate Python integration
    Py_Finalize();

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}
