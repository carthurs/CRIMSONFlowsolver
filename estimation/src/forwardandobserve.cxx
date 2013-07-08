#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE
#define VERDANDI_WITH_PETSC
#define VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
#define VERDANDI_OBSERVATION_ERROR_SPARSE

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS
//#define VERDANDI_WITH_MPI
//#define VERDANDI_WITH_OMP

//#define VERDANDI_LOGGING_LEVEL -7

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

#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"
#include "seldon/vector/PetscVector.cxx"
#include "seldon/matrix/PetscMatrix.cxx"

#include "SimvascularVerdandiModel.cxx"
#include "SimvascularObservationManager.cxx"
#include "method/ReducedOrderUnscentedKalmanFilter.cxx"

#include "ROUKFModified.cxx"

using namespace Verdandi;
using namespace std;

static char help[] = "ROUKF driver.\n\n";

int main(int argc, char** argv)
{

	TRY;

	if (argc < 2) {
		string mesg = "Usage:\n";
		mesg +=
				string("  ") + argv[0]
						+ " Nproc [assimilation method configuration file]";
		cout << mesg << endl;
		return 1;
	}
	else if (argc > 2) {
		static volatile int debuggerPresent = 0;
		while (!debuggerPresent)
			; // assign debuggerPresent=1
	}

	PetscInitialize(&argc, &argv, (char *)0, help);

    std::string current_directory = getenv("PWD");
    std::string input_filename = argv[1];

    cout << input_filename.substr(input_filename.find_last_of("\\/")+1,input_filename.length()-1) << endl;

    ROUKFModified<double,SimvascularVerdandiModel, SimvascularObservationManager> driver;

    driver.Initialize(argv[1], true);

    driver.GetObservationManager().SetTime(driver.GetModel(),driver.GetModel().GetTime());
    driver.GetObservationManager().SaveObservationSingleLocal(driver.GetModel().GetState());

    while (!driver.HasFinished())
    {
        driver.InitializeStep();



    	driver.GetModel().Forward();



        driver.GetObservationManager().SetTime(driver.GetModel(),driver.GetModel().GetTime());
        driver.GetObservationManager().SaveObservationSingleLocal(driver.GetModel().GetState());

    }

    driver.GetModel().Finalize();

    END;

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}