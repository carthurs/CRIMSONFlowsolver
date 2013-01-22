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

//    Verdandi::ReducedOrderUnscentedKalmanFilter<double,
//        SimvascularVerdandiModel,
//        SimvascularObservationManager > driver;

	ROUKFModified<double,SimvascularVerdandiModel, SimvascularObservationManager> driver;

	SimvascularVerdandiModel::reduced_state_error_variance Pred;

    driver.Initialize(argv[1], true);

    // debugging output

    std::ostringstream ostr; //output string stream
    std::string filename = "Pred_";
    ostr << driver.GetModel().GetRank();
    filename = filename+ostr.str();
    ofstream outfile;

    if (driver.GetModel().GetRank() == 0)
    	outfile.open(filename.c_str());


    while (!driver.HasFinished())
    {
        driver.InitializeStep();

//    	driver.GetObservationManager().SetTime(driver.GetModel(),driver.GetModel().GetTime());
//    	driver.GetObservationManager().SaveObservationSingleLocal(driver.GetModel().GetState());
//        driver.GetModel().Forward();

        driver.Forward();
        driver.Analyze();
        driver.GetModel().ForwardFinalize();

        driver.GetReducedStateErrorVariance(Pred);

        // write out the covariance matrix for the reduced state estimate
        if (driver.GetModel().GetRank() == 0) {
        	outfile << driver.GetModel().GetTime() << endl;
        	Pred.WriteText(outfile);
        }

    }

    driver.GetModel().Finalize();

    if (driver.GetModel().GetRank() == 0)
    	outfile.close();

    END;

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}
