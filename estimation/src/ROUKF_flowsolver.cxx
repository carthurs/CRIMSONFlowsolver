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

#if defined(SELDON_WITH_MKL)
	mkl_cbwr_set(MKL_CBWR_AUTO);
#endif

	PetscInitialize(&argc, &argv, (char *)0, help);

	SimvascularVerdandiModel::reduced_state_error_variance Pred, obsGram;

    std::string current_directory = getenv("PWD");
    std::string input_filename = argv[1];

    cout << input_filename.substr(input_filename.find_last_of("\\/")+1,input_filename.length()-1) << endl;

    ROUKFModified<double,SimvascularVerdandiModel, SimvascularObservationManager> driver;

    driver.Initialize(argv[1], true);

    // debugging output
    //std::ostringstream ostr; //output string stream
    std::string Pfilename = "Pred",filename_ext = ".dat";
    std::string Gfilename = "Gram";
    //ostr << driver.GetModel().GetRank();
    Pfilename = Pfilename+filename_ext;
    Gfilename = Gfilename+filename_ext;
    ofstream Poutfile,Goutfile;

    if (driver.GetModel().GetRank() == driver.GetModel().GetNumProcs() - 1) {
    	Poutfile.open(Pfilename.c_str());
    	Goutfile.open(Gfilename.c_str());
    }

    driver.GetObservationManager().SetTime(driver.GetModel(),driver.GetModel().GetTime());
    driver.GetObservationManager().SaveObservationSingleLocal(driver.GetModel().GetState());


    while (!driver.HasFinished())
    {
        driver.InitializeStep();

        driver.Forward();
        driver.Analyze();

        driver.GetModel().ForwardFinalize();

        driver.GetReducedStateErrorVariance(Pred);

        driver.GetObservabilityGramian(obsGram);

        // write out the covariance matrix for the reduced state estimate
        if (driver.GetModel().GetRank() == driver.GetModel().GetNumProcs() - 1) {
        	Poutfile << driver.GetModel().GetTime() << endl;
        	Pred.WriteText(Poutfile);
        	Goutfile << driver.GetModel().GetTime() << endl;
        	obsGram.WriteText(Goutfile);
        }

        driver.GetObservationManager().SetTime(driver.GetModel(),driver.GetModel().GetTime());
        driver.GetObservationManager().SaveObservationSingleLocal(driver.GetModel().GetState());

    }

    driver.GetModel().Finalize();

    if (driver.GetModel().GetRank() == 0) {
    	Poutfile.close();
    	Goutfile.close();
    }
    END;

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}
