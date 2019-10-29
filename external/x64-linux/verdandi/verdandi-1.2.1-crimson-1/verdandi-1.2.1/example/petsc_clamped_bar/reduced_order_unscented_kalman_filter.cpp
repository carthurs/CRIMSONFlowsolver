#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE

#define VERDANDI_WITH_PETSC

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS
//#define VERDANDI_WITH_MPI
//#define VERDANDI_WITH_OMP

//#define VERDANDI_LOGGING_LEVEL -7

#if defined(VERDANDI_WITH_MPI)
#include <mpi.h>
#endif

#if defined(VERDANDI_WITH_OMP)
#include <omp.h>
#endif

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"


#include "seldon/vector/PetscVector.cxx"
#include "seldon/matrix/PetscMatrix.cxx"

#include "model/PetscClampedBar.cxx"
#include "observation_manager/PetscLinearObservationManager.cxx"
#include "method/ReducedOrderUnscentedKalmanFilter.cxx"

static char help[] = "ROUKF driver.\n\n";

using namespace std;

int main(int argc, char** argv)
{

    TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        std::cout << mesg << std::endl;
        return 1;
    }

    PetscInitialize(&argc, &argv, (char *)0, help);

    typedef double real;

    Verdandi::ReducedOrderUnscentedKalmanFilter<real,
        Verdandi::PetscClampedBar<real>,
        Verdandi::PetscLinearObservationManager<real> > driver;
    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.Analyze();
    }

    driver.Finalize();

    END;

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}
