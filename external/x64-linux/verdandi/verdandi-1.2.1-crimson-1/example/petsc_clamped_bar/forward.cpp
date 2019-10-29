#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE

#define VERDANDI_WITH_PETSC

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"

#include "seldon/vector/PetscVector.cxx"
#include "seldon/matrix/PetscMatrix.cxx"

#include "model/PetscClampedBar.cxx"
#include "method/ForwardDriver.cxx"

using namespace std;

static char help[] = "Forward driver.\n\n";

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

    Verdandi::ForwardDriver<Verdandi::PetscClampedBar<double> >
        driver;

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    driver.Finalize();

    END;

    int ierr;
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;

}
