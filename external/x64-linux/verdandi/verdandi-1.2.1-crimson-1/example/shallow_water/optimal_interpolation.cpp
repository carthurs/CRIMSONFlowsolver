// You may use 'LinearObservationManager' or 'GridToNetworkObservationManager'
// as observation operator.
#ifndef OBSERVATION_OPERATOR
#define OBSERVATION_OPERATOR LinearObservationManager
#endif


#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU

#define VERDANDI_DENSE

#include "Verdandi.hxx"

#include "seldon/SeldonSolver.hxx"

#include "model/ShallowWater.cxx"
#include "method/OptimalInterpolation.cxx"

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#include QUOTE(OBSERVATION_OPERATOR.cxx)


int main(int argc, char** argv)
{

    VERDANDI_TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        std::cout << mesg << std::endl;
        return 1;
    }

    typedef double real;

    Verdandi::OptimalInterpolation<real, Verdandi::ShallowWater<real>,
        Verdandi::OBSERVATION_OPERATOR<real> > driver;

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.Analyze();
        driver.FinalizeStep();
    }

    driver.Finalize();

    VERDANDI_END;

    return 0;

}
