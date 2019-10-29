#define VERDANDI_DEBUG_LEVEL_4

#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"

#include "model/Lorenz.cxx"
#include "method/ForwardDriver.cxx"


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

    Verdandi::ForwardDriver<Verdandi::Lorenz<double> > driver;

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    VERDANDI_END;

    return 0;
}
