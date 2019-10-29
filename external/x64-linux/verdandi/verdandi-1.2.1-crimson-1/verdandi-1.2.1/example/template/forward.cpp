#define VERDANDI_DEBUG_LEVEL_4

#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"
using namespace Verdandi;

#include "model/ModelTemplate.cxx"
#include "method/ForwardDriver.cxx"


int main(int argc, char** argv)
{
    VERDANDI_TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        cout << mesg << endl;
        return 1;
    }

    ForwardDriver<ModelTemplate> driver;

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.Finalize();
    }

    driver.Finalize();

    VERDANDI_END;

    return 0;
}
