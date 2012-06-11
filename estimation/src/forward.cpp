#define VERDANDI_DEBUG_LEVEL_4

#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"
using namespace Verdandi;

#include "SimvascularVerdandiModel.cxx"
#include "method/ForwardDriver.cxx"

int main(int argc, char** argv)
{
    TRY;

    if (argc < 3)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + "Nproc [model configuration file] [assimilation method configuration file]";
        cout << mesg << endl;
		return 1;
	} /*else {
		static volatile int debuggerPresent = 0;
		while (!debuggerPresent)
			; // assign debuggerPresent=1
	}*/

    ForwardDriver<SimvascularVerdandiModel> driver;

    driver.GetModel().Initialize(argc, argv);
    driver.Initialize(argv[2], false);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    END;

    return 0;
}
