#include "CRIMSONPython.hxx"
#include <stdexcept>
#include "boost/filesystem.hpp"
#include <sstream>

void initialisePython()
{
  // char pySearchPath[] = "/usr/lib/Python2.7";
   char pySearchPath[] = "/usr";
   // const char* pySearchPath = std::getenv("PYTHONHOME");
   Py_Initialize();
   Py_SetPythonHome(pySearchPath);
   PyRun_SimpleString("import sys");
   
   char* crimsonFlowsolverHome;
   crimsonFlowsolverHome = getenv("CRIMSON_FLOWSOLVER_HOME");
   if (crimsonFlowsolverHome == NULL)
   {
     throw std::runtime_error("EE: Please set environmental variable CRIMSON_FLOWSOLVER_HOME to the root of the CRIMSON flowsolver source tree\n");
   }
   boost::filesystem::path crimsonFlowsolverHomePath(crimsonFlowsolverHome);
   if (!boost::filesystem::exists(crimsonFlowsolverHomePath))
   {
     throw std::runtime_error("EE: Error relating to environmental variable CRIMSON_FLOWSOLVER_HOME. Please check it is correctly set.\n");
   }

   // Construct a relative path with the location of the python flow control script we need:
   boost::filesystem::path pathOfCRIMONPythonLibraryRelativeToCrimsonFlowsolverHome("basicControlScripts/lib/");
   // Append to crimsonFlowsolverHomePath to get to the location of the python script we need:
   boost::filesystem::path pathToCRIMSONPythonLibrary = crimsonFlowsolverHomePath;
   pathToCRIMSONPythonLibrary /= pathOfCRIMONPythonLibraryRelativeToCrimsonFlowsolverHome;

   std::stringstream pythonCRIMSONimportString;
   pythonCRIMSONimportString << "sys.path.append('" << pathToCRIMSONPythonLibrary.string() << "')";
   PyRun_SimpleString(pythonCRIMSONimportString.str().c_str());
}

void safe_Py_DECREF(PyObject* objectToDecref)
{
  if (objectToDecref!=NULL)
  {
    Py_DECREF(objectToDecref);
  }
}