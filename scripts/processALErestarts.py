#
# TO DO: UPDATE THIS TO WORK WITH THE LATEST PYTHON CONFIGS !!
#


###########################################################################################################################################
import sys
# sys.path.append(r'C:\Users\kl12\Documents\Work\Code\CRIMSON\PythonModules') 
sys.path.append(r'/home/klau/dev/crimson/PythonModules') 
import vtk
from vtk.util import numpy_support
import PythonQtMock as PythonQt
sys.modules["PythonQt"] = PythonQt
import CRIMSONSolver
from CRIMSONSolver.SolverStudies import PhastaSolverIO, PhastaConfig
import numpy as np
from os import listdir
from re import findall
###########################################################################################################################################

print "Using VTK Version "+str(vtk.VTK_MAJOR_VERSION)+"."+str(vtk.VTK_MINOR_VERSION)+"."+str(vtk.VTK_BUILD_VERSION)

# set step number from arguments 
stepNumber = sys.argv[1]

# read geombc
geombcFilename = 'geombc.dat.1'
geombcConfig = PhastaConfig.geombcConfig
geombcReader = PhastaSolverIO.PhastaFileReader(geombcFilename, geombcConfig)

# read restart
restartFilename = 'restart.'+stepNumber+'.0'   
print restartFilename
restartConfig = PhastaConfig.restartConfig
restartReader = PhastaSolverIO.PhastaFileReader(restartFilename, restartConfig)

# generate geombc dictionary
geombcDictionary = {}
for i in xrange(geombcReader.solutionStorage.getNArrays()):
    geombcDictionary[geombcReader.solutionStorage.getArrayName(i)] = i

# generate restart dictionary
restartDictionary = {}
for i in xrange(restartReader.solutionStorage.getNArrays()):
    restartDictionary[restartReader.solutionStorage.getArrayName(i)] = i

# load coordinates from restart and convert to vtk points
coordinates = restartReader.solutionStorage.getArrayData(restartDictionary['updated mesh coordinates']).copy()
numberOfNodes = coordinates.shape[0]/3
coordinates = np.reshape(coordinates, (numberOfNodes, 3), order='C')
points = vtk.vtkPoints()    

# load interior connectivity and convert to 0->n-1 indexing
interiorConnectivity = geombcReader.solutionStorage.getArrayData(geombcDictionary['connectivity interior linear tetrahedron']).copy()
numberOfInteriorElements = interiorConnectivity.shape[0]/4
interiorConnectivity = np.reshape(interiorConnectivity, (numberOfInteriorElements, 4), order='C')
interiorConnectivity -= 1 

# load pressure
pressure = restartReader.solutionStorage.getArrayData(restartDictionary['pressure']).copy()
pressure = np.reshape(pressure, (numberOfNodes, 1), order='C')
vtkPressure = vtk.vtkDoubleArray()
vtkPressure.SetNumberOfComponents(1)
vtkPressure.SetName("pressure")
vtkPressure = numpy_support.numpy_to_vtk(pressure)

# set coordinates
points.SetData(numpy_support.numpy_to_vtk(coordinates))

# set unstructured grid
unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.SetPoints(points)
unstructuredGrid.Allocate(numberOfInteriorElements, 1)
for i, tetrahedron in enumerate(interiorConnectivity):
  tetra = vtk.vtkTetra()
  for j in range(0,4):
    tetra.GetPointIds().SetId(j, tetrahedron[j])
  unstructuredGrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

        

# if (unstructuredGrid.GetPointData().GetArray("pressure") != None):
#   unstructuredGrid.GetPointData().GetArray("pressure").Reset()
      
unstructuredGrid.GetPointData().AddArray(vtkPressure)
unstructuredGrid.GetPointData().GetArray(0).SetName("pressure")



# write out to file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("restart."+stepNumber+".vtu")
if vtk.VTK_MAJOR_VERSION < 6:
    writer.SetInput(unstructuredGrid)            
writer.Write()

