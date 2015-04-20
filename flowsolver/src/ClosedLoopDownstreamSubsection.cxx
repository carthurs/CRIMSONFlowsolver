#include "ClosedLoopDownstreamSubsection.hxx"
#include <utility>

bool ClosedLoopDownstreamSubsection::boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const
{
	return (m_setOfAttachedBoundaryConditionIndices.count(boundaryConditionIndex) == 1);
}

void ClosedLoopDownstreamSubsection::setPointerToNeighbouringBoundaryConditionCircuit(boost::shared_ptr<NetlistCircuit> upstreamBCCircuit)
{
    m_upstreamBoundaryConditionCircuits.push_back(upstreamBCCircuit);
    // m_numberOfUpstreamCircuits gets updated every time we add an upstreamBCCircuit
    // (although its value won't be used - at the time of writing this comment - until all upstreamBCCircuits have been set.)
    m_numberOfUpstreamCircuits = m_upstreamBoundaryConditionCircuits.size();
}

void ClosedLoopDownstreamSubsection::initialiseModel()
{
    // Get the input data
    mp_NetlistCircuit->createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    // mp_NetlistCircuit->identifyAtomicSubcircuits();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    mp_NetlistCircuit->closeAllDiodes();

    mp_NetlistCircuit->initialiseCircuit();

    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    //
}

void ClosedLoopDownstreamSubsection::buildAndSolveLinearSystemIfNotYetDone()
{
    // Check whether the linear system still needs to be built and solved; if not, do nothing.
    if (!m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep)
    {
        terminatePetscArrays(); // This does nothing if the arrays don't exist. \todo consider refactoring the matrix creation/termination.
        // Call the upstream boundary conditions to ask for their contributions to the (closed loop)-type
        // linear system:
        for (auto upstreamBCCircuit = m_upstreamBoundaryConditionCircuits.begin(); upstreamBCCircuit != m_upstreamBoundaryConditionCircuits.end(); upstreamBCCircuit++)
        {
            Mat matrixContribution;
            Vec rhsContribuiton;

            boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBCCircuit);
            downcastCircuit->getMatrixContribution(matrixContribution);
            downcastCircuit->getRHSContribution(rhsContribuiton);
            m_matrixContributionsFromUpstreamBoundaryConditions.push(matrixContribution);
            m_rhsContributionsFromUpstreamBoundaryConditions.push(rhsContribuiton);

            // Get the offsets which tell us where we will find 
            downcastCircuit->getLocationsToLookForImplicitCoefficientInfoInLinearSystem();
            //\todo actually get and store these....

            m_systemSize += downcastCircuit->getNumberOfDegreesOfFreedom();
        }

        assert(m_systemSize > 0); // defensive
        
        // Get the final system size by adding in the number of degrees of freedom in the downstream closed loop subsection circuit.
        m_systemSize += mp_NetlistCircuit->getNumberOfDegreesOfFreedom();

        createVectorsAndMatricesForCircuitLinearSystem();


        PetscErrorCode errFlag;
        // Tile the matrices to make the full closed loop system matrix
        {
            // these variables will be used during tiling to mark where the top-left corner of
            // the next matrix that we tile into the system goes:
            int m_nextBlankSystemMatrixRow = 0;
            int m_nextBlankSystemMatrixColumn = 0;

            assert(m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.size() == 0);

            for (int upstreamCircuit = 0; upstreamCircuit < m_numberOfUpstreamCircuits; upstreamCircuit++)
            {
                Mat nextMatrixToAddToSystem = m_matrixContributionsFromUpstreamBoundaryConditions.pop();
                PetscInt numberOfRows;
                PetscInt numberOfColumns;
                errFlag = MatGetSize(nextMatrixToAddToSystem, &numberOfRows, &numberOfColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                // Extract the actual data array in the petsc matrix nextMatrixToAddToSystem, so we can pass it to MatSetValuesBlocked, for inclusion in m_closedLoopSystemMatrix:
                PetscScalar* rawDataInNextMatrixToAddToSystem;
                errFlag = MatGetArray(nextMatrixToAddToSystem, &rawDataInNextMatrixToAddToSystem); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                // Create location data for the rows and columns where we will place nextMatrixToAddToSystem in m_closedLoopSystemMatrix:
                PetscInt globalRowIndices[numberOfRows];
                createContiguousIntegerRange(m_nextBlankSystemMatrixRow, numberOfRows, globalRowIndices);
                PetscInt globalColumnIndices[numberOfColumns];
                createContiguousIntegerRange(m_nextBlankSystemMatrixColumn, numberOfColumns, globalColumnIndices);

                errFlag = MatSetValuesBlocked(m_closedLoopSystemMatrix, numberOfRows, globalRowIndices, numberOfColumns, globalColumnIndices, rawDataInNextMatrixToAddToSystem, INSERT_VALUES);

                // I'm not convinced this call is necessary given that we're about to destroy nextMatrixToAddToSystem anyway,
                // but the Petsc documentation says I MUST call it after MatGetArray once access to an array is no longer needed,
                // and who am I to argue with a block-capital imperative?
                errFlag = MatRestoreArray(nextMatrixToAddToSystem, rawDataInNextMatrixToAddToSystem); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                m_nextBlankSystemMatrixRow += numberOfRows;

                m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.push_back(m_nextBlankSystemMatrixColumn);
                m_nextBlankSystemMatrixColumn += numberOfColumns;
            }

            // This adds the location of the first column of the downstream closed loop circuit in m_closedLoopSystemMatrix.
            m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.push_back(m_nextBlankSystemMatrixColumn);

            // Scoping unit containing just the addition of the closed loop downstream circuit matrix to the full system matrix, m_closedLoopSystemMatrix:
            {
                // Add the closed loop downstream circuit's matrix:
                Mat matrixContribution;
                mp_NetlistCircuit->getMatrixContribution(matrixContribution);

                PetscInt numberOfRows;
                PetscInt numberOfColumns;
                errFlag = MatGetSize(matrixContribution, &numberOfRows, &numberOfColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                // Extract the actual data array in the petsc matrix matrixContribution, so we can pass it to MatSetValuesBlocked, for inclusion in m_closedLoopSystemMatrix:
                PetscScalar* rawDataInMatrix;
                errFlag = MatGetArray(matrixContribution, &rawDataInMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                // Create location data for the rows and columns where we will place matrixContribution in m_closedLoopSystemMatrix:
                PetscInt globalRowIndices[numberOfRows];
                createContiguousIntegerRange(m_nextBlankSystemMatrixRow, numberOfRows, globalRowIndices);
                PetscInt globalColumnIndices[numberOfColumns];
                createContiguousIntegerRange(m_nextBlankSystemMatrixColumn, numberOfColumns, globalColumnIndices);

                errFlag = MatSetValuesBlocked(m_closedLoopSystemMatrix, numberOfRows, globalRowIndices, numberOfColumns, globalColumnIndices, rawDataInMatrix, INSERT_VALUES);

                // as noted above, I'm not sure this is needed... but just in case...
                errFlag = MatRestoreArray(matrixContribution, rawDataInMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                m_nextBlankSystemMatrixRow += numberOfRows;
                m_nextBlankSystemMatrixColumn += numberOfColumns;
            }
        }

        // By this stage, the matrices for all the subcircuits of the closed loop
        // are tiled into the big matrix, m_closedLoopSystemMatrix. If this assert
        // fails, it means there are no rows left to actually add equations which
        // couple the circuit together.
        assert(m_nextBlankSystemMatrixRow < m_systemSize);

        assert(m_nextBlankSystemMatrixColumn < m_systemSize);

        // Add the Kirchoff laws for the connecting nodes
        appendKirchoffLawsAtInterfacesBetweenCircuits();

        // At the interface between an "upstream" boundary condition and the "downstream" 
        // closed loop circuit, the interfacing nodes are duplicated (as they belong to both
        // the upstream and the downstream circuit). We must enforce equality of pressure
        // at such copies of the same node. Do that now.
        enforcePressureEqualityBetweenDuplicatedNodes();

        // Finalise the matrix construction
        errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
        // std::cout << "System matrix for closed loop " << m_index << ":" << std::endl;
        //  errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

        // LU factor m_closedLoopSystemMatrix
        PetscErrorCode errFlag = MatLUFactor(m_closedLoopSystemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);


        // Tile the RHS contributions into our closed loop RHS:
        {
            // To track where the start of the next rhs contribution goes in m_closedLoopRHS:
            int m_nextBlankRhsRow = 0;

            // Zero out the RHS ready for the construction on this timestep:
            errFlag = VecZeroEntries(m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

            for (int upstreamCircuit = 0; upstreamCircuit < m_numberOfUpstreamCircuits; upstreamCircuit++)
            {
                Vec nextVectorToAddToClosedLoopRHS = m_rhsContributionsFromUpstreamBoundaryConditions.pop();
                PetscInt numberOfRows;
                errFlag = VecGetSize(nextVectorToAddToClosedLoopRHS, &numberOfRows); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                // Get the raw data to add to the RHS:
                PetscScalar* rawDataInNextVectorToAddToClosedLoopRHS;
                errFlag = VecGetArray(nextVectorToAddToClosedLoopRHS, &rawDataInNextVectorToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);
                
                PetscInt globalRowIndices[numberOfRows];
                createContiguousIntegerRange(m_nextBlankRhsRow, numberOfRows, globalRowIndices);
                errFlag = VecSetValues(m_closedLoopRHS, numberOfRows, globalRowIndices, rawDataInNextVectorToAddToClosedLoopRHS, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                errFlag = VecRestoreArray(nextVectorToAddToClosedLoopRHS, rawDataInNextVectorToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                m_nextBlankRhsRow += numberOfRows;
            }

            // Scoping unit just to add the downstream closed loop subsection to m_closedLoopRHS:
            {
                // Add the closed loop downstream circuit's RHS:
                Vec rhsContribuiton;
                mp_NetlistCircuit->getRHSContribution(rhsContribuiton);

                PetscInt numberOfRows;
                errFlag = VecGetSize(rhsContribuiton, &numberOfRows); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                // Get the raw data to add to the RHS:
                PetscScalar* rawDataToAddToClosedLoopRHS;
                errFlag = VecGetArray(rhsContribuiton, &rawDataToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                PetscInt globalRowIndices[numberOfRows];
                createContiguousIntegerRange(m_nextBlankRhsRow, numberOfRows, globalRowIndices);
                errFlag = VecSetValues(m_closedLoopRHS, numberOfRows, globalRowIndices, rawDataToAddToClosedLoopRHS, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                errFlag = VecRestoreArray(rhsContribuiton, rawDataToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

                m_nextBlankRhsRow += numberOfRows;
            }
        }

        // By this stage, the matrices for all the subcircuits of the closed loop
        // are tiled into the big matrix, m_closedLoopSystemMatrix. If this assert
        // fails, it means there are no rows left to actually add equations which
        // couple the circuit together.
        assert(m_nextBlankRhsRow < m_systemSize);

        // get the inverse of the system matrix:
        errFlag = MatMatSolve(m_systemMatrix,m_identityMatrixForPetscInversionHack,m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // Release the m_systemMatrix so we can edit it again on the next iteration (we only need the just-computed m_inverseOfClosedLoopMatrix for computations on this step now.)
        errFlag = MatSetUnfactored(m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        
        // Solve the system
        errFlag = MatMult(m_inverseOfClosedLoopMatrix,m_closedLoopRHS,m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);


        // Set the "done for this timestep" flag. Reset this with ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain() (from BC manager for now)
        m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = true;
    }
}

void ClosedLoopDownstreamSubsection::createVectorsAndMatricesForCircuitLinearSystem()
{
    // Create a vector to hold the solution
    errFlag = VecCreate(PETSC_COMM_SELF,&m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecSetType(m_solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
    errFlag = VecSetSizes(m_solutionVector,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecZeroEntries(m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyBegin(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create a vector to hold the RHS (m_closedLoopRHS)
    errFlag = VecCreate(PETSC_COMM_SELF,&m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);m_closedLoopRHS
    errFlag = VecSetType(m_closedLoopRHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
    errFlag = VecSetSizes(m_closedLoopRHS,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecZeroEntries(m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyBegin(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create the vector to hold the system matrix
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create a matrix to store the inverse, m_inverseOfClosedLoopMatrix:
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create an identity matrix for use when inverting the system matrix:
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // Fill the diagonal with ones:
    for (int ii=0; ii<m_systemSize; ii++)
    {
        errFlag = MatSetValue(m_identityMatrixForPetscInversionHack,ii,ii,1.0,INSERT_VALUES);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    errFlag = MatAssemblyBegin(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatAssemblyEnd(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void ClosedLoopDownstreamSubsection::initialisePetscArrayNames()
{
    m_closedLoopRHS = PETSC_NULL;
    m_solutionVector = PETSC_NULL;
    m_closedLoopSystemMatrix = PETSC_NULL;
    m_inverseOfClosedLoopMatrix = PETSC_NULL;
    m_identityMatrixForPetscInversionHack = PETSC_NULL;
}

void ClosedLoopDownstreamSubsection::terminatePetscArrays()
{
    PetscErrorCode errFlag;
    if (m_closedLoopRHS)
    {
        errFlag = VecDestroy(&m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_solutionVector)
    {
        errFlag = VecDestroy(&m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_closedLoopSystemMatrix)
    {
        errFlag = MatDestroy(&m_closedLoopSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_inverseOfClosedLoopMatrix)
    {
        errFlag = MatDestroy(&m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_identityMatrixForPetscInversionHack)
    {
        errFlag = MatDestroy(&m_identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
}

void ClosedLoopDownstreamSubsection::appendKirchoffLawsAtInterfacesBetweenCircuits()
{
    // Loop over the upstream boundary conditions which form part of this closed loop circuit:
    {
        auto upstreamBCCircuit = m_upstreamBoundaryConditionCircuits.begin();
        int upstreamCircuitIndex = 0;
        while (upstreamBCCircuit != m_upstreamBoundaryConditionCircuits.end())
        {
            boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBCCircuit);
            
            // Get the circuit data for this boundary (so we can check the sign for
            // the Kirchoff equation by seeing if a node is the start node or end
            // node of a component):
            boost::shared_ptr<CircuitData> circuitData = downcastCircuit->getCircuitDescription();

            // Get the deferred (i.e. interfacing-with-closed-loop-downstream-subsection) multiple incident current nodes for this boundary
            std::vector<int> multipleIncidentCurrentNodes = downcastCircuit->getNodesWithDeferredKirchoffEquations();

            // Get the column where the data block for the current upstream boundarycondition circuit
            // begins in the big matrix, m_closedLoopSystemMatrix
            int columnOffsetOfCurrentUpstreamCircuit = m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.at(upstreamCircuitIndex);
            // Loop the multiple incident current nodes:
            for (auto multipleIncidentCurrentNode = multipleIncidentCurrentNodes.begin(); multipleIncidentCurrentNode != multipleIncidentCurrentNodes.end(); multipleIncidentCurrentNode++)
            {
                int numberOfHistoryPressures_upstreamCircuit = downcastCircuit->getNumberOfHistoryPressures();

                // Write the part of the Kirchoff equation for this node (multipleIncidentCurrentNode)
                // which corresponds to the components incident at multipleIndidentCurrentNode in
                // the upstream boundary condition:
                writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(circuitData, multipleIncidentCurrentNode, m_nextBlankSystemMatrixRow, numberOfHistoryPressures_upstreamCircuit, columnOffsetOfCurrentUpstreamCircuit);

                // Write the part of the Kirchoff equation for this node (multipleIncidentCurrentNode)
                // which corresponds to the components incident at multipleIndidentCurrentNode in
                // the downstream closed loop subsection circuit:
                int numberOfHistoryPressures_downstreamCircuit
                int upstreamIndexOfMultipleIncidentCurrentNode = mp_NetlistCircuit->convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(multipleIncidentCurrentNode);
                int columnOffsetOfDownstreamClosedLoopCircuit = m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.back();
                writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(mp_NetlistCircuit, upstreamIndexOfMultipleIncidentCurrentNode, m_nextBlankSystemMatrixRow, numberOfHistoryPressures_downstreamCircuit, columnOffsetOfDownstreamClosedLoopCircuit);
                
                // Move to the next available row in the matrix that isn't used for an equation yet:
                m_nextBlankSystemMatrixRow++;
            }

            upstreamCircuitIndex++;
            // Loop control code:
            upstreamBCCircuit++;
        }
    }
}

void ClosedLoopDownstreamSubsection::writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(const boost::shared_ptr<const CircuitData> circuitData, const int multipleIncidentCurrentNode, const int row, const int numberOfHistoryPressures, const int columnOffset)
{
    // Do the equations for the nodes with multiple incident currents (just the part for the upstream boundary condition circuit first...
    // we'll append the currents for the components in the downstream closed loop circuit in a moment...)
    for (int ll=0; ll < circuitData->numberOfComponents; ll++)
    {
      int column = ll + circuitData->numberOfPressureNodes + numberOfHistoryPressures + columnOffset;
      bool foundMultipleIncidentCurrentsForEndNode = (circuitData->components.at(ll)->endNode->getIndex() == multipleIncidentCurrentNode); 
      if (foundMultipleIncidentCurrentsForEndNode)
      {
        errFlag = MatSetValue(m_closedLoopSystemMatrix,row,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      }

      bool foundMultipleIncidentCurrentsForStartNode = (circuitData->components.at(ll)->startNode->getIndex() == multipleIncidentCurrentNode);
      if (foundMultipleIncidentCurrentsForStartNode)
      {
        errFlag = MatSetValue(m_closedLoopSystemMatrix,row,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      }
    }
}

int ClosedLoopDownstreamSubsection::getCircuitIndexFromSurfaceIndex(const int upstreamSurfaceIndex) const
{
    for (auto upstreamBC = m_upstreamBoundaryConditionCircuits.begin(); upstreamBC != m_upstreamBoundaryConditionCircuits.end(); upstreamBC++)
    {
        if ((*upstreamBC)->surfaceIndexMatches(upstreamSurfaceIndex))
        {
            boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBC);
            return downcastCircuit->getCircuitIndex();
        }
    }
}

void ClosedLoopDownstreamSubsection::enforcePressureEqualityBetweenDuplicatedNodes()
{
    // Vectors to hold the pass-by-reference return from getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices
    std::vector<int> downstreamNodeIndices;
    std::vector<int> upstreamNodeIndices;
    std::vector<int> upstreamSurfaceIndices;
    mp_NetlistCircuit->getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(downstreamNodeIndices, upstreamNodeIndices, upstreamSurfaceIndices);

    // The upstreamSurfaceIndices (which are solver.inp indices for the surfaces of the 3D model) need to be converted to 
    // the indices of the circuits themselves (i.e. 0th, 1st, 2nd,... circuit; whereas the surfaces may have non-consecutive 
    // arbitrary numbering)
    std::vector<int> upstreamCircuitIndices;
    for (auto upstreamSurfaceIndex = upstreamSurfaceIndices.begin(); upstreamSurfaceIndex != upstreamSurfaceIndices.end(); upstreamCircuitIndex++)
    {
        int upstreamCircuitIndex = getCircuitIndexFromSurfaceIndex(*upstreamSurfaceIndex);
        upstreamCircuitIndices.push_back(upstreamCircuitIndex);
    }

    for (int sharedNodeIndex = 0; sharedNodeIndex < downstreamNodeIndices.size(); sharedNodeIndex++)
    {
        // The entry for the upstream copy of the pressure node:
        {
            int currentUpstreamCircuitIndex = upstreamCircuitIndices.at(sharedNodeIndex);
            int upstreamNodeIndex = upstreamNodeIndices.at(sharedNodeIndex);
            int column = m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.at(currentUpstreamCircuitIndex) + upstreamNodeIndex;
            errFlag = MatSetValue(m_closedLoopSystemMatrix,m_nextBlankSystemMatrixRow,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }

        // The entry for the downstream (i.e. local) copy of the pressure node:
        {
            int downstreamNodeIndex = downstreamNodeIndices.at(sharedNodeIndex);
            int column = m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.back() + downstreamNodeIndex;
            errFlag = MatSetValue(m_closedLoopSystemMatrix,m_nextBlankSystemMatrixRow,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }

        m_nextBlankSystemMatrixRow++;
    } 
}

void ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain()
{
    m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = false;
}

std::pair<double,double> ClosedLoopDownstreamSubsection::getImplicitCoefficients(const int boundaryConditionIndex) const
{
    // assert the linear system has been solved:
    assert(m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep);

    // The linear system is solved, so we can just extract the necessary values from the resulting solution
    // vector and inverted system matrix:
}

void ClosedLoopDownstreamSubsection::createContiguousIntegerRange(const int startingInteger, const int numberOfIntegers, PetscInt* const arrayToFill)
{
    for (int ii = 0; ii < numberOfIntegers; ii++)
    {
        arrayToFill[ii] = startingInteger + ii;
    }
}
