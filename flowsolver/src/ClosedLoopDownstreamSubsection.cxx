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


        PetscErrorCode errFlag;
        // Tile the matrices to make the full closed loop system matrix
        {
            errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatZeroEntries(m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // these variables will be used during tiling to mark where the top-left corner of
            // the next matrix that we tile into the system goes:
            int m_nextBlankSystemMatrixRow = 0;
            int m_nextBlankSystemMatrixColumn = 0;

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
                m_nextBlankSystemMatrixColumn += numberOfColumns;
            }

            // Scopign unit containing just the addition of the closed loop circuit matrix to the full system matrix, m_closedLoopSystemMatrix:
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

        // Add the Kirchoff laws for the connecting nodes, and also set their pressures to be equal
        appendKirchoffLawsToClosedLoopLinearSystem();

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

            // Create a vector to hold the RHS (m_closedLoopRHS)
            errFlag = VecCreate(PETSC_COMM_SELF,&m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);m_closedLoopRHS
            errFlag = VecSetType(m_closedLoopRHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
            errFlag = VecSetSizes(m_closedLoopRHS,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = VecZeroEntries(m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = VecAssemblyBegin(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = VecAssemblyEnd(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            assert(m_indicesOfFirstRowOfEachUpstreamContributionInClosedLoopMatrix.size() == 0);

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

                m_indicesOfFirstRowOfEachUpstreamContributionInClosedLoopMatrix.push_back(m_nextBlankRhsRow);
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
        // First, create a matrix to store the inverse, m_inverseOfClosedLoopMatrix:
        errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = MatZeroEntries(m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

        errFlag = MatMatSolve(m_systemMatrix,m_identityMatrixForPetscInversionHack,m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // Release the m_systemMatrix so we can edit it again on the next iteration (we only need the just-computed m_inverseOfClosedLoopMatrix for computations on this step now.)
        errFlag = MatSetUnfactored(m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        
        // Solve the system
        // Create a vector to hold the solution
        errFlag = VecCreate(PETSC_COMM_SELF,&m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = VecSetType(m_solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
        errFlag = VecSetSizes(m_solutionVector,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = VecZeroEntries(m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = VecAssemblyBegin(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = VecAssemblyEnd(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        
        errFlag = MatMult(m_inverseOfClosedLoopMatrix,m_closedLoopRHS,m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

        // Clean up
        errFlag = MatDestroy(&m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = MatDestroy(&m_closedLoopSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        errFlag = VecDestroy(&m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

        // Set the "done for this timestep" flag. Reset this with ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain() (from BC manager for now)
        m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = true;
    }
}

void ClosedLoopDownstreamSubsection::appendKirchoffLawsToClosedLoopLinearSystem()
{
    // Loop over the upstream boundary conditions which form part of this closed loop circuit:
    {
        auto upstreamBCCircuit = m_upstreamBoundaryConditionCircuits.begin();
        int upstreamCircuitIndex = 0;
        while (upstreamBCCircuit != m_upstreamBoundaryConditionCircuits.end())
        {
            // Get the multiple incident current nodes for this boundary
            boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBCCircuit);
            downcastCircuit->

            // Get the circuit data for this boundary (so we can check the sign for
            // the Kirchoff equation by seeing if a node is the start node or end
            // node of a component):
            boost::shared_ptr<CircuitData> circuitData = downcastCircuit->getCircuitDescription();

            // Get the 
            int columnOffsetOfCurrentUpstreamCircuit = m_indicesOfFirstRowOfEachUpstreamContributionInClosedLoopMatrix.at(upstreamCircuitIndex);
            // Loop the multiple incident current nodes:
            for ...
            {
                int numberOfHistoryPressures = downcastCircuit->getNumberOfHistoryPressures();

                // Write the Kirchoff equation
                writeKirchoffEquationIntoClosedLoopSysteMatrix(numberOfComponents, multipleIncidentCurrentNode, m_nextBlankSystemMatrixRow, numberOfHistoryPressures, columnOffsetOfCurrentUpstreamCircuit);
                

                // Move to the next available row in the matrix that isn't used for an equation yet:
                m_nextBlankSystemMatrixRow++;
            }

            upstreamCircuitIndex++;
            // Loop control code:
            upstreamBCCircuit++;
        }
    }
    // Finally, do the Kirchoff equations for multiple incident current nodes
    // which are internal to the closed loop downstream subsection (i.e. those
    // which do not lie at an interface with an upstream boundary condition circuit)
    m_closedLoopSystemMatrix
    
}
// std::vector<int>& closedLoopNodesWithMultipleIncidentCurrents
void ClosedLoopDownstreamSubsection::writeKirchoffEquationIntoClosedLoopSysteMatrix(const boost::shared_ptr<const CircuitData> circuitData, const int multipleIncidentCurrentNode, const int row, const int numberOfHistoryPressures, const int columnOffset)
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

    // Check whether the multiple incident current node is also
    // an interface node with the closed loop circuit. If so,
    // extend the Kirchoff equation to include this in the balance
    
}

void ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain()
{
    m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = false;
}

std::pair<double,double> ClosedLoopDownstreamSubsection::getImplicitCoefficients(const int boundaryConditionIndex) const
{
    // assert the linear system has been solved:
    assert(m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep);

    // Use the known offsets in the matrix of the boundary condition with index boundaryConditionIndex
    // to extract the implicit coefficients
}

void ClosedLoopDownstreamSubsection::createContiguousIntegerRange(const int startingInteger, const int numberOfIntegers, PetscInt* const arrayToFill)
{
    for (int ii = 0; ii < numberOfIntegers; ii++)
    {
        arrayToFill[ii] = startingInteger + ii;
    }
}
