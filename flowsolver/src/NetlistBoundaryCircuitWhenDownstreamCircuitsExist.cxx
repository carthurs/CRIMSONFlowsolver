#include "NetlistBoundaryCircuitWhenDownstreamCircuitsExist.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::initialiseCircuit()
{
    // Discover which pressure nodes don't need Kirchoff law applications yet, because such laws will be applied later
    // once this circuit is combined with the downstreamCircuit to make a (closed loop)-type boundary circuit.
    m_pressureNodesWhichConnectToDownstreamCircuits = NetlistDownstreamCircuitReader::Instance()->getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(m_surfaceIndex);
    m_numberOfNodesConnectingToAnotherCircuit = m_pressureNodesWhichConnectToDownstreamCircuits.size();

    // This function exists just so we can modify what initialiseCircuit does in subclasses without repeating code.
    initialiseCircuit_common();

    m_numberOfSystemRows = m_numberOfSystemColumns - m_numberOfNodesConnectingToAnotherCircuit;

    createVectorsAndMatricesForCircuitLinearSystem();
}

bool NetlistBoundaryCircuitWhenDownstreamCircuitsExist::kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const
{
    bool nodeInterfacesWithDownstreamCircuit = false;
    if (m_pressureNodesWhichConnectToDownstreamCircuits.count(nodeIndex) == 1)
    {
        nodeInterfacesWithDownstreamCircuit = true;
    }
    return nodeInterfacesWithDownstreamCircuit;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::updateLPN(const int timestepNumber)
{
	// Call the downstream circuit(s?) and tell them to contact each NetlistBoundaryCondition to get
    // their contributions to the (closed loop)-type matrix, build the full matrix and solve it
    // 
    // Internally, that call will only do anything if this is the first boundary condition to try to computeImplicitCoefficients
    // this time. Otherwise, the system state has already been computed, and we don't need to do anything
	for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        downstreamSubcircuit->lock()->buildAndSolveLinearSystemForUpdateIfNotYetDone(timestepNumber, m_delt);
    }

    // Extract & get the m_solutionVector for just this circuit back from the big downstream
    // loop:
    assert(m_netlistDownstreamLoopClosingSubcircuits.size() <= 1); // \todo Actually, we can't have a given boundary circuit connecting to more than one downstream loop-closing circuit - otherwise which one do we ask for the solution vector here? Therefore, need to remove support for having more than one downstream loop-closing circuit per boundary.
    std::vector<PetscScalar> dataForSolutionVectorForThisSurface = m_netlistDownstreamLoopClosingSubcircuits.at(0).lock()->getSolutionVectorEntriesCorrespondingToSurface(m_surfaceIndex);
    // Put the raw data we just received into m_solutionVector, so that
    // the function calls below can find it where they expect it to be:
    PetscErrorCode errFlag;
    for (int vectorEntry = 0; vectorEntry < dataForSolutionVectorForThisSurface.size(); vectorEntry++)
    {
    	PetscScalar valueToInsert = dataForSolutionVectorForThisSurface.at(vectorEntry);
    	errFlag = VecSetValue(m_solutionVector,vectorEntry,valueToInsert,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }

    // Get the updated nodal pressures:
    giveNodesTheirPressuresFromSolutionVector();

    // Get the updated component flows:
    giveComponentsTheirFlowsFromSolutionVector();

    // Get the updated volumes:
    giveComponentsTheirVolumesFromSolutionVector();

    // Call the downstream subsections, and tell them to get their own parts of the full closed loop
    // solution vector, and extract the information to give to the components of the circuit
    // (pressure, flows, volumes)
    //
    // We don't need to do a solve in this case, because the solve has just been done above in this function.
    // \todo this is .at(0) because there should only ever be one of these (see the comment with the assert guardin this, above).
	m_netlistDownstreamLoopClosingSubcircuits.at(0).lock()->giveNodesAndComponentsTheirUpdatedValues();

}

std::pair<double,double> NetlistBoundaryCircuitWhenDownstreamCircuitsExist::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    // Call the downstream circuit(s?) and tell them to contact each NetlistBoundaryCondition to get
    // their contributions to the (closed loop)-type matrix, build the full matrix and solve it
    // 
    // Internally, that call will only do anything if this is the first boundary condition to try to computeImplicitCoefficients
    // this time. Otherwise, the system state has already been computed, and we don't need to do anything
    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        downstreamSubcircuit->lock()->buildAndSolveLinearSystemIfNotYetDone(timestepNumber, alfi_delt);
    }

    // Call the downstream circuits to get the implicit coefficients that they computed for this
    // boundary
    std::pair<double,double> returnValue;
    int counterToDetectErrors = 0;
    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        if (downstreamSubcircuit->lock()->boundaryConditionCircuitConnectsToThisDownstreamSubsection(m_surfaceIndex))
        {
            returnValue = downstreamSubcircuit->lock()->getImplicitCoefficients(m_IndexOfThisNetlistLPNInInputFile);
            counterToDetectErrors++;
        }
    }

    assert(counterToDetectErrors == 1);

    return returnValue;
}

std::pair<boundary_data_t,double> NetlistBoundaryCircuitWhenDownstreamCircuitsExist::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
	// Call the downstream circuit(s?) and tell them to contact each NetlistBoundaryCondition to get
    // their contributions to the (closed loop)-type matrix, build the full matrix and solve it
    // 
    // Internally, that call will only do anything if this is the first boundary condition to try to computeImplicitCoefficients
    // this time. Otherwise, the system state has already been computed, and we don't need to do anything
    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        downstreamSubcircuit->lock()->buildAndSolveLinearSystemIfNotYetDone(timestepNumber, m_delt);
    }

    // Call the downstream circuits to get the pressure of flow that they computed for this
	// boundary interface:
	std::pair<boundary_data_t,double> returnValue;
	if (mp_circuitData->hasPrescribedFlowAcrossInterface()) // This boundary condition is receiving flow and returning pressure (Neumann mode if NetlistSubcircuit is a boundary condition)
	{
	    int counterToDetectErrors = 0;
	    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
	    {
	        if (downstreamSubcircuit->lock()->boundaryConditionCircuitConnectsToThisDownstreamSubsection(m_surfaceIndex))
	        {
	            m_interfacePressure = downstreamSubcircuit->lock()->getComputedInterfacePressure(m_IndexOfThisNetlistLPNInInputFile);
	            counterToDetectErrors++;
	        }
	    }
	    assert(counterToDetectErrors == 1);

		returnValue = std::make_pair(Boundary_Pressure,m_interfacePressure);
	}
	else if (mp_circuitData->hasPrescribedPressureAcrossInterface()) // This boundary condition is receiving pressure and returning flow (Dirichlet mode if NetlistSubcircuit is a boundary condition)
	{
		int counterToDetectErrors = 0;
	    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
	    {
	        if (downstreamSubcircuit->lock()->boundaryConditionCircuitConnectsToThisDownstreamSubsection(m_surfaceIndex))
	        {
	            m_interfaceFlow = downstreamSubcircuit->lock()->getComputedInterfaceFlow(m_IndexOfThisNetlistLPNInInputFile);
	            counterToDetectErrors++;
	        }
	    }
	    assert(counterToDetectErrors == 1);
	    
		returnValue = std::make_pair(Boundary_Flow,m_interfaceFlow);
	}
	else
	{
		std::stringstream errorMessage;
		errorMessage << "EE: Internal error when computing flow or pressure to pass to zero-D domain replacement." << std::endl;
		throw std::logic_error(errorMessage.str());
	}

	return returnValue;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getMatrixContribution(const double alfi_delt, Mat& matrixFromThisBoundary)
{
    generateLinearSystemWithoutFactorisation(alfi_delt);
    matrixFromThisBoundary = m_systemMatrix;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getRHSContribution(const int timestepNumber, Vec& rhsFromThisBoundary)
{
    assembleRHS(timestepNumber);
    rhsFromThisBoundary = m_RHS;
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getLocationOf3DInterfaceComputedFlowInSolutionVector() const
{
    return m_locationOf3DInterfaceComputedFlowInSolutionVector.at(0);
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getLocationOf3DInterfaceComputedPressureInSolutionVector() const
{
	return m_locationOf3DInterfaceComputedPressureInSolutionVector.at(0);
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getColumnOf3DInterfacePrescribedPressureInLinearSystem() const
{
    return m_columnOf3DInterfacePrescribedPressureInLinearSystem.at(0);
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getColumnOf3DInterfacePrescribedFlowInLinearSystem() const
{
    return m_columnOf3DInterfacePrescribedFlowInLinearSystem.at(0);
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getCircuitIndex() const
{
    return m_IndexOfThisNetlistLPNInInputFile;
}