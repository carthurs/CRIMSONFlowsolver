#ifndef DATATYPESINCPP_HXX_
#define DATATYPESINCPP_HXX_

#include <vector>
#include <utility>

enum boundary_condition_t {BoundaryCondition_Null, BoundaryCondition_RCR, BoundaryCondition_ControlledCoronary, BoundaryCondition_Netlist, BoundaryCondition_Impedance};

enum circuit_component_t {Component_Null, Component_Resistor, Component_Capacitor, Component_Inductor, Component_Diode, Component_MonopolePressureNode, Component_VolumeTracking, Component_VolumeTrackingPressureChamber};

enum circuit_nodal_pressure_prescription_t {Pressure_Null, Pressure_NotPrescribed, Pressure_Fixed, Pressure_LeftVentricular, Pressure_3DInterface};

enum circuit_component_flow_prescription_t {Flow_Null, Flow_NotPrescribed, Flow_MonopoleSoUndefined, Flow_Fixed, Flow_3DInterface};

enum circuit_diode_node_t {Node_Null, Node_ConnectsCircuit, Node_IsMonopolar};

enum parameter_controller_t {Controller_Null, Controller_LeftVentricularElastance, Controller_BleedResistance, Controller_BleedCompliance, Controller_CustomPythonComponentParameter, Controller_CustomPythonComponentFlowFile, Controller_CustomPythonNode, Controller_CustomPythonNodePressureFile, Controller_CustomPythonComponentUnstressedVolume};

enum boundary_data_t {Boundary_Pressure, Boundary_Flow};

enum circuit_item_t {Circuit_Component, Circuit_Node};

typedef std::vector<std::pair<double,double>> TimeValuePairVector;

#endif