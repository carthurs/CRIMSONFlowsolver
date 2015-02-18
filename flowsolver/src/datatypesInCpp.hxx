#ifndef DATATYPESINCPP_HXX_
#define DATATYPESINCPP_HXX_

enum circuit_component_t {Component_Null, Component_Resistor, Component_Capacitor, Component_Inductor, Component_Diode, Component_MonopolePressureNode};

enum circuit_nodal_pressure_prescription_t {Pressure_Null, Pressure_NotPrescribed, Pressure_Fixed, Pressure_LeftVentricular, Pressure_3DInterface, Pressure_VolumeDependent};

enum circuit_component_flow_prescription_t {Flow_Null, Flow_NotPrescribed, Flow_MonopoleSoUndefined, Flow_Fixed, Flow_3DInterface};

enum circuit_diode_node_t {Node_Null, Node_ConnectsCircuit, Node_IsMonopolar};

#endif