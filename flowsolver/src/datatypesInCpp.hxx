#ifndef DATATYPESINCPP_HXX_
#define DATATYPESINCPP_HXX_

enum circuit_component_t {Component_Resistor, Component_Capacitor, Component_Inductor};

enum circuit_nodal_pressure_prescription_t {Pressure_Fixed, Pressure_LeftVentricular};

enum circuit_component_flow_prescription_t {Flow_Fixed, Flow_3DInterface};

#endif