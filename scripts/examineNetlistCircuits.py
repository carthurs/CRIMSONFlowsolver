#!/apps/anaconda/bin/python

import xml.etree.ElementTree as ElementTree


def is_resistor(circuit_component):
    if circuit_component.find('type').text.lower() == 'resistor':
        return True
    else:
        return False

def is_inductor(circuit_component):
    if circuit_component.find('type').text.lower() == 'inductor':
        return True
    else:
        return False

def is_diode(circuit_component):
    if circuit_component.find('type').text.lower() == 'diode':
        return True
    else:
        return False

def is_capacitor(circuit_component):
    if circuit_component.find('type').text.lower() == 'capacitor':
        return True
    else:
        return False

def is_volume_tracking_ressure_chamber(circuit_component):
    if circuit_component.find('type').text.lower() == 'volumetrackingpressurechamber':
        return True
    else:
        return False


def a_or_an(next_word):
    if next_word[0] in ['a','e','i','o','u']:
        return "an " + next_word
    else:
        return "a " + next_word


if __name__ == '__main__':
    netlist_surfaces_xml_tree = ElementTree.parse('netlist_surfaces.xml')
    all_circuits = netlist_surfaces_xml_tree.getroot()

#    print all_circuits.tag, all_circuits.attrib
    total_reciprocal_resistance_all_boundaries = 0.0
    for circuit in all_circuits:
#        print circuit.tag, circuit.attrib
        print "    ==== Circuit number", circuit.find('circuitIndex').text + " ===="
        components = circuit.find('components')

        total_circuit_resistance = 0.0

        for component in components:
            control_info = component.find('control')
            if control_info is not None:
                control_info_string = " and control of type " + control_info.find('type').text
                control_source = control_info.find('source')
                if control_source is not None:
                    control_info_string += " with source name " + control_source.text
            else:
                control_info_string = ""

            get_xml_attribute = lambda xml_tag_name : component.find(xml_tag_name).text
            print "Component", get_xml_attribute('index'), "is", a_or_an(get_xml_attribute('type')), "with parameter", get_xml_attribute('parameterValue') + control_info_string

            if is_resistor(component):
                total_circuit_resistance += float(get_xml_attribute('parameterValue'))

        if total_circuit_resistance != 0.0:
            total_reciprocal_resistance_all_boundaries += 1.0/total_circuit_resistance
        print "\nSum of circuit resistance (resistors only; ignoring diodes - note that this is a straight sum, NOT an equivalent resistance):", total_circuit_resistance
        print "============================================"

    print "\n\nPotentially Useless Resistance Metric: (NOT total peripheral resistance - unless all Netlists are outflows and all resistances in outflows are sequential)", 1.0/total_reciprocal_resistance_all_boundaries
    print "\n"

