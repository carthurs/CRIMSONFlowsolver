#!/apps/anaconda/bin/python

def a_or_an(next_word):
    if next_word[0] in ['a','e','i','o','u']:
        return "an " + next_word
    else:
        return "a " + next_word


if __name__ == '__main__':

    with open('rcrt.dat', 'r') as file:
        rcrt_data = file.readlines()

    # strip whitespace
    rcrt_data = [line.strip() for line in rcrt_data]

    # Strip initial file metadata:
    del rcrt_data[0]

    # Break up into single-circuit sublists:
    all_circuit_data_blocks = []
    current_line_index = 0

    while current_line_index < len(rcrt_data):
         number_of_time_varying_data_lines = int(rcrt_data[current_line_index])
         single_circuit_block_data = []
         for block_line_index in range(4 + number_of_time_varying_data_lines):
             single_circuit_block_data.append(rcrt_data[current_line_index])
             current_line_index += 1
         all_circuit_data_blocks.append(single_circuit_block_data)
         
         

    circuit_number = 1
    total_resistance_all_rcrs = 0.0
    for circuit in all_circuit_data_blocks:
        print "    ==== Circuit number", str(circuit_number) + " ===="
        print "Resistor Rp has resistance", float(circuit[1])
        print "Capacitor has compliance", float(circuit[2])
        print "Resistor Rd has resistance", float(circuit[3])
        total_circuit_resistance = float(circuit[1]) + float(circuit[3])

        print "Total resistance is", total_circuit_resistance
        total_resistance_all_rcrs += total_circuit_resistance
        circuit_number += 1
        print "============================================"

    print "\n\nTotal resistance of RCRs contained in rcrt.dat: (This is only equal to total peripheral resistance if there are only rcrt.dat RCRs at outflows)", total_resistance_all_rcrs
    print "\n"

