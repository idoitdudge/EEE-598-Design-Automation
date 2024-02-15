import sys, os
import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description='Parser takes in arguments depending on the files the user wants to parse')

parser.add_argument('--read_ckt', default="c17.bench", type=str, help='File name of the circuit bench file')
parser.add_argument('--slews', action="store_true", help='Output slew data only from LUT')
parser.add_argument('--delays', action="store_true", help='Output delay data only from LUT')
parser.add_argument('--read_nldm', default="sample_NLDM.lib", type=str, help='Liberty file name to instantiate LUT')

args = parser.parse_args() # argument parser

class Node:
    def __init__(self): # data structure to save each node
        self.name = "" # number only
        self.outname = "" # gate type
        self.Cload = 0.0
        self.inputs = [] # gate number for fan-in
        self.fan_in = [] # full gate name for fan-in
        self.outputs = [] # gate number for fan-out
        self.fan_out = [] # full gate name for fan-out
        self.Tau_in = []
        self.inp_arrival = []
        self.outp_arrival = []
        self.max_out_arrival = 0.0
        self.Tau_out = 0.0

if args.read_ckt:
    with open(args.read_ckt, "r") as f:  # open file and read lines
        content = f.readlines()

    gates = {} # dict to save gates
    outputs = {} # dict to save outputs
    gates_only = [] # array to save only gates, no inputs
    input_num = "" # number of inputs
    output_num = "" # number of outputs
    total_gates = [] # total number of gates

    for line in content:
        line.strip() # remove trailing white spaces

        if ('inputs' in line): # parse first input line
            input_num = line.split()[1]
            
        elif ('outputs' in line): # parse second output line
            output_num = line.split()[1]
        
        elif ('gates' in line): # parse total gates line
            gate_type = re.search('(.*?) \(\s*(.*?)\)', line, flags=re.DOTALL).group(2)
            total_gates.append(gate_type.split(' + '))

        if (line.startswith ("INPUT")): # parse input gates
            gate_name = re.search('INPUT\((.*?)\)', line, flags=re.DOTALL).group(1)
            gates[gate_name] = Node()
            gates[gate_name].name = str(gate_name)
            gates[gate_name].outname = "INPUT-" + str(gate_name)
        
        elif (line.startswith ("OUTPUT")): # parse input gates
            gate_name = re.search('OUTPUT\((.*?)\)', line, flags=re.DOTALL).group(1)
            outputs[gate_name] = Node()
            outputs[gate_name].name = str(gate_name)
            outputs[gate_name].outname = "OUTPUT-" + str(gate_name)
        
        elif ('=' in line) and (',' in line): # parse gates with more than one inputs
            gate_name = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(1)
            gate_type = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(2)
            gates[gate_name] = Node()
            gates[gate_name].name = gate_name
            gates[gate_name].outname = (gate_type + "-" + gate_name)
            fan_in = line.split('(')[1].split(')')[0].split(', ')
            gates_only.append(gate_name)
           
            for i in fan_in: # get fanin for current gate
                gates[gate_name].inputs.append(i)

        elif ('NOT' in line) or ('BUFF' in line): # parse NOT and BUFFER gates
            gate_name = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(1)
            gate_type = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(2)
            gates[gate_name] = Node()
            gates[gate_name].name = gate_name
            gates[gate_name].outname = (gate_type + "-" + gate_name)
            fan_in = line.split('(')[1].split(')')[0].split(', ')
            gates_only.append(gate_name)

            for i in fan_in: # get fanin for current gate
                gates[gate_name].inputs.append(i)

    for gate in gates:
        for i in gates[gate].inputs:
            gates[gate].fan_in.append(gates[i].outname)

    for gate in gates: # get fanouts for each gates
        if gates[gate].inputs:
            for i in gates[gate].inputs:
                gates[i].outputs.append(gate)
                gates[i].fan_out.append(gates[gate].outname)

    for output in outputs: # match OUTPUT- gates to corresponding final output gates
        gates[output].outputs.append(output)
        gates[output].fan_out.append(outputs[output].outname)

    # test = '431'

    # print(gates[test].outname)
    # print(gates[test].fan_in)
    # print(gates[test].fan_out)

    with open ('ckt_details.txt', 'w') as ckt_detail: # write outputs to ckt_details.txt
        ckt_detail.write(input_num + " primary inputs\n")
        ckt_detail.write(output_num + " primary outputs\n")
        for i in total_gates:
            ckt_detail.write(', '.join(i) + "\n")

        ckt_detail.write("Fanout...\n")
        for gate_num in gates_only:
            ckt_detail.write(gates[gate_num].outname + ": " + ', '.join(gates[gate_num].fan_out))
            ckt_detail.write("\n")

        ckt_detail.write("\nFanin...\n")
        for gate_num in gates_only:
            ckt_detail.write(gates[gate_num].outname + ": " + ', '.join(gates[gate_num].fan_in))
            ckt_detail.write("\n")

    ckt_detail.close()



class LUT:
    def __init__(self):
        self.Allgate_name =  np.array([]) #name of cell
        self.All_delays = np.empty((0, 7, 7)) #2D numpy array delay LUTs for each cell
        self.All_slews = np.empty((0, 7, 7)) #2D numpy array to store output slew LUTs for each cell
        self.Cload_vals_slew = np.empty((0, 7)) #1D numpy array corresponds to the 2nd index in the LUT
        self.Tau_in_vals_slew = np.empty((0, 7)) #1D numpy array corresponds to the 1st index in the LUT
        self.Cload_vals_delay = np.empty((0, 7))
        self.Tau_in_vals_delay = np.empty((0, 7))
        self.Cvals = np.array([])

    def assign_arrays (self, NLDM_file):
        with open (NLDM_file, 'r') as lib_file:
            content = lib_file.read()

        for cell_name in re.finditer('cell \((.*?)\)', content):
            #print(cell_name.group(1))
            self.Allgate_name = np.append(self.Allgate_name, cell_name.group(1))
        
        for cap_value in re.finditer('capacitance\s*:\s*(.*);', content):
            #print(cap_value.group(1))
            self.Cvals = np.append(self.Cvals, cap_value.group(1))
        
        for cell_name in self.Allgate_name:
            current_cell = re.search('cell \(' + cell_name + '\)(.*?)(cell \(|$)', content, flags=re.DOTALL).group(1)
            for delay_slew in ['cell_delay', 'output_slew']:
                data_type = re.search(delay_slew + '(.*?){(.*?)}', current_cell, flags=re.DOTALL).group(2)

                if (delay_slew == 'cell_delay'):
                    input_slew_vals = re.search('index_1 \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1)
                    self.Tau_in_vals_delay = np.vstack([self.Tau_in_vals_delay, input_slew_vals.split(',')])
                    output_cload_vals = re.search('index_2 \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1)
                    self.Cload_vals_delay = np.vstack([self.Cload_vals_delay, output_cload_vals.split(',')])
                    proxy_array = np.empty([0,7])
                    delay_values = re.search('values \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1).replace("\n", "").replace("\n", "").replace("\t", "").replace(" ", "")
                    processed_delay_vals = delay_values.split('",\\"')
                    for i in processed_delay_vals:
                        delay_arr = i.split(',')
                        proxy_array = np.vstack([proxy_array, delay_arr])
                    self.All_delays = np.vstack([self.All_delays, [proxy_array]])
                
                else:
                    input_slew_vals = re.search('index_1 \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1)
                    self.Tau_in_vals_slew = np.vstack([self.Tau_in_vals_slew, input_slew_vals.split(',')])
                    output_cload_vals = re.search('index_2 \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1)
                    self.Cload_vals_slew = np.vstack([self.Cload_vals_slew, output_cload_vals.split(',')])
                    proxy_array = np.empty([0,7])
                    slew_values = re.search('values \(\"(.*?)\"\)', data_type, flags=re.DOTALL).group(1).replace("\n", "").replace("\n", "").replace("\t", "").replace(" ", "")
                    processed_slew_vals = slew_values.split('",\\"')
                    for i in processed_slew_vals:
                        slew_arr = i.split(',')
                        proxy_array = np.vstack([proxy_array, slew_arr])
                    self.All_slews = np.vstack([self.All_slews, [proxy_array]])
                    #print(self.All_slews)

        if(args.slews):
            with open ('slew_LUT.txt', 'w') as slew_LUT: # write outputs to slew_LUT.txt
                for i in range(0, len(self.Allgate_name)):
                    slew_LUT.write("Cell name: " + self.Allgate_name[i] + "\n")

                    slew_LUT.write("Input slews: " + ', '.join(map(str, self.Tau_in_vals_slew[i])) + "\n")
                    slew_LUT.write("Load capacitance: " + ', '.join(map(str, self.Cload_vals_slew[i])) + "\n")

                    slew_LUT.write("Slews:")
                    for j in range(0, len(self.Tau_in_vals_slew)):
                        slew_LUT.write(', '.join(map(str, self.All_slews[i][j])) + ';\n')
                
                    slew_LUT.write("\n")
        
        elif(args.delays):
            with open ('delay_LUT.txt', 'w') as delay_LUT: # write outputs to slew_LUT.txt
                for i in range(0, len(self.Allgate_name)):
                    delay_LUT.write("Cell name: " + self.Allgate_name[i] + "\n")

                    delay_LUT.write("Input slews: " + ', '.join(map(str, self.Tau_in_vals_delay[i])) + "\n")
                    delay_LUT.write("Load capacitance: " + ', '.join(map(str, self.Cload_vals_delay[i])) + "\n")

                    delay_LUT.write("Delays:")
                    for j in range(0, len(self.Tau_in_vals_delay)):
                        delay_LUT.write(', '.join(map(str, self.All_delays[i][j])) + ';\n')
                    
                    delay_LUT.write("\n")

lut_instance = LUT()
lut_instance.assign_arrays(args.read_nldm)
