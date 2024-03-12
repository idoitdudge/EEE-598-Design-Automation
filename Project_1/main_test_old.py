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
        self.outname = "" # full gate name
        self.gatetype = "" # gate type
        self.Cload = 0.0
        self.inputs = [] # gate number for fan-in
        self.fan_in = [] # full gate name for fan-in
        self.outputs = [] # gate number for fan-out
        self.fan_out = [] # full gate name for fan-out
        self.delays = []
        self.Tau_in = []
        self.Tau_calcs = []
        self.inp_arrival = []
        self.outp_arrival = []
        self.marked = []
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
    num_of_gates = {}

    for line in content:
        line.strip() # remove trailing white spaces

        if ('inputs' in line): # parse first input line
            input_num = line.split()[1]
            
        elif ('outputs' in line): # parse second output line
            output_num = line.split()[1]
        
        # elif ('gates' in line): # parse total gates line
        #     gate_type = re.search('(.*?) \(\s*(.*?)\)', line, flags=re.DOTALL).group(2)
        #     total_gates.append(gate_type.split(' + '))

        if (line.startswith ("INPUT")): # parse input gates
            gate_name = re.search('INPUT\((.*?)\)', line, flags=re.DOTALL).group(1)
            gates[gate_name] = Node()
            gates[gate_name].name = str(gate_name)
            gates[gate_name].gatetype = str("INPUT")
            gates[gate_name].outname = "INPUT-" + str(gate_name)
        
        elif (line.startswith ("OUTPUT")): # parse input gates
            gate_name = re.search('OUTPUT\((.*?)\)', line, flags=re.DOTALL).group(1)
            outputs[gate_name] = Node()
            outputs[gate_name].name = str(gate_name)
            outputs[gate_name].gatetype = str("OUTPUT")
            outputs[gate_name].outname = "OUTPUT-" + str(gate_name)
        
        elif ('=' in line) and (',' in line): # parse gates with more than one inputs
            gate_name = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(1)
            gate_type = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(2)
            
            if gate_type in num_of_gates:
                num_of_gates[gate_type] += 1
            
            else:
                num_of_gates[gate_type] = 1

            gates[gate_name] = Node()
            gates[gate_name].name = gate_name
            gates[gate_name].gatetype = gate_type
            gates[gate_name].outname = (gate_type + "-" + gate_name)
            fan_in = line.split('(')[1].split(')')[0].split(', ')
            gates_only.append(gate_name)
           
            for i in fan_in: # get fanin for current gate
                gates[gate_name].inputs.append(i)

        elif ('NOT' in line) or ('BUFF' in line): # parse NOT and BUFFER gates
            gate_name = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(1)
            gate_type = re.search('(.*?) = (.*?)\((.*?)\)', line, flags=re.DOTALL).group(2)

            if gate_type in num_of_gates:
                num_of_gates[gate_type] += 1
            
            else:
                num_of_gates[gate_type] = 1

    
            gates[gate_name] = Node()
            gates[gate_name].name = gate_name
            gates[gate_name].gatetype = gate_type
            gates[gate_name].outname = (gate_type + "-" + gate_name)
            fan_in = line.split('(')[1].split(')')[0].split(', ')
            gates_only.append(gate_name)

            for i in fan_in: # get fanin for current gate
                gates[gate_name].inputs.append(i)

    for gate in gates:
        if (gate in outputs) and gates[gate].gatetype == 'INPUT':
            gates_only.append(gate)
           # print(gate)

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

class LUT:
    def __init__(self):
        self.Allgate_name =  np.array([]) #name of cell
        self.Allgate_name_clean = np.array([]) # clean name of cell
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

        self.Allgate_name_clean = [gate.replace ('2_X1', '') for gate in self.Allgate_name]
        self.Allgate_name_clean = [gate.replace ('_X1', '') for gate in self.Allgate_name_clean]

        
        for i, cell_name in enumerate (self.Allgate_name_clean):
            if 'BUF' in cell_name:
                self.Allgate_name_clean[i] = cell_name.replace('BUF', 'BUFF')
            
            elif 'INV' in cell_name:
                self.Allgate_name_clean[i] = cell_name.replace('INV', 'NOT')
        
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

def initialize (gate_queue):
    for gate in gates:
        for i in range (len(gates[gate].fan_in)):
            gates[gate].marked.append(-1)

    for gate in gates:
        if gates[gate].gatetype == 'INPUT':
            gates[gate].Tau_out = 0.002
            gates[gate].marked.append(1)
            gates[gate].max_out_arrival = 0
            queue_insertion_check(gate_queue, gate)
        
def check_marked (gate):
    for i, item in enumerate(gates[gate].marked):
        if item == -1:
            gates[gate].marked[i] = 1
            break

def queue_insertion_check (gate_queue, gate):
    if all(element == 1 for element in gates[gate].marked):
        gate_queue.append(gate)

def find_indices (list, val):
    index1, index2 = 0, 0

    for i, value in enumerate(list):
        if value == val:
            return i

        elif value < val:
            index1 = i
        
        elif value > val:
            index2 = i
            break
    
    return index1, index2

def interpolation (v11, v12, v21, v22, c, c1, c2, tau, tau1, tau2):
    u = v11 * (c2 - c) * (tau2 - tau)
    v = v12 * (c - c1) * (tau2 - tau)
    w = v21 * (c2 - c) * (tau - tau1)
    x = v22 * (c - c1) * (tau - tau1)
    y = (c2 - c1) * (tau2 - tau1) # normalize factor
    return (u + v + w + x) / y

def delay_search (lut_instance, name, tau_in, cload, n_value):
    gate_index = lut_instance.Allgate_name_clean.index(name)
    
    tau_values = list(map(float, lut_instance.Tau_in_vals_delay[gate_index]))
    cload_values = list(map(float, lut_instance.Cload_vals_delay[gate_index]))
    
    slew_index = find_indices (tau_values, tau_in)
    cap_index = find_indices (cload_values, cload)

    if isinstance(slew_index, tuple) or isinstance(cap_index, tuple):
        if isinstance(slew_index, int):
            if n_value <= 2:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1])
            
            else:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1]) * (n_value / 2)
        
        elif isinstance(cap_index, int):
            if n_value <= 2:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]])

            else:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]]) * (n_value / 2)
            
        else:
            if n_value <= 2:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]])

            else:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]]) * (n_value / 2)


    else:
        if n_value <= 2:
            return float(lut_instance.All_delays[slew_index][cap_index])
        
        else:
            return float(lut_instance.All_delays[slew_index][cap_index]) * (n_value / 2)


def slew_search (lut_instance, name, tau_in, cload, n_value):
    gate_index = lut_instance.Allgate_name_clean.index(name)
    
    tau_values = list(map(float, lut_instance.Tau_in_vals_slew[gate_index]))
    cload_values = list(map(float, lut_instance.Cload_vals_slew[gate_index]))
    
    slew_index = find_indices (tau_values, tau_in)
    cap_index = find_indices (cload_values, cload)

    if isinstance(slew_index, tuple) or isinstance(cap_index, tuple):
        if isinstance(slew_index, int):
            if n_value <= 2:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1])
            
            else:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1]) * (n_value / 2)
        
        elif isinstance(cap_index, int):
            if n_value <= 2:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]])
            
            else:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]]) * (n_value / 2)

        else:
            if n_value <= 2:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]])
            
            else:
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]]) * (n_value / 2)

    else:
        if n_value <= 2:
            return float(lut_instance.All_delays[slew_index][cap_index])
        
        else:
            return float(lut_instance.All_delays[slew_index][cap_index]) * (n_value / 2)

def main ():
    lut_instance = LUT()
    lut_instance.assign_arrays(args.read_nldm)

    gate_queue = []

    initialize(gate_queue)

    while gate_queue:
        curr_gate = gate_queue.pop(0)

        for fan_outs in gates[curr_gate].outputs:
            if fan_outs != gates[curr_gate].name:
                check_marked (fan_outs)
                
                cap_index = lut_instance.Allgate_name_clean.index(gates[fan_outs].gatetype)

                gates[curr_gate].Cload = gates[curr_gate].Cload + float(lut_instance.Cvals[cap_index])

                # if fan_outs != gates[curr_gate].name:
                queue_insertion_check (gate_queue, fan_outs)

            else:
                gates[curr_gate].Cload = 4 * float(lut_instance.Cvals[lut_instance.Allgate_name_clean.index('NOT')])
        
        if (gates[curr_gate].gatetype == 'INPUT'):
            for fan_outs in gates[curr_gate].outputs:
                gates[fan_outs].Tau_in.append(gates[curr_gate].Tau_out)
                gates[fan_outs].inp_arrival.append(gates[curr_gate].max_out_arrival)
        
        else:
            for tau_in in gates[curr_gate].Tau_in:
                n_value = len(gates[curr_gate].Tau_in)
                gates[curr_gate].delays.append(delay_search (lut_instance, gates[curr_gate].gatetype, tau_in, gates[curr_gate].Cload, n_value))
                gates[curr_gate].Tau_calcs.append(slew_search (lut_instance, gates[curr_gate].gatetype, tau_in, gates[curr_gate].Cload, n_value))

            for i, arrival_time in enumerate(gates[curr_gate].inp_arrival):
                gates[curr_gate].outp_arrival.append(arrival_time + gates[curr_gate].delays[i])
            
            gates[curr_gate].max_out_arrival = max(gates[curr_gate].outp_arrival)
            max_index = gates[curr_gate].outp_arrival.index(gates[curr_gate].max_out_arrival)
            gates[curr_gate].Tau_out = gates[curr_gate].Tau_calcs[max_index]

            for fan_outs in gates[curr_gate].outputs:
                gates[fan_outs].inp_arrival.append(gates[curr_gate].max_out_arrival)
                gates[fan_outs].Tau_in.append(gates[curr_gate].Tau_out)

    circuit_delay = float('-inf')

    for gate in gates:
        for outs in outputs:
            if gates[gate].name == outputs[outs].name:
                if gates[gate].max_out_arrival > circuit_delay:
                    #print(gates[gate].max_out_arrival)
                    circuit_delay = gates[gate].max_out_arrival
    
    print("Circuit delay: ", circuit_delay * pow(10, 3), "picoseconds")

