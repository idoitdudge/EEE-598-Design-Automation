import sys, os
import argparse
import numpy as np
import re
from collections import deque

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
        self.delays = [] # delay values
        self.Tau_in = [] # tau in values
        self.Tau_calcs = [] # tau values during interpolation
        self.inp_arrival = [] # input arrival time
        self.outp_arrival = [] # output arrival time
        self.marked = [] # node marker for traversal
        self.max_out_arrival = 0.0 # max arrival time
        self.Tau_out = 0.0 # output slew
 
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
            #print(gate)

    for gate in gates:
        for i in gates[gate].inputs:
            gates[gate].fan_in.append(gates[i].outname)

    for gate in gates: # get fanouts for each gates
        if gates[gate].inputs:
            for i in gates[gate].inputs:
                gates[i].outputs.append(gate)
                gates[i].fan_out.append(gates[gate].outname)

    for output in outputs: # match OUTPUT- gates to corresponding final output gates
        if (gates[output].gatetype == 'INPUT'):
            gates[output].inputs.append(output)
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

    def assign_arrays (self, NLDM_file): # Parse LUT file
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

        
        for i, cell_name in enumerate (self.Allgate_name_clean): # Clean buffer and inverter name to a new array
            if 'BUF' in cell_name:
                self.Allgate_name_clean[i] = cell_name.replace('BUF', 'BUFF')
            
            elif 'INV' in cell_name:
                self.Allgate_name_clean[i] = cell_name.replace('INV', 'NOT')
        
        for cell_name in self.Allgate_name: # Put index1, index2, slew and delay values in arrays 
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

def initialize (gate_queue): # Initialize queue and other attributes for forward traversal
    for gate in gates:
        for i in range (len(gates[gate].inputs)):
            gates[gate].marked.append(-1)
            gates[gate].Tau_in.append(-1)
            gates[gate].inp_arrival.append(-1)

    for gate in gates:
        if gates[gate].gatetype == 'INPUT':
            gates[gate].Tau_out = 0.002
            gates[gate].Tau_in = []
            gates[gate].Tau_in.append(-1)
            gates[gate].marked = []
            gates[gate].marked.append(1)
            gates[gate].max_out_arrival = 0
            queue_insertion_check(gate_queue, gate)
        
def check_marked (gate): # Mark Node.marked attribute for traversal
    for i, item in enumerate(gates[gate].marked):
        if item == -1:
            gates[gate].marked[i] = 1
            break

def queue_insertion_check (gate_queue, gate): # Queue insertion checker. Only enters queue if marked is all 1
    if all(element == 1 for element in gates[gate].marked):
        gate_queue.append(gate)

def find_indices (list, val): # Find indices for index1 and index2. If they don't belong, return a tuple of the high and low index
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

def interpolation (v11, v12, v21, v22, c, c1, c2, tau, tau1, tau2): # 2D bilinear interpolation function
    u = v11 * (c2 - c) * (tau2 - tau)
    v = v12 * (c - c1) * (tau2 - tau)
    w = v21 * (c2 - c) * (tau - tau1)
    x = v22 * (c - c1) * (tau - tau1)
    y = (c2 - c1) * (tau2 - tau1) # normalize factor
    return (u + v + w + x) / y

def delay_search (lut_instance, name, tau_in, cload, n_value): # Delay search function. Takes in Cload and Tau in and performs 2D Interpolation
    gate_index = lut_instance.Allgate_name_clean.index(name)
    tau_values = list(map(float, lut_instance.Tau_in_vals_delay[gate_index]))
    cload_values = list(map(float, lut_instance.Cload_vals_delay[gate_index]))
    
    slew_index = find_indices (tau_values, tau_in)
    cap_index = find_indices (cload_values, cload)

    if (tau_in > tau_values[len(tau_values) - 1]): # If values exceed the array, tuple becomes the last 2 indexes
        slew_index = (len(tau_values) - 2, len(tau_values) - 1)

    if (cload > cload_values[len(cload_values) - 1]):
        cap_index = (len(cload_values) - 2, len(cload_values) - 1) # Same as above

    if isinstance(slew_index, tuple) or isinstance(cap_index, tuple):
        if isinstance(slew_index, int):
            #print("1 int")
            if n_value <= 2:
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1])
            
            else:
                #print("1 int")
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1]) * (n_value / 2)
        
        elif isinstance(cap_index, int):
            if n_value <= 2:
                #print("1 int")
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]])

            else:
                #print("1 int")
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]]) * (n_value / 2)
            
        else:
            if n_value <= 2:
                #print("interp 2")
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]])

            else:
                #print("interp 2")
                return interpolation(float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_delays[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]]) * (n_value / 2)


    else:
        if n_value <= 2:
            #print("2 int")
            return float(lut_instance.All_delays[slew_index][cap_index])
        
        else:
            #print("2 int")
            return float(lut_instance.All_delays[slew_index][cap_index]) * (n_value / 2)


def slew_search (lut_instance, name, tau_in, cload, n_value): # Slew search function. Takes in Cload and Tau in and performs 2D Interpolation
    gate_index = lut_instance.Allgate_name_clean.index(name)
    
    tau_values = list(map(float, lut_instance.Tau_in_vals_slew[gate_index]))
    cload_values = list(map(float, lut_instance.Cload_vals_slew[gate_index]))



    slew_index = find_indices (tau_values, tau_in)
    cap_index = find_indices (cload_values, cload)

    if (tau_in > tau_values[len(tau_values) - 1]): # Same as the delay version. Last 2 indices if exceed final element's value
        slew_index = (len(tau_values) - 2, len(tau_values) - 1)

    if (cload > cload_values[len(cload_values) - 1]):
        cap_index = (len(cload_values) - 2, len(cload_values) - 1)

    if isinstance(slew_index, tuple) or isinstance(cap_index, tuple):
        if isinstance(slew_index, int):
            if n_value <= 2:
                #print("1 int")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1])
            
            else:
                #print("1 int")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index][cap_index[1]]), 
                                        float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index + 1][cap_index[1]]),
                                        cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                        tau_in, tau_values[slew_index], tau_values[slew_index + 1]) * (n_value / 2)
        
        elif isinstance(cap_index, int):
            if n_value <= 2:
                #print("1 int")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]])
            
            else:
                #print("1 int")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index + 1]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index + 1]),
                                cload, cload_values[cap_index], cload_values[gate_index][cap_index + 1],
                                tau_in, tau_values[slew_index[0]], tau_values[gate_index][slew_index[1]]) * (n_value / 2)

        else:
            if n_value <= 2:
                #print("interp 2")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]])
            
            else:
                #print("interp 2")
                return interpolation(float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[0]][cap_index[1]]), 
                                float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[0]]), float(lut_instance.All_slews[gate_index][slew_index[1]][cap_index[1]]),
                                cload, cload_values[cap_index[0]], cload_values[cap_index[1]],
                                tau_in, tau_values[slew_index[0]], tau_values[slew_index[1]]) * (n_value / 2)

    else:
        if n_value <= 2:
            #print("2 int")
            return float(lut_instance.All_delays[slew_index][cap_index])
        
        else:
            #print("2 int")
            return float(lut_instance.All_delays[slew_index][cap_index]) * (n_value / 2)

#Purpose of function is to sort gates in topological order from output to input
def sort(gates):
    visited = {}
    order = []

    def visit(gate): #checks if gates are visited, being visited, or are not visited
        if gate in visited:
            return
        
        visited[gate] = "visiting" 
        
        for in_gate in gate.inputs: 
            input_gate = gates[in_gate]
            if input_gate not in visited:
                visit(input_gate)
        
        visited[gate] = "visited"
        order.append(gate.name) #adds gate to queue 

    for gate in gates.values():
        if gate not in visited:
            visit(gate)
            
    return order 

#Purpose of this function is to calculate the required arrival time for each input, output, and gate
def arrival(gates, delay):
    max_time = float('inf') #set original time for everything to infinity (arbitraty value)
    RAT = {gate.name: max_time for gate in gates.values()} #Requried arrival time
    init_time = 1.1 * delay #sets intial time for everything to 1.1 * circuit_delay
    
    for gate in gates.values(): #Makes sure that the gates, inputs, outputs are created as node
        for input in gate.inputs:
            if input not in gates:
                gates.input = Node(input)
        for output in gate.outputs:
            if output not in gates:
                gates.output = Node(output)
    
    for gate in gates.values():
        if all(output not in gate2.inputs for output in gate.outputs for gate2 in gates.values() if gate2 != gate): #checks for primary outputs
            RAT[gate.name] = init_time
    
    processed_gates = []
    processed_order = sort(gates)
    
    for gate_name in reversed(processed_order): #Traverses circuit backwards using sort function to find required arrival time
        
        current_gate = gates[gate_name]  
        
        #checks inputs of gates, outputs of gates, and primary inputs to calculate required arrival
        if all(fanout in processed_gates for fanout in current_gate.outputs) or all(any(output not in gate2.inputs for gate2 in gates.values()) for output in current_gate.outputs):
            for i, in_gate in enumerate(current_gate.inputs):
                input_gate = gates[in_gate]
                
                if not input_gate.outputs:
                    input_gate.delays = [0]  # Check if delays is not empty
                
                if current_gate.delays:
                        rt_input_gate =  RAT[current_gate.name] - (current_gate.delays[i]) #calculating required arrival time using arcs of gates
                        RAT[in_gate] = min(RAT[in_gate], rt_input_gate) 
        
        
        # print("Gate: ", {current_gate.outname})
        # print(f"Inputs: {current_gate.inputs}")
        # print(f"Current gate delay of gate: ", current_gate.outname, (current_gate.delays))
        # print("Required arrival time: ", RAT[current_gate.name] * 1000, "ps") 
    
    return RAT

#Purpose of function is to traverse forward through the circuit to calculate forward arrival time to use in slack
def forwards(gates, outputs):
    FATimes = {}
    
    for gate in gates.values():
        FATimes[gate.name] = gate.max_out_arrival 
    
    return FATimes


#Purpose of this function is to find the critical path of the graph (path with greatest delay)
def critpath(gates, outputs):
    
    smol_slack_out = None
    
    for output_name in outputs:
        output = gates[output_name]
        if smol_slack_out is None or output.slack < smol_slack_out.slack:
            smol_slack_out = output

    crit = [smol_slack_out]
    current_gate = smol_slack_out
    
    while current_gate.inputs:
        
        smol_slack_in = None
        
        for input_name in current_gate.inputs:
            input_gate = gates[input_name]
            
            if smol_slack_in is None or input_gate.slack < gates[smol_slack_in].slack:
                smol_slack_in = input_name
        
        current_gate = gates[smol_slack_in]
        crit.insert(0, current_gate)
    
    crit_path = []
    last_gate = ""
    for gate in crit:
        crit_path.append(gate.outname)
        last_gate = gate.name
    crit_path.append("OUTPUT-" + last_gate)
    return crit_path

#Purpose of this function finds the output with the minimum slack
def slack(gates, outputs, delay, RAT, FATimes):        
    if isinstance(RAT, dict) and isinstance(FATimes, dict):  
        for gate in gates.values():
            if gate.name in FATimes and gate.name in RAT:
                forward_arrival_time = FATimes[gate.name]    
                required_arrival_time = RAT[gate.name]
                slack = required_arrival_time - forward_arrival_time
                gate.slack = slack
                    
def main (): # Main function for forward and backwards traversal
    lut_instance = LUT()
    lut_instance.assign_arrays(args.read_nldm)

    gate_queue = []

    initialize(gate_queue) # Initialize gate queue for forward traversal
    while gate_queue:
        curr_gate = gate_queue.pop(0)
        for fan_outs in gates[curr_gate].outputs: # Add values to queue
            if fan_outs != gates[curr_gate].name: # Edge case checker for if input doesn't go straight to output
                check_marked (fan_outs)
                
                cap_index = lut_instance.Allgate_name_clean.index(gates[fan_outs].gatetype)

                gates[curr_gate].Cload = gates[curr_gate].Cload + float(lut_instance.Cvals[cap_index])
            
                queue_insertion_check (gate_queue, fan_outs)

            elif (fan_outs == gates[curr_gate].name): # Edge case. When traversal arrives at primary outputs
                gates[curr_gate].Cload = gates[curr_gate].Cload + 4 * float(lut_instance.Cvals[lut_instance.Allgate_name_clean.index('NOT')])
        
        if (gates[curr_gate].gatetype == 'INPUT'): # Initialize Primary Input slew and max arrival to 0.002 and 0
            for fan_outs in gates[curr_gate].outputs:

                index = gates[fan_outs].inputs.index(gates[curr_gate].name)

                gates[fan_outs].Tau_in[index] = gates[curr_gate].Tau_out

                gates[fan_outs].inp_arrival[index] = gates[curr_gate].max_out_arrival
        else: # Normal delay and slew calculations if not primary 
            for tau_in in gates[curr_gate].Tau_in:
                n_value = len(gates[curr_gate].Tau_in)
                gates[curr_gate].delays.append(delay_search (lut_instance, gates[curr_gate].gatetype, tau_in, gates[curr_gate].Cload, n_value)) # Append to delay array
                gates[curr_gate].Tau_calcs.append(slew_search (lut_instance, gates[curr_gate].gatetype, tau_in, gates[curr_gate].Cload, n_value)) # Append to slew array
            
            for i, arrival_time in enumerate(gates[curr_gate].inp_arrival): # Calculate arrival times
                gates[curr_gate].outp_arrival.append(arrival_time + gates[curr_gate].delays[i])

            gates[curr_gate].max_out_arrival = max(gates[curr_gate].outp_arrival) # Find max arrival time out of all arrival times
            max_index = gates[curr_gate].outp_arrival.index(gates[curr_gate].max_out_arrival)
            gates[curr_gate].Tau_out = gates[curr_gate].Tau_calcs[max_index] # Find argmax slew of max arrival time

            for fan_outs in gates[curr_gate].outputs: # Add input slew as output slew of previous gate
                if (fan_outs == gates[curr_gate].name): # For edge case (primary output)
                    gates[fan_outs].inp_arrival.append(gates[curr_gate].max_out_arrival)
                    gates[fan_outs].Tau_in.append(gates[curr_gate].Tau_out)
                
                else: # Everything else
                    index = gates[fan_outs].inputs.index(gates[curr_gate].name)
                    gates[fan_outs].inp_arrival[index] = gates[curr_gate].max_out_arrival
                    gates[fan_outs].Tau_in[index] = (gates[curr_gate].Tau_out)
    delay = float('-inf')

    for gate in gates:
        for outs in outputs:
            if gates[gate].name == outputs[outs].name:
                if gates[gate].max_out_arrival > delay: 
                    delay = gates[gate].max_out_arrival # Output delay
    
    #print("Circuit delay: ", delay * pow(10, 3), "picoseconds") 
    
    sort(gates)


    RAT = arrival(gates, delay)
   
    
    FATimes = forwards(gates, outputs)
   
    slack(gates, outputs, delay, RAT, FATimes)

    cry_path = []
    cry_path = critpath(gates, outputs)
    final_path = ",".join(cry_path)
    #print("Critical path: ", cry_path)

    delay = delay * pow(10, 3)
    delay = str('{:.6g}'.format(delay))

    with open ('ckt_traversal.txt', 'w') as sta:

        sta.write(f"Circuit delay: {delay} picoseconds\n\n")
        sta.write(f"Gate slacks:\n")
        for gate in gates.values():
            sta.write(f"{gate.outname}: {str('{:.6g}'.format(gate.slack * pow(10, 3)))} picoseconds\n")
            for output in outputs:
                if gate.name == outputs[output].name:
                    sta.write(f"OUTPUT-{gate.name}: {str('{:.6g}'.format(gate.slack * pow(10, 3)))} picoseconds\n")
        
        sta.write(f"\nCritical path:\n")
        sta.write(f"{final_path} ")

main()