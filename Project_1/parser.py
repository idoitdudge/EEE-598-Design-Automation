import sys, os
import argparse

class Node:
    # '''Configuration and Argument parser to generate information on inquired data'''

    # def __init__(self, args):
    #     self.parser.add_argument('--read_ckt', default=1, type=str, help='read and output bench file information')

    def __init__(self): # data structure to save each node
        self.name = "" # number only
        self.outname = "" # gate type + gate number
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
    
parser = argparse.ArgumentParser(description='Parser takes in arguments depending on the files the user wants to parse')

parser.add_argument('--read_ckt', type=str, help='File name of the circuit bench file')

args = parser.parse_args() # argument parser

with open(args.read_ckt, "r") as f:  # open file and read lines
    content = f.readlines()

gates = {} # array to save gates
outputs = {} # array to save outputs
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
        gate_type = line.split('( ')[1].rstrip(' )\n')
        total_gates.append(gate_type.split(' + '))
        total_gates = [[s.rstrip('s') for s in sublist] for sublist in total_gates]

    elif (line.startswith ("INPUT")): # parse input gates
        num = line.split('(')[1].rstrip(')\n')
        gates[num] = Node()
        gates[num].name = str(num)
        gates[num].outname = "INPUT-" + str(num)
    
    elif (line.startswith ("OUTPUT")): # parse input gates
        num = line.split('(')[1].rstrip(')\n')
        outputs[num] = Node()
        outputs[num].name = str(num)
        outputs[num].outname = "OUTPUT-" + str(num)
    
    elif (',' in line): # parse gates with more than one inputs
        gate_type = line.split('= ')[1].split('(')[0].strip()
        num = line.split(' =')[0].strip()
        fan_in = line.split('(')[1].split(')')[0].split(', ')
        gates_only.append(num)
        gates[num] = Node()
        gates[num].name = str(num)
        gates[num].outname = str(gate_type) + "-" + str(num)
        for i in fan_in: # get fanin for current gate
            gates[num].inputs.append(i)
            gates[num].fan_in.append(gates[i].outname)

    elif ('NOT' in line) or ('BUFF' in line): # parse NOT and BUFFER gates
        gate_type = line.split('= ')[1].split('(')[0].strip()
        num = line.split(' =')[0].strip()
        fan_in = line.split('(')[1].split(')')[0].split(', ')
        gates_only.append(num)
        gates[num] = Node()
        gates[num].name = str(num)
        gates[num].outname = str(gate_type) + "-" + str(num)
        for i in fan_in: # get fanin for current gate
            gates[num].inputs.append(i)
            gates[num].fan_in.append(gates[i].outname)

    
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
        ckt_detail.write(', '.join(i) + " gates\n")

    ckt_detail.write("Fanout...\n")
    for gate_num in gates_only:
        ckt_detail.write(gates[gate_num].outname + ": " + ', '.join(gates[gate_num].fan_out))
        ckt_detail.write("\n")

    ckt_detail.write("Fanin...\n")
    for gate_num in gates_only:
        ckt_detail.write(gates[gate_num].outname + ": " + ', '.join(gates[gate_num].fan_in))
        ckt_detail.write("\n")
