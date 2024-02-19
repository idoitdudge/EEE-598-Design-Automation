This folder consists of all the files relating to Project 1: Developing a topological traversal algorithm to perform Static Timing Analysis for circuits <br />

# Commands

To use the parser, run the command to get circuit information:
'''

python3.7 parser.py --read_ckt file_name.bench

'''

To get delay data, run command for delay LUT:
'''

python3.7 parser.py --delays --read_nldm sample_NLDM.lib 

'''

To get slew data, run command for delay LUT:
'''

python3.7 parser.py --slews --read_nldm sample_NLDM.lib

'''