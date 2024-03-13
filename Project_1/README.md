# Static Timing Analysis Topological Traversal Project
This folder consists of all the files relating to Project 1: Developing a topological traversal algorithm to perform Static Timing Analysis for circuits. 

## About This Project
This project utilizes topological traversal to perform static timing analysis on various different bench files and look-up tables. To find the various timing elements of each gate, we perform 2D bilinear interpolation in order to find the values that are not listed in the LUT.

### Files
This project consists of 3 major files. The commands to use each file is listed below.
The first file is ***requirements.txt***, where the required files to run this repository is listed within.
The second file is ***parser.py***, which lists out the circuit details, delay, and slew details of the LUT.
The third file is ***main_sta.py***, which contains the main functionalities of the STA traversal function.
#### Requirements
Requirements.txt contains the dependencies needed to run the program. To load all of them, use the command:
```

pip3 install -r requirements.txt

```

##### Commands

To use the parser, run the command to get circuit information:
```

python3.7 parser.py --read_ckt file_name.bench

```

To get delay data, run command for delay LUT:
```

python3.7 parser.py --delays --read_nldm sample_NLDM.lib 

```

To get slew data, run command for delay LUT:
```

python3.7 parser.py --slews --read_nldm sample_NLDM.lib

```

To un the main program to perform Static Timing Analysis, run the following command:
```

python3.7 main_sta.py --read_ckt file_name.bench --read_nldm sample_NLDM.lib

```



