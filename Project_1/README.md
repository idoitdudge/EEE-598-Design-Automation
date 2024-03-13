This folder consists of all the files relating to Project 1: Developing a topological traversal algorithm to perform Static Timing Analysis for circuits <br />

# Requirements
Requirements.txt contains the dependencies needed to run the program. To load all of them, use the command:
```

pip3 install -r requirements.txt

```

## Commands

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

To get circuit delay, critical path, and slack values, run command for c17 bench:
```

python3.7 main_sta.py --read_ckt c17.bench --read_nldm sample_NLDM.lib

```

To get circuit delay, critical path, and slack values, run command for c432 bench:
```

python3.7 main_sta.py --read_ckt c432.bench --read_nldm sample_NLDM.lib

```


To get circuit delay, critical path, and slack values, run command for c6288 bench:
```

python3.7 main_sta.py --read_ckt c6288.bench --read_nldm sample_NLDM.lib

```


To get circuit delay, critical path, and slack values, run command for c7522 bench:
```

python3.7 main_sta.py --read_ckt c7522.bench --read_nldm sample_NLDM.lib

```


To get circuit delay, critical path, and slack values, run command for b01_C bench:
```

python3.7 main_sta.py --read_ckt b01_C.bench --read_nldm sample_NLDM.lib

```


To get circuit delay, critical path, and slack values, run command for b02_C bench:
```

python3.7 main_sta.py --read_ckt b02_C.bench --read_nldm sample_NLDM.lib

```

To get circuit delay, critical path, and slack values, run command for b15_C bench:
```
python3.7 main_sta.py --read_ckt b15_C.bench --read_nldm sample_NLDM.lib

```

To get circuit delay, critical path, and slack values, run command for b20_C bench:
```
python3.7 main_sta.py --read_ckt b20_C.bench --read_nldm sample_NLDM.lib

```
