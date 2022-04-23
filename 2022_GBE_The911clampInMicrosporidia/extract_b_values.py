
#!/usr/bin/python

name = "extract_b_values.py"
version = "0.2.8"
updated = "2021-10-18"

import pickle
from sys import argv

## Load in .pkl elements into a dictionary
with open(argv[1],"rb") as pickled_pdb:
	unpickled_content = pickle.load(pickled_pdb)

## Grab b-values from the pickled information
b_values = unpickled_content['plddt']

## Write out the b-values for the given model
with open(f"{argv[1]}.b_values","w") as b_value_file:
	for val in b_values:
		print(val,end='\n',file=b_value_file)