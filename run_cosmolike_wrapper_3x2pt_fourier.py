import sys
import os.path
import argparse
import yaml
import numpy as np
from run_cosmolike_3x2pt_fourier import *

########## main #######
parser = argparse.ArgumentParser(description='call run_cosmolike_mpp outside the pipeline')
parser.add_argument("parameter_file", help="YAML configuration file")
args = parser.parse_args()
try:
    param_file = args.parameter_file
    print(param_file)
except SystemExit:
	sys.exit(1)
load_yaml = 1
try:
	params = yaml.load(open(param_file), Loader=yaml.FullLoader)
	for sc in params.get("scale_cuts", []):
		params.update(yaml.load(open(sc)))
except:
	load_yaml = 0
if (not load_yaml):
	try:
		params = yaml.full_load(open(param_file))
		for sc in params.get("scale_cuts", []):
			params.update(yaml.full_load(open(sc)))
	except:
		print("Could not read yaml file\n")
		exit(1)

print("Alive")
run_cosmolike(params)
