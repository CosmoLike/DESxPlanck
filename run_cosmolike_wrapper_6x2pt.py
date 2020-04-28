import sys
import os.path
import argparse
import yaml
import numpy as np
from run_cosmolike_6x2pt import *

def make_mask(params, MASK):
	nl = params['lbins']
	ell = np.zeros(nl)
	l_min = params['lbounds'][0]
	l_max = params['lbounds'][1]
	dlgl = (np.log(l_max) - np.log(l_min))/nl
	for i in range(0,nl):
		lmin = np.exp(np.log(l_min)+i*dlgl)
		lmax = np.exp(np.log(l_min)+(i+1)*dlgl)
		ell[i] = 2./3.*(lmax*lmax*lmax-lmin*lmin*lmin)/(lmax*lmax-lmin*lmin)
	ntomo_source = params['ntomo_source']
	ntomo_lens = params['ntomo_lens']
	ndata = 0
	mask = []
	if "xip" in params['twoptnames']:
		for i in range(0,ntomo_source):
			for j in range(i,ntomo_source):
				param = "angle_range_xip_{}_{}".format(i+1,j+1)
				tmin, tmax = params[param]
				for k in range(0,ntheta):
					if ((theta[k] < tmin) | (theta[k] > tmax)):
						mask.append(0.)
					else:
						mask.append(1.)
						ndata += 1
	else:
		for i in range(0,ntomo_source):
			for j in range(i,ntomo_source):
				for k in range(0,ntheta):
		 			mask.append(0.)
	if "xim" in params['twoptnames']:
		for i in range(0,ntomo_source):
			for j in range(i,ntomo_source):
				param = "angle_range_xim_{}_{}".format(i+1,j+1)
				tmin, tmax = params[param]
				for k in range(0,ntheta):
					if ((theta[k] < tmin) | (theta[k] > tmax)):
						mask.append(0.)
					else:
						mask.append(1.)
						ndata += 1
	else:
		for i in range(0,ntomo_source):
			for j in range(i,ntomo_source):
				for k in range(0,ntheta):
		 			mask.append(0.)
	if "gammat" in params['twoptnames']:
		for i in range(0,ntomo_lens):
			for j in range(0,ntomo_source):
#				if ggltomo[i*ntomo_source+j,2]:
					param = "angle_range_gammat_{}_{}".format(i+1,j+1)
					tmin, tmax = params[param]
					for k in range(0,ntheta):
						if ((theta[k] < tmin) | (theta[k] > tmax)):
							mask.append(0.)
						else:
							mask.append(1.)
							ndata += 1
	else:
		for i in range(0,ntomo_lens):
			for j in range(0,ntomo_source):
				#if ggltomo[i*ntomo_source+j,2]:
					for k in range(0,ntheta):
			 			mask.append(0.)
	if "wtheta" in params['twoptnames']:
		for i in range(0,ntomo_lens):
			param = "angle_range_wtheta_{}_{}".format(i+1,i+1)
			tmin, tmax = params[param]
			for k in range(0,ntheta):
				if ((theta[k] < tmin) | (theta[k] > tmax)):
					mask.append(0.)
				else:
					mask.append(1.)
					ndata += 1
	else:
		for i in range(0,ntomo_lens):
			for k in range(0,ntheta):
					mask.append(0.)
	params['mask_checksum'] = ndata
	if MASK:
		f = open(params['mask_file'],"w")
		for i in range(0,len(mask)):
			f.write("%d %.1f\n"%(i, mask[i]))
		f.close()
		print("Wrote make_file %s\nEXIT\n"%(params['mask_file']))
		exit(1)
	#params['new_mask_file'] =maskfile
########## main #######
parser = argparse.ArgumentParser(description='call run_cosmolike_mpp outside the pipeline')
parser.add_argument("parameter_file", help="YAML configuration file")
parser.add_argument('-mask',dest = "mask", action='store_true',help ="create mask file only")
args = parser.parse_args()
try:
    param_file = args.parameter_file
except SystemExit:
	sys.exit(1)
load_yaml = 1
try:
	params = yaml.load(open(param_file))
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

if ((not os.path.isfile(params['mask_file'])) and (not args.mask)):
	print("\nmask_file %s doesn't exist!\nRerun \"run_cosmolike_wrapper.py YAML_FILE -mask\" to create mask_file from scale cuts\n"%(params['mask_file']))
	exit(1)
# make_mask(params, args.mask)
run_cosmolike(params)
