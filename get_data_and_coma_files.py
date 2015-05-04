import yaml
import sys
import os
from create_mdep_data_files import write_files



def make_data_files(card, path):
	"""
	Creates data and coma files according to the definitions in card
	@param card: Name of a YAML file
	@type card: str
	@param path: Location for the files
	@type path: str
	"""
	with open(card) as inin:
		definitions = yaml.load(inin)
	method = definitions["method"]
	waves = []
	upper = []
	lower = []
	for wave in definitions["waves"]:
		waves.append(wave["name"])
		upper.append(float(wave["mmax"]))
		lower.append(float(wave["mmin"]))
	name = definitions["fit_name"]
	if len(definitions["t_binning"]) != len(definitions["mass_independent_fit"]):
		raise IndexError # Not enough sources for all t' bins given

	try:
		if len(definitions["data_files"]) == len(definitions["t_binning"]) and len(definitions["coma_files"]) == len(definitions["t_binning"]):
			print "Files already seem to exists, do not create them"
			return
	except KeyError:
		pass
	dataNames, comaNames = write_files(path, definitions["mass_independent_fit"], waves, upper, lower, method, "data_"+name+"_", "coma_"+name+"_")
	
	with open(card,'a') as inout:
		inout.write("data_files:\n")
		for dataFile in dataNames:
			inout.write("- "+dataFile+"\n")
		inout.write("coma_files:\n")
		for comaFile in comaNames:
			inout.write("- "+comaFile+"\n")
	

if __name__ == "__main__":

# old_method works, other modes still to be tested

	make_data_files(sys.argv[1],'/nfs/mds/user/fkrinner/data_mdep_fit/testFiles')

