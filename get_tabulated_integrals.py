import sys
import os
"""
Writes a file tabulated_integrals.h
This file can be included to give the integrals wanted for the mass dependent fitter
This is at the moment not very elegant, but not avoidable at the moment
"""
waves = [
		'1-(0-+)0+ f0(980) pi S', 
		'1-(1++)0+ f0(980) pi P',
		'1-(1++)0+ rho pi D',
		'1-(1++)0+ rho pi S',
		'1-(2++)1+ f2 pi P',
		'1-(2++)1+ rho pi D',
		'1-(2++)2+ rho pi D',
		'1-(2-+)0+ f2 pi D',
		'1-(2-+)0+ f2 pi S',
		'1-(2-+)0+ rho pi F',
		'1-(2-+)1+ f2 pi S',
		'1-(4++)1+ f2 pi F',
		'1-(4++)1+ rho pi G',
		'1-(1-+)1+ rho pi P',
		'1-(3++)1+ rho pi G',
		'1-(3++)0+ rho3 pi S',
					]

numbers={
	"1-(0-+)0+ f0(980) pi S":0,
	"1-(4++)1+ rho pi G":1,
	"1-(1++)0+ rho pi S":2,
	"1-(1++)0+ f0(980) pi P":3,
	"1-(2-+)0+ f2 pi S":4,
	"1-(2++)1+ rho pi D":5,
	"1-(4++)1+ f2 pi F":6,
	"1-(1++)0+ rho pi D":7,
	"1-(2++)1+ f2 pi P":8,
	"1-(2++)2+ rho pi D":9,
	"1-(2-+)1+ f2 pi S":10,
	"1-(2-+)0+ rho pi F":11,
	"1-(2-+)0+ f2 pi D":12,
	"1-(1-+)1+ rho pi P":13,
	"1-(3++)1+ rho pi G":14,
	"1-(3++)0+ rho3 pi S":15
}
nwaves=str(len(waves))
numbers_inverse={}
for wave in numbers.iterkeys():
	number = numbers[wave]
	numbers_inverse[number]=wave



int_file = '/nfs/hicran/project/compass/analysis/sschmeing/PWA-massdependend/Integrals/NewIntegrals/work/intdiag_txt.dat'


ints={}
for wave in waves:
	vals = []
	integrals = open(int_file,'r')
	for line in integrals.readlines():
		if wave in line:
			chunks=line.split()
			vals.append([float(chunks[-3]),chunks[-2]])
	vals.sort()
	ints[wave]=vals
integrals.close()

npoints=str(len(ints[waves[0]]))

out_int = open("/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared_retry/tabulated_integrals.h",'w')
#out_int = open("tabulated_integrals.h",'w')
dot_template = open("tabulated_integrals.template",'r')

out_int.write("// This file was created by a python script:\n")
out_int.write("// "+os.path.realpath(__file__)+"\n")


for line in dot_template.readlines():
	if "<<INSERT_MAP>>" in line:
		for wave in waves:
			out_int.write("//"+wave+"\t"+str(numbers[wave])+"\n")
	elif "<<INSERT_NMAX>>" in line:
		out_int.write("\tint nnnMax = "+str(len(waves))+";\n")
	elif "<<INSERT_TABLES>>" in line:
		out_int.write("\t double ints["+nwaves+"]["+npoints+"] = {\n")
		for i in range(len(waves)):
			strs = [ppp[1] for ppp in ints[numbers_inverse[i]]]
			out_int.write("{"+", ".join(strs)+'}')
			if i<len(waves)-1:
				out_int.write(",\n")
			else:
				out_int.write("};\n")
	else:
		out_int.write(line)

dot_template.close()
out_int.close()	


