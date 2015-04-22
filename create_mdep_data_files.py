import numpy as np
import numpy.linalg as la
import os
import sys
from convert_text_output import readTextFile, getBestFits, getComaData

def write_anchor(target_dir,
		direct, 		# directory of the fit
		tbin_identifier,
		waves,			# wave names to be used
		upLims,			# mass limits for the waves
		loLims,
		name_data = 'data_',
		name_coma = 'coma_',
		CONJUGATE=True ):	# Conjugate the fit result.
	"""
	Write text files for the anchor_t method
	@param target_dir: Dorectory to write the files in
	@type: target_dir: str
	@param direct: Directory of the mass independent fit
	@type direct: str	
	@param tbin_identifier: Some name to identify the differnt t' bins. 
	@type tbin_identifier: str
	@param waves: List of wave names
	@type waves: list
	@param upLims: Upper mass limits for the waves, list of floats
	@type upLims: list
	@param loLims: Lower mass limits for the waves, list of floats
	@type loLims: list
	@param name_data: Some name base for the data files
	@type name_data: str
	@param name_coma: Some name base for the covariance matrix files
	@type name_coma: str
	@param CONJUGATE: Conjugae the fit result yes/no (switched some time, can also be done later in the fitter)
	@type CONJUGATE: bool
	@return: Names of the created files
	@rtype: str
	"""
	nWaves = len(waves)
	for i in range(1,len(upLims)):
		if upLims[i] > upLims[0] or loLims[i] < loLims[0]:
			print "First wave has to be active everywhere"
			raise ValueError

	data = getComaData(waves,upLims,loLims,direct, flagg='ANCHOR_FIRST', eps=1.E-3, CONJUGATE=CONJUGATE)
	pointsSorted=[]
	indices=[0]
	for i in range(1,nWaves):
		indices.append(i*nWaves) 
		indices.append(i)

	reduced_dat = []
	reduced_coma= []
	for bin in range(len(data[0])):
		dats = data[0][bin]
		dats_red=[]
		for i in indices:
			dats_red.append(dats[i])
		reduced_dat.append(dats_red)
		coma = data[1][bin]
		coma_red = []
		for i in indices:
			coma_line = []
			for j in indices:				
				coma_line.append(coma[i][j])
			coma_red.append(coma_line)
		reduced_coma.append(coma_red)
	dataName = target_dir+os.sep+name_data+tbin_identifier
	comaName = target_dir+os.sep+name_coma+tbin_identifier
	datFile = open(dataName,'w')
	comaFile= open(comaName,'w')
	for bin in range(len(reduced_dat)):
		dat = reduced_dat[bin]
		for ptt in dat:
			datFile.write(str(ptt)+'   ')
		coma = reduced_coma[bin]
		for line in coma:
			for ptt in line:
				comaFile.write(str(ptt)+'   ')
		comaFile.write('\n')
		datFile.write('\n')
	comaFile.close()
	datFile.close()
	return dataName, comaName

def get_raw_data_coma(
			fileName, 
			full_waveset):
	"""
	Reads a textoutput-file and returns the data and coma correspondign to the waveset
	@param fileName: File name of the fot result file
	@type fileName: str
	@param full_waveset: List of waves to use (list of strings)
	@type full_waveset: list
	@return: Data and covariance from tyhe file for the given waveset
	"""
	waveSet = []
	for wave in full_waveset:
		waveSet.append(wave.strip())
	rawData = readTextFile(fileName)
	nEvents = rawData[0]['nevents']
	rawRawData = []
	for wave in waveSet:
		for key in rawData[0].iterkeys():
			if wave == key.strip():
				rawRawData.append(rawData[0][key])		

	coma = []
	for i in range(len(waveSet)): # Factor 2 from re/im
		coma.append([0.]*(2*len(waveSet)))
		coma.append([0.]*(2*len(waveSet)))
	data = [0.]*(2*len(waveSet))

	for i in range(len(waveSet)):
		point = rawRawData[i]
		data[2*i  ] = point[0]*(float(nEvents)**.5)
		data[2*i+1] = point[1]*(float(nEvents)**.5)
		iii = point[2]
		for j in range(len(rawRawData)):
			jjj = rawRawData[j][2]
			coma[2*i  ][2*j  ] = rawData[1][2*iii  ][2*jjj  ]*nEvents
			coma[2*i+1][2*j  ] = rawData[1][2*iii+1][2*jjj  ]*nEvents
			coma[2*i  ][2*j+1] = rawData[1][2*iii  ][2*jjj+1]*nEvents
			coma[2*i+1][2*j+1] = rawData[1][2*iii+1][2*jjj+1]*nEvents
	return [data,coma]

def own_pinv(matrix, minev =0.):
	"""
	Own implementation of the pseudo inverse for symmetric matrices
	@param matrix: Matrix to invert (list of lists)
	@type matrix: list
	@param minev: Minimal absolute value for an eigenvalue that is not assumed to be zero
	@type minev: float
	@return: Pseudo-inverte matrix (list of lists)
	@rtype: list
	"""
	val,vec = la.eig(matrix)
	dim = len(val)
	new = np.zeros((dim,dim))
	for i in range(dim):
		if abs(val[i]) > minev:
			new[i,i] = 1/val[i]
	ret = np.dot(vec,np.dot(new,np.transpose(vec)))
	return ret


def invert_coma(matrix,method = 'pinv'):
	"""
	Does 'inversion' of covariance matrices
	@param matrix: Matrix to invert (list of lists)
	@type matrix: list
	@param method: Method to be used
		pinv: Pesudo inverse
		diag: Matrix_{ij}^{-1} = 1/Matrix_{ij} if Matrix_{ij} != 0. else Matrix_{ij}^{-1} = 0.
	@type method: str
	@return: Inverted matrix (list of lists)
	@rtype: list
	"""
	dim = len(matrix)
	if method == 'pinv':
		dim = len(matrix)
		matrix_np = np.zeros((dim,dim))
		for i in range(dim):
			for j in range(dim):
				matrix_np[i,j] = matrix[i][j]
#		matrix_pinv = la.pinv(matrix)
		matrix_pinv = own_pinv(matrix_np,1.E-6)
		print matrix_pinv.shape
		ret = []
		for i in range(dim):
			retLine = []
			for j in range(dim):
				retLine.append(matrix_pinv[i,j])
			ret.append(retLine)
		return ret
	if method == 'diag':
		ret = []
		for i in range(dim):
			ret.append([0.]*dim)
		for i in range(dim):
			if not matrix[i][i] == 0.:
				ret[i][i] = 1./matrix[i][i]
		return ret
	raise ValueError # Unknown method

def get_SDM_data_coma(raw_data_coma):
	"""
	Creates SpinDensityMatrix data and coma from the amplitudes
	@param raw_data_coma: Raw data and covariance matrix, output of get_raw_data_coma
	@return: Transformed data and covariance matrix
	"""
	raw_data = raw_data_coma[0]
	raw_coma = raw_data_coma[1]
	nWaves = len(raw_data)/2
	SDM = [0.]*nWaves**2
	jacobian = []
	for i in range(2*nWaves):
		jacobian.append([0.]*nWaves**2)

	for i in range(nWaves):
		for j in range(nWaves):
			one_index = i*nWaves+j
			if i==j: # Intensities
				SDM[one_index] = raw_data[2*i  ]**2 + raw_data[2*i+1]**2
				jacobian[2*i  ][one_index] = 2*raw_data[2*i  ]
				jacobian[2*i+1][one_index] = 2*raw_data[2*i+1]
			if i<j: # Real part
				SDM[one_index] = raw_data[2*i  ]*raw_data[2*j  ] + raw_data[2*i+1]*raw_data[2*j+1]
				jacobian[2*i  ][one_index] =   raw_data[2*j  ]
				jacobian[2*j  ][one_index] =   raw_data[2*i  ]
				jacobian[2*i+1][one_index] =   raw_data[2*j+1]
				jacobian[2*j+1][one_index] =   raw_data[2*i+1]
			if i>j: # Imag Part
				SDM[one_index] = - raw_data[2*j  ]*raw_data[2*i+1] + raw_data[2*i  ]*raw_data[2*j+1]
				jacobian[2*i  ][one_index] =   raw_data[2*j+1]
				jacobian[2*j  ][one_index] = - raw_data[2*i+1]
				jacobian[2*i+1][one_index] = - raw_data[2*j  ]
				jacobian[2*j+1][one_index] =   raw_data[2*i  ]
	SDM_coma = []
	for i in range(nWaves**2):
		SDM_coma.append([0.]*nWaves**2)

	for i in range(nWaves**2):
		for j in range(2*nWaves):
			for k in range(2*nWaves):
				for l in range(nWaves**2):
					SDM_coma[i][l]+=jacobian[j][i]*jacobian[k][l]*raw_coma[j][k]
	return [SDM,SDM_coma]



def write_files(
			target_dir, 
			source_dirs, 
			waveset, 
			upper_lims, 
			lower_lims, 
			Method = 'pinv',
			name_base_data = 'data_', 
			name_base_coma = 'coma_'):
	"""
	Writes data and coma files for fit results in the source_dirs (one for each t' bin) for waveset to  target_dir (Selects right method)
	@param target_dir: Dorectory to write the files in
	@type: target_dir: str
	@param souce_dir: Directory of the mass independent fit
	@type source_dir: str
	@param waveset: List of wave names
	@type waveset: list
	@param upper_lims: Upper mass limits for the waves, list of floats
	@type upper_lims: list
	@param lower_lims: Lower mass limits for the waves, list of floats
	@type lower_lims: list
	@param method: Select the method
		full_covariance: Use the pseudo inverse of the covariance matrix
		anchor_t: Use just interferences with the first wave: Anchor wave
		old_method: Use just the diagonal errors
		full_diagonal: Use just the diagonal errors but with the full covariance matrix (for test only)
	@type method: str
	@param name_base_data: Some name base for the data files
	@type name_base_data: str
	@param name_base_coma: Some name base for the covariance matrix files
	@type name_base_coma: str
	@return: Names of the created files
	@rtype: str
	"""
	i=0
##################################################################
# Chose method by only one flag
##################################################################
	if not Method in ["old_method","full_covariance","full_diagonal","anchor_t"]:
		print "Method not defined"
		raise ValueError
	if Method == "old_method":
		method = "diag"
		FULL_COMA = False
	elif Method == "full_diagonal":
		method = "diag"
		FULL_COMA = True
	elif Method == "full_covariance":
		method = "pinv"
		FULL_COMA = True
	elif Method == "anchor_t":
		method = "anch"
##################################################################
	dataNames=[]
	comaNames=[]
	for source_dir in source_dirs:
		if method == 'anch':
			dataName, comaName = write_anchor(
						target_dir,
						source_dir, 		# directory of the fit
						str(i),
						waveset,			# wave names to be used
						upper_lims,			# mass limits for the waves
						lower_lims,
						name_base_data,
						name_base_coma)
		else:
			dataName,comaName = write_files_tbin(
						str(i),
						target_dir, 
						source_dir, 
						waveset, 
						upper_lims,
						lower_lims, 
						method = method, 
						name_base_data = name_base_data, 
						name_base_coma = name_base_coma,
						FULL_COMA = FULL_COMA)
		dataNames.append(dataName)
		comaNames.append(comaName)
		i+=1
	return dataNames, comaNames

def write_files_tbin(
			indentifier,
			target_dir,
			source_dir,
			waveset,
			upper_lims,
			lower_lims,
			method = 'pinv',
			name_base_data = 'data_', 
			name_base_coma = 'coma_',
			FULL_COMA = True):
	"""
	Writes data and coma files for one t' bin
	@param identifier: Some name to identify the differnt t' bins. 
	@type identifier: str
	@param target_dir: Dorectory to write the files in
	@type: target_dir: str
	@param souce_dir: Directory of the mass independent fit
	@type source_dir: str
	@param waveset: List of wave names
	@type waveset: list
	@param upper_lims: Upper mass limits for the waves, list of floats
	@type upper_lims: list
	@param lower_lims: Lower mass limits for the waves, list of floats
	@type lower_lims: list
	@param method: Method to use:
		pinv: Use paseudo inverse of the covariance matrix
		diag: Invert only th diagonal errors
	@type method: str
	@param name_base_data: Some name base for the data files
	@type name_base_data: str
	@param name_base_coma: Some name base for the covariance matrix files
	@type name_base_coma: str
	@param FULL_COMA: False: Write only diagonal elements of the inverted coma
	@type FULL_COMA: bool
	@return: Names of the files created
	@rtype: str
	"""
	file_list = getBestFits(source_dir)
	file_list.sort(key = lambda x:x[1])
	nWaves = len(waveset)	
	nFile = 0
	nWrite= 0
	dataName = target_dir+os.sep+name_base_data+indentifier
	comaName = target_dir+os.sep+name_base_coma+indentifier
	with open(dataName,'w') as out_data:
		with open(comaName,'w') as out_coma:
			for fn in file_list:
				nFile+=1
				m_low = fn[1]
				m_up  = fn[2]
				m = (m_up+m_low)/2.
				SDM_data_coma = get_SDM_data_coma(get_raw_data_coma(fn[0], waveset))
				for i in range(nWaves): # Apply mass limits
					for j in range(nWaves):
						index = i*nWaves+j
						if upper_lims[i] < m or upper_lims[j] < m or lower_lims[i] > m or lower_lims[j] > m:
							SDM_data_coma[0][index] = 0.
							for index2 in range(nWaves**2):
								SDM_data_coma[1][index ][index2] = 0.
								SDM_data_coma[1][index2][index ] = 0.
				for val in SDM_data_coma[0]:
					out_data.write(str(val)+' ')
				out_data.write('\n')
				inverted = invert_coma(SDM_data_coma[1],method)
				written = 0
				for i in range(nWaves**2):
					if FULL_COMA:
						for j in range(nWaves**2):
							written +=1
							if not inverted[i][j].imag == 0.:
								print i,j,inverted[i][j]
								raise ValueError # Non Real Coma
							out_coma.write(str(inverted[i][j].real)+' ')
							nWrite+=1
					else:
						written+=1
						if not inverted[i][i].imag == 0.:
							print i,i,inverted[i][i]
							raise ValueError # Non Real Coma
						out_coma.write(str(inverted[i][i].real)+' ')
						nWrite+=1
				out_coma.write('\n')
	return dataName, comaName

def is_tprime_dir(string):
	"""
	Checks if a foldername cooresponds to a t' bin as Dima's program does
	@param string: Name of the folder to check
	@type string: str
	@rtype: bool
	"""
	try: 
		float(string.split('-')[0])
	except ValueError:
		return False
	try:
		string.split('-')[1]
		try:
			float(string.split('-')[1])
		except ValueError:
			return False
	except IndexError:
		return False
	return True


if __name__ == "__main__":

	sources = []
#	sourcedir = '/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_22tbins/'
	sourcedir = '/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/'
#	sourcedir = '/nfs/mds/user/fkrinner/massIndepententFits/fits/deck_fabi_thresh/fit/'
	for fn in os.listdir(sourcedir):
		if is_tprime_dir(fn):
			sources.append(sourcedir+os.sep+fn)
	sources.sort()



#	waves = [	"1-(1++)0+ rho pi S",
#			"1-(1-+)1+ rho pi P",
#			"1-(1++)0+ rho pi D",
#			"1-(0-+)0+ f0(980) pi S",
#			"1-(2++)1+ f2 pi P",
#			"1-(2++)1+ rho pi D"]#,
#			"1-(2++)2+ rho pi D",
#			"1-(2-+)0+ f2 pi D",
#			"1-(2-+)0+ f2 pi S",
#			"1-(2-+)0+ rho pi F",
#			"1-(2-+)1+ f2 pi S",
#			"1-(4++)1+ f2 pi F",
#			"1-(4++)1+ rho pi G"]

	waves = [	"1-(3++)1+ rho pi G",
			"1-(3++)0+ rho3 pi S",
			"1-(4++)1+ f2 pi F"]



#	uppers_13 = [2.3,2.0,2.1,2.3,2.0,2.0,2.0,2.3,2.3,2.1,2.3,2.3,2.3 ]
#	lowers_13 = [0.9,0.9,0.9,1.2,1.0,0.9,1.0,1.6,1.4,1.2,1.4,1.4,1.25]

#	uppers = [uppers_13[0],uppers_13[5]]
#	lowers = [lowers_13[0],lowers_13[5]]

	uppers=[2.4,2.4,2.3]
	lowers=[1.0,1.5,1.25]

	write_files(	"/nfs/mds/user/fkrinner/data_mdep_fit/3pp_fits",
			sources,
			waves,
			uppers,
			lowers,
			"anchor_t",
			"data_testulon_",
			"coma_testulon_")
	


	
