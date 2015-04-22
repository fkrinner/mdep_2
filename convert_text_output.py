import os
import numpy
import numpy.linalg as la
#------------------------------------------------------------------------------------------------------------------------------------
def readTextFile(inFile): 
	"""
	Reads the text-output 'inFile' and returns [{wave: ...},[[COMA]]]. Works for arbitrary rank.
	@param inFile: name of the input file
	@type inFile: str
	@return: Data and covariance matrix from the input file
	"""
	waves={}
	print "  OPEN: "+inFile
	data=open(inFile,'r')
	print 'Reading data file: \n   - '+inFile
	enterTheMatrix=False
	nextisNevents=False
	nextisMasses=False
	nextisTprime=False
	nextisLike=False
	covarianceMatrix=[]
	prevWaves=[]
	prevRank=1
        FLAT_count = 0 # Handle the FLAT special, since it may contain DECK amplitudes
	for line in data.readlines():
		if line.startswith("'"):
			wave=line[1:61]
			if "FLAT" in wave:
				wave = get_flat_name(inFile,FLAT_count) # determines, whether the FLAT is really a Deck...
				FLAT_count+=1
			values = get_values(line.split("'")[-1])
			re=values[0][0]
			im=values[0][1]
			if prevRank == len(values)-1:# Add omitted entries for ranks != 1 (0.000,0.000)
				rankstr ='R0'+str(len(values))
				for i in range(1,len(values)):
					prevWave = prevWaves[-i]
					waves[prevWave[:60-len(rankstr)]+rankstr] = [0.,0.]	
				if prevRank == 1: # Remove the rankless entry if necessary
					waves[prevWave[:57]+'R01'] = waves[prevWave]	
					waves.pop(prevWave)
			if len(values)== 1:	
				waves[wave]=[re,im]		
			else:									
				for i in range(len(values)):
					rankstr='R0'+str(i+1)
					waves[wave[:60-len(rankstr)]+rankstr] = [values[i][0],values[i][1]]
			prevRank = len(values)
			prevWaves.append(wave)
		if enterTheMatrix:
			modLine=line.replace('-',' -') # Separate enries, where the '-' takes the whitespace
			modLine=modLine.replace('E -','E-') # undo the separation from the line before, if an exponent was affected
			modLine=modLine.strip()
			covarianceMatrixLine = [float(chunk) for chunk in modLine.split()]
			covarianceMatrix.append(covarianceMatrixLine)
		if nextisTprime:
			nextisTprime=False
			waves['tprime']=[float(ttt) for ttt in line.split(';')]
		if nextisMasses:
			nextisMasses=False
			waves['m3Pi']=[float(mmm) for mmm in line.split(';')]
		if nextisNevents:
			nextisNevents=False
			waves['nevents']=int(line)
		if nextisLike:
			nextisLike=  False
			waves['likelihood']=float(line)
		if 'missing rank entries' in line:
			enterTheMatrix=True
		if "t' bin" in line:
			nextisTprime=True
		if 'Mass bin ' in line:
			nextisMasses=True
		if 'Number of events' in line:
			nextisNevents=True
		if 'log(likelihood)' in line:
			nextisLike=True
	nWave = 0
	for wave in prevWaves: #Count the waves..
		try:
			waves[wave].append(nWave)
			nWave+=1
		except KeyError:
			rank =1
			countWTR = 0
			while True:
				countWTR +=1
				try:
					rankstr='R0'+str(rank)
					waves[wave[:60-len(rankstr)]+rankstr].append(nWave)
					nWave+=1
					rank+=1
				except KeyError:
					break
				if countWTR > 1000:
					raise Exception("'while True:' loop seems to be stuck. Rank > 1000 seems not right")
	return [waves,covarianceMatrix]	
#------------------------------------------------------------------------------------------------------------------------------------
def getBestFits(inDirect=os.curdir+os.sep): 
	"""
	Returns a list of files with the best log-likelihoods
	@param inDirect: Derctory of the mass-independent fit
	@type inDirect: str
	@return: list of bes likelihood files
	@rtype: list
	"""
	direct=inDirect
	while True: # Remove // from path so it can be found in the 'bestFits.txt' even if //, /// or ... is given in the path
		directOld=direct
		direct=direct.replace(os.sep+os.sep,os.sep)
		if direct == directOld:
			break
	fits=[]
	foundDirect=False
	foundDirectTot=False
	if not os.path.isfile('bestFits.txt'): #Store results, so if the best likelihoods are determined once, they can be found more easyly
		open('bestFits.txt', 'a').close()
	read=open('bestFits.txt','r')
	for line in read.readlines():
#		print line
		chunks=line.split()
		if chunks[0]=='DIRECTORY:' and foundDirect:
			foundDirect=False
		if foundDirect:
			fits.append([chunks[0],float(chunks[1]),float(chunks[2])])
		if chunks[0]=='DIRECTORY:':
			if chunks[1]==direct:
				print "Found directory in 'bestFits.txt': Do not scan the folders, but take the results from the file."
				foundDirect=True	
				foundDirectTot=True	
	if not foundDirectTot:
		for fn in os.listdir(direct):
			if 'text_fit_' in fn:
				chunks=fn.split('_')
				m3PiMin=float(chunks[5])/1000
				m3PiMax=float(chunks[6])/1000	
				bestFile=getBestLike(direct+os.sep+fn)
				fits.append([direct+os.sep+fn+os.sep+bestFile,m3PiMin,m3PiMax])
		store=open('bestFits.txt','a')
		store.write('DIRECTORY: '+direct+'\n')
		for i in range(0,len(fits)):
			store.write(fits[i][0]+'  '+str(fits[i][1])+'  '+str(fits[i][2])+'\n')
		store.close()
	return fits
#------------------------------------------------------------------------------------------------------------------------------------	
def get_values(string):
	"""
	Extracts the values from '(...,...)(...,...)(...,...) shaped string, where every ... is a float
	@param string: Inpuot string to extract values from
	@type string: str
	@return: List of values
	@rtype: list
	"""
	valst = string.strip()[1:-1] # remove (...) at the end
	chunks = valst.split(')(')
	values = []
	for chunk in chunks:
		values.append([float(bit) for bit in chunk.split(',')])
	return values
#------------------------------------------------------------------------------------------------------------------------------------
def get_flat_name(
			inFile, 
			FLAT_number,	
			search_type = "FIT_RESULT"	):
	"""
	Returns the type of a 'FLAT' wave, since e.g. Deck amplitudes are also named Deck
	@param inFile: Name of the input file
	@type inFile: str
	@param FLAT_number: Number of the appearing FLAT wave in the file
	@type FLAT_number: int
	@param search_type: Selects type of input file: (integral or fit-result file)
	@type search_type: str
	@return: Name of the FLAT wave
	@rtype: str
	"""
	if search_type == "FIT_RESULT":		
		direct = os.path.dirname(inFile)+os.sep+".."+os.sep
	if search_type == "INTEGRALS":
		direct = os.path.dirname(inFile)+os.sep
	card_name=''
	for fn in os.listdir(direct):
		if "card_" in fn or fn == "card.dat":
			card_name = direct+fn
			break

	if card_name=='':
		return 'unnamed_FLAT(no_card)                                        '

	card = open(card_name,'r')
	FLAT_count = 0
	was_FLAT = False
	unnamed_FLAT = False
	for line in card.readlines():
		if was_FLAT:
			if FLAT_count == FLAT_number:
				if line.startswith("*IWAVENAM_NAMDEP"):
					name = line.split("'")[1]
				else:
					name = "FLAT"
					if not unnamed_FLAT:
						name+=str(FLAT_count)
						unnamed_FLAT = True
				while len(name) < 60:
					name+=' '
				return name
			else:
				FLAT_count+=1
		if line.startswith("*IWAVENAM") and "FLAT" in line and not line[0] =="C":
			was_FLAT = True
		else:
			was_FLAT = False
	return 'unnamed_FLAT                                                '
#------------------------------------------------------------------------------------------------------------------------------------	
def getComaData(	waves,		# List of waves
			up,		# Upper limits
			low,		# Lower Limirs for the 'waves'
			direct, 	# Directory
			flagg ='PINV', 	# Mode of matrix-inversion (default: PSEUDO-INVERSE)
			eps=1.E-3, 	# Not needed at the moment
			CONJUGATE = True):
	"""
	Prepares the data to be used by the mass dependent fitter
	@param waves: List of wave names to use
	@type waves: list
	@param up: List of upper mass limits for the waves
	@type up: list
	@param low: List of lower mass limits for the waves
	@type low: list
	@param direct: Directory to get the fit from
	@type direct: str
	@param flagg: Defines the inversion method (default: PINV = Pseudo-inverse)
	@type flagg: str
	@param eps: Numerical limit
	@type eps: float
	@param CONJUGATE:
	@type CONJUGATE: bool
	@return: Data sets and inverted covariance matrices for the waveset
	"""
	raw_data = getRelevantData(waves,direct)
	nWaves = len(waves)
	nBins  = len(raw_data)
	final_comas_inv=[]
	data_points=[]
	if CONJUGATE:
		conj=-1
	else:
		conj=1
	for point in raw_data:
		masss = (point[0]+point[1])/2
		wave_on=[]						# Check, which waves are active in the current bin
		for i in range(len(up)):				# Set all other elements to zero
			if masss > low[i] and masss < up[i]:		# So that the wrong correlations won't be taken into account
				wave_on.append(1.)
			else:
				wave_on.append(0.)
		for i in range(len(wave_on)):
			point[2][2*i  ]*=wave_on[i]
			point[2][2*i+1]*=wave_on[i]	
			for j in range(len(wave_on)):
				point[3][2*i  ][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i  ][2*j+1]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j+1]*=wave_on[i]*wave_on[j]		
#		prettyPrint(point[3])
		jacobian = []
		data_point=[]
		raw_coma=numpy.matrix(point[3])
		for i in range(nWaves**2):
			jacLine =[0.]*2*nWaves
			ii = int(i/nWaves)
			jj = i - nWaves*ii
			if ii == jj: # Intensities
				addpoint = point[2][2*ii]**2 + point[2][2*ii+1]**2
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = 2*point[2][2*ii]   # dInt_i / dRe_i = 2*Re_i
					jacLine[2*ii+1]= 2*point[2][2*ii+1] # dInt_i / dIm_i = 2*Im_i
			elif ii > jj: # Real Part
				addpoint = point[2][2*ii]*point[2][2*jj] + point[2][2*ii+1]*point[2][2*jj+1]
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = point[2][2*jj]	    # dRe_ij / dRe_i = Re_j			
					jacLine[2*jj]  = point[2][2*ii]	    # dRe_ij / dRe_j = Re_i
					jacLine[2*ii+1]= point[2][2*jj+1]   # dRe_ij / dIm_i = Im_j
					jacLine[2*jj+1]= point[2][2*ii+1]   # dRe_ij / dIm_j = Im_i
			else: # Imaginary Part		# POTENTIAL MINUS SIGN HERE ... INVESTIGATE
			#			R_j      *   I_i           -     R_i       *     I_j
				addpoint = point[2][2*jj]*point[2][2*ii+1] - point[2][2*ii]*point[2][2*jj+1]
				addpoint*=conj
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  =-point[2][2*jj+1] *conj  # dIm_ij / dRe_i =-Im_j				
					jacLine[2*jj]  = point[2][2*ii+1] *conj  # dIm_ij / dRe_j = Im_i
					jacLine[2*ii+1]= point[2][2*jj]	  *conj  # dIm_ij / dIm_i = Re_j
					jacLine[2*jj+1]=-point[2][2*ii]   *conj  # dIm_ij / dIm_j =-Re_i
			jacobian.append(jacLine)
		jacobian = numpy.matrix(jacobian)
		final_coma = jacobian*raw_coma*jacobian.T
#		print data_point
#		final_coma_inv = invert_right_submatrix(final_coma,flagg, eps)
		if not flagg == "ANCHOR_FIRST":
			raise ValueError # Other methods not supported at the moment
		final_coma_inv = la.pinv(final_coma).tolist()
		final_comas_inv.append(final_coma_inv)
		data_points.append(data_point)
	return [data_points,final_comas_inv]
#------------------------------------------------------------------------------------------------------------------------------------
def getRelevantData(waves,direct):
	"""
	Gives the T and covariance matrix for all 'waves[i]' in 'direct'
	@param waves: List of waves to use
	@type waves: list
	@param direct: Directory of the mass independent fit-result
	@type direct: str
	@return: Amplitudes of the waves and their covariances
	"""
	points =[]
	wavesStrip=[]
	for wave in waves:
		wavesStrip.append(wave.strip())
	fileList=getBestFits(direct)
	waveNumbers={}
	dataZero = readTextFile(fileList[0][0])
	for wave in dataZero[0].iterkeys():
		if wave.strip() in wavesStrip:
			waveNumbers[wave.strip()] = [wave,dataZero[0][wave][2]]
	for fil in fileList:
		fitData = readTextFile(fil[0])
		mMin = fil[1]
		mMax = fil[2]
		Ts=[]
		coma =[]
		for wave in wavesStrip:
			nevents = fitData[0]['nevents']
			re = fitData[0][waveNumbers[wave][0]][0]*nevents**.5
			im = fitData[0][waveNumbers[wave][0]][1]*nevents**.5
			nn = waveNumbers[wave][1]
			comaLine1=[]
			comaLine2=[]
			for wave2 in wavesStrip:
				mm = waveNumbers[wave2][1]
				comaLine1.append(fitData[1][2*nn][2*mm]*nevents)
				comaLine1.append(fitData[1][2*nn][2*mm+1]*nevents)
				comaLine2.append(fitData[1][2*nn+1][2*mm]*nevents)
				comaLine2.append(fitData[1][2*nn+1][2*mm+1]*nevents)
			Ts.append(re)
			Ts.append(im)
			coma.append(comaLine1)
			coma.append(comaLine2)
		points.append([mMin,mMax,Ts,coma])
	points.sort()
	return points
