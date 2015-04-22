import sys
from math import pi, atan2
from rootpy.plotting import HistStack, Hist, Hist2D, Canvas
from rootpy.io import root_open
from rootpy import ROOT, asrootpy
from convert_text_output import getBestFits, readTextFile

class spinDensityMatrixBin:
	"""
	Class holding a single-bin Spin Density Matrix
	"""
	def __init__(self,wave_names,fit_file,mMin,mMax):
		"""
		@param: wave_names: List of wave names as in the mass-independet fit
		@type wave_names: list
		@param: fit_file: Fit-result-file to generate the SDM from
		@type fit_file: str
		@param mMin: Lower mass limit for the fit-result
		@type mMin: float
		@param: mMax: Upper mass limit for the fit-result
		@type: float
		"""
		self.mMin = mMin
		self.mMax = mMax
		nWaves = len(wave_names)
		self.phases = [[0.]*nWaves for i in range(nWaves)]
		self.errors_phase = [[0.]*nWaves for i in range(nWaves)]

		self.re   = [[0.]*nWaves for i in range(nWaves)]
		self.errors_re = [[0.]*nWaves for i in range(nWaves)]

		self.im   = [[0.]*nWaves for i in range(nWaves)]
		self.errors_im = [[0.]*nWaves for i in range(nWaves)]


		self.intens = [0.]*nWaves
		self.errors_intens = [0.]*nWaves
		
		result = readTextFile(fit_file)
		nEvents= result[0]['nevents']
		amplitudes = [None]*nWaves
		for key in result[0].iterkeys():
			for iWave, wave in enumerate(wave_names):	
				if wave in key:
					amplitudes[iWave] = result[0][key]
		for i in range(nWaves):
			for j in range(nWaves):
				if i==j:
					self.intens[i] = (amplitudes[i][0]**2 + amplitudes[j][1]**2)*nEvents
					index = amplitudes[i][2]
					jac = [2*amplitudes[i][0],2*amplitudes[i][1]]
					coma= [ [result[1][2*index  ][2*index  ],result[1][2*index+1][2*index  ]],
						[result[1][2*index  ][2*index+1],result[1][2*index+1][2*index+1]]]	
					error = 0.
					for ii in range(2):
						for jj in range(2):
							error+=jac[ii]*jac[jj]*coma[ii][jj]


					self.errors_intens[i] = error**.5*nEvents
				
				index1 = amplitudes[i][2]
				index2 = amplitudes[j][2]
				coma = [[result[1][2*index1  ][2*index1  ],result[1][2*index1  ][2*index1+1],result[1][2*index1  ][2*index2  ],result[1][2*index1  ][2*index2+1]],
					[result[1][2*index1+1][2*index1  ],result[1][2*index1+1][2*index1+1],result[1][2*index1+1][2*index2  ],result[1][2*index1+1][2*index2+1]],
					[result[1][2*index2  ][2*index1  ],result[1][2*index2  ][2*index1+1],result[1][2*index2  ][2*index2  ],result[1][2*index2  ][2*index2+1]],
					[result[1][2*index2+1][2*index1  ],result[1][2*index2+1][2*index1+1],result[1][2*index2+1][2*index2  ],result[1][2*index2+1][2*index2+1]]]
				re1 = amplitudes[i][0]
				im1 = amplitudes[i][1]
				re2 = amplitudes[j][0]
				im2 = amplitudes[j][1]
				re=re1*re2+im1*im2
				im=im1*re2-im2*re1
				phase=atan2(im,re)

				reJac=[re2,im2,re1,im1]
				imJac=[-im2,re2,im1,-re1]
				phaseJac=[0.,0.,0.,0.]
				for ii in range(0,4):
					try:
						phaseJac[ii]=1./(1.+im**2./re**2.)*(1./re*imJac[ii]-im/re**2.*reJac[ii])
					except ZeroDivisionError:
						phaseJac[ii]=0.
				phaseErr = 0.
				reErr = 0.
				imErr = 0.
				for iii in range(4):
					for jjj in range(4):
						phaseErr+= phaseJac[iii]*phaseJac[jjj]*coma[iii][jjj]
						reErr+=reJac[iii]*reJac[jjj]*coma[iii][jjj]
						imErr+=imJac[iii]*imJac[jjj]*coma[iii][jjj]
				self.phases[i][j] = phase
				self.im[i][j] = im*nEvents
				self.re[i][j] = re*nEvents
				try:
					if not self.phases[i][j] ==0.:
						self.errors_phase[i][j] = phaseErr**.5
				except ValueError:
					print "Waning, negative varince:",self.errors_phase[i][j]
				try:
					if not self.re[i][j] == 0.:
						self.errors_re[i][j] = reErr**.5*nEvents
				except ValueError:
					print "Waning, negative varince:",self.errors_re[i][j]
				try:
					if not self.im[i][j] ==0.:
						self.errors_im[i][j] = imErr**.5*nEvents
				except ValueError:
					print "Waning, negative varince:",self.errors_im[i][j]

	def center(self):
		"""
		Cets the bin-center of the described mass bin
		@return: Mass-bin center
		@rtype: float
		"""
		return (self.mMin+self.mMax)/2.



class spinDensityMatrix:
	"""
	Class describing a binned spin-density matrix from a mass-independent fit
	"""
	def __init__(self,wave_names,fit_folder):
		"""
		@param wave_names: List of wave names as in the mass-independet fit
		@type wave_names: list
		@param fit_folder: Folder, where the mass-independent fit is located
		@type fit_fodler: str
		"""
		self.wave_numbers = {}
		for i_wave,wave in enumerate(wave_names):
			self.wave_numbers[wave] = i_wave
			
		self.bins = []
		self.nWaves = len(wave_names)
		best_files = getBestFits(fit_folder)
		for best_file in best_files:
			self.bins.append(spinDensityMatrixBin(wave_names,best_file[0],best_file[1],best_file[2]))

	def get_mmin_mmax_nbin(self):
		"""
		Gets the binning of the spin-density matrix
		@return: (minimum mass, maximum mass, number of bins)
		@rtype mmin: float
		@rtype mmax: float
		@rtype nbins: int
		"""
		mmin = float('inf')
		mmax = float('-inf')
		for bin in self.bins:
			mmin = min(mmin,bin.mMin)
			mmax = max(mmax,bin.mMax)
		nbins = len(self.bins)
		return mmin,mmax,nbins


	def get_hist(self,typ,wave1, wave2 = None):
		"""
		Creates a ROOT histogram for the gibten waves and type
		@param typ: Type of the histogram:
			intensity: Intensity of wave1 (No wave2 needed)
			real: Real part of the intereference of wave1 and wave2
			imag: Imaginary part of the intereference of wave1 and wave2
			phase: Complex phase of the intereference of wave1 and wave2
		@type typ: str
		@param wave1: Number of the first wave
		@type wave1: int
		@param wave2: Number of the second wave
		@type wave2: int
		@return: Chosen ROOT histogram		
		@rtype: Hist
		"""
		mmin,mmax,nbin = self.get_mmin_mmax_nbin()
		hist = Hist(nbin,mmin,mmax)
		for bin in self.bins:
			mass = bin.center()
			n_bin = hist.FindBin(mass)
			if typ == "intensity":
				hist.SetBinContent(n_bin,bin.intens[wave1])
				hist.SetBinError(n_bin,bin.errors_intens[wave1])
			elif typ == "real":
				hist.SetBinContent(n_bin,bin.re[wave1][wave2])
				hist.SetBinError(n_bin,bin.errors_re[wave1][wave2])
			elif typ == "imag":
				hist.SetBinContent(n_bin,bin.im[wave1][wave2])
				hist.SetBinError(n_bin,bin.errors_im[wave1][wave2])
			elif typ == "phase":
				hist.SetBinContent(n_bin,bin.phases[wave1][wave2]*180./pi)
				hist.SetBinError(n_bin,bin.errors_phase[wave1][wave2]*180./pi)
			else:
				raise KeyError # Type not defined
		hist.SetTitle("mass-independent")
		hist.SetName("mass-independent")
		return hist

if __name__ == "__main__":
	"""
	Setup for test
	"""
	waves = [	"1-(4++)1+ f2 pi F",
			"1-(3++)1+ rho pi G",
			"1-(3++)0+ rho3 pi S"]


	SDM = spinDensityMatrix(waves,'/nfs/mds/user/fkrinner/massIndepententFits/fits/deck_fabi4/fit/0.100000-0.112853')

	hist = SDM.get_hist("phase",0,1)
	hist.Draw()
	raw_input()

