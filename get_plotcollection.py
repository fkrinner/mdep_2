import os
import sys
import yaml
import math
sys.path.append("build")
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared_retry/plotting/pycpwa/src/')

from libchi2py import chi2
from rootpy.plotting import HistStack, Hist, Hist2D, Canvas
from rootpy.io import root_open
from rootpy import ROOT, asrootpy
from pycpwa.spindensityplots.plotcollection import PlotCollection, TBinPlots
from mass_independent_result import spinDensityMatrix
from pycpwa.definitions.wave import Wave
from ROOT import TH1D
import numpy as np


def get_fit_histogram(chi2,typ,tbin,wave1, wave2 = None, mmin = .5, mmax = 2.5, nbin = 500,component = -1, parameters = None):
	"""
	Gets a single histogram from a chi2 object
	@param chi2: Chi2 from where the histogrm is built
	@type chi2: chi2
	@param typ: type of the histogram:
		intensity: Intensity for wave1
		amp_re: Real part of wave1
		amp_im: Imag part of wave1
		amp_phase: Phase of wave1
		real: Real part of the interference of wave1 and wave2
		imag: Imag part of the interference of wave1 and wave2
		phase: Phase of the interference of wave1 and wave2
	@type param: str
	@param wave1: Number of the first wave
	@type wave1: int
	@param wave2: Number of the first wave (only for interferences, ignored otherwise)
	@type wave2: int
	@param mmin: Minimum mass for the histogram
	@type mmin: float
	@param mmax: Maximum mass for the histogram
	@type mmax: float
	@param nbin: Number of bins for the histogram
	@type nbin: int
	@param component: Number of active component, all active if component == -1 (default)
	@type component: int
	@return Specified histogram
	@rtyp Hist
	"""
	nWaves = chi2.nWaves()
	nTbin= chi2.nTbin()
	perT = chi2.nBrCpl()
	nCpl = perT*nTbin


	hist = Hist(nbin,mmin,mmax)
	if typ.startswith("all_waves"):
		if not component == -1:
			raise ValueError # Not possible with specified component
		hist = []
		for i in range(nWaves):
			histLine = []
			for j in range(nWaves):
				histLine.append(Hist(nbin,mmin,mmax))
			hist.append(histLine)

	step = (mmax - mmin)/nbin
	if not parameters:
		parameters = chi2.fullParameters()
	print '-----------------------------------------------------'
	print '-----------------------------------------------------'
	print tbin,wave1,component
	print '-----------------------------------------------------'
	print parameters[:2*nTbin*perT]
	if not component == -1:
		params_to_keep = chi2.getFuncParameters(component)
		for i in range(2*nCpl):
			if not i in params_to_keep:
				parameters[i] = 0.
	print parameters[:2*nTbin*perT]
	for i in range(nTbin):
		if not i == tbin:
			for j in range(perT):
				parameters[2*(j+i*perT)  ] = 0. # Set eveything in the wrong tbin to zero, just to be shure
				parameters[2*(j+i*perT)+1] = 0.
	print '-----------------------------------------------------'
	print parameters[:2*nTbin*perT]

	for i in range(1,nbin+1):
		m = mmin+(i-.5)*step
		amps = chi2.Amplitudes(m,tbin, parameters)

		if typ == "intensity":
			hist.SetBinContent(i,abs(amps[wave1])**2)
		elif typ == "ampl_real":
			hist.SetBinContent(i,amps[wave1].real)
		elif typ == "ampl_imag":
			hist.SetBinContent(i,amps[wave1].imag)
		elif typ == "ampl_phase":
			hist.SetBinContent(i,math.atan2(amps[wave1].imag,amps[wave1].real)*180./math.pi)
		elif typ.startswith("all_waves"):
			for k in range(nWaves):
				for l in range(nWaves):
					interference = amps[k]*amps[l].conjugate()
					if typ.endswith("real"):
						hist[k][l].SetBinContent(i,interference.real)
					if typ.endswith("imag"):
						hist[k][l].SetBinContent(i,interference.imag)
					if typ.endswith("phase"):
						hist[k][l].SetBinConetnt(i,math.atan2(interference.imag,interference.real)*180./math.pi)
		else:
			interference = amps[wave1]*amps[wave2].conjugate()
			if typ =="real":
				hist.SetBinContent(i,interference.real)
			elif typ == "imag":
				hist.SetBinContent(i,interference.imag)
			elif typ == "phase":
				phase=math.atan2(interference.imag,interference.real)*180./math.pi
				hist.SetBinContent(i,phase)
	if component == -1:
		name = "mass-dependent"
	else:
		name = chi2.get_component_name(component)
	if not typ.startswith("all_waves"):
		hist.SetTitle(name)
		hist.SetName(name)	
	else:
		for i in range(nWaves):
			for j in range(nWaves):
				hist[i][j].SetTitle(name)
				hist[i][j].SetName(name)	
	return hist

def get_mcmc_value(chi2, tbin, wave, component=-1,mmin=.5,mmax=2.5,nPoints = 1000, nEval = 1000, nStepPerEval = 100):
	"""
	Gets the mcmc integral value with error for the specified component and t' bin
	@param chi2: Chi2 to be evaluated
	@type chi2: chi2
	@param tbin: t' bin to be used
	@type tbin: int
	@param wave: Number of the wave to be used
	@type wave: int
	@param component: Component to be used
	@type component: int
	@param mmin: Lower integral limit
	@type mmin: float
	@param mmax: Upper integral limit
	@type mmax: float
	@param nPoints: Number of points used for integration
	@type nPoints: int
	@param nEval: Number of evaluations of the Integral
	@type nEval: int
	@param nStepPerEval: Number of mcmc steps between the integral evaluations
	@type nStepPerEval:int
	"""
	nCpl = chi2.nCpl()
	nTbin = chi2.nTbin()
	for t in range(nTbin):
		if not tbin == t:
			chi2.setEvalTbin(t,False)
	

	params_to_vary = None # Still to write
	parameters = Chi2.parameters()
	startpar = []
	for i in params_to_vary:
		startpar.append(parameters[i])
	def EvalFunc(par):
		for i,j in enumerate(params_to_vary):
			parameters[j] = par[i]
		return chi2.Eval(parameters)
		
	markov = markov(startpar, EvalFunc)
	ints = []
	for i in range(nEval):
		for j in range(nStepPerEval):
			markov.step()
		ints.append(get_integral_value(chi2,tbin,wave,component,mmin,mmax,nPoints,parameters))
	mean = 0.
	for val in ints:
		mean+=val/nEval
	err = 0.
	for val in ints:
		err+=(val-mean)**2/(nEval-1)
	err**=.5
	for t in range(mTbin):
		chi2.setEvalTbin(t,True)
	return mean,err



def get_integral_value(chi2,tbin,wave,component=-1,mmin=.5,mmax=2.5,nPoints=1000, parameters=None):
	"""
	Returns integral for one specified component and t' bin
	@param chi2: Chi2 to be evaluated
	@type chi2: chi2
	@param tbin: t' bin to be used
	@type tbin: int
	@param wave: Number of the wave to be used
	@type wave: int
	@param component: Component to be used
	@type component: int
	@param mmin: Lower integral limit
	@type mmin: float
	@param mmax: Upper integral limit
	@type mmax: float
	@param nPoints: Number of points used for integration
	@type nPoints: int
	@param parameters: List pf parameters to bes used
	@type paramters: list
	@return: Integrated intensity of the components
	@rtype: float
	"""
	oneHist = get_fit_histogram(chi2,'intensity',tbin,wave,mmin=mmin,mmax=mmax,nbin=nPoints,component = component)
	value = 0.
	name = oneHist.GetTitle()
	for bin in range(1,oneHist.GetNbinsX()+1):
		value+=oneHist.GetBinContent(bin)
	return value

def get_fit_t_dependece(chi2, t_binning, wave, component = -1, mmin = .5, mmax = 2.5, nPoints = 1000):
	"""
	Returns t' dependence histogram specified
	@param chi2: Chi2 to be evaluated
	@type chi2: chi2
	@param tbin: t_binning bin to be used
	@type tbin: list
	@param wave: Number of the wave to be used
	@type wave: int
	@param component: Component to be used
	@type component: int
	@param mmin: Lower integral limit
	@type mmin: float
	@param mmax: Upper integral limit
	@type mmax: float
	@param nPoints: Number of points used for integration
	@type nPoints: int
	@return: t' dependence histogram
	@rtype: Hist
	"""
	nTbin = len(t_binning)-1
	if not nTbin == chi2.nTbin():
		raise IndexError # Number of tBins does not match
	values = []
	name = ''
	for tbin in range(nTbin):
		values.append(get_integral_value(chi2,tbin,wave,component,mmin,mmax,nPoints))
	hist = TH1D(name,name,nTbin,np.asarray(t_binning,dtype = np.float64))
	hist.SetTitle(name)
	hist.SetName(name)
	for i in range(nTbin):
		hist.SetBinContent(i+1,values[i]/(t_binning[i+1]-t_binning[i]))
	with root_open("samuel.root","RECREATE"):
		hist.Write()
	return hist


def get_plotcollection(chi2, mmin = .5, mmax = 2.5, nbins = 500):
	"""
	Creates the PlotCollection for a given chi2
	@param chi2: Chi2 from which the collection is produced
	@type chi2: chi2
	@param mmin: Minimum mass for the histograms
	@type mmin: float
	@param mmax: Maximum mass for the histograms
	@type mmax: float
	@param nbin: Number of bins for the histograms
	@type nbin: int
	@return: PlotCollection for chi2
	@rtype PlotCollection
	"""
	collection = PlotCollection()
	wave_names = chi2.waveNames()
	nTbin = chi2.nTbin()
	nWaves = chi2.nWaves()
	YAML_file = chi2.YAML_file()
	borders = chi2.borders_waves()
	LaTeX_namemap = {}
	limits_map = {}
	component_types = {"mass-dependent":1,"mass-independent":2}
	with open(YAML_file,'r') as inin:
		card = yaml.load(inin)
	parametrization_file = card["parametrization_file"]
	if "<card_path>" in parametrization_file:
		parametrization_file = parametrization_file.replace("<card_path>",os.path.dirname(YAML_file)+os.sep)
	with open(parametrization_file, 'r') as inin:
		parametrizations = yaml.load(inin)
	tbin_borders = []
	for i_tbin in range(nTbin):
		fit_folder = card["mass_independent_fit"][i_tbin]
		result = spinDensityMatrix(wave_names,fit_folder)
		tstring = ''
		i_tstring = 1
		while tstring == '':
			tstring = fit_folder.split('/')[-i_tstring]
			i_tstring+=1
		tbin_lower = float(tstring.split('-')[0])
		tbin_upper = float(tstring.split('-')[1])

		if not tbin_lower in tbin_borders:
			tbin_borders.append(tbin_lower)
		if not tbin_upper in tbin_borders:
			tbin_borders.append(tbin_upper)

		tbin_plots = TBinPlots(tbin_lower,tbin_upper)
		hists_re = get_fit_histogram(chi2,"all_waves_real",i_tbin,None,None,mmin=mmin,mmax=mmax,nbin=nbins)
		hists_im = get_fit_histogram(chi2,"all_waves_real",i_tbin,None,None,mmin=mmin,mmax=mmax,nbin=nbins)
		hists_ph = get_fit_histogram(chi2,"all_waves_real",i_tbin,None,None,mmin=mmin,mmax=mmax,nbin=nbins)

		for i_wave, wave in enumerate(wave_names):
			print wave
			hist = result.get_hist("intensity",i_wave)
			act_wave = Wave(wave)
			LaTeX_namemap[wave] = act_wave.toRootLatex()
			limits_map[wave]= (chi2.lowerLims()[i_wave],chi2.upperLims()[i_wave])
			tbin_plots.intensity_[act_wave].Add(hist)
			hist = get_fit_histogram(chi2,"intensity",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins)
			tbin_plots.intensity_[act_wave].Add(hist)
			hist = get_fit_histogram(chi2,"ampl_real",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins)
			tbin_plots.ampl_real_[act_wave].Add(hist)
			hist = get_fit_histogram(chi2,"ampl_imag",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins)
			tbin_plots.ampl_imag_[act_wave].Add(hist)
			if i_wave == 0:
				begin = 0
			else:
				begin = borders[i_wave-1]
			stop  = borders[i_wave]
			for component in range(begin,stop):
				hist = get_fit_histogram(chi2,"intensity",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins,component = component)
				param_name = hist.GetTitle()
				component_types[param_name] = 0
				if parametrizations[param_name].has_key("type"):
						if parametrizations[param_name]["type"] == "resonant":
							component_types[param_name] = 3
						elif parametrizations[param_name]["type"] == "non-resonant":
							component_types[param_name] = 4
						else:
							raise ValueError # Unknown component type
				tbin_plots.intensity_[act_wave].Add(hist)
				hist = get_fit_histogram(chi2,"ampl_real",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins,component = component)
				tbin_plots.ampl_real_[act_wave].Add(hist)
				hist = get_fit_histogram(chi2,"ampl_imag",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins,component = component)
				tbin_plots.ampl_imag_[act_wave].Add(hist)
			for j_wave, wave2 in enumerate(wave_names):
				act_wave2 = Wave(wave2)
				print wave2
				# Real
				hist = result.get_hist("real",i_wave,j_wave)
				tbin_plots.real_[act_wave][act_wave2].Add(hist)
				hist = hists_re[i_wave][j_wave]
				tbin_plots.real_[act_wave][act_wave2].Add(hist)
				# Imag
				hist = result.get_hist("imag",i_wave,j_wave)
				tbin_plots.imag_[act_wave][act_wave2].Add(hist)
				hist = hists_im[i_wave][j_wave]
				tbin_plots.imag_[act_wave][act_wave2].Add(hist)
				# Phase
				hist = result.get_hist("phase",i_wave,j_wave)
				tbin_plots.phase_[act_wave][act_wave2].Add(hist)
				hist = hists_ph[i_wave][j_wave]
				tbin_plots.phase_[act_wave][act_wave2].Add(hist)

		collection._tbins.append(tbin_plots)
	tbin_borders.sort()
	collection.setWavesFromIntensities()
	collection.setLatexnames(LaTeX_namemap)
	collection.setMassdependentFitRanges(limits_map)
	collection.setComponentType(component_types)
	return collection

