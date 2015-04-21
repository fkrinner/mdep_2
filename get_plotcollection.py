import sys
import yaml
sys.path.append("build")
from libchi2py import chi2
from rootpy.plotting import HistStack, Hist, Hist2D, Canvas
from rootpy.io import root_open
from rootpy import ROOT, asrootpy
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared_retry/plotting/pycpwa/src/')
from pycpwa.spindensityplots.plotcollection import PlotCollection, TBinPlots
from spinDensityMatrix import spinDensityMatrix
from pycpwa.definitions.wave import Wave

import math
def get_fit_histogram(chi2,typ,tbin,wave1, wave2 = None, mmin = .5, mmax = 2.5, nbin = 500,component = -1):
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
	@param wave1: Number of the forst wave
	@type wave1: int
	@param wave2: Number of the forst wave (only for interferences, ignored otherwise)
	@type wave2: int
	@param mmin: Minimum mass for the histogram
	@type mmin: double
	@param mmax: Maximum mass for the histogram
	@type mmax: double
	@param nbin: Number of bins for the histogram
	@type nbin: int
	@param component: Number of active component, all active if component == -1 (default)
	@type component: int
	@return Specified histogram
	@rtyp Hist
	"""
	hist = Hist(nbin,mmin,mmax)
	step = (mmax - mmin)/nbin
	parameters = chi2.parameters()
	nTbin= chi2.nTbin()
	nCpl = chi2.nCpl()/nTbin
	if not component == -1:
		for tBin in range(nTbin):
			pass
			for cpl in range(nCpl):
				if not component == cpl:
					parameters[2*nCpl*tBin + 2*cpl  ] = 0.
					parameters[2*nCpl*tBin + 2*cpl+1] = 0.
	for i in range(1,nbin+1):
		m = mmin+(i-.5)*step
		amps = chi2.Amplitudes(m,tbin, parameters)
		if typ == "intensity":
			hist.SetBinContent(i,abs(amps[wave1])**2)
		elif typ == "amp_re":
			hist.SetBinContent(i,amps[wave1].real)
		elif typ == "amp_im":
			hist.SetBinContent(i,amps[wave1].imag)
		elif typ == "amp_phase":
			hist.SetBinContent(i,math.atan2(amps[wave1].imag,amps[wave1].real))
		else:
			interference = amps[wave1]*amps[wave2].conjugate()
			if typ =="real":
				hist.SetBinContent(i,interference.real)
			elif typ == "imag":
				hist.SetBinContent(i,interference.imag)
			elif typ == "phase":
				phase=math.atan2(interference.imag,interference.real)
				hist.SetBinContent(i,phase)
	if component == -1:
		name = "mass-dependent"
	else:
		name = chi2.get_component_name(component)
	hist.SetTitle(name)
	hist.SetName(name)	
	return hist



def get_plotcollection(chi2, mmin = .5, mmax = 2.5, nbins = 500):
	"""
	Creates the PlotCollection for a given chi2
	@param chi2: Chi2 from which the collection is produced
	@type chi2: chi2
	@param mmin: Minimum mass for the histograms
	@type mmin: double
	@param mmax: Maximum mass for the histograms
	@type mmax: double
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
	for i_tbin in range(nTbin):
		result = spinDensityMatrix(wave_names,card["mass_independent_fit"][i_tbin])
		tbin_plots = TBinPlots()
		for i_wave, wave in enumerate(wave_names):
			hist = result.get_hist("intensity",i_wave)
			act_wave = Wave(wave)
			LaTeX_namemap[wave] = act_wave.toRootLatex()
			limits_map[wave]= (chi2.lowerLims()[i_wave],chi2.upperLims()[i_wave])
			tbin_plots.intensity_[act_wave].Add(hist)
			hist = get_fit_histogram(chi2,"intensity",i_tbin,i_wave,mmin=mmin,mmax=mmax,nbin=nbins)
			tbin_plots.intensity_[act_wave].Add(hist)
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
			for j_wave, wave2 in enumerate(wave_names):
				act_wave2 = Wave(wave2)
				# Real
				hist = result.get_hist("real",i_wave,j_wave)
				tbin_plots.real_[act_wave][act_wave2].Add(hist)
				hist = get_fit_histogram(chi2,"real",i_tbin,i_wave,j_wave,mmin=mmin,mmax=mmax,nbin=nbins)
				tbin_plots.real_[act_wave][act_wave2].Add(hist)
				# Imag
				hist = result.get_hist("imag",i_wave,j_wave)
				tbin_plots.imag_[act_wave][act_wave2].Add(hist)
				hist = get_fit_histogram(chi2,"imag",i_tbin,i_wave,j_wave,mmin=mmin,mmax=mmax,nbin=nbins)
				tbin_plots.imag_[act_wave][act_wave2].Add(hist)
				# Phase
				hist = result.get_hist("phase",i_wave,j_wave)
				tbin_plots.phase_[act_wave][act_wave2].Add(hist)
				hist = get_fit_histogram(chi2,"phase",i_tbin,i_wave,j_wave,mmin=mmin,mmax=mmax,nbin=nbins)
				tbin_plots.phase_[act_wave][act_wave2].Add(hist)

		collection._tbins.append(tbin_plots)
	collection.setWavesFromIntensities()
	collection.setLatexnames(LaTeX_namemap)
	collection.setMassdependentFitRanges(limits_map)
	collection.setComponentType(component_types)
	print LaTeX_namemap
	print
	print limits_map
	print
	print component_types
	return collection

