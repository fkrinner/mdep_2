#include<iostream>
#include"waveset.h"
#include"phaseSpace.h"

#include"matrix_utilities.h"
#include"string_utilities.h"


double UPPER_MASS_LIMIT = 2.5;
double LOWER_MASS_LIMIT = 0.5;
int DEFAULT_L =0;

waveset::waveset():
	_nWaves(0),
	_nFuncs(0),
	_globalPs(0),
	_nPar(0),
	_nFtw(0),
	_nPoints(0),
	_nIso(0),
	_maxBinIso(0),
	_nBranch(0),
	_nTbin(1),
	_nBins(1),
	_write_out(false),
	_has_isobars(false),
	_binning(std::vector<double>(2)){
	_YAML_file = "";
	setTbinning(std::vector<std::vector<double> >(2,std::vector<double>(1,0.)));
};
//########################################################################################################################################################
#ifdef USE_YAML
///Constructror from YAML files
waveset::waveset(
							std::string 						card):

	_nWaves(0),
	_nFuncs(0),
	_globalPs(0),
	_nPar(0),
	_nFtw(0),
	_nPoints(0),
	_nIso(0),
	_maxBinIso(0),
	_nBranch(0),
	_nTbin(1),
	_nBins(1),
	_write_out(false),
	_has_isobars(false),
	_binning(std::vector<double>(2)){

	_YAML_file = card;
	setTbinning(std::vector<std::vector<double> >(2,std::vector<double>(1,0.)));
	YAML::Node Ycard   = YAML::LoadFile(card);
	std::string parametrizations = get_relative_path(Ycard["parametrization_file"].as<std::string>(),_YAML_file);
	std::string waves = get_relative_path(Ycard["wave_file"].as<std::string>(),_YAML_file);
	YAML::Node Yparam  = YAML::LoadFile(parametrizations);
	YAML::Node Ywaves  = YAML::LoadFile(waves);
	std::cout<<"waveset::waveset(...): Take parametrizations from: "<<parametrizations<<std::endl;
	std::cout<<"waveset::waveset(...): Take wave-definitions from: "<<waves<<std::endl;
	std::cout<<"waveset::waveset(...): Initialize waveset from YAML file: "<<card<<std::endl;
	std::cout<<"waveset::waveset(...): Load phase space"<<std::endl;
	loadGlobalPhaseSpace(Ycard);
	std::cout<<"waveset::waveset(...): Phase space loaded\nLoad functions"<<std::endl;
	std::map<std::string,int> fMap = loadFunctions(Ycard, Yparam);
	std::cout<<"waveset::waveset(...): Functions loaded\nLoad waves"<<std::endl;
	loadWaves(Ycard, Ywaves);
	std::cout<<"waveset::waveset(...): Waves loaded\nSet parametrizations"<<std::endl;
	loadFtw(Ycard, fMap);
	std::cout<<"waveset::waveset(...): Parametrizations set\nSetBranchings"<<std::endl;
	loadBranchings(Ycard);
	std::cout<<"waveset::waveset(...): Branchings set\nSet binnings"<<std::endl;
	loadBinnings(Ycard);
	std::cout<<"waveset::waveset(...): Binning set"<<std::endl;
	if(_has_isobars){
		std::cout<<"waveset::waveset(...): Try to set isobar binnings"<<std::endl;
		std::string binning_file;
		if (Ycard["isobar_binnings"]){
			binning_file = Ycard["isobar_binnings"].as<std::string>();
		}else{
			std::cerr<<"waveset::waveset(...): Error: No isobar binnings given"<<std::endl;
		};
		std::cout<<"waveset::waveset(...): Open: "<<binning_file<<std::endl;
		YAML::Node Yisob = YAML::LoadFile(binning_file);
		loadIsoBinnings(Ycard,Yisob);
		std::cout<<"waveset::waveset(...): Isobar binnings set"<<std::endl;
	};
	std::cout<<"waveset::waveset(...): Waveset initialized from YAML file"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Gives the amplitudes for all waves and isobar-mass bins
template<typename xdouble>
std::vector<std::complex<xdouble> > waveset::amps(
							const double 						*m,
							const std::complex<xdouble>	 			*cpl,
							const xdouble	 					*par,
							std::vector<std::vector<std::complex<xdouble> > > 	&funcEvals2pi,
							bool 							ignore_limits)				const{

	std::vector<std::complex<xdouble> > funcEval = funcs(m,par,ignore_limits); // Evalulated BW-functions
	std::vector<std::complex<xdouble> > ampl = std::vector<std::complex<xdouble> >(_nPoints,std::complex<xdouble>(0.,0.)); // Vector of final amplitudes
	std::vector<double> ps = phase_space(m);
	std::complex<xdouble> amp; // Actual wave amplitude
	int upBor=0; // Upper limit for function number
	int loBor=0; // Lower limit for function number
	int amplcount=0;
	for (size_t wave =0;wave<_nWaves;wave++){
		int n_iso_bin = _wave_binning_pts[wave];
		loBor = upBor;
		upBor = _borders_waves[wave];
		if ((m[0] >= _lowerLims[wave] and m[0]<_upperLims[wave]) or ignore_limits){
			if(-1==n_iso_bin){
				std::complex<xdouble> amp(0.,0.);
				for (int nFunc = loBor; nFunc<upBor; nFunc++){
					int func = _funcs_to_waves[nFunc]; // Number of function contributing to the actual wave
					amp+=std::complex<xdouble>(ps[wave],0.)*cpl[nFunc]*funcEval[func];
				};
				ampl[amplcount]=amp;
				amplcount++;
			}else{
				for(int bin2=0;bin2<n_iso_bin;bin2++){
					std::complex<xdouble> amp(0.,0.);
					for (int nFunc = loBor; nFunc<upBor;nFunc++){
						int func =  _funcs_to_waves[nFunc];
						int iso  =  _iso_to_waves[nFunc];
						amp+=std::complex<xdouble>(ps[wave],0.)*cpl[nFunc]*funcEval[func]*funcEvals2pi[iso][bin2];
					};
					ampl[amplcount]=amp;
					amplcount++;
				};
			};
		}else{
				amplcount+=abs(n_iso_bin);
		};
	};
	return ampl;
};
template std::vector<std::complex<double> > waveset::amps(const double *m, const std::complex<double> *cpl,const double *par, std::vector<std::vector<std::complex<double> > > &funcEvals2pi, bool ignore_limits) const;
//########################################################################################################################################################
///Returns the function values at m3pi = m and shape parameters par
template<typename xdouble>
std::vector<std::complex<xdouble> > waveset::funcs(
							const double 						*m, // Contains also other possible variables
							const xdouble	 					*par,
							bool							ignore_limits)				const{

	std::vector<std::complex<xdouble> > f = std::vector<std::complex<xdouble> >(_nFuncs);
	int upPar=0;
	int loPar=0;
	for (size_t func=0; func<_nFuncs; func++){
		loPar = upPar;
		upPar = _borders_par[func];
		if ((m[0] >= _funcLowerLims[func] and m[0] < _funcUpperLims[func]) or ignore_limits){ // Only calculate needed functions
			f[func]=_amp_funcs[func]->Eval(m,par+loPar);
		}else{
			f[func]=std::complex<xdouble>(0.,0.);
		};

	};
	return f;
};
template std::vector<std::complex<double> > waveset::funcs(const double *m,const double *par, bool ignore_limits) const;
//########################################################################################################################################################
///Evaluates the isobar paramterizations at ALL masses, so they do not have to be recalculated each time
template<typename xdouble>
std::vector<std::vector<std::complex<xdouble> > > waveset::iso_funcs(
							const xdouble	  					*par		)			const{

	if (_has_isobars){
		std::vector<std::vector<std::complex<xdouble> > > f = std::vector<std::vector<std::complex<xdouble> > >(_nIso,std::vector<std::complex<xdouble> >(_maxBinIso,std::complex<xdouble>(0.,0.)));
		int upPar=0;
		int loPar=0;
		for (size_t iso=0;iso<_nIso;iso++){
			loPar = upPar;
			upPar = _iso_borders_par[iso];
			int nBins = _iso_binning_pts[iso];
			int nBinning = _iso_n_binning[iso];
			for(int bin=0;bin<nBins;bin++){
				double m = (_iso_binnings[nBinning][bin]+_iso_binnings[nBinning][bin+1])/2;
				f[iso][bin]=_amp_isos[iso]->Eval(&m,par+loPar);
			};
		};
		return f;
	};
	return std::vector<std::vector<std::complex<xdouble> > >();
};
template std::vector<std::vector<std::complex<double> > > waveset::iso_funcs(const double *par) const;
//########################################################################################################################################################
///Gives a vector with phase space factors for each wave at m3pi = m
std::vector<double> waveset::phase_space(
								const double 					*m)					const{

	double global_ps = phaseSpace(m[0],_globalPs,0,0.);
	std::vector<double> ps = std::vector<double>(_nWaves);
	for(size_t wave=0;wave<_nWaves;wave++){
		ps[wave]=global_ps*_wave_phase_space[wave](m);
	};
	return ps;
};
//########################################################################################################################################################
///Enable Auto-diff, if needed
#ifdef ADOL_ON
template std::vector<std::complex<adtl::adouble> > waveset::amps(const double* m, const std::complex<adtl::adouble> *cpl,const adtl::adouble *par, std::vector<std::vector<std::complex<adtl::adouble> > > &funcEvals2pi, bool ignore_limits) const;
template std::vector<std::complex<adtl::adouble> > waveset::funcs(const double* m,const adtl::adouble *par, bool ignore_limits) const;
template std::vector<std::vector<std::complex<adtl::adouble> > > waveset::iso_funcs(const adtl::adouble *par) const;
#endif//ADOL_ON
//########################################################################################################################################################
///Adds a wave, increases the internal definitions
size_t waveset::add_wave(){

	_nWaves+=1;
	size_t nBor = _borders_waves.size();
	if (nBor==0){
		_borders_waves.push_back(0);
	}else{
		_borders_waves.push_back(_borders_waves[nBor-1]);
	};
	_L.push_back(DEFAULT_L);
	_upperLims.push_back(UPPER_MASS_LIMIT);
	_lowerLims.push_back(LOWER_MASS_LIMIT);
	_wave_phase_space.push_back(phase_space_func());
	_waveNames.push_back("unnamed_wave");
	_L_iso.push_back(DEFAULT_L);
	_wave_n_binning.push_back(-1);
	_wave_binning_pts.push_back(-1);
	return _nWaves;
};
//########################################################################################################################################################
///Adds a function - Sets internal definitions accordingly
size_t waveset::add_func(
							int 							i){ 	// # of function

	_nFuncs+=1;
	amplitude* new_amp = get_amplitude(i);
	new_amp->setL(DEFAULT_L);
	int nParNon=new_amp->nPar();
	_amp_funcs.push_back(new_amp);
	size_t nBor = _borders_par.size();
	if (nBor == 0){
		_borders_par.push_back(nParNon);
	}else{
		int nUp = _borders_par[nBor-1];
		_borders_par.push_back(nUp+nParNon);
	};
	_funcUpperLims.push_back(UPPER_MASS_LIMIT);
	_funcLowerLims.push_back(LOWER_MASS_LIMIT);
	return _amp_funcs.size();
};
//########################################################################################################################################################
///Adds a new isobar, sets the internal definitions accordingly
size_t waveset::add_iso(
							int 							i){ 	// # of function for isobar
	_nIso+=1;
	amplitude* new_amp = get_amplitude(i);
	new_amp->setL(DEFAULT_L);
	int nParNon = new_amp->nPar();
//	int nPar = nParNon + new_amp->nCon();
	_amp_isos.push_back(new_amp);
	size_t nBor = _iso_borders_par.size();
	if (nBor == 0){
		_iso_borders_par.push_back(nParNon);
	}else{
		int nUp = _iso_borders_par[nBor-1];
//		std::cout<<nUp<<"::::"<<nParNon<<std::endl;
		_iso_borders_par.push_back(nUp+nParNon);
	};
	_iso_binning_pts.push_back(-1);
	_iso_n_binning.push_back(-1);
	return _amp_isos.size();
};
//########################################################################################################################################################
///Adds a function without an isobar to a wave
void waveset::add_func_to_wave(
							int 							wave,
							int 							func){

	size_t border = _borders_waves[wave];
	std::vector<int> new_relations;
	std::vector<int> new_relations_iso;
	size_t nRel = _funcs_to_waves.size();
	for (size_t i=0; i<nRel;i++){
		new_relations.push_back(_funcs_to_waves[i]);
		new_relations_iso.push_back(_iso_to_waves[i]);
		if (i==border-1){
			new_relations.push_back(func);
			new_relations_iso.push_back(-1);
		};
	};
	if (0==new_relations.size()){
		new_relations.push_back(func);
		new_relations_iso.push_back(-1);
	};
	_funcs_to_waves=new_relations;
	_iso_to_waves=new_relations_iso;
	for (size_t i=wave;i<_nWaves;i++){
		_borders_waves[i]+=1;
	};
	updateFuncLims();
	updateFuncSpin();
	updateIsobar();
	updateNftw();
	handle_branchings(wave,func);
};
//########################################################################################################################################################
///Adds a function-isobar pair to a wave
void waveset::add_funcs_to_wave(
							int 							wave,
							int 							func,
							int 							func2){ // func2 == isobar

	size_t border = _borders_waves[wave];
	std::vector<int> new_relations;
	std::vector<int> new_relations_iso;
	size_t nRel = _funcs_to_waves.size();
	for (size_t i=0; i<nRel;i++){
		new_relations.push_back(_funcs_to_waves[i]);
		new_relations_iso.push_back(_iso_to_waves[i]);
		if (i==border-1){
			new_relations.push_back(func);
			new_relations_iso.push_back(func2);
			_has_isobars=true;
		};
	};
	if (0==new_relations.size()){
		new_relations.push_back(func);
		new_relations_iso.push_back(func2);
		_has_isobars=true; // Set falg to indicate isobar parametrization
	};
	_funcs_to_waves=new_relations;
	_iso_to_waves=new_relations_iso;
	for (size_t i=wave;i<_nWaves;i++){
		_borders_waves[i]+=1;
	};
	updateFuncLims();
	updateFuncSpin();
	updateIsobar();
	updateNftw();
	handle_branchings(wave,func);
};
//########################################################################################################################################################
///Couples two functions: {C1(t), C2(t)} -> {B1*C(t), B2*C(t)}
void waveset::couple_funcs(
							int 						i1,	// # of wave 1
							int 						i2){ 	// # of wave 2

	if (i1 == i2){
		std::cerr<<"waveset::couple_funcs(...): Error: Trying to couple two times the same functions"<<std::endl;
		return;
	};
	int branch1 = _coupled[i1];
	int branch2 = _coupled[i2];
	if (branch1 == branch2 and branch1 != -1){
		return;
	};
	if (-1==branch1 and -1==branch2){
		int last_branch = -1;
		int nCpl = _coupled.size();
		for (int i = 0;i<nCpl;i++){
			if (_coupled[i] > last_branch){
				last_branch = _coupled[i];
			};
		};
		last_branch++;
		_coupled[i1] = last_branch;
		_coupled[i2] = last_branch;
	}else if(-1==branch1){
		_coupled[i1] = branch2;
	}else if(-1==branch2){
		_coupled[i2] = branch1;
	}else{
		std::cerr<<"waveset::couple_funcs(...): Error: function["<<i1<<"] and function["<<i2<<"] are already in different branchings"<<std::endl;
	};
	update_n_cpls();
	update_n_branch();
};
//########################################################################################################################################################
///Simple setter for the mass limits of the single waves
void waveset::setWaveLimits(
							int 							i, 	// # of wave
							double 							lower,
							double 							upper){

	_upperLims[i] = upper;
	_lowerLims[i] = lower;
	updateFuncLims();
	update_min_max_bin();
};
//########################################################################################################################################################
///Simple setter for the wave spins
void waveset::setWaveSpin(
							int 							i, 	// # of wave
							int 							L){
	_L[i]=L;
	updateFuncSpin();
};
//########################################################################################################################################################
///Sets the isobar spin for a wave
void waveset::setWaveIsobarSpin(
							int 							wave,
							int 							L){

	_L_iso[wave]=L;
	updateIsobar();
};
//########################################################################################################################################################
///Simple setter for the global phase space factor
void waveset::setGlobalPhaseSpace(
							int 							i){	// # of phase space function

	_globalPs=i;
};
//########################################################################################################################################################
///Simple setter for the wave-specific phase space factors
void waveset::setWavePhaseSpace(
							int 							i,	// # of wave
							std::string						name,	// Wavename
							std::string						integralFile){	//IntegralFile

	_wave_phase_space[i] = phase_space_func(name,integralFile);
};
//########################################################################################################################################################
///Sets the isobar binning for a certain wave
void waveset::setWaveIsobarBinning(
							int 							wave,
							int 							binning){

	_wave_n_binning[wave] = binning;
	int npoints = _iso_binnings[binning].size()-1;
	_wave_binning_pts[wave] = npoints;
	updateIsobar();
};
//########################################################################################################################################################
///Simple setter for the constants used in one function
void waveset::setConst(
							int 							i, 	// # of const
							double 							con){

	int nConst =0;
	for (size_t j =0;j<_nFuncs;j++){
		nConst+=_amp_funcs[j]->nCon();
		if(nConst > i){
			nConst-=_amp_funcs[j]->nCon();
			_amp_funcs[j]->setCon(i-nConst,con);
			break;
		};
	};
};
//########################################################################################################################################################
///Sets constants in isobar paramterizations
void waveset::set_iso_const(
							int 							con, 	// # of const
							double 							value){

	int cCount =0;
	for (size_t i =0;i<_nIso;i++){
		cCount+=_amp_isos[i]->nCon();
		if(cCount>con){
			cCount-=_amp_isos[i]->nCon();
			_amp_isos[i]->setCon(con-cCount,value);
			break;
		};
	};

};

//########################################################################################################################################################
///Adds a new isobar binning
void waveset::add_isobar_binning(
							std::vector<double> 					binning){

	_iso_binnings.push_back(binning);
	size_t nnew = binning.size()-1;
	if (_maxBinIso<nnew){
		_maxBinIso=nnew;
	};
};
//########################################################################################################################################################
///Simple setter for wave names
void waveset::setWaveName(
							int 							i, 	// # of wave
							std::string 						name){

	_waveNames[i] = name;
};
//########################################################################################################################################################
///Simple setter for function name
void waveset::setFunctionName(
							int 							i, 	// # of function
							std::string 						name){

	_amp_funcs[i]->setFunctionName(name);
};
//########################################################################################################################################################
///Simple setter for paremeter name
void waveset::setParameterName(
							int 							i, 	// # of parameter
							std::string 						name){

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_funcs[f]->nPar();
		if (pCount>i){
			pCount-=_amp_funcs[f]->nPar();			
			_amp_funcs[f]->setParName(i-pCount,name);
			break;
		};
	};
};
//########################################################################################################################################################
///Simple setter for paremeter 
bool waveset::setPar(
							int 							i, 	// # of parameter
							double 							val){

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_funcs[f]->nPar();
		if (pCount>i){
			pCount-=_amp_funcs[f]->nPar();			
			return _amp_funcs[f]->setPar(i-pCount,val);
		};
	};
	return false;
};
//########################################################################################################################################################
///Simple getter for paremeter 
double waveset::getPar(
							int 							i)		const{ 	// # of parameter

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_funcs[f]->nPar();
		if (pCount>i){
			pCount-=_amp_funcs[f]->nPar();			
			return _amp_funcs[f]->getParameter(i-pCount);
		};
	};
	return std::numeric_limits<double>::quiet_NaN();
};
//########################################################################################################################################################
///Simple getter for constant 
double waveset::getCon(
							int 							i)		const{ 	// # of const

	int cCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		cCount+=_amp_funcs[f]->nCon();
		if (cCount>i){
			cCount-=_amp_funcs[f]->nCon();			
			return _amp_funcs[f]->getConstant(i-cCount);
		};
	};
	return std::numeric_limits<double>::quiet_NaN();
};
//########################################################################################################################################################
///Simple setter for paremeter 
bool waveset::setIsoPar(
							int 							i, 	// # of parameter
							double 							val){

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_isos[f]->nPar();
		if (pCount>i){
			pCount-=_amp_isos[f]->nPar();			
			return _amp_isos[f]->setPar(i-pCount,val);
		};
	};
	return false;
};
//########################################################################################################################################################
///Simple getter for paremeter 
double waveset::getIsoPar(
							int 							i)		const{ 	// # of parameter

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_isos[f]->nPar();
		if (pCount>i){
			pCount-=_amp_isos[f]->nPar();			
			return _amp_isos[f]->getParameter(i-pCount);
		};
	};
	return std::numeric_limits<double>::quiet_NaN();
};
//########################################################################################################################################################
///Simple getter for constant 
double waveset::getIsoCon(
							int 							i)		const{ 	// # of const

	int cCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		cCount+=_amp_isos[f]->nCon();
		if (cCount>i){
			cCount-=_amp_isos[f]->nCon();			
			return _amp_isos[f]->getConstant(i-cCount);
		};
	};
	return std::numeric_limits<double>::quiet_NaN();
};
//########################################################################################################################################################
///Simple setter for constant name
void waveset::setConstantName(
							int 							i, 	// # of constant
							std::string 						name){

	int nConst =0;
	for (size_t j =0;j<_nFuncs;j++){
		nConst+=_amp_funcs[j]->nCon();
		if(nConst > i){
			nConst-=_amp_funcs[j]->nCon();
			_amp_funcs[j]->setConName(i-nConst,name);
			break;
		};
	};
};
//########################################################################################################################################################
///Simple setter for isobar names
void waveset::setIsobarName(
							int 							i,	// # of isobar
							std::string 						name){

	_amp_isos[i]->setFunctionName(name);
};
//########################################################################################################################################################
///Simple setter for isopar parameter names
void waveset::setIsoParName(
							int 							i,	// # of isobar parameter
							std::string 						name){

	int pCount =0;
	for(size_t iso=0;iso<_nIso;iso++){
		pCount+=_amp_isos[iso]->nPar();
		if(pCount>i){
			pCount-=_amp_isos[iso]->nPar();
			_amp_isos[iso]->setParName(i-pCount,name);
			break;
		};
	};
};
//########################################################################################################################################################
///Simple setter for isobar constant names
void waveset::setIsoConstName(
							int 							i,	// # of isobar constant
							std::string 						name){
	int cCount=0;
	for (size_t iso=0;iso<_nIso;iso++){
		cCount+=_amp_isos[iso]->nCon();
		if (cCount>i){
			cCount-=_amp_isos[iso]->nCon();
			_amp_isos[iso]->setConName(i-cCount,name);
			break;
		};
	};
};
//########################################################################################################################################################
///Sets the binning in m3pi
void waveset::setBinning(
							std::vector<double> 					binning){

	_binning = binning;
	_nBins = binning.size()-1;
	_mMin = binning[0];
	_mMax = binning[_nBins];
	update_min_max_bin();
};
//########################################################################################################################################################
///Sets the binning in t'
void waveset::setTbinning(
							std::vector<std::vector<double> >			binning){

	_t_binning = binning;
	_nTbin = binning.size();
	if (binning.size() >0){
		_nTvar = binning[0].size();
	}else{
		_nTvar = 0;
	};
	if (_eval_tbin.size() != _nTbin){
		_eval_tbin = std::vector<bool>(_nTbin,true);
	};
};
//########################################################################################################################################################
///Switches single t' bins on/off
void waveset::setEvalTbin(
							int 							i,	// # of t' bin
							bool 							flag){
	_eval_tbin[i] = flag;
};
//########################################################################################################################################################
///Sets minimum bin evaluated my hand
void waveset::setMinBin(
							int 							in){

_minBin = in;
};
//########################################################################################################################################################
///Sets maximum bin evaluated my hand
void waveset::setMaxBin(
							int 							in){

_maxBin = in;
};
//########################################################################################################################################################
#ifdef USE_YAML
///Loads the phase space from a YAML file
bool waveset::loadGlobalPhaseSpace(
							YAML::Node 						&waveset){

	bool sgps = false;
	if (waveset["standard_phase_space"]){
		if (waveset["standard_phase_space"].as<bool>()){
			setGlobalPhaseSpace(10);
			sgps=true;
		};
	};
	int nPS = waveset["global_phase_space"].size();
	for (int i=0;i<nPS;i++){
		setGlobalPhaseSpace(waveset["global_phase_space"][i].as<int>());
		sgps=true;
	};
	if (not sgps){
		std::cout<<"waveset::loadGlobalPhaseSpace(...): Warning: No global phase space set" << std::endl;
	};
	return sgps;
};
//########################################################################################################################################################
///Loads the functions from a YAML file
std::map<std::string,int> waveset::loadFunctions(
							YAML::Node						&waveset,
							YAML::Node						&param){

	std::map<std::string,int> fMap;
	int fCount = 1; // Starts as one, since map[nonExistingKey] = 0 -> use this to determine, whether a function exists.
	int iCount = 1;
	int pCount = 0;
	int ipCount= 0;
	int cCount = 0;
	int icCount= 0;
	int nWaves = waveset["waves"].size();
	std::vector<double> pVals;
	std::vector<double> pSteps;
	std::vector<bool> pRels;
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int j=0;j<nFunc;j++){
			std::string fName;
			std::string iName = "none";
			bool isisobar = false;
			if (waveset["waves"][i]["parametrizations"][0].size() != waveset["waves"][i]["parametrizations"][j].size()){
				std::cerr<<"waveset::loadFunctions(...): Error: Isobarred and de-isobarred parametrizations in the same wave"<<std::endl;
			};
			if (0==waveset["waves"][i]["parametrizations"][j].size()){
				fName = waveset["waves"][i]["parametrizations"][j].as<std::string>();
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(2==waveset["waves"][i]["parametrizations"][j].size()) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][j][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][j][1].as<std::string>();
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
			if(param[fName]){
				if (0==fMap[fName]){
					fMap[fName]=fCount;
					int nParam = param[fName]["parametrization"].as<int>();
					size_t newf = add_func(nParam)-1;
					int nPar = _amp_funcs[newf]->nPar();
					int nCon = _amp_funcs[newf]->nCon();
					_amp_funcs[newf]->print();
					setFunctionName(fCount-1,fName);
					fCount++;
					int nParDef = param[fName]["parameters"].size();
					if (nParDef == nPar){
						for(int par=0;par<nPar;par++){
							std::string pName = param[fName]["parameters"][par]["name"].as<std::string>();
							setParameterName(pCount,pName);
							double value = param[fName]["parameters"][par]["value"].as<double>();
							setPar(pCount,value);
							pCount++;
						};
					}else{
						std::cerr<<"waveset::loadFunctions(...): Error: Number of defined parameters ("<<nParDef<<") does not match required number for "<<fName<<" ("<<nPar<<")"<<std::endl;
					};
					int nConDef = param[fName]["constants"].size();
					if(nConDef == nCon){
						for (int con=0;con<nCon;con++){
							std::string cName = param[fName]["constants"][con]["name"].as<std::string>();
							setConstantName(cCount ,cName);
							double value  = param[fName]["constants"][con]["value"].as<double>();
							setConst(cCount,value);
							cCount++;
						};
					}else{
						std::cerr<<"waveset::loadFunctions(...): Error: Number of defined constants ("<<nConDef<<") does not match required number for "<<fName<<" ("<<nCon<<")"<<std::endl;
					};
				};
			}else{
				std::cerr << "waveset::loadFunctions(...): Error: '"<<fName<<"' not defined in the parametrization-file"<<std::endl;
			};
			if (isisobar){
				if(param[iName]){
					if (0==fMap[iName]){ // Write isobars also in the fMap
						fMap[iName]=iCount;
						int nParam = param[iName]["parametrization"].as<int>();
						int nPar = _amp_isos[iCount-1]->nPar();
						int nCon = _amp_isos[iCount-1]->nCon();
						add_iso(nParam);
						setIsobarName(iCount-1,iName);
						iCount++;
						int nParDef = param[iName]["parameters"].size();
						if (nParDef == nPar){
							for (int par=0;par<nPar;par++){
								std::string pName = param[iName]["parameters"][par]["name"].as<std::string>();
								setIsoParName(ipCount,pName);
								ipCount++;
							};
						}else{
							std::cerr<<"waveset::loadFunctions(...): Error: Number of defined parameters does not match required number for "<<iName<<std::endl;
						};
						int nConDef = param[iName]["constants"].size();
						if(nConDef == nCon){
							for (int con=0;con<nCon;con++){
								std::string cName = param[iName]["constants"][con]["name"].as<std::string>();
								setIsoConstName(icCount,cName);
								double value = param[iName]["constants"][con]["value"].as<double>();
								set_iso_const(icCount,value);
								icCount++;
							};
						}else{
							std::cerr<<"waveset::loadFunctions(...): Error: Number of defined constants does not match required number for "<<iName<<std::endl;
						};
					};
				}else{
					std::cerr << "waveset::loadFunctions(...): Error: '"<<iName<<"' not defined in the parametrization-file"<<std::endl;
				};
			};
		};
	};
	return fMap;
};
//########################################################################################################################################################
///Loads the waves from a YAML file
bool waveset::loadWaves(
							YAML::Node						&waveset,
							YAML::Node						&defs){

	bool ookk = true;
	int nWaves = waveset["waves"].size();
//	double mmin = waveset["waves"][0]["mmin"].as<double>();
//	double mmax = waveset["waves"][0]["mmax"].as<double>();
	for (int i=0;i<nWaves;i++){
		add_wave();
		std::string wName = waveset["waves"][i]["name"].as<std::string>();
		setWaveName(i,wName);
		if (not defs[wName]){
			ookk=false;
			std::cerr<<"waveset::loadWaves(...): Error: "<< wName<<" not defined in the wave-definitions file"<<std::endl;
		};
		setWaveSpin(i,defs[wName]["spin"].as<int>());
		if(defs[wName]["isobar_spin"]){
			setWaveIsobarSpin(i,defs[wName]["isobar_spin"].as<int>());
		};
		double mmin_act = waveset["waves"][i]["mmin"].as<double>();
		double mmax_act = waveset["waves"][i]["mmax"].as<double>();
/*		if (mmin_act<mmin){
			ookk=false;
			std::cerr << "Error: Lower limit of '"<<waveset["waves"][i]["name"]<<"' below anchor mass limit"<<std::endl;
		};
		if (mmax_act>mmax){
			ookk=false;
			std::cerr << "Error: Upper limit of '"<<waveset["waves"][i]["name"]<<"' above anchor mass limit"<<std::endl;
		};
*/// Comment this check, since this is not required anymore for some methods
		setWaveLimits(i,mmin_act,mmax_act);
		bool set_phase_space = true; 
		if (waveset["divide_integrals"]){
			if (waveset["divide_integrals"].as<bool>()){
				std::cout<<"'divide_integrals' flag found. Do not set phase_space"<<std::endl;
				set_phase_space = false; // If the phase space is already divided out, do not apply it.
			};
		};
		if (set_phase_space){
			setWavePhaseSpace(i,wName,waveset["integralFile"].as<std::string>());
		};
	};
	return ookk;
};
//########################################################################################################################################################
///Puts funtions and waves together // Does not support isobars at the moment
void waveset::loadFtw(
							YAML::Node						&waveset,
							std::map<std::string,int>				&fMap){

	int nWaves = waveset["waves"].size();
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int func=0;func<nFunc;func++){
			std::string fName;
			std::string iName;
			bool isisobar = false;
			if (0==waveset["waves"][i]["parametrizations"][func].size()){
				fName = waveset["waves"][i]["parametrizations"][func].as<std::string>();
				int nf = fMap[fName]-1;
				add_func_to_wave(i,nf);
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(2==waveset["waves"][i]["parametrizations"][func].size()) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][func][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][func][1].as<std::string>();
				int nf = fMap[fName]-1;
				int ni = fMap[iName]-1;
				add_funcs_to_wave(i,nf,ni);
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
		};
	};
};
//########################################################################################################################################################
///Loads branchings from a YAML file
void waveset::loadBranchings(
							YAML::Node						&waveset){

	int nWaves = waveset["waves"].size();
	std::vector<int> fixed_branchings;
//	int nBr=0;
	int nBranch = waveset["branchings"].size();
	for (int i=0;i<nBranch;i++){
		std::vector<int> wave_numbers;
		int nnn = waveset["branchings"][i].size();
		for (int j=0;j<nnn;j++){
			int func_count = 0;
			std::string wnam_b = waveset["branchings"][i][j][0].as<std::string>();
			std::string fnam_b = waveset["branchings"][i][j][1].as<std::string>();
			for (int wave=0;wave<nWaves;wave++){
				std::string wnam_d = waveset["waves"][wave]["name"].as<std::string>();
				int nfunc = waveset["waves"][wave]["parametrizations"].size();
				for (int func=0;func<nfunc;func++){
					std::string fnam_d;
					if (0==waveset["waves"][wave]["parametrizations"][func].size()){
						fnam_d = waveset["waves"][wave]["parametrizations"][func].as<std::string>();
					}else if(2==waveset["waves"][wave]["parametrizations"][func].size()){
						fnam_d = waveset["waves"][wave]["parametrizations"][func][0].as<std::string>(); // When a coupled function appear with different isobars in one wave, all of them will be "branched"
					};
					if (fnam_d == fnam_b and wnam_d == wnam_b){
						wave_numbers.push_back(func_count);
						};
					func_count++;
				};
			};
		};
		if (wave_numbers.size()<2){
			std::cerr<<"waveset::loadBranchings(...): Error: Branching definition with less than two waves encountered."<<std::endl;
		}else{
			for (size_t cplto = 1; cplto < wave_numbers.size(); cplto++){
				couple_funcs(wave_numbers[0],wave_numbers[cplto]);
			};
		};
	};
};
//########################################################################################################################################################
///Loads binnings (m3Pi and t') from YAML file
void waveset::loadBinnings(
							YAML::Node						&waveset){

	double M = waveset["mmin"].as<double>();
	double m_max = waveset["mmax"].as<double>();
	double binwidth = waveset["binwidth"].as<double>();
	std::vector<double> binning;
	binning.push_back(M);
	while(M<m_max){
		M+=binwidth;
		binning.push_back(M);
	};
	setBinning(binning);
	std::vector<std::vector<double> > tbinning;
	if (waveset["t_binning"]){
		int nTbin = waveset["t_binning"].size();
		for (int i=0; i<nTbin;i++){
			int nvaract = waveset["t_binning"][i].size();
			std::vector<double> act_var;
			for (int jj=0;jj<nvaract;++jj){
				act_var.push_back(waveset["t_binning"][i][jj].as<double>());
			};
			tbinning.push_back(act_var);
		};
		setTbinning(tbinning);
	};
	update_min_max_bin();
};
//########################################################################################################################################################
///Load isobar binnings from YAML file
void waveset::loadIsoBinnings(
							YAML::Node 						&waveset,
							YAML::Node						&binnings){

	std::map<std::string,int> binnings_map;
	int nWaves = waveset["waves"].size();
	int nBinn = 1;// Start with 1 for int maps, to determine, if entry exists
	for (int i=0;i<nWaves;i++){ // Create the binnigns map
		std::string wName = waveset["waves"][i]["name"].as<std::string>();
		int nPar = waveset["waves"][i]["parametrizations"].size();
		if (nPar>0){
			int nPer = waveset["waves"][i]["parametrizations"][0].size();
			if (2==nPer){
				if(waveset["waves"][i]["isobar_binning"]){
					std::string binning_name = waveset["waves"][i]["isobar_binning"].as<std::string>();
					if (0==binnings_map[binning_name]){
						binnings_map[binning_name] = nBinn;
						nBinn++;
						std::vector<double> binning;
						int binning_size = binnings[binning_name].size();
						for (int bin=0;bin<binning_size;bin++){
							binning.push_back(binnings[binning_name][bin].as<double>());
						};
						add_isobar_binning(binning);
					};
					setWaveIsobarBinning(i,binnings_map[binning_name]-1);
				}else{
					std::cerr<<"waveset::loadIsoBinnings(...): Error: De-isobarred wave '"<<wName<<"' does not have an isobar binning"<<std::endl;
				};
			};
		};
	};
	updateNpoints();
};
#endif//USE_YAML
//########################################################################################################################################################
///Returns the name of the funtion for component i
std::string waveset::get_component_name(		size_t							i)	const{	// # of component

	size_t fcn = _funcs_to_waves[i];
	return _amp_funcs[fcn]->name();
};
//########################################################################################################################################################
///Returns the number of parameters (without branchings or couplings)
size_t waveset::getNpar()												const{

	if (0==_borders_par.size()){
		return 0;
	};
	return _borders_par[_borders_par.size() -1];
};
//########################################################################################################################################################
///Returns the number of parameters for isobars
size_t waveset::getNiso()												const{

	size_t siz = _iso_borders_par.size();
	if (0==siz){
		return 0;
	};
	return _iso_borders_par[siz-1];
};
//########################################################################################################################################################
///Returns the name of waves
std::string waveset::getWaveName(
							int 							i)	const{ 	// # of wave

	return _waveNames[i];
};
//########################################################################################################################################################
/// Gets all functions used by wave
std::vector<int> waveset::get_wave_functions(
							int 							wave)	const{

	std::vector<int> wafe;
	if(_borders_waves.size()==0){
		return wafe;
	};
	int upper = _borders_waves[wave];
	int lower=0;
	if (wave>0){
		lower = _borders_waves[wave-1];
	};
	for (int nnn = lower;nnn<upper;nnn++){
		int coil = _funcs_to_waves[nnn];
		bool already = false;
		for (size_t aaa=0;aaa<wafe.size();aaa++){
			if (wafe[aaa] == coil){
				already = true;
				break;
			};
		};
		if (not already){
			wafe.push_back(coil);
		};
	};
	return wafe;
};
//########################################################################################################################################################
///Gets the parameters numbers for wave
std::vector<int> waveset::get_wave_pars(
							int 							wave)	const{

	std::vector<int> pars;
	std::vector<int> funcs = get_wave_functions(wave);
	size_t nFu = funcs.size();
	for (size_t i=0;i<nFu;i++){
		std::vector<int> funcPars = get_function_pars(funcs[i]);
		pars.insert(pars.end(),funcPars.begin(),funcPars.end());
	};
	return pars;
};
//#######################################################################################################################################################
///Gets the constants numbers for wave
std::vector<int> waveset::get_wave_const(
							int 							wave)	const{

	std::vector<int> consts;
	std::vector<int> funcs = get_wave_functions(wave);
	size_t nFu = funcs.size();
	for (size_t i=0;i<nFu;i++){
		std::vector<int> funcConst = get_function_const(funcs[i]);
		consts.insert(consts.end(),funcConst.begin(),funcConst.end());
	};
	return consts;
};
//########################################################################################################################################################
///Gets all isobars used by wave
std::vector<int> waveset::get_wave_isobars(
							int 							wave)	const{

	int upper = _borders_waves[wave];
	int lower=0;
	if (wave>0){
		lower = _borders_waves[wave-1];
	};
	std::vector<int> isos;
	for (int i=lower;i<upper;i++){
		bool already = false;
		for (size_t aaa=0;aaa<isos.size();aaa++){
			if (isos[aaa]==_iso_to_waves[i]){
				already = false;
				break;
			};
		};
		if (_iso_to_waves[i] >-1 and not already){
			isos.push_back(_iso_to_waves[i]);
		};
	};
	return isos;
};
//########################################################################################################################################################
///Gets the isobar paramters for wave
std::vector<int> waveset::get_wave_iso_pars(
							int 							wave)	const{

	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> pars;
	for (size_t i=0;i<isos.size();i++){
		std::vector<int> isopars = get_isobar_pars(isos[i]);
		pars.insert(pars.end(),isopars.begin(),isopars.end());
	};
	return pars;
};
//########################################################################################################################################################
///Gets the isobar constants for wave
std::vector<int> waveset::get_wave_iso_const(
							int 							wave)	const{

	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> consts;
	for (size_t i=0;i<isos.size();i++){
		std::vector<int> isoconst = get_isobar_const(isos[i]);
		consts.insert(consts.end(),isoconst.begin(),isoconst.end());
	};
	return consts;
};
//########################################################################################################################################################
/// Simple getter for function names
std::string waveset::getFunctionName(
							int 							i)	const{

	return _amp_funcs[i]->name();
};
//########################################################################################################################################################
///Gets the numbers of parameters for each func
std::vector<int> waveset::get_nPars()											const{

	std::vector<int> nPars;
	nPars.push_back(_borders_par[0]);
	size_t i_max = _borders_par.size();
	for (size_t i=1;i<i_max;i++){
		int diff = _borders_par[i]-_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
//########################################################################################################################################################
///Gets the numbers of constants for each func
std::vector<int> waveset::get_nConst()											const{

	std::vector<int> nConst;
	for (size_t f=0;f<_nFuncs;f++){
		int nnn = _amp_funcs[f]->nCon();
		nConst.push_back(nnn);
	};
	return nConst;
};
//########################################################################################################################################################
///Gets the paramters numbers for func
std::vector<int> waveset::get_function_pars(
							int 							func)	const{

	std::vector<int> pars;
	if (_borders_par.size()==0){
		return pars;
	};
	int upper = _borders_par[func];
	int lower = 0;
	if (func >0){
		lower = _borders_par[func-1];
	};
	for (int nnn = lower;nnn<upper;nnn++){
		pars.push_back(nnn);
	};
	return pars;
};
//########################################################################################################################################################
///Get the constant numbers for func
std::vector<int> waveset::get_function_const(
							int 							func)	const{

	int before = 0;
	for (int i =0;i<func-1;i++){
		before+=_amp_funcs[i]->nCon();
	};
	std::vector<int> ret;
	for (size_t i=0; i<_amp_funcs[func]->nCon();i++){
		ret.push_back(before+i);
	};
	return ret;
};
//########################################################################################################################################################
///Get the waves that employ the function number 'func'
std::vector<int> waveset::get_function_waves(
							int 							func)	const{

	std::vector<int> waves;
	for (size_t wave=0;wave<_nWaves;wave++){
		std::vector<int> funcs = get_wave_functions(wave);
		size_t i_max = funcs.size();
		for (size_t i=0;i<i_max;i++){
			if (funcs[i] == func){
				waves.push_back(wave);
			};
		};
	};
	return waves;
};
//########################################################################################################################################################
///Simple getter for isobar names
std::string waveset::getIsobarName(
							int 							i)	const{

	return _amp_isos[i]->name();

};
//########################################################################################################################################################
///Gets the numbers of paramters for each isobar
std::vector<int> waveset::get_nParsIso()										const{

	std::vector<int> nPars;
	nPars.push_back(_iso_borders_par[0]);
	size_t i_max = _iso_borders_par.size();
	for (size_t i=1;i<i_max;i++){
		int diff = _iso_borders_par[i]-_iso_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
//########################################################################################################################################################
///Gets the numbers of constants for each isobar
std::vector<int> waveset::get_nConstIso()										const{

	std::vector<int> nConst;
	for (size_t iso=0;iso<_nIso;iso++){
		nConst.push_back(_amp_isos[iso]->nCon());
	};
	return nConst;
};
//########################################################################################################################################################
///Gets the numbers of paramters for the 'func' isobar paramterization
std::vector<int> waveset::get_isobar_pars(
							int 							func)	const{

	std::vector<int> pars;
	if(0==_iso_borders_par.size()){
		return pars;
	};
	int upper = _iso_borders_par[func];
	int lower = 0;
	if (func>0){
		lower = _iso_borders_par[func-1];
	};
	for (int i=lower;i<upper;i++){
		pars.push_back(i);
	};
	return pars;
};
//########################################################################################################################################################
///Gets the numbers of constants for the 'func' isobar paramterization
std::vector<int> waveset::get_isobar_const(
							int 							func)	const{

	std::vector<int> consts;
	int cCount=0;
	for(int iso=0;iso<func-1;iso++){
		cCount+=_amp_isos[iso]->nCon();
	};
	for(size_t con=0;con<_amp_isos[func]->nCon();con++){
		consts.push_back(cCount+con);
	};
	return consts;
};
//########################################################################################################################################################
///Get the waves that employ a certain isobar
std::vector<int> waveset::get_isobar_waves(
							int 							iso)	const{

	std::vector<int> waves;
	for (size_t i=0;i<_nWaves;i++){
		std::vector<int> isos = get_wave_isobars(i);
		for (size_t ff=0; ff<isos.size();ff++){
			if(iso==isos[ff]){
				waves.push_back(i);
				break;
			};
		};
	};
	return waves;
};
//########################################################################################################################################################
///Simple getter for parameter names
std::string waveset::getParameterName(
							int 							i)	const{ 	// # of paramter

	int pCount = 0;
	for (size_t f=0;f<_nFuncs;f++){
		pCount+=_amp_funcs[f]->nPar();
		if (pCount>i){
			pCount-=_amp_funcs[f]->nPar();			
			return _amp_funcs[f]->getParName(i-pCount);
		};
	};
	std::cerr<<"waveset::getParameterName(...): Error: Parameter #"<<i<<" does not exist"<<std::endl;
	return "not_used";

};
//########################################################################################################################################################
///Simple getter for constant names
std::string waveset::getConstantName(
							int 							i)	const{ 	// # of constant

	int nConst =0;
	for (size_t j =0;j<_nFuncs;j++){
		nConst+=_amp_funcs[j]->nCon();
		if(nConst > i){
			nConst-=_amp_funcs[j]->nCon();
			return _amp_funcs[j]->getConName(i-nConst);
		};
	};
	std::cerr<<"waveset::getConstantName(...): Error: Constant #"<<i<<" does not exist"<<std::endl;
	return "not_used";
};
//########################################################################################################################################################
///Simple getter for isobar parameter names
std::string waveset::getIsoParName(
							int 							i)	const{ 	// # of isobar paramter

	int pCount =0;
	for(size_t iso=0;iso<_nIso;iso++){
		pCount+=_amp_isos[iso]->nPar();
		if(pCount>i){
			pCount-=_amp_isos[iso]->nPar();
			return _amp_isos[iso]->getParName(i-pCount);
		};
	};
	std::cerr<<"waveset::getIsoParName(...): Error: Isobar Parameter #"<<i<<" does not exist"<<std::endl;
	return "not_used";
};
//########################################################################################################################################################
///Simple getter for isobar constant names
std::string waveset::getIsoConstName(
							int 							i)	const{ 	// # of isobar constant
	int cCount=0;
	for (size_t iso=0;iso<_nIso;iso++){
		cCount+=_amp_isos[iso]->nCon();
		if (cCount>i){
			cCount-=_amp_isos[iso]->nCon();
			return _amp_isos[iso]->getConName(i-cCount);
		};
	};
	std::cerr<<"waveset::getIsoConstName(...): Error: Isobar Constant #"<<i<<" does not exist"<<std::endl;
	return "not_used";
};
//########################################################################################################################################################
///Gets all parameter indices used for function i
std::vector<size_t> waveset::getFuncParameters(
							size_t							ftw)	const{

	std::vector<size_t> ret;
	size_t nCpl = _n_cpls.at(ftw);
	for (size_t tbin = 0;tbin<_nTbin;++tbin){
		ret.push_back(2*tbin*_nBrCpl+2*nCpl  );
		ret.push_back(2*tbin*_nBrCpl+2*nCpl+1);
	};
	std::vector<int> params = get_function_pars(_funcs_to_waves.at(ftw));
	for (size_t ii=0; ii<params.size();++ii){
		ret.push_back((size_t)(2*_nBrCpl*_nTbin+params.at(ii)));
	};
	int nBra = _n_branch.at(ftw);
	std::vector<int> lastpars = get_function_pars(_funcs_to_waves.at(_nFtw-1));
	size_t nPar = lastpars.at(lastpars.size()-1)+1;
	if (nBra > -1){
		ret.push_back((size_t)(2*_nBrCpl*_nTbin+nPar+2*nBra));
		ret.push_back((size_t)(2*_nBrCpl*_nTbin+nPar+2*nBra+1));
	};
	int iso = _iso_to_waves.at(ftw);
	if (iso > -1){
		std::vector<int> iso_params = get_isobar_pars(iso);
		for (size_t ii=0;ii<iso_params.size();++ii){
			ret.push_back(2*_nBrCpl*_nTbin+nPar+2*_nBranch+iso_params.at(ii));
		};
	};
	return ret;
};
//########################################################################################################################################################
///Gets the mass bin for a certain m3Pi
int waveset::get_bin(
							double 							mass)	const{

	for(size_t i=0;i<_nBins;i++){
		if (_binning[i] <= mass and _binning[i+1]>mass){
			return i;
		};
	};
	std::cerr<<"waveset.cxx::get_bin(): Error: Mass not in range"<<std::endl;
	return 0;
};
//########################################################################################################################################################
///Get the first function for each branching-set
std::vector<int> waveset::getFirstBranch()										const{

	std::vector<int> firstBranch;
	int maxBranch = -1;
	for (size_t i=0;i<_coupled.size();i++){
		if (_coupled[i] > maxBranch){
			maxBranch = _n_branch[i];
		};
	};
	for (int i=0;i<maxBranch+1;i++){
		for (size_t j=0;j<_coupled.size();j++){
			if (_coupled[j] == i){
				firstBranch.push_back(_n_branch[j]);
				break;
			};
		};
	};
	return firstBranch;
};
//########################################################################################################################################################
///Updates the number of function-to-wave couplings
void waveset::updateNftw(){

	_nFtw = _funcs_to_waves.size();
};
//########################################################################################################################################################
///Updates the number of employed points
void waveset::updateNpoints(){

	int nnn = 0;
	std::vector<int> point_to_wave;
	_point_borders_wave=std::vector<size_t>(0);
	for (size_t wave=0;wave<_nWaves;wave++){
		int plus = abs(_wave_binning_pts[wave]);
		nnn+=plus;
		_point_borders_wave.push_back(nnn);
		for (int i=0;i<plus;i++){
			point_to_wave.push_back(wave);
		};
	};
	_nPoints = nnn;
	_point_to_wave = point_to_wave;
};
//########################################################################################################################################################
///Updates the mass limits for the funtions from the mass limits for the waves
void waveset::updateFuncLims(){

	for (size_t func=0;func<_nFuncs;func++){
			double upper =  LOWER_MASS_LIMIT;// This is exchanged (upper <-> LOWER), since
			double lower =  UPPER_MASS_LIMIT;// the smallest mass region is wanted
			std::vector<int> waves = get_function_waves(func);
			size_t i_max = waves.size();
			for (size_t i=0;i<i_max;i++){
				int wave = waves[i];
				double waveUpper = _upperLims[wave];
				double waveLower = _lowerLims[wave];
				if (waveUpper > upper){
					upper = waveUpper;
				};
				if (waveLower < lower){
					lower = waveLower;
				};
			};
		_funcUpperLims[func] = upper;
		_funcLowerLims[func] = lower;
	};
};
//########################################################################################################################################################
//Updates the function spin according to their corresponding wave spin (if there is a mistake (different spins for the same function), nothing will happen)
void waveset::updateFuncSpin(){

	for (size_t func=0;func<_nFuncs;func++){
		std::vector<int> waves = get_function_waves(func);
		int i_max = waves.size();
		if (i_max >0){
			int L = _L[waves[0]];
			_amp_funcs[func]->setL(L);
		};
	};
};
//########################################################################################################################################################
///Updates the internal definitions for the isobar paramterizations
void waveset::updateIsobar(){

	for (size_t wave =0;wave<_nWaves;wave++){
		std::vector<int> isos = get_wave_isobars(wave);
		for (size_t niso=0;niso<isos.size();niso++){
			int iso = isos[niso];
			_amp_isos[iso]->setL(_L_iso[wave]);
			_iso_binning_pts[iso]=_wave_binning_pts[wave];
			_iso_n_binning[iso]  =_wave_n_binning[wave];
		};
	};
	updateNpoints();
};
//########################################################################################################################################################
///Updates the internal ranges for evaluation
void waveset::update_min_max_bin(){

	_minBin = 0;
//	std::cout<<"Call update_min_max_bin()"<<std::endl;
	_maxBin = _binning.size();
	if (_lowerLims.size() == 0 or _upperLims.size() == 0){
		std::cout<< "waveset::update_min_max_bin(): Warning: Wave limits not set, omitting internal setting of _minBin and _maxBin"<<std::endl;
		std::cout<< "waveset::update_min_max_bin(): Warning: Setting _minBin = 0; _maxBin = _binning.size() = "<<_binning.size()<<std::endl;
		_minBin =0;
		_maxBin =_binning.size();
		return;
	};
	double ancMin = _lowerLims[0];
	double ancMax = _upperLims[0];
	for (size_t i=0;i<_lowerLims.size();i++){
		if(_lowerLims[i] < ancMin){
			ancMin = _lowerLims[i];
		};
		if(_upperLims[i] > ancMax){
			ancMax = _upperLims[i];
		};
	};
//	std::cout<<ancMin<<"-"<<ancMax<<std::endl;
	if (_binning.size() == 0){
		std::cout<< "waveset::update_min_max_bin(): Warning: No binning set, cannot determine _minBin, _maxBin" <<std::endl;
	};
	for (size_t i=0;i<_binning.size()-1;i++){
		double up = _binning[i+1];
		double low= _binning[i];
		if (ancMin >= low and ancMin <up){
			_minBin = i;
		};
		if (ancMax > low and ancMax <=up){
			_maxBin = i+1;
		};
	};
//	std::cout<<"_minBin: "<<_minBin<<std::endl;
//	std::cout<<"_maxBin: "<<_maxBin<<std::endl;
};
//########################################################################################################################################################
///Handles the branchings after a function/wave coupling was added
void waveset::handle_branchings(
							int 							wave,
							int 							func){

	size_t border = _borders_waves[wave]-1;
	std::vector<int> coupled;
	std::vector<int> n_branch;
	size_t nRel = _coupled.size();
	for (size_t i=0; i<nRel;i++){
		coupled.push_back(_coupled[i]);
		n_branch.push_back(_n_branch[i]);
		if (i==border-1){
			n_branch.push_back(-1);
			coupled.push_back(-1);
		};
	};
	if (coupled.size() ==0){
		n_branch.push_back(-1);
		coupled.push_back(-1);
	};
	_coupled=coupled;
	_n_branch = n_branch;
	update_n_cpls();
	update_n_branch();
};
//########################################################################################################################################################
///Updates the number of couplings
void waveset::update_n_cpls(){

	std::vector<size_t> new_cpl;
	int max_cpl = 0;
	std::vector<int> brs;
	std::vector<int> their_cpls;
	for (size_t i=0;i<_coupled.size();i++){
		if (_coupled[i] ==-1){
			new_cpl.push_back(max_cpl);
			max_cpl++;
		}else{
			bool already = false;
			for (size_t j=0;j<brs.size();j++){
				if(brs[j] == _coupled[i]){
					new_cpl.push_back(their_cpls[j]);
					already = true;
					break;
				};
			};
			if (not already){
				new_cpl.push_back(max_cpl);
				brs.push_back(_coupled[i]);
				their_cpls.push_back(max_cpl);
				max_cpl++;
			};
		};
	};
	_nBrCpl = max_cpl;
	_n_cpls = new_cpl;
};
//########################################################################################################################################################
///Updates the branchings
void waveset::update_n_branch(){

	std::vector<int> n_branch;
	int act_branch = 0;
	for (size_t i=0;i<_coupled.size();i++){
		if (_coupled[i] == -1){
			n_branch.push_back(-1);
		}else{
			n_branch.push_back(act_branch);
			act_branch++;
		};
	};
	_n_branch = n_branch;
	_nBranch = act_branch;
};
//########################################################################################################################################################
///Sets the value for t' fot the corresponding t' bin [tbin] for each constant that is t'
std::vector<double> waveset::getVar(
							double							m,
							int 							tbin)		const{

	std::vector<double> vec(_nTvar+1);
	vec[0] = m;
//	vec[1] = (_t_binning[tbin][0] * 0.7 + _t_binning[tbin+1][0] * 0.3);
	for (size_t i=0;i<_nTvar;i++){
		vec[i+1] = _t_binning[tbin][i];
	};
	return vec;
};
//########################################################################################################################################################
///Does some internal consistency checks.
bool waveset::checkConsistency()											const{

	int nErr =0;
	if(_nWaves != _waveNames.size()){
		nErr+=1;
		std::cout<< "waveset::checkConsistency(): Inconsistency found: _nWaves does not match _waveNames.size()"<<std::endl;
	};
	if(_nFuncs != _amp_funcs.size()){
		nErr+=1;
		std::cout << "waveset::checkConsistency(): Inconsistency found: _nFuncs does not match _amp_funcs.size()"<<std::endl;
	};
	if( _nFuncs != _amp_funcs.size()){
		nErr+=1;
		std::cout << "waveset::checkConsistency(): Inconsistecy found: _nFuncs = "<<_nFuncs<<" does not match _amp_funcs.size() = "<<_amp_funcs.size()<<std::endl;
	};
	for (size_t i=1;i<_borders_waves.size();i++){
		if (_borders_waves[i] < _borders_waves[i-1]){
			nErr+=1;
			std::cout << "waveset::checkConsistency(): Inconsistency found: _borders_waves[] is not sorted"<<std::endl;
		};
	};
	for (size_t i=1;i<_borders_par.size();i++){
		if (_borders_par[i] < _borders_par[i-1]){
			nErr+=1;
			std::cout << "waveset::checkConsistency(): Inconsistency found: _borders_par[] is not sorted"<<std::endl;
		};
	};
	for (size_t i=0;i<_funcs_to_waves.size();i++){
		if (_funcs_to_waves[i] > (int)_nFuncs-1){
			nErr+=1;
			std::cout << "waveset::checkConsistency(): Inconsistency found: Function with number "<<_funcs_to_waves[i]<< " does not exist"<<std::endl;
		};
	};
	if (_wave_phase_space.size() != _nWaves){
		nErr+=1;
		std::cout << "waveset::checkConsistency(): Inconsistency found: _wave_phase_space.size() does not match _nWaves" << std::endl;
	};
	if(_nWaves != _L.size()){
		nErr+=1;
		std::cout<< "waveset::checkConsistency(): Inconsistency found: _nWaves does not match _L.size()"<<std::endl;
	};
	if(_nWaves != _upperLims.size()){
		nErr+=1;
		std::cout<< "waveset::checkConsistency(): Inconsistency found: _nWaves does not match _upperLims.size()"<<std::endl;
	};
	if(_nWaves != _lowerLims.size()){
		nErr+=1;
		std::cout<< "waveset::checkConsistency(): Inconsistency found: _nWaves does not match _lowerLims.size()"<<std::endl;
	};
	if (_iso_to_waves.size()>0){
		if(_iso_to_waves[0] != -1){
			nErr++;
			std::cout<<"waveset::checkConsistency(): Inconsistency found: First wave is de-isobarred"<<std::endl;
		};
	};
	for (size_t iso=0;iso<_nIso;iso++){
		std::vector<int> waves = get_isobar_waves(iso);
		if(waves.size()>0){
			int binn = _wave_n_binning[waves[0]];
			for (size_t wave=1;wave<waves.size();wave++){
				if (_wave_n_binning[waves[wave]] != binn){
					nErr++;
					std::cout<<"waveset::checkConsistency(): Inconsistency found: Different binnings for the same isobar"<<std::endl;
				};
			};
		};
	};
	if (nErr>0){
		return false;
	};
	return true;
	if (nErr>0){
		return false;
	};
	return true;
};
//########################################################################################################################################################
///Prints the internal status of the class instance
void waveset::printStatus()												const{

	std::cout <<"Internal status of waveset:"<<std::endl<<std::endl;
	std::cout<<"WAVES:"<<std::endl;
	std::cout << "_nWaves: " <<_nWaves<<std::endl;
	std::cout << std::endl;
	std::cout << "_waveNames" << std::endl;
	print_vector(_waveNames);
	std::cout<<std::endl<<"_nPoints: "<<_nPoints<<std::endl;
	std::cout << std::endl;
	std::cout<<"_nFtw (number of funtions to waves): "<<_nFtw<<std::endl;
	std::cout<<std::endl;
	std::cout << "_upperLims" << std::endl;
	print_vector(_upperLims);
	std::cout << "_lowerLims" << std::endl;
	print_vector(_lowerLims);
	std::cout << std::endl;
	std::cout << "_L" << std::endl;
	print_vector(_L);
	std::cout<<std::endl<<"_L_iso"<<std::endl;
	print_vector(_L_iso);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_wave_binning_pts"<<std::endl;
	print_vector(_wave_binning_pts);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_point_to_wave"<<std::endl;
	print_vector(_point_to_wave);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_wave_n_binning"<<std::endl;
	print_vector(_wave_n_binning);
	std::cout<<std::endl<<std::endl<<"FUNCTIONS: "<<std::endl;
	std::cout << "_funcs_to_waves" << std::endl;
	print_vector(_funcs_to_waves);
	std::cout<<std::endl<<"_n_cpls"<<std::endl;
	print_vector(_n_cpls);
	std::cout << std::endl;
	std::cout<<"_nPar: "<<_nPar<<std::endl;
	std::cout << std::endl;
	std::cout << "_nFuncs: " << _nFuncs<<std::endl;
	std::cout << std::endl;
	std::cout<<"FUNCTIONS:"<<std::endl;
	for (size_t i=0;i<_nFuncs;i++){
		_amp_funcs[i]->print();
		std::cout<<std::endl;
	};
	std::cout << std::endl;
	std::cout << "_borders_waves" << std::endl;
	print_vector(_borders_waves);
	std::cout<<std::endl;
	std::cout << "_funcUpperLims" << std::endl;
	print_vector(_funcUpperLims);
	std::cout << "_funcLowerLims" << std::endl;
	print_vector(_funcLowerLims);
	std::cout << std::endl;
	std::cout<<std::endl<<std::endl<<"ISOBARS: "<<std::endl;
	std::cout<<std::endl<<"_iso_to_waves"<<std::endl;
	print_vector(_iso_to_waves);
	std::cout << std::endl;
	std::cout<<std::endl<<"_nIso: "<<_nIso<<std::endl;
	std::cout << std::endl;
	std::cout<<std::endl<<"_isos"<<std::endl;
	std::cout<<"ISOBARS:"<<std::endl;
	for(size_t i=0;i<_nIso;i++){
		_amp_isos[i]->print();
	};
	std::cout << std::endl;
	std::cout<<std::endl<<"_iso_n_binning"<<std::endl;
	print_vector(_iso_n_binning);
	std::cout<<std::endl<<"_iso_binnings"<<std::endl;
	for (size_t i=0;i<_iso_binnings.size();i++){
		print_vector(_iso_binnings[i]);
	};
	std::cout << std::endl;
	std::cout<<std::endl<<"_point_borders_wave"<<std::endl;
	print_vector(_point_borders_wave);
	std::cout << std::endl;
	std::cout<<std::endl<<std::endl<<"PARAMETERS & CONSTANTS: "<<std::endl;
	std::cout << std::endl;
	std::cout << "_borders_par" << std::endl;
	print_vector(_borders_par);
	std::cout<<std::endl<<"_iso_borders_par"<<std::endl;
	print_vector(_iso_borders_par);
	std::cout<<std::endl<<std::endl<<"PHASE SPACE: "<<std::endl;
	std::cout<<"_globalPs: "<<_globalPs<<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout<<std::endl<<std::endl<<"BINNING: "<<std::endl;
	std::cout<<"_nBins: "<<_nBins<<std::endl<<std::endl<<"_binning"<<std::endl;
	print_vector(_binning);
	std::cout<<std::endl<<"_minBin: "<<_minBin<<"; _maxBin: "<<_maxBin<<std::endl;
	std::cout<<std::endl<<"_mMin: "<<_mMin<<"; _mMax: "<<_mMax<<std::endl;
	std::cout<<"_nTbin: "<<_nTbin<<std::endl<<std::endl<<"_t_binning"<<std::endl;
	for (size_t i=0;i<_t_binning.size();++i){
		print_vector(_t_binning[i]);
	};
	std::cout<<"_nTvar: "<<_nTvar<<std::endl<<std::endl;
	std::cout<<std::endl;
	print_vector(_eval_tbin);
	std::cout<<std::endl<<std::endl<<"BRANCHING: "<<std::endl;
	std::cout<<"_nBranch: "<<_nBranch<<"; _nBrCpl: "<<_nBrCpl<<std::endl<<"_coupled"<<std::endl;
	print_vector(_coupled);
	std::cout<<std::endl<<"_n_branch"<<std::endl;
	print_vector(_n_branch);
	std::cout<<std::endl<<std::endl<<"INTERNAL: "<<std::endl;
	std::cout<<"_write_out: "<<_write_out<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl<<"_has_isobars: "<<_has_isobars<<std::endl;
};
//########################################################################################################################################################
///Prints the paramters in a 'nice' way
void waveset::printParameters()												const{

	int func_count = 0;
	int iso_count =0;
	for (size_t i=0;i<_nWaves;i++){
		std::cout<<"────────────────────────────────────────────────"<<std::endl;
		std::cout<< _waveNames[i]<<std::endl;
		std::cout<<"┬───────────────────────────────────────────────"<<std::endl;
		std::vector<int> funcs = get_wave_functions(i);
		std::vector<int> isos  = get_wave_isobars(i);
		size_t j_max = funcs.size();
		size_t is_max= isos.size();
		for (size_t j=0;j<j_max;j++){
			int func = funcs[j];
			if (j==funcs.size()-1 and 0==is_max){
				std::cout << "└─";
			}else{
				std::cout << "├─";
			};
			std::cout << _amp_funcs[func]->name()<<" #"<< func_count  << std::endl;
			func_count++;
			std::vector<int> pars = get_function_pars(func);
			std::vector<int> cons = get_function_const(func);
			size_t k_max = pars.size();
			for(size_t k=0;k<k_max;k++){
				int par = pars[k];
				if (j==funcs.size()-1 and 0==is_max){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==pars.size()-1 and cons.size()==0 ){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< getParameterName(par)<<": "<<getPar(par)<<std::endl;
			};
			k_max = cons.size();
			for (size_t k=0;k<k_max;k++){
				int con = cons[k];
				if (j==funcs.size()-1 and 0==is_max){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==cons.size()-1 ){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< getConstantName(con)<<": "<<getCon(con)<<std::endl;
			};

		};
		if (is_max > 0){
			std::cout<<"│ISOBARS:"<<std::endl;
			for (size_t is=0;is<is_max;is++){
				int iso = isos[is];
				if ( isos.size()-1==is){
					std::cout<<"└─";
				}else{
					std::cout << "├─";
				};
				std::cout << _amp_isos[iso]->name()<<" #"<< iso_count  << std::endl;
				iso_count++;
				std::vector<int> pars = get_isobar_pars(iso);
				std::vector<int> cons = get_isobar_const(iso);
				size_t k_max = pars.size();
				for(size_t k=0;k<k_max;k++){
					int par = pars[k];
					if (is==isos.size()-1){
						std::cout<<"  ";
					}else{
						std::cout <<"│ ";
					};
					if (k==pars.size()-1 and cons.size()==0){
						std::cout <<"└─";
					}else{
						std::cout <<"├─";
					};
					std::cout<< getIsoParName(par)<<std::endl;
				};
				k_max = cons.size();
				for (size_t k=0;k<k_max;k++){
					int con = cons[k];
					if (is==isos.size()-1){
						std::cout<<"  ";
					}else{
						std::cout <<"│ ";
					};
					if (k==cons.size()-1){
						std::cout <<"└─";
					}else{
						std::cout <<"├─";
					};
					std::cout<< getIsoConstName(con)<<std::endl;
				};

			};
		};
	std::cout<<std::endl;
	};
};
//########################################################################################################################################################
///Initializes the text output. Opens a file with name
void waveset::open_output(
							std::string 							name){

	_outStream = new std::ofstream;
	_outStream->open(name.c_str());
	_write_out = true;
};
//########################################################################################################################################################
///Closes the text output
void waveset::close_output(){

	_outStream->close();
	_write_out = false;
};
//########################################################################################################################################################
