#include "minimize.h"

#include<vector>
#include<complex>
#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<stdexcept>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#ifdef ADOL_ON
#include "adolc/adolc.h"
#endif//ADOL_ON

#include"matrix_utilities.h"

minimize::minimize(): 
	_init(false),
	_randRange(100.), 
	_maxFunctionCalls(1000000),
	_maxIterations(100000),
	_tolerance(1.),
	_minStepSize(0.0001){
	_method = new METHOD();
};
#ifdef USE_YAML
//########################################################################################################################################################
///Constructor from YAML file
minimize::minimize(
							std::string 					card):
	_init(false),
	_randRange(100.), 
	_maxFunctionCalls(1000000),
	_maxIterations(100000),
	_tolerance(1.),
	_minStepSize(0.0001){
	YAML::Node Ycard   = YAML::LoadFile(card);
	size_t method_type = get_method(Ycard);
	_method_type = method_type;
	if (method_type == 0){
		_method = new anchor_t(card);	
	}else if(method_type == 1){
		_method = new full_covariance(card);	
	}else if(method_type==2){
		_method = new old_method(card);	
	};
//	int nmeth = 0;
//	if (nmeth==0){
//		_method = anchor_t(card);
//	}else if(nmeth==1){
//		_method = full_covariance(card);
//	}else if(nmeth==2){
//		_method = old_method(card);
//	};

//	_method = anchor_t(card);


	std::cout<<"Load fitter definitions"<<std::endl;
	loadFitterDefinitions(Ycard);
	std::cout<<"Fitter definitions loaded\nFinish setup"<<std::endl;
	finish_setUp();
	std::cout<<"Setup finished"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Actual call for fitter. At the moment the instance is copied and fitted, this might be improved...
double minimize::fit(){ 

	reload_par_definitions();
	print_vector(_released);
	if (_method_type == 0){
		_f=ROOT::Math::Functor(*((anchor_t*)_method),_method->nTot());
	}else if(_method_type == 1){
		_f=ROOT::Math::Functor(*((full_covariance*)_method),_method->nTot());
	}else if(_method_type ==2){
		_f=ROOT::Math::Functor(*((old_method*)_method),_method->nTot());
	};
	if(_init){
		_min->Minimize();
		const double *xs = _min->X();
		std::vector<double> best_par(_method->nTot());
		for (size_t i=0;i<_method->nTot();i++){
			_method->setParameter(i,xs[i]);
		};
		return (*_method)(xs);
	}else{
		std::cerr<<"Error: Fitter not initialized"<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//########################################################################################################################################################
///Finds start values for couplings and branchings // works only for anchot_t, for full_covariance not
void minimize::initCouplings(
							size_t 						nSeeds){ 


	std::cout<<"Initialize couplings"<<std::endl;

	size_t nTbin = _method->Waveset()->nTbin();
	size_t cpls = _method->nCpl()/nTbin;
	std::vector<double> bestChi2Tbin(nTbin);
	std::vector<std::vector<double> >bestCplTbin(nTbin,std::vector<double>(2*cpls));
	size_t nBra = _method->nBra();
	std::vector<std::vector<double> > bestBras = std::vector<std::vector<double> >(nTbin,std::vector<double >(2*nBra));

	if (not	_method->setUseBranch(false)){
		for (size_t i=2*_method->nCpl() + _method->nPar();i<2*_method->nCpl() + _method->nPar()+2*_method->nBra();++i){
			relPar(i);
		};
	};
	for (size_t seed=0;seed<nSeeds;++seed){
		std::cout<<"Seed #"<<seed<<std::endl;
		setRandomCpl(); // Set random couplings
		if (not _method->setUseBranch(false)){
			setRandomBra();
		};
		for(size_t tbin=0; tbin<_method->Waveset()->nTbin();tbin++){ // Switch off all t' bins
			_method->Waveset()->setEvalTbin(tbin,false);
		};

		for(size_t tbin=0; tbin<_method->Waveset()->nTbin();tbin++){ // Find cpl for each t' bin
			_method->Waveset()->setEvalTbin(tbin,true);
			std::cout<<"tBin #"<<tbin<<std::endl;
			for (size_t i =0;i<2*_method->nCpl();i++){
//				std::cout<<"fix "<<i<<std::endl;
				fixPar(i);
			};
			for (size_t i=0;i<2*cpls;i++){
//				std::cout<<"rel "<<i<<std::endl;
				relPar(2*cpls*tbin+i);
			};
			setRandomCpl(); // Set here again, otherwise different seed make no sense...
			double onetbinchi2 = fit();
			if (seed==0 or onetbinchi2 < bestChi2Tbin[tbin]){
				std::vector<double> actpar = _method->parameters();
				bestChi2Tbin[tbin] = onetbinchi2;
				for (size_t bestpar=0;bestpar<2*cpls;++bestpar){
					bestCplTbin[tbin][bestpar] = actpar[2*tbin*cpls+bestpar];
				};
				for(size_t bestbra=0;bestbra<2*nBra;++bestbra){
					bestBras[tbin][bestbra] = actpar[2*_method->nCpl()+_method->nPar()+bestbra];
				};
			};
			std::cout <<"... Chi2 = "<<onetbinchi2<<std::endl;
			_method->Waveset()->setEvalTbin(tbin,false);
		};
	};
	for (size_t tbin=0;tbin<nTbin;++tbin){
		std::cout<<"best Chi2 for tbin #"<<tbin<<": "<<bestChi2Tbin[tbin]<<std::endl;

		for (size_t cpl=0;cpl<2*cpls;++cpl){
			setParameter(2*cpls*tbin+cpl,bestCplTbin[tbin][cpl]);
		};
	};

	for(size_t tbin=0;tbin<_method->Waveset()->nTbin();tbin++){ // Switch on all t' bins
		_method->Waveset()->setEvalTbin(tbin,true);
	};
	std::vector<std::complex<double> > couplings(_method->nCpl());
	std::vector<double> par(_method->nPar());
	for (size_t i=0;i<_method->nCpl();i++){
		couplings[i] = std::complex<double>(_method->parameters()[2*i],_method->parameters()[2*i+1]);
	};
	for (size_t tbin=0;tbin<nTbin;++tbin){
		std::complex<double> firstcpl = couplings[cpls*tbin];
		std::complex<double> phase = conj(firstcpl)/abs(firstcpl);
		for (size_t cpl =0; cpl<cpls; ++cpl){
			size_t nn = tbin*cpls+cpl;
			couplings[nn]*=phase;
			setParameter(2*nn  ,couplings[nn].real());
			setParameter(2*nn+1,couplings[nn].imag());
		};
		for (size_t bra=0;bra<nBra;++bra){
			std::complex<double> cc(bestBras[tbin][2*bra],bestBras[tbin][2*bra+1]);
			cc*=phase;
			bestBras[tbin][2*bra  ] = cc.real();
			bestBras[tbin][2*bra+1] = cc.imag();
		};


	};


	for (size_t i=0;i<_method->nPar();i++){
		par[i] = _method->parameters()[2*_method->nCpl()+i];
	};
	std::vector<double> iso_par(_method->nIso());
	for (size_t i=0;i<_method->nIso();i++){
		iso_par[i] = _method->parameters()[2*_method->nCpl()+_method->nPar()+2*_method->nBra()+i];
	};
//	std::cout << "Total with _method->EvalAutoCpl() (For consistency check): "<< _method->EvalAutoCpl(&couplings[0],&par[0],&iso_par[0])<<std::endl; Removed do to convertability
	_method->setUseBranch(true);
	if (_method->nBra()>0){
		std::vector<double> branchings(2*nBra,0.);
		for (size_t bra = 0;bra<2*nBra;++bra){
			for (size_t tbin=0;tbin<nTbin;++tbin){
				branchings[bra]+=bestBras[tbin][bra];
			};
			branchings[bra]/=nTbin;
			setParameter(2*_method->nCpl()+_method->nPar()+bra,branchings[bra]);
		};

		// branchCouplingsToOne(); // Set all coupled couplings to one, since all should be in the branchings right now // Somehow Chi2 is better, when this is not done
//std::cout << "With the found branchings, Chi2(...)="<< _method->EvalAutoCplBranch(&bra[0],&couplings[0],&par[0],&iso_par[0])<<" ('_method->EvalAutoCplBranch(...)')"<<std::endl; //[0]//
		for (size_t i =0;i<2*_method->nCpl();i++){ // Fix couplings
			fixPar(i);
		};
		for (size_t i=0;i<2*_method->nBra();i++){ // Rel Branchings
			relPar(2*_method->nCpl()+_method->nPar()+i);
		};
		fit();
		std::cout<<"Couplings and branchings"<<std::endl;
		for (size_t i =0;i<2*_method->nCpl();i++){ // Rel couplings
			relPar(i);
		};
		fit();
	}else{
		for(size_t i=0;i<_method->nCpl();i++){
			relPar(2*i);
			relPar(2*i+1);
		};
	};
	std::cout<<"Total: "<<(*_method)(_min->X())<<std::endl;
	std::cout<<"Couplings and branchings found"<<std::endl;
	std::cout<<"Setting automatic limits for couplings and branchings"<<std::endl;
	for (size_t i=0;i<_method->nCpl();i++){
		double val = std::max(_method->parameters()[2*i]*_method->parameters()[2*i],_method->parameters()[2*i+1]*_method->parameters()[2*i+1]);
		val = pow(val,.5);
		_method->setParLimits(2*i  ,3*val,-3*val);
		_method->setParLimits(2*i+1,3*val,-3*val);
	};
	int par_bef = 2*_method->nCpl() +_method->nPar();
	for (size_t i=0;i<_method->nBra();i++){
		double val = std::max(_method->parameters()[par_bef+2*i]*_method->parameters()[par_bef+2*i],_method->parameters()[par_bef+2*i+1]*_method->parameters()[par_bef+2*i+1]);
		val = pow(val,.5);
		_method->setParLimits(par_bef+2*i  ,3*val,-3*val);
		_method->setParLimits(par_bef+2*i+1,3*val,-3*val);
	};
};
//########################################################################################################################################################
///Call lower 'setParamter' and sets internal definitions
void minimize::setParameter(
							int 						i, 	// # of parameter
							double 						par){

	_method->setParameter(i,par);
	reload_par_definitions(i);
};
//########################################################################################################################################################
///Sets a parameter by name
void minimize::setParameter(
							std::string 					name, 
							double 						par){

	int number = _method->getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_method->nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setParameter(number,par);
	};
};
//########################################################################################################################################################
///Set Parameter limits by number
void minimize::setParLimits(				int 						i,
							double 						upper,
							double 						lower					){

	_method->setParLimits(i,upper,lower);
	reload_par_definitions(i);
};
//########################################################################################################################################################
/// Set Parameter limits by name
void minimize::setParLimits(				std::string 					name,
							double						upper,
							double						lower){

	int number = _method->getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_method->nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setParLimits(number,upper,lower);
	};
};
//########################################################################################################################################################
///Sets all parameters at once
void minimize::setParameters(
							std::vector<double> 				pars){

	_method->setParameters(pars);
	reload_par_definitions();
};
//########################################################################################################################################################
///Sets parameter step size
void minimize::setStepSize(
							int 						i, 	// # of parameter
							double 						step){

	if (_step_sizes.size() < _method->nTot()){
		_step_sizes = std::vector<double>(_method->nTot(),_minStepSize);
		reload_par_definitions();
	};
	_step_sizes[i] = step;
	reload_par_definitions(i);
};

//########################################################################################################################################################
///Sets the step size by name
void minimize::setStepSize(
							std::string 					name, 
							double 						par){

	int number = _method->getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_method->nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setStepSize(number,par);
	};
};
//########################################################################################################################################################
///Sets all step sizes
void minimize::setStepSizes(
							std::vector<double> 				steps){

	if (steps.size() >= _method->nTot()){
		_step_sizes=steps;
		reload_par_definitions();
	}else{
		std::cerr<<"Error: Input steps.size() too small"<<std::endl;
	};
};
//########################################################################################################################################################
///Sets the range for generating random couplings and branchings
void minimize::setRandRange(
							double 						range){

	_randRange=range;
};
//########################################################################################################################################################
///Getter for parameters (for simplicty to avoid minimize::method()->parameters()[]
double minimize::getParameter(				size_t						i)						const{

	size_t max = _method->nTot();
	if (i>max){
		throw std::invalid_argument("Parameter i does not exist");
	};
	return _method->parameters()[i];
};
//########################################################################################################################################################
///Release parameter by number
void minimize::relPar(
							int 						i){	// # of parameter

	if (_released.size() == _method->nTot()){
		_released[i] = true;
		reload_par_definitions(i);
	}else{
		std::cout<<_method->nTot()<<" "<<_method->parameters().size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};
//########################################################################################################################################################
///Fix parameter by number
void minimize::fixPar(
							int 						i){	// # of parameter

	if (_released.size() == _method->nTot()){
		_released[i] = false;
		reload_par_definitions(i);
	}else{
		std::cout<<_method->nTot()<<" "<<_method->parameters().size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};
//########################################################################################################################################################
///Release parameter by name
void minimize::relPar(
							std::string 					name){

	int number = _method->getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_method->nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		relPar(number);
	};
};
//########################################################################################################################################################
///Fix parameter by name
void minimize::fixPar(
							std::string 					name){

	int number = _method->getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_method->nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		fixPar(number);
	};
};
//########################################################################################################################################################
///Prints the internal status
void minimize::printStatus(){ 

	_method->printStatus();
	std::cout<<std::endl<<"_best_par"<<std::endl;
	print_vector(_best_par);
	std::cout<<std::endl<<"_randRange: "<<_randRange<<std::endl;
	std::cout<<std::endl<<"_minStepSize: "<<_minStepSize<<std::endl;
	std::cout<<std::endl<<"_init: "<<_init<<std::endl;
	std::cout<<std::endl<<"_maxFunctionCalls: "<<_maxFunctionCalls<<std::endl;
	std::cout<<std::endl<<"_maxIterations: "<<_maxIterations<<std::endl;
	std::cout<<std::endl<<"_tolerance: "<<_tolerance<<std::endl;
	std::cout<<std::endl<<"_step_sizes"<<std::endl;
	print_vector(_step_sizes);
	std::cout<<std::endl<<"_released:"<<std::endl;
	print_vector(_released);
};
//########################################################################################################################################################
///Updates internal definitions
void minimize::update_definitions(){ 

	_method->update_definitions();
	std::vector<bool> rels;
	for (size_t i=0;i<_method->nCpl();i++){
		rels.push_back(true);
		rels.push_back(true);
	};
	for (size_t i=0;i<_method->nPar();i++){
		rels.push_back(false);
	};
	for (size_t i=0;i<_method->nBra();i++){
		rels.push_back(false);
		rels.push_back(false);
	};
	for (size_t i=0;i<_method->nIso();i++){
		rels.push_back(false);
	};
	_released = rels;
	if(_init){
		_min->SetTolerance(_tolerance);
		_min->SetMaxFunctionCalls(_maxFunctionCalls);
		_min->SetMaxIterations(_maxIterations);
	};
};
//########################################################################################################################################################
///Update the parameter definitions for parameter number mara_peter, if mara_peter ==-1, then all are updated
void minimize::reload_par_definitions(
							int 						mara_peter){ 

	int uLim = 0;
	int oLim = _method->nTot();
	if (mara_peter > -1){
		uLim = mara_peter;
		oLim = mara_peter+1;
	};
	if((*_method->lower_parameter_limits()).size() != _method->parameters().size()){
		_method->init_lower_limits(_method->parameters().size());
	};
	if((*_method->upper_parameter_limits()).size() != _method->parameters().size()){
		_method->init_upper_limits(_method->parameters().size());
	};
	if(_init){
		for(int i=uLim;i<oLim;i++){
			if(_released[i]){
				if((*_method->lower_parameter_limits())[i] < (*_method->upper_parameter_limits())[i]){
					_min->SetLimitedVariable(i,(*_method->parNames())[i],_method->parameters()[i],_step_sizes[i],(*_method->lower_parameter_limits())[i],(*_method->upper_parameter_limits())[i]);
//					std::cout<<(*_method->parNames())[i]<<" limited to "<<(*_method->lower_parameter_limits())[i]<<"-"<<(*_method->upper_parameter_limits())[i]<<std::endl;
				}else{
					_min->SetVariable(i,(*_method->parNames())[i],_method->parameters()[i],_step_sizes[i]);
				};
			}else{
				_min->SetFixedVariable(i,(*_method->parNames())[i],_method->parameters()[i]);
			};
		};
	};
};
//########################################################################################################################################################
///Initializes the fitter
bool minimize::initialize(std::string s1, std::string s2){ 

	if (_method->parameters().size()<_method->nTot()){
		std::cout<<"_method->parameters().size() < _method->nTot(). Abort initialize(initialization."<<std::endl;
		return false;
	};
	if ((*_method->parNames()).size()<_method->nTot()){
		std::cout<<"(*_method->parNames()).size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_step_sizes.size()<_method->nTot()){
		std::cout<<"_step_sizes.size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_released.size()<_method->nTot()){
		std::cout<<"_released.size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	_min = ROOT::Math::Factory::CreateMinimizer(s1,s2);
	_min->SetMaxFunctionCalls(_maxFunctionCalls);
	_min->SetMaxIterations(_maxIterations);
	_min->SetTolerance(_tolerance);
	if (_method_type == 0){
		_f=ROOT::Math::Functor(*((anchor_t*)_method),_method->nTot());
	}else if(_method_type == 1){
		_f=ROOT::Math::Functor(*((full_covariance*)_method),_method->nTot());
	}else if(_method_type ==2){
		_f=ROOT::Math::Functor(*((old_method*)_method),_method->nTot());
	};

	_min->SetFunction(_f);
	_init = true;
	update_definitions();
	reload_par_definitions();
	return true;
};
//########################################################################################################################################################
///Randomize couplings within _randRange
void minimize::setRandomCpl(){

	for(size_t i=0; i<2*_method->nCpl();i++){
		if (_released[i]){
			_method->setParameter(i,2*((double)rand()/RAND_MAX-0.5)*_randRange);
		};
	};
	reload_par_definitions();
};
//########################################################################################################################################################
///Randomize branchings within _randRange
void minimize::setRandomBra(){ 

	for (size_t i=0;i<2*_method->nBra();i++){
		if (_released[2*_method->nCpl()+_method->nPar()+i]){
			_method->setParameter(2*_method->nCpl()+_method->nPar()+i, 2*((double)rand()/RAND_MAX-0.5)*_randRange);
		};
	};
	reload_par_definitions();
};
//########################################################################################################################################################
/// Find random range for the couplings
void minimize::findRandRange(){

	double basis = 2.;
	double minChi2 = std::numeric_limits<double>::max();
	size_t minDim = 0;
	for (size_t i=0;i<20;++i){
		setRandRange(pow(basis,i));
		for (size_t j=0;j<10;++j){
			setRandomCpl();
			double actChi2 = (*_method)(&(_method->parameters())[0]);
//			std::cout<<"findRandRange: "<<basis<<"^"<<i<<": "<<actChi2<<std::endl;
			if (actChi2<minChi2){
				minChi2 = actChi2;
				minDim = i;
			};
		};
	};
	std::cout<<"Found best _randRange to be: "<<basis<<"^"<<minDim<<" with Chi2 approx.: "<<minChi2<<std::endl; 
	setRandRange(pow(basis,minDim));
	setRandomCpl();
};
//########################################################################################################################################################
///Sets some internal vaiables accordingly
void minimize::finish_setUp(){

	std::vector<int> first_branch = _method->Waveset()->getFirstBranch();
	if (_released.size() != _method->parameters().size()){
	//		std::cout<<"Warning: No paramter status set, releasing couplings, fixing all others."<<std::endl;
		std::vector<bool> std_rel(_method->parameters().size(),false);
		for (size_t i=0;i<2*_method->nCpl();i++){
			std_rel[i]=true;
		};
		_released = std_rel;
	};
	for (size_t i=0; i<first_branch.size();i++){
		setParameter(2*_method->nCpl()+_method->nPar()+2*first_branch[i]  ,1.); // Re(Br)
		setParameter(2*_method->nCpl()+_method->nPar()+2*first_branch[i]+1,0.); // Im(Br)
		fixPar(2*_method->nCpl()+_method->nPar()+2*first_branch[i]  );
		fixPar(2*_method->nCpl()+_method->nPar()+2*first_branch[i]+1);
	};
	setRandomCpl();
	setRandomBra();
	initialize();
};
#ifdef USE_YAML
//########################################################################################################################################################
///Loads some fitter definitions from a YAML file
void minimize::loadFitterDefinitions(
							YAML::Node 					&waveset){

	if (waveset["min_step_size"]){
		_minStepSize = waveset["min_step_size"].as<double>();
	};
	int nPar = _method->Waveset()->getNpar();
	int nCpl = _method->nCpl();
	for (int par = 0;par<nPar;par++){
		setStepSize(2*nCpl+par,std::max(_minStepSize, 0.0001*fabs(_method->parameters()[2*nCpl+par])));
	};
	if (waveset["real_anchor_cpl"]){ /// Also include this here, so no extra method is necessary
		if(waveset["real_anchor_cpl"].as<bool>()){
			setParameter(1,0.);
			fixPar(1);
		};
	};
};
//########################################################################################################################################################
size_t minimize::get_method(
							YAML::Node					&card)						const{

	if (not card["method"]){
		throw std::runtime_error("Method not defined in the 'card'");
	};
	if (card["method"].as<std::string>() == "anchor_t"){
		return 0;
	};
	if (card["method"].as<std::string>() == "full_covariance"){
		return 1;
	};
	if (card["method"].as<std::string>() == "old_method"){
		return 2;
	};
	throw runtime_error("Invalid method ginven in the 'card'");
	return 1111111;
};
#endif//USE_YAML
//########################################################################################################################################################
///Cube required by the MultiNest package
void minimize::cube(					double						*in)						const{
	
	for (size_t i=0;i<_method->nTot();i++){
		in[i] = (1-in[i])*(*_method->lower_parameter_limits())[i]+ in[i]*(*_method->upper_parameter_limits())[i];
	};
};
//########################################################################################################################################################
