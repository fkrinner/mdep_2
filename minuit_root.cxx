#include"minuit_root.h"

//########################################################################################################################################################
minuit_root::minuit_root():
	minimize(),
	_s1("Minuit2"),
	_s2("Migrad"){

	initialize();
};
//########################################################################################################################################################
minuit_root::minuit_root(std::string card):
	minimize(card),
	_s1("Minuit2"),
	_s2("Migrad"){

	initialize();
};
//########################################################################################################################################################
///Actual call for fitter. At the moment the instance is copied and fitted, this might be improved...
double minuit_root::fit(){ 

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
		std::cerr<<"minuit_root::fit(): Error: Fitter not initialized"<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//########################################################################################################################################################
///Update the parameter definitions in the fitter
void minuit_root::reload_par_definitions_fitter(
							int 						uLim,
							int						oLim){ 

	for(int i=uLim;i<oLim;i++){
		if(_released[i]){
			if((*_method->lower_parameter_limits())[i] < (*_method->upper_parameter_limits())[i]){
				_min->SetLimitedVariable(i,(*_method->parNames())[i],_method->parameters()[i],_step_sizes[i],(*_method->lower_parameter_limits())[i],(*_method->upper_parameter_limits())[i]);
			}else{
				_min->SetVariable(i,(*_method->parNames())[i],_method->parameters()[i],_step_sizes[i]);
			};
		}else{
			_min->SetFixedVariable(i,(*_method->parNames())[i],_method->parameters()[i]);
		};
	};
};
//########################################################################################################################################################
///Initializes the fitter
bool minuit_root::initialize(){ 

	if (_method->parameters().size()<_method->nTot()){
		std::cerr<<"minuit_root::initialize(...): Error: _method->parameters().size() < _method->nTot(). Abort initialize(initialization."<<std::endl;
		return false;
	};
	if ((*_method->parNames()).size()<_method->nTot()){
		std::cerr<<"minuit_root::initialize(...): Error: (*_method->parNames()).size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_step_sizes.size()<_method->nTot()){
		std::cerr<<"minuit_root::initialize(...): Error: _step_sizes.size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_released.size()<_method->nTot()){
		std::cerr<<"minuit_root::initialize(...): Error: _released.size() < _method->nTot(). Abort initialization."<<std::endl;
		return false;
	};
	_min = ROOT::Math::Factory::CreateMinimizer(_s1,_s2);
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
///Updates internal definitions
void minuit_root::update_definitions_fitter(){
	_min->SetTolerance(_tolerance);
	_min->SetMaxFunctionCalls(_maxFunctionCalls);
	_min->SetMaxIterations(_maxIterations);
};
//########################################################################################################################################################
///Sets minuit_root specific stuff
void	minuit_root::setMinimizerSpecifications(			int 						spec_int, 
									double 						spec_double, 
									std::string 					spec_string){
	
	std::cout<<"minuit_root::setMinimizerSpecifications(...): spec_double = "<<spec_double<<" will be ignored"<<std::endl;
	if (spec_int == 1){
		_s1 = spec_string;
		std::cout<<"minuit_root::setMinimizerSpecifications(...): Minuit defining string '_s1' set to: "<<_s1<<std::endl;
	}else if (spec_int ==2){
		_s2 = spec_string;
		std::cout<<"minuit_root::setMinimizerSpecifications(...): Minuit defining string '_s2' set to: "<<_s2<<std::endl;
	}else{
		std::cerr<<"minuit_root::setMinimizerSpecifications(...): Error: Configured specification not defined"<<std::endl;
	};
};
//########################################################################################################################################################
