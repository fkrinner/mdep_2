#include"nlopt_class.h"
#include<vector>

//double nlopt_function(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
//
//	nlopt_class* nlomin = new nlopt_class();//reinterpret_cast<nlopt_class*>(f_data); // Not beautiful, but workaround for nlopt
//	return instance->nlopt_function(x,grad);
//
//};
double ff_nlopt(const std::vector<double> &x, std::vector<double> &grad, void* f_data){

	nlopt_class* inst = reinterpret_cast<nlopt_class*>(f_data);
	double ret = inst->nlopt_function(x,grad);
	return ret;
};

//########################################################################################################################################################
nlopt_class::nlopt_class():
	minimize(),
	_algorithm(nlopt::LN_COBYLA){

	initialize();
};
//########################################################################################################################################################
nlopt_class::nlopt_class(std::string card):
	minimize(card),
	_algorithm(nlopt::LN_COBYLA){
	if (not initialize_nlopt(card)){
		std::cerr<<"nlopt_class::nlopt_class(std::string): Error: initialize_nlopt(std::string) failed"<<std::endl;
		throw;
	};
	initialize();
};
//########################################################################################################################################################
///Fuction to be called by the fitter
double nlopt_class::nlopt_function(const std::vector<double> &x, std::vector<double> &grad){

	size_t nTot = _method->nTot();
	size_t count = 0;
	for (size_t i=0;i<nTot;++i){
		if(_released[i]){
			_fitter_parameters[i] = x[count];
			++count;
		};
	};
	if (grad.size() != 0){
		std::vector<double> diff = _method->Diff(_fitter_parameters);
		count = 0;
		for (size_t i=0;i<nTot;++i){
			if(_released[i]){
				grad[count] = diff[i];
				++count;
			};
		};
	};
	return (*_method)(_fitter_parameters);
};
//########################################################################################################################################################
///Actual call for fitter. At the moment the instance is copied and fitted, this might be improved...
double nlopt_class::fit(){ 

	reload_par_definitions();
	print_vector(_released);
	_init = initialize();
	if(_init){
		double bestval;
		std::vector<double> parameters = _method->parameters();
		_opt.optimize(parameters, bestval);
		for (size_t i=0;i<_method->nTot();i++){
			_method->setParameter(i,parameters[i]);
		};
		return bestval;
	}else{
		std::cerr<<"nlopt_class::fit(): Error: Fitter not initialized"<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//########################################################################################################################################################
/// Initialize method
bool nlopt_class::initialize(){
	
	std::vector<double> par = _method->parameters();
	if (_fitter_parameters.size() != par.size()){
		_fitter_parameters = par;
	};
	_opt = nlopt::opt(_algorithm, _method->getNtot());
	_opt.set_ftol_abs(1.E-2);
	_opt.set_min_objective(ff_nlopt, this);
	return true;
};
//########################################################################################################################################################
///Method to set minimizerSpecifications
void nlopt_class::setMinimizerSpecifications(int spec_int, double spec_double, std::string spec_string){
	

/*Much more specification stuff can be written here:
 - More algorithms better implemented
 - Stopping criteria
 - ...
*/
	if(spec_int == 1){ // Specify algorithm
		if (spec_string == "NELDERMEAD"){
			_algorithm = nlopt::LN_NELDERMEAD;
		}else if (spec_string == "SBPLX"){
			_algorithm = nlopt::LN_SBPLX;
		}else if (spec_string == "LBFGS"){
			_algorithm = nlopt::LD_LBFGS;
		}else{
			std::cerr<<"nlopt_class::setMinimizerSpecifications(...): Error: Unknown algirithm: "<<spec_string<<std::endl;
			throw;
		};
		std::cout<<"nlopt_class::setMinimizerSpecifications(...): Set algorithm to "<<nlopt::algorithm_name(_algorithm)<<std::endl;
	};
};
//########################################################################################################################################################
///Load nlopt_definitions from card
bool nlopt_class::initialize_nlopt(std::string card){
	YAML::Node Ycard = YAML::LoadFile(card);
	if (Ycard["nlopt"]){

		return true;
	}else if(Ycard["nlopt_definition_file"]){
		return initialize_nlopt(Ycard["nlopt_definition_file"].as<std::string>());
		if (Ycard["nlopt"]["algorithm"]){
			setMinimizerSpecifications(1,0.,Ycard["nlopt"]["algorithm"].as<std::string>());
		};
	};
	std::cout<< "nlopt_class::initilize_nlopt(std::string): Warning: No definitions given in "<<card<<std::endl;
	return true;
};
