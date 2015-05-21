#include"chi2container.h"


chi2container::chi2container(){
	_minimizeChi2 = new minimize();
};

chi2container::chi2container(std::string card){

	YAML::Node Ycard   = YAML::LoadFile(card);
	if (not Ycard["minimizer"]){
		std::cerr<<"chi2container::chi2container(std::string): Error: minimizer not defined in the card ('"<<card<<"')"<<std::endl;
	};
	std::string minimizerType = Ycard["minimizer"].as<std::string>();
	if (minimizerType == "minuit_root"){
		_minimizeChi2 = new minuit_root(card);
	}else if( minimizerType == "nlopt_class"){
		_minimizeChi2 = new nlopt_class(card);
	}else{
		std::cerr<<"chi2container::chi2container(std::string): Error: Invalid minimizer type: "<<minimizerType<<std::endl;
		throw;
	};
};
