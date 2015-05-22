#ifndef NLOPT_CLASS
#define NLOPT_CLASS
#include<string>
#include"minimize.h"
#include <nlopt.hpp>
class nlopt_class:public minimize{
	public:
		nlopt_class();
		nlopt_class(std::string card);

		double          fit                             ();
		void            update_definitions_fitter       (){};
		void            reload_par_definitions_fitter   (int ulim, int olim){};
		bool            initialize                      ();

		bool            initialize_nlopt                (std::string card);

		double nlopt_function(const std::vector<double> &x, std::vector<double> &grad);
		const double *minimizerParameters               () const {return &_fitter_parameters[0];};

		void setMinimizerSpecifications (int spec_int = 0, double spec_double = 0., std::string spec_string = "");
	protected:
		nlopt::opt             _opt;                // nlopt optimizer
		nlopt::algorithm       _algorithm;
		int                    _opt_method;         // nlopt method to be used
		std::vector<double>    _fitter_parameters;  // Paramters, to hold the constant ones
};
#endif//NLOPT_CLASS
