#ifndef MINIMIZE_MINI_MICE
#define MINIMIZE_MINI_MICE
#include<vector>
#include<complex>
#include<string>

#include"method.h"
#include"anchor_t.h"
#include"full_covariance.h"
#include"old_method.h"

typedef method METHOD; // Otherwise, the class method 'method()' and the type 'method' would have the same name --> typedef method METHOD

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML
double MIN_STEP_SIZE = 0.00001;

class minimize{
	public:
		minimize();
		~minimize(){delete _method;};
#ifdef USE_YAML
		minimize(std::string card);
#endif//USE_YAML
																			// In python
		double 			operator()			()				{return (*_method)(&_best_par[0]);};	
		double 			operator()			(std::vector<double>&xx)	{return (*_method)(xx);};
		double 			operator()			(const double*xx)		{return (*_method)(xx);};

		// Fitting routines
		void 			initCouplings			(size_t nSeeds = 1, int ntbin = -1);

	// Setters and getters
		METHOD*			method()							{return _method;};
		void 			setParameter			(int i, double par		);
		void 			setParameter			(std::string name, double par	);
		void			setParLimits			(int i, double upper, double lower);
		void			setParLimits			(std::string name, double upper, double lower);
		void 			setParameters			(std::vector<double> pars	);
		void 			setStepSize			(int i, double step		);
		void 			setStepSize			(std::string name, double par	);
		void 			setStepSizes			(std::vector<double> steps	);
		void 			setRandRange			(double range			);
		
		double			getParameter			(size_t i)	const;			

	// Fixing and releasing
		void 			relPar(int i);
		void 			fixPar(int i);
		void 			relPar(std::string name);
		void 			fixPar(std::string name);
		const std::vector<bool>*getReleased			()		const		{return &_released;};


	// Print routines
		std::string 		className			()		const		{return "minimize";};				// X
		void 			printStatus();													// X

	// Internal handlers
		void			update_definitions		();
		void			reload_par_definitions		(int mara_peter = -1);

		virtual double 		fit				()				{std::cout<<"minimize.h: fit() not overwritten"<<std::endl;				throw;return 0.;};
		virtual void		update_definitions_fitter	()				{std::cout<<"minimize.h: update_definitions_fitter() not overwritten"<<std::endl;	throw;};
		virtual void		reload_par_definitions_fitter	(int ulim, int olim)		{std::cout<<"minimize.h: reload_par_definitions_fitter() not overwritten"<<std::endl;	throw;};
		virtual bool		initialize			()				{std::cout<<"minimize.h: initialize() not overwritten"<<std::endl;			throw;return false;};
		virtual const double *	minimizerParameters		()		const 		{std::cout<<"minimize.h: minimizerParameters() not overwritten"<<std::endl;		throw;return 0;};

// Genuine method to pass integers, doubles or strings to the minimizer
		virtual void		setMinimizerSpecifications	(int spec_int = 0, double spec_double = 0., std::string spec_string = ""){std::cout<<"minimize.h: setMinimizerSpecifications() not overwritten"<<std::endl;throw;}; 


		void 			setRandomCpl			();										// X
		void 			setRandomBra			();										// X
		void			findRandRange			();
		void 			finish_setUp			();
		void			setMaxCalls			(size_t nCalls);

	// MULTINEST


#ifdef USE_YAML	
		void 			loadFitterDefinitions(YAML::Node &waveset);
		size_t			get_method(YAML::Node &card)			const;
#endif//USE_YAML
	protected:
	//METHOD
		METHOD*			_method;						// The method used (at the moment anchor_t)

		size_t			_method_type;
	// OWN STUFF
		std::vector<double> 	_best_par; 						// Best paramters
		double 			_randRange; 						// Range for random paramters (couplings and branchings)
		double 			_minStepSize;						// Minimal step size

	// MINIMIZER STUFF
		bool 			_init; 							// Flag for the initialization of the minimizer
		size_t 			_maxFunctionCalls;					// Miminizer definition
		size_t 			_maxIterations;						// Miminizer definition
		double 			_tolerance;						// Miminizer definition
		std::vector<double> 	_step_sizes;						// Step Size for each paramter	
		std::vector<bool> 	_released; 						// Status of each paramters
};
#endif//MINIMIZE_MINI_MICE
