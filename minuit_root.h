#ifndef ROOT_MINUIT_MINIMIZE
#define ROOT_MINUIT_MINIMIZE
#include<string>
#include"minimize.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
class minuit_root:public minimize{
	public:
		minuit_root();
		minuit_root(std::string card);

		double 		fit				();
		void		update_definitions_fitter	();
		void		reload_par_definitions_fitter	(int ulim, int olim);
		bool		initialize			();
		const double *	minimizerParameters		()		const 		{return _min->X();};

		void		setMinimizerSpecifications	(int spec_int = 0, double spec_double = 0., std::string spec_string = "");

	protected:
		ROOT::Math::Minimizer* 	_min;							// ROOT Minimizer
		ROOT::Math::Functor 	_f;							// ROOT Functor object
		std::string		_s1;							// Minuit definition1 std::string s1="Minuit2"
		std::string 		_s2;							// Method definition  std::string s2="Migrad"



};
#endif//ROOT_MINUIT_MINIMIZE
