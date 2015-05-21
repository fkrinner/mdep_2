#ifndef CHI2CONTAINTER
#define CHI2CONTAINTER
#include"minimize.h"
#include"minuit_root.h"
#include"nlopt_class.h"
#include<string>

class chi2container {
	public:
		chi2container();
		chi2container(std::string card);

		const minimize* minimizeChi2(){return _minimizeChi2;};

		std::string className(){return "chi2container";};
		void printStatus(){_minimizeChi2->printStatus();}; // Kepp the option to add somethig here

		void setMinimizerSpecifications(int spec_int = 0, double spec_double = 0., std::string spec_string = ""){
			_minimizeChi2->setMinimizerSpecifications(spec_int, spec_double, spec_string);
		};

	protected:
		minimize*	_minimizeChi2;

};
#endif//CHI2CONTAINTER
