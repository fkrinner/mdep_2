#ifndef ___FULLL___COMA___
#define ___FULLL___COMA___
#include"waveset.h"
#include"full_SDM.h"
#include<complex>
#include<vector>
#include<string>

#include"matrix_utilities.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class full_covariance : public full_SDM{
	public:
	// CONSTRUCTOR
		full_covariance();
#ifdef USE_YAML
		full_covariance(
								std::string 						card);
#endif//USE_YAML
	// EVALUATION METHODS
		double operator()();

		double operator()(
								std::vector<double> 					&xx);

		double operator()(
								const double						*xx);

		template<typename xdouble>
		xdouble EvalCP(
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalBranch(
								const std::complex<xdouble>				*branch,
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		double EvalTbin(
								int 							tbin,
								const std::complex<double> 				*cpl,
								const double	 					*par,
								const double						*iso_par)			const;
#ifdef ADOL_ON
		adouble EvalTbin(
								int 							tbin,
								const std::complex<adouble> 				*cpl,
								const adouble	 					*par,
								const adouble						*iso_par)			const;
#endif//ADOL_ON


		template<typename xdouble>
		xdouble EvalBin(
								int 							tbin,
								int 							bin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval)			const;

		template<typename xdouble>
		std::vector<xdouble> delta(
								int 							tbin,
								int 							bin,
								double 							mass,
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval)			const;

	// DERIVATIVES
#ifdef ADOL_ON
		std::vector<double> 			Diff(std::vector<double> &xx)								const;
		std::vector<double> 			Diff(const double* xx)									const;
#endif//ADOL_ON
	// DATA HANDLING
		bool 					set_data(int tbin, int bin, std::vector<double> data);
		bool 					set_coma(int tbin, int bin, std::vector<std::vector<double> > coma);
		void 					loadComa(int tbin, const char* comaFile);
		void 					conjugate();
	// OTHER SETTERS & GETTERS

		const std::vector<std::vector<double> >*get_coma(int tbin, int bin)const	{return &_coma[tbin][bin];};

		std::vector<double>			getParameters		()const;
	// OTHER METHODS
		std::string 				className		()const		{return "full_covariance";};
		void 					setTbinning(std::vector<std::vector<double> > binning);
	// PLOTTING
		using					method::write_plots;
		void					write_plots(std::string filename, int tbin,const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;

	protected:
		// PARAMETERS AND DATA
		std::vector<std::vector<std::vector<std::vector<double> > > > 		_coma; 			// Covariance matrix
};


#endif//___FULLL___COMA___
