#ifndef OLDA___METHODA
#define OLDA___METHODA
#include"waveset.h"
#include"full_SDM.h"
#include<complex>
#include<vector>
#include<string>

#include"matrix_utilities.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class old_method : public full_SDM{
	public:
	// CONSTRUCTOR
		old_method();
#ifdef USE_YAML
		old_method(
								std::string 						card);
#endif//USE_YAML
	// EVALUATION METHODS
		double mainEval(
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
		adtl::adouble EvalTbin(
								int 							tbin,
								const std::complex<adtl::adouble> 			*cpl,
								const adtl::adouble	 				*par,
								const adtl::adouble					*iso_par)			const;
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
		std::vector<double>                     DiffAuto(
		                                                const std::vector<double>                               &xx)                            const;
#endif//ADOL_ON
		std::vector<double>                     Diff(
		                                                const std::vector<double>                               &xx)                            const;
                                                               
		std::vector<double>                     DiffTbin(
		                                                size_t                                                  tbin,
                                                                const std::vector<double>                               &xx)                            const;

		std::vector<double>                    DiffBin(
		                                                size_t                                                  bin,
		                                                size_t                                                  tbin, 
		                                                const std::vector<double>                               &xx, 
		                                                bool                                                    ignore_limits)                  const;
	// DATA HANDLING
		bool 					set_coma(int tbin, int bin, std::vector<double> coma);
		void 					loadComa(int tbin, const char* comaFile);
		void 					conjugate();
	// OTHER SETTERS & GETTERS
		const std::vector<double>*		get_coma(int tbin, int bin)const	{return &_coma[tbin][bin];};

		std::vector<double>			getParameters		()const;
	// OTHER METHODS
		std::string 				className		()const		{return "old_method";};
		void 					setTbinning(std::vector<std::vector<double> > binning);
	// PLOTTING
		using 					method::write_plots;
		void					write_plots(std::string filename, int tbin, const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;

	protected:
		// PARAMETERS AND DATA
		std::vector<std::vector<std::vector<double> > > 			_coma; 			// Just the errors in this case
};


#endif//OLDA___METHODA
