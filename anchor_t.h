#ifndef ANCHOR_ROHCNA
#define ANCHOR_ROHCNA
#include"waveset.h"
#include"method.h"
#include<complex>
#include<vector>
#include<string>

#include"matrix_utilities.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class anchor_t : public method{
	public:
	// CONSTRUCTOR
		anchor_t();
#ifdef USE_YAML
		anchor_t(
								std::string 						card);
#endif//USE_YAML
	// EVALUATION METHODS
		double 			mainEval						(const double *xx);

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

		template<typename xdouble>
		xdouble EvalAutoCpl(
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalAutoCplBranch(
								const std::complex<xdouble>				*bra,
								const std::complex<xdouble>				*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalAutoCplTbin(
								int 							tbin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalTbin(
								int 							tbin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								const xdouble						*iso_par)			const;

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

	// AUTO CPL METHODS
		template<typename xdouble>
		AandB<xdouble> get_AB( // For de-isobarred
								int 							tbin,
								const std::complex<xdouble>	 			*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCpl(
								int 							tbin,
								const std::complex<xdouble> 				*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCplBra(
								int 							tbin,
								const std::complex<xdouble>	 			*branch,
								const std::complex<xdouble> 				*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;
	// DERIVATIVES
#ifdef ADOL_ON
		std::vector<double> 			Diff			(const double* xx)							const;
#endif//ADOL_ON
	// PARAMETER HANDLING

		std::vector<std::complex<double> > 	get_branchings		(const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<double> &iso_par) const;
		std::vector<std::complex<double> >	getAllCouplings		(int tbin,const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;
		void 					branchCouplingsToOne	();



	// DATA HANDLING
		bool 					set_data		(int tbin, int bin, std::vector<double> data);
		bool 					set_coma		(int tbin, int bin, std::vector<std::vector<double> > coma);
		void 					loadData		(int tbin, const char* dataFile);
		void 					loadComa		(int tbin, const char* comaFile);
		void 					nullify			();
		void 					conjugate		();
	// OTHER SETTERS & GETTERS
		bool					useBranch		()							const			{return _useBranch;};
		const std::vector<std::vector<double> >*get_coma		(int tbin, int bin)					const			{return &_coma[tbin][bin];};

		std::vector<double>			getParameters		()							const;

		size_t 					getNtotAnc		()							const;
		size_t 					getNanc			()							const;
		size_t 					getNtot			()							const			{return getNtotAnc();};
		size_t					getNcpl			()							const			{return getNanc();};
		
		bool 					setUseBranch		(bool in);
		const std::vector<double>		fullParameters		()							const;


	// OTHER METHODS
		std::string 				className		()							const			{return "anchor_t";};
		void 					printStatus		()							const;
		void 					set_is_ampl		(bool is_ampl);
		void 					setTbinning		(std::vector<std::vector<double> > binning);
		void 					update_n_cpls		();
		void					update_is_active	();

	// PLOTTING
		using					method::write_plots;
		void					write_plots		(std::string filename, int tbin,const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;
		std::vector<std::complex<double> > 	amplitudes		(double mass, int tbin, const std::vector<double> &param) const;


	protected:
		size_t 									_nBrCplAnc; // Number of couplings with branchings in the anchor wave
		// PARAMETERS AND DATA
		std::vector<std::vector<std::vector<std::vector<double> > > > 		_coma; 			// Covariance matrix
		std::vector<std::vector<std::vector<bool> > > _is_active; // Flag, which point is actually active
		// OTHER MEMBERS
		bool 									_is_ampl; 		// Enable Amplitunde fitting, needs to be tested.
		bool 									_useBranch; 		// Switches the usage of branchings on/off (only in the operator() method
};


#endif//ANCHOR_ROHCNA
