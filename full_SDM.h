#ifndef FULL_SDDMMM
#define FULL_SDDMMM
#include"method.h"
#include<string>

class full_SDM : public method{
	public:
		full_SDM():method(){};
		full_SDM(std::string card) : method(card){};

		std::vector<std::complex<double> > 	get_branchings(const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<double> &iso_par) const;
		std::vector<std::complex<double> >	getAllCouplings(int tbin,const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;

		virtual std::string			className() const {return "undefined_full_SDM_method";};
		void 					branchCouplingsToOne();
		bool 					set_data(int tbin, int bin, std::vector<double> data);
		void 					loadData(int tbin, const char* dataFile);

		std::vector<std::vector<std::complex<double> > > full_to_br_cpl(std::vector<std::complex<double> > &cpl);
		void 					update_n_cpls();
		void 					nullify();
		size_t 					getNtot()		const;
		size_t 					getNcpl()		const;
		void 					printStatus()		const;
		void					update_is_active();
	protected:
		std::vector<std::vector<bool> >						_is_active;		// Flag, which point is actually active
};
#endif//FULL_SDDMMM
