#include<string>
#include<iostream>
#include<fstream>


#include"minimize.h"
#include"method.h"
#include"old_method.h"
#include"anchor_t.h"

#include"../chi_squared/breitWigners.h"
#include"yaml-cpp/yaml.h"
#include"string_utilities.h"
#include"matrix_utilities.h"

std::vector<double> load_file(std::string fileName){
	std::vector<double> ret	;
	std::fstream data(fileName.c_str(),std::ios_base::in);
	double val;
	while (data>>val){
		ret.push_back(val);
	};
	return ret;
};


int main(int argc, char* argv[]){

	std::string card = std::string(argv[1]);
	if (argc > 2){
		int seed = atoi(argv[2]);
		std::srand(seed); 
	};

	std::cout<<"Load definitions from:"<<std::endl;
	std::cout<<card<<std::endl;
	std::cout<<"--------start------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	minimize Chi2(card);
	std::cout<<"got out?"<<std::endl;
//	Chi2.finish_setUp();
	//int n = Chi2.getParNumber("G_a1(1260)");
	//Chi2.setParLimits(n,.13,13.);
	//Chi2.printStatus();
	//return 0;

	size_t nPar = Chi2.method()->Waveset()->getNpar();
	size_t nCpl = Chi2.method()->nCpl();
	size_t nBra = Chi2.method()->Waveset()->nBranch();
	size_t nTot = Chi2.method()->nTot();
	size_t nIso = Chi2.method()->Waveset()->getNiso();
	size_t nTbin= Chi2.method()->Waveset()->nTbin();
	size_t nFtw = Chi2.method()->Waveset()->nFtw();

	size_t dummy = nTbin*nFtw;

	dummy+=dummy;

	Chi2.printStatus();
	Chi2.method()->Waveset()->printParameters();

	std::cout<< "nPar: "<<nPar<<"; nCpl: "<<nCpl<<"; nBra: "<<nBra<<"; nIso: "<<nIso<<" => nTot: "<<nTot<<std::endl;

	std::vector<double> parrrr = load_file("./erte");
	std::cout<<"Found file with "<<parrrr.size()<<" paramters"<<std::endl;

	std::cout<<(*Chi2.method())()<<std::endl;
	std::cout<<Chi2(&Chi2.method()->parameters()[0])<<std::endl;
	std::cout<<Chi2.method()->className()<<std::endl;
//	return 0;


//	for (size_t i=0;i<2*nCpl/nTbin;++i){
//		Chi2.setParameter(i,parrrr[i]);
//	};


	std::cout<<"--------init-------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;

	Chi2.initCouplings(5);	
//	Chi2.relPar("M_pi1(1600)");
//	Chi2.relPar("G_pi1(1600)");
//	Chi2.fit();
//	Chi2.fit();

	std::cout<<"-------inited------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	std::cout<<"best couplings found:"<<std::endl;

	int mmm = Chi2.method()->getParNumber("M_pi1(1600)");
	int ggg = Chi2.method()->getParNumber("G_pi1(1600)");

	for (size_t i=0;i<parrrr.size();++i){
//		Chi2.setParameter(i,parrrr[i]);
	};
	std::cout<<"M: "<< Chi2.method()->getParameter(mmm)<<std::endl;
	std::cout<<"G: "<< Chi2.method()->getParameter(ggg)<<std::endl;

//	Chi2.fit();
	std::vector<double> bestpars = Chi2.method()->parameters();
	ofstream olo;
	olo.open("erte");
	for(size_t i=0;i<bestpars.size();++i){
		olo<<bestpars[i]<<" ";
	};
	olo.close();
//
//	std::cout<<"isolreit? "<<Chi2()<<std::endl;
	Chi2.method()->write_plots("plots0.txt",0);
	return 0;


/*
// From init couplings without branchs
	std::vector<std::complex<double> > actCpl(nFtw);
	for (size_t tbin=0;tbin<nTbin;++tbin){
		Chi2.method()->Waveset()->setEvalTbin(tbin,false);
	};
	std::vector<std::complex<double> > bestBra(nBra,std::complex<double>(0.,0.));
	for (size_t tbin=0;tbin<nTbin;++tbin){
		Chi2.method()->Waveset()->setEvalTbin(tbin,true);
		for (size_t cpl=0;cpl<nFtw;++cpl){
			size_t iii = 2*tbin*nFtw+2*cpl;
			actCpl[cpl] = std::complex<double>(cplll[iii],cplll[iii+1]);
		};
		std::vector<std::vector<std::complex<double> > > cplBra = Chi2.method()->full_to_br_cpl(actCpl);
		for (size_t cpl=0;cpl<nCpl/nTbin;++cpl){
			Chi2.setParameter(2*nCpl/nTbin*tbin+2*cpl  ,cplBra[0][cpl].real());
			Chi2.setParameter(2*nCpl/nTbin*tbin+2*cpl+1,cplBra[0][cpl].imag());
			Chi2.relPar(2*nCpl/nTbin*tbin+2*cpl  );
			Chi2.relPar(2*nCpl/nTbin*tbin+2*cpl+1);
		};
		for (size_t bra=0;bra<nBra;++bra){
			bestBra[bra]+=cplBra[1][bra]/std::complex<double>(nTbin,0.);
			Chi2.setParameter(2*nCpl+nPar+2*bra  ,cplBra[1][bra].real());
			Chi2.setParameter(2*nCpl+nPar+2*bra+1,cplBra[1][bra].imag());
		};

		std::cout<<tbin<<":::;:::"<<Chi2()<<std::endl;
		Chi2.method()->Waveset()->setEvalTbin(tbin,false);
	};
	for (size_t bra=0;bra<nBra;++bra){
		Chi2.setParameter(2*nCpl+nPar+2*bra  ,bestBra[bra].real());
		Chi2.setParameter(2*nCpl+nPar+2*bra+1,bestBra[bra].imag());
		Chi2.relPar(2*nCpl+nPar+2*bra  );
		Chi2.relPar(2*nCpl+nPar+2*bra+1);

	};


	for (size_t tbin=0;tbin<nTbin;++tbin){
		Chi2.method()->Waveset()->setEvalTbin(tbin,true);
		
	};
	std::cout<<"!!!!!!!!"<<Chi2()<<std::endl;
	
	Chi2.fit();
//Fritresult from brenched
	for (size_t i=0;i<650;++i){
		Chi2.setParameter(i,ppar[i]);
	};
	std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
//Unbranched fit result with straing values above



	for (size_t i=0;i<770;++i){
		Chi2.setParameter(i,cplolo[i]);
	};
	for (size_t tbin =0;tbin<nTbin;++tbin){
		stringstream nam;
		nam<<"whatever"<<tbin;
		Chi2.method()->write_plots(nam.str(),tbin);
	};
	for(size_t i=1;i<11;++i){
		Chi2.method()->Waveset()->setEvalTbin(i,false);
	};


	Chi2.method()->Waveset()->open_output("le_out_off_se");
	std::cout<<Chi2()<<std::endl;
	Chi2.method()->Waveset()->close_output();
//	std::cout<<Chi2.fit()<<std::endl;
//	print_vector(Chi2.method()->parameters());
*/
};

