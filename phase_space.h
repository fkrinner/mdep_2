#ifndef PHASOV_SPACOV
#define PHASOV_SPACOV
#include <fstream>
#include <sstream>
#include<string>
#include<complex>
#include<iostream>
#include<limits>
#include"amplitude_functions.h"



class phase_space_func{
	public:
		phase_space_func();
		phase_space_func(std::string waveName, std::string integralFile);

		bool get_binning(std::string integralFile);
		bool fill_integral(std::string waveName, std::string integralFile);
		std::string className()const{return "phase_space_func";};
		double operator()(const double* m) const;
	private:
		bool _isempty;
		std::string _wave_name;
		size_t _nStep;
		double _mmin;
		double _mmax;
		double _step;
		std::vector<double> _integrals;
};

phase_space_func::phase_space_func(){
	_isempty = true;
};


phase_space_func::phase_space_func(std::string waveName, std::string integralFile){
 
	_wave_name = waveName;
	_isempty = !get_binning(integralFile);
	fill_integral(waveName, integralFile);	
};

bool phase_space_func::get_binning(std::string integralFile){

	std::vector<double> masses;
	std::ifstream infile(integralFile.c_str());
	std::string line;
	while(std::getline(infile, line)){
		std::string buffer;
		std::stringstream strstr(line);
		bool isnext = false;
		while(strstr >> buffer){
			if (isnext){
				std::istringstream massstr(buffer);
				double mass;
				massstr >> mass;
				if (masses.size() == 0){
					masses.push_back(mass);
				};
				if (masses[masses.size()-1] != mass){
					masses.push_back(mass);
				};
			};
			if (buffer == "'"){
				isnext = true;
			}else{
				isnext = false;
			};
		};
	};
	if (masses.size() <2){
		std::cerr<< "phase_space_func::get_binning(...): Error: Less than two masses found in binning."<<std::endl;
		return false;
	};

	double binwidth0 = masses[1]-masses[0];
	for (size_t i=0;i<masses.size()-1;++i){
		if (masses[i+1] <= masses[i]){
			std::cerr<<"phase_space_func::get_binning(...): Error: Found non-ordered binning"<<std::endl;
		};
	};
	_mmin = masses[0]; // At the moment assume equidistant binning of the integrals
	_mmax = masses[masses.size()-1];
	_step = binwidth0;
	_nStep= masses.size();
	_integrals = std::vector<double>(_nStep);
	return true;
};

bool phase_space_func::fill_integral(std::string waveName, std::string integralFile){
	std::ifstream infile(integralFile.c_str());
	std::string line;
	size_t int_count =0;
	while(std::getline(infile, line)){
		if (line.find(waveName) != std::string::npos){
			std::string buffer;
			std::stringstream strstr(line);
			int countdown =-1;
			while(strstr >> buffer){
				if (countdown ==0){
					std::istringstream integralstr(buffer);
					integralstr >> _integrals[int_count];
					int_count++;
				};
				if (buffer == "'"){
					countdown = 2;
				};
				countdown--;
			};
		};
	};
	if (int_count != _nStep){
		std::cerr<<"phase_space::fill_integral(...): Error: Not all integrals found."<<std::endl;
	};
	return true;
};


double phase_space_func::operator()(const double* m)const{
	double m3Pi = m[0];
	if (_isempty){
		return 1.;
	};
	if (m3Pi < _mmin or m3Pi > _mmax){
		std::cerr << "phase_space_func::operator()(...): Error: Phase space not defined at m3Pi = " << m3Pi << std::endl;
		return 0.;
	};
	int pos = (m3Pi - _mmin)/_step + 1.; // The +1. is due to historical reasons at position (*)
	double frac = (m3Pi - _mmin-(pos-1)*_step)/_step;
	return pow(_integrals[pos-1]*(1-frac) + _integrals[pos]*frac,.5);
};




#endif//PHASOV_SPACOV
