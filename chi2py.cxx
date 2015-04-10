#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>
#include"minimize.h"
#include<string>
namespace bp = boost::python;

template<class T>
bp::list std_vector_to_py_list(const std::vector<T>& v){
	bp::object get_iter = bp::iterator<std::vector<T> >();
	bp::object iter = get_iter(v);
	bp::list l(iter);
	return l;
};


template< typename T >
std::vector<T> to_std_vector(const bp::object& iterable){
    return std::vector<T>(bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>( ));
};

template< typename T >
std::vector<std::vector<T> > to_std_vector_vector(const bp::object& iterable){
	std::vector<bp::object> intermediate = to_std_vector<bp::object>(iterable);
	std::vector<T> line;
	std::vector<std::vector<T> > final;
	for (unsigned int i=0;i<intermediate.size();i++){
		line=to_std_vector<T>(intermediate[i]);
		final.push_back(line);
	};
	return final;	
};


struct chi2py:public minimize{

	chi2py(std::string card):
		minimize(card){};

	size_t nTot(){return _method->nTot();};
	size_t nTbin(){return _method->Waveset()->nTbin();};
	size_t nCpl(){return _method->nCpl();};
	size_t nPar(){return _method->Waveset()->getNpar();};
	size_t nBra(){return _method->Waveset()->nBranch();};
	size_t nFtw(){return _method->Waveset()->nFtw();};
	size_t nWaves(){return _method->Waveset()->nWaves();};

	void relPar(std::string inin){
		minimize::relPar(inin);
	};
	void fixPar(std::string inin){
		minimize::fixPar(inin);
	};
	void printParameters(){
		_method->Waveset()->printParameters();
	};
	std::string method(){
		return _method->className();
	};

	double Eval(bp::object par){
		std::vector<double> xxx = to_std_vector<double>(par);
		return (*_method)(xxx);
	};
	double EvalSelf(){
		return (*_method)(&_method->parameters()[0]);
	};

	double getParameter(std::string name){
		int i = _method->getParNumber(name);
		if (i<0){
			std::cout<<"Error: Parameter '"<<name<<"' not definded"<<std::endl;
	 		return std::numeric_limits<double>::quiet_NaN();
		};
		return minimize::getParameter(i);
	};
	void setParameter(std::string name, double val){
		minimize::setParameter(name,val);
	};
	void write_plots(std::string file_name, size_t tbin){
		_method->write_plots(file_name,tbin);
	};

	bp::list parameterNames(){
		std::vector<std::string> names = *_method->parNames();
		return std_vector_to_py_list(names);
	};
	bp::list parameters(){
		std::vector<double> xxx = _method->parameters();
		return std_vector_to_py_list(xxx);
	};
	bp::list Amplitudes(double mass, int tbin, bp::object param_in){
		std::vector<double> par = to_std_vector<double>(param_in);
		std::vector<std::complex<double> > amps = _method->amplitudes(mass, tbin, par);
		return std_vector_to_py_list(amps);
	};
	bp::list AmplitudesSelf(double mass, int tbin){
		std::vector<std::complex<double> > ampl = _method->amplitudes(mass, tbin, _method->parameters());
		return std_vector_to_py_list(ampl);
	};
	bp::list waveNames(){
		std::vector<std::string> names = *_method->Waveset()->waveNames();
		return std_vector_to_py_list(names);
	};

	void setParLimits(std::string name, double upper, double lower){
		minimize::setParLimits(name,upper,lower);
	};
	void setParameterFile(std::string fileName){
		_method->setParameterFile(fileName);
	};
	void setNout(size_t n){
		_method->setNoutFile(n);
	};
	void writeParameters(std::string fileName){
		_method->writeParameters(fileName);
	};
	void readParameters(std::string fileName){
		_method->readParameters(fileName);
	};
	void initCouplings(size_t nseeds){
		minimize::initCouplings(nseeds,-1);
	};
	void initSingleTbin(size_t nseeds, size_t tbin){
		minimize::initCouplings(nseeds,tbin);
	};
	std::string YAML_file(){
		_method->Waveset()->YAML_file();
	};
};

struct chi2py_wrap: chi2py, bp::wrapper<chi2py>{
	
	chi2py_wrap(const chi2py &inin) : chi2py(inin), bp::wrapper<chi2py>(){};
	chi2py_wrap(std::string inin)   : chi2py(inin), bp::wrapper<chi2py>(){};

};


int wolen(){
	return 0;
};

BOOST_PYTHON_MODULE(libchi2py){
	// Expose the std::vectors<...> to python
	
	bp::class_<std::vector<double> >("vector_double")
		.def(bp::vector_indexing_suite<std::vector<double> >())
	;
	bp::class_<std::vector<bool> >("vector_bool")
		.def(bp::vector_indexing_suite<std::vector<bool> >())
	;
	bp::class_<std::vector<std::string> >("vector_string")
		.def(bp::vector_indexing_suite<std::vector<std::string> >())
	;
	bp::class_<std::vector<int> >("vector_int")
		.def(bp::vector_indexing_suite<std::vector<int> >())
	;

	bp::class_<std::vector<std::complex<double> > >("vector_complex")
		.def(bp::vector_indexing_suite<std::vector<std::complex<double> > >())
	;

	bp::class_<std::vector<std::vector<double> > >("vector_array_double")
		.def(bp::vector_indexing_suite<std::vector<std::vector<double> > >())
	;


	bp::class_<std::vector<std::vector<std::vector<double> > > >("vector_array_array_double") // To be able to have a vector of covariance matrices.
		.def(bp::vector_indexing_suite<std::vector<std::vector<std::vector<double> > > >())
	;

	bp::class_<chi2py_wrap> chi2 = bp::class_<chi2py_wrap>("chi2",bp::init<std::string>());


	chi2.def("Eval",				&chi2py::Eval					);
	chi2.def("EvalSelf",				&chi2py::EvalSelf				);

	chi2.def("Amplitudes",				&chi2py::Amplitudes				);
	chi2.def("AmplitudesSelf",			&chi2py::AmplitudesSelf				);

	chi2.def("nTot",				&chi2py::nTot					);
	chi2.def("nCpl",				&chi2py::nCpl					);
	chi2.def("nPar",				&chi2py::nPar					);
	chi2.def("nBra",				&chi2py::nBra					);
	chi2.def("nFtw",				&chi2py::nFtw					);
	chi2.def("nTbin",				&chi2py::nTbin					);
	chi2.def("nWaves",				&chi2py::nWaves					);

	chi2.def("method",				&chi2py::method					);
	chi2.def("className",				&chi2py::className				);
	chi2.def("printStatus",				&chi2py::printStatus				);
	chi2.def("printParameters",			&chi2py::printParameters			);
	chi2.def("setParameter",			&chi2py::setParameter				);
	chi2.def("getParameter",			&chi2py::getParameter				);
	
	chi2.def("parameterNames",			&chi2py::parameterNames				);
	chi2.def("setRandomCpl",			&chi2py::setRandomCpl				);
	chi2.def("setRandomBra",			&chi2py::setRandomBra				);

	chi2.def("fit",					&chi2py::fit					);
	chi2.def("initCouplings",			&chi2py::initCouplings				);
	chi2.def("initSingleTbin",			&chi2py::initSingleTbin			);
	chi2.def("relPar",				&chi2py::relPar					);
	chi2.def("fixPar",				&chi2py::fixPar					);
	chi2.def("parameters",				&chi2py::parameters				);

	chi2.def("write_plots",				&chi2py::write_plots				);
	chi2.def("setParLimits",			&chi2py::setParLimits				);
	chi2.def("setParameterFile",			&chi2py::setParameterFile			);
	chi2.def("setNout",				&chi2py::setNout				);
	chi2.def("writeParameters",			&chi2py::writeParameters			);
	chi2.def("readParameters",			&chi2py::readParameters				);
	chi2.def("setMaxCalls",				&chi2py::setMaxCalls				);
	chi2.def("YAML_file",				&chi2py::YAML_file				);
	chi2.def("waveNames",				&chi2py::waveNames				);
};























