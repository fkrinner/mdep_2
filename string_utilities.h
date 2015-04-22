#ifndef STOLEN_DATETIME_FUNCTION
#define STOLEN_DATETIME_FUNCTION
#include<ctime>
const std::string currentDateTime(){
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
};

std::string dirname(const std::string inFile, const char sep = '/'){
	bool found = false;
	size_t last_sep = 0;
	size_t last_first_sep = 0;
	for (size_t i=0;i<inFile.size();++i){
		if (inFile.at(i) == sep){
			if (i>last_sep+1){
				last_first_sep = i;
			};			
			last_sep = i;
			found = true;
		};
	};

	if (found){
		return inFile.substr(0,last_first_sep+1);
	}else{
		return inFile;
	};
};


std::string get_relative_path(const std::string in_path, const std::string card_file, const std::string marker = std::string("<card_path>") ,const char sep = '/'){
	size_t position = in_path.find(marker);
	if (position == std::string::npos){
		return in_path;
	};
	if (position != 0){
		std::cerr<<"string_utilities.h::get_relative_path(...): Error: Marker found at wrong position > 0"<<std::endl;
		throw;
	};
	if (in_path.size() == marker.size()){ // Complicated to avoid double separatprs
		return dirname(card_file)+in_path.substr(marker.size(),in_path.size()-1);
	};
	if (in_path.at(marker.size()) == sep){
		return dirname(card_file)+in_path.substr(marker.size()+1,in_path.size()-1);
	}else{
		return dirname(card_file)+in_path.substr(marker.size(),in_path.size()-1);
	};
};
#endif//STOLEN_DATETIME_FUNCTION
