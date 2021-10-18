#include <memory>
#include <random>

inline void pauser(){
    /// a portable way to pause a program
    std::string dummy;
    std::cout << "Press enter to continue...";
    std::getline(std::cin, dummy);
}

inline int64_t getFilesize(const std::string& fname){
    /// return with the size of a file in bytes, return -1 for error
    std::ifstream f;
    f.open(fname, std::ios::in | std::ios::binary); //open the file and check for error
    if (!f.is_open()){
       std::cout<< "error: open file for size check failed!" <<std::endl;
       std::cout<< "Cannot open file: " << fname << std::endl;
       return -1;
    }
    f.seekg (0, f.end); //getting the size of the file
    return f.tellg(); //no need to close f, the destructor handles it
}

inline bool FileExists(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {return false;}
}

inline int CallBash(const std::string& cmd){
  const int ret = std::system(cmd.c_str());
  #pragma omp critical
  std::cout << "Return value of " << cmd << " : " << ret << std::endl;
  return ret;
}

inline int CallBash_quiet(const std::string& cmd){
	const int ret = std::system(cmd.c_str());
	#pragma omp critical
	if (ret!=0) std::cout << "Return value of " << cmd << " : " << ret << std::endl;
	return ret;
}

std::string GetCommandOutput(const std::string& cmdstr) {
    ///based on https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
    const char * const cmd = cmdstr.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe){
        std::cout<<"popen() failed!"<<std::endl;
        assert(pipe);
    }
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

template <typename T>
std::string to_string_prec(const T a_value, const int n = 8){
    std::ostringstream out;
    out << std::setprecision(n) << std::fixed << a_value;
    return out.str();
}

template <typename T>
std::string getSpacer(const T value, const uint32_t maxWidth){
	std::string spacer;
	for(uint32_t i=1; i<maxWidth; i++){
		if (value < std::pow(10, i)) spacer += ' ';
	}
	return spacer;
}

std::mt19937_64 initPRNG(){
  //prepare a uniform PRNG
  std::random_device rd;
  std::mt19937_64 PRNG(rd()); //seed PRNG using /dev/urandom or similar OS provided RNG
  std::uniform_int_distribution<> tmpdist{0, 255};
  //make sure the internal state of the PRNG is properly mixed by generating 10M random numbers
  //PRNGs often have unreliable distribution uniformity and other statistical properties before their internal state is sufficiently mixed
  int tmprngint;
  for (uint32_t i=0;i<10000000;i++) tmprngint=tmpdist(PRNG);
  return PRNG;
}

inline bool fpeq(const double a, const double b){
	if(std::abs(a-b) < 1.0E-14) return true;
	return false;
}

inline bool fpwhole(const double a){
	return fpeq(std::round(a), a);
}

inline std::string prettybool(const bool a){
	if (a){return "YES";}else{return "NO";}
}

std::vector<std::string> getCLIargs(const int argc, const char * const argv[]){
	std::vector<std::string> cliArgs;
	for (int i=1; i<argc; i++) cliArgs.push_back(argv[i]);
	return cliArgs;
}

void saveUint(const std::string fname, const uint32_t i){
	std::ofstream outfile;
	outfile.open(fname);
	assert(outfile.is_open());
	outfile<< i<< std::endl;
	assert(outfile.good());
}

uint32_t loadUint(const std::string fname){
	uint32_t tmp;
	std::ifstream infile;
	infile.open(fname);
	assert(infile.is_open());
	infile >> tmp;
	return tmp;
}
