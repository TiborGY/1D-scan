class atom{
    /// stores the element symbol of an atom and its XYZ coordinates
    public:
    char symbol[3];
    double xcoord, ycoord, zcoord;
	//constructors
	inline atom() : symbol{0,0,0}{}
	inline atom(const char symbol_in_1, const char symbol_in_2, const double xcoord_in, const double ycoord_in, const double zcoord_in)
		 : symbol{symbol_in_1, symbol_in_2, 0}, xcoord(xcoord_in), ycoord(ycoord_in), zcoord(zcoord_in){}
	inline atom(const char* const symbol_in, const double xcoord_in, const double ycoord_in, const double zcoord_in)
		 : xcoord(xcoord_in), ycoord(ycoord_in), zcoord(zcoord_in)
	{
		symbol[0] = symbol_in[0];
		symbol[1] = symbol_in[1];
		symbol[2] = 0;
	}
};

class atom_vel{
    public:
    double xvel, yvel, zvel;
};

class xyzfile{
    /// stores the entire contents of a single XYZ-like file
    public:
    uint32_t natoms;
    double energy;
    std::string extraline;
    std::vector<atom> atomvec;
	
	template <typename... Ts>
	xyzfile getXYZ([[maybe_unused]] Ts... dummyPack) const {
		///this provides compatibility with the getXYZ function in multiXYZ
		return *this;
	}
};

class gemfile : public xyzfile{
    /// stores the entire contents of a single XYZ-like file, with trajectory number and step
    /// inherits all of the members of xyzfile
    public:
    uint32_t trajNum, trajStep;
};

class dosfile : public gemfile{
    /// stores the entire contents of a single XYZ-like file, with velocity vectors, trajectory numbers and step numbers
    /// inherits all of the members of gemfile
    public:
    std::vector<atom_vel> velVec;
};

class multiXYZ{
	///stores a collection of XYZ files as a structure of arrays
	///this can only store XYZ files that have the same number of atoms
	public:
	bool isValidated;
	uint32_t natoms;
	std::vector<double> energyVec;
	std::vector<std::string> extralineVec;
	std::vector<char> symbolVec;
	std::vector<double> xcoordVec, ycoordVec, zcoordVec;

	//CONSTRUCTORS
	//the default constructor constructs an empty multiXYZ object
	multiXYZ() : isValidated(false), natoms(12345){}
	
	//construct a multiXYZ object from a single, multi-atom geometry
	// multiXYZ(const uint32_t natoms_in, const double energy_in, const std::string& extraline_in, const char* const symbols_in,
			 // const double* const xcoord_in, const double* const ycoord_in, const double* const zcoord_in) : isValidated(false), natoms(natoms_in)
	// {
		// energyVec.push_back(energy_in);
		// extralineVec.push_back(extraline_in);
		// symbolVec.resize(3*natoms_in);
		// xcoordVec.resize(natoms_in);
		// ycoordVec.resize(natoms_in);
		// zcoordVec.resize(natoms_in);
		// for(uint32_t i=0; i<natoms_in; i++){
			// symbolVec[i*3] = symbols_in[i*3];
			// symbolVec[i*3 + 1] = symbols_in[i*3 + 1];
			// symbolVec[i*3 + 2] = symbols_in[i*3 + 2];
			// xcoordVec[i] = xcoord_in[i];
			// ycoordVec[i] = ycoord_in[i];
			// zcoordVec[i] = zcoord_in[i];
		// }
	// }
	
	//construct a multiXYZ object from a single, multi-atom geometry (STL vector version)
	// multiXYZ(const double energy_in, const std::string& extraline_in, const std::vector<char>& symbols_in,
			 // const std::vector<double>& xcoord_in, const std::vector<double>& ycoord_in,
			 // const std::vector<double>& zcoord_in) : isValidated(false), natoms(xcoord_in.size())
	// {
		// energyVec.push_back(energy_in);
		// extralineVec.push_back(extraline_in);
		// symbolVec = symbols_in;
		// xcoordVec = xcoord_in;
		// ycoordVec = ycoord_in;
		// zcoordVec = zcoord_in;
	// }
	
	//construct a multiXYZ object from a single, multi-atom geometry (STL array version)
	// template<std::size_t SIZE>
	// multiXYZ(const double energy_in, const std::string& extraline_in, const std::array<char,SIZE*3>& symbols_in,
			 // const std::array<double,SIZE>& xcoord_in, const std::array<double,SIZE>& ycoord_in, const std::array<double,SIZE>& zcoord_in)
			 // : isValidated(false), natoms(SIZE), symbolVec(symbols_in.cbegin(),symbols_in.cend()), xcoordVec(xcoord_in.cbegin(),xcoord_in.cend()),
			 // ycoordVec(ycoord_in.cbegin(),ycoord_in.cend()), zcoordVec(zcoord_in.cbegin(),zcoord_in.cend())
	// {
		// energyVec.push_back(energy_in);
		// extralineVec.push_back(extraline_in);
	// }
	
	//construct a multiXYZ object from a single xyzfile object
	// multiXYZ(const xyzfile& xyz_in): isValidated(false), natoms(xyz_in.natoms){
		// energyVec.push_back(xyz_in.energy);
		// extralineVec.push_back(xyz_in.extraline);
		// symbolVec.resize(3*xyz_in.natoms);
		// xcoordVec.resize(xyz_in.natoms);
		// ycoordVec.resize(xyz_in.natoms);
		// zcoordVec.resize(xyz_in.natoms);
		// for(uint32_t i=0; i<xyz_in.natoms; i++){
			// symbolVec[i*3] = xyz_in.atomvec[i].symbol[0];
			// symbolVec[i*3 + 1] = xyz_in.atomvec[i].symbol[1];
			// symbolVec[i*3 + 2] = xyz_in.atomvec[i].symbol[2];
			// xcoordVec[i] = xyz_in.atomvec[i].xcoord;
			// ycoordVec[i] = xyz_in.atomvec[i].ycoord;
			// zcoordVec[i] = xyz_in.atomvec[i].zcoord;
		// }
	// }
	
	//construct a multiXYZ object from an STL vector of xyzfile objects
	// multiXYZ(const std::vector<xyzfile>& xyzvec_in): isValidated(false), natoms(xyzvec_in[0].natoms){
		// const auto Natom = xyzvec_in[0].natoms;
		// const auto Ngeom = xyzvec_in.size();
		// energyVec.resize(Ngeom);
		// extralineVec.resize(Ngeom);
		// symbolVec.resize(3*Ngeom*Natom);
		// xcoordVec.resize(Ngeom*Natom);
		// ycoordVec.resize(Ngeom*Natom);
		// zcoordVec.resize(Ngeom*Natom);
		// for(uint32_t j=0; j<Ngeom; j++){
			// assert(xyzvec_in[j].natoms == Natom);
			// assert(xyzvec_in[j].atomvec.size() == Natom);
			// energyVec[j] = xyzvec_in[j].energy;
			// extralineVec[j] = xyzvec_in[j].extraline;
			// for(uint32_t i=0; i<Natom; i++){
				// symbolVec[j*Natom*3 + i*3] = xyzvec_in[j].atomvec[i].symbol[0];
				// symbolVec[j*Natom*3 + i*3 + 1] = xyzvec_in[j].atomvec[i].symbol[1];
				// symbolVec[j*Natom*3 + i*3 + 2] = xyzvec_in[j].atomvec[i].symbol[2];
				// xcoordVec[i] = xyzvec_in[j].atomvec[i].xcoord;
				// ycoordVec[i] = xyzvec_in[j].atomvec[i].ycoord;
				// zcoordVec[i] = xyzvec_in[j].atomvec[i].zcoord;
			// }
		// }
	// }

	//DATA ACCESS FUNCTIONS
	inline uint32_t size() const {
		if(isValidated){
			return energyVec.size();
		}else{
			std::cout<<"multiXYZ-like object is not validated. size() should not be used unless validated"<<std::endl;
			assert(false);
		}
	}
	inline double& accessX (const uint32_t geomNum, const uint32_t atomNum){
		return xcoordVec[geomNum*natoms + atomNum];
	}
	inline double& accessY (const uint32_t geomNum, const uint32_t atomNum){
		return ycoordVec[geomNum*natoms + atomNum];
	}
	inline double& accessZ (const uint32_t geomNum, const uint32_t atomNum){
		return zcoordVec[geomNum*natoms + atomNum];
	}
	inline char* accessSymbol (const uint32_t geomNum, const uint32_t atomNum){
		isValidated = false; //invalidate object, since changing a symbol can render it inconsistent
		return &symbolVec[geomNum*natoms*3 + atomNum*3];
	}
	inline char* accessSymbol_noninvalidating (const uint32_t geomNum, const uint32_t atomNum){
		///allows write access to the symbols without invalidation, USE WITH CARE
		return &symbolVec[geomNum*natoms*3 + atomNum*3];
	}
	
	inline const double& getX (const uint32_t geomNum, const uint32_t atomNum) const {
		return xcoordVec[geomNum*natoms + atomNum];
	}
	inline const double& getY (const uint32_t geomNum, const uint32_t atomNum) const {
		return ycoordVec[geomNum*natoms + atomNum];
	}
	inline const double& getZ (const uint32_t geomNum, const uint32_t atomNum) const {
		return zcoordVec[geomNum*natoms + atomNum];
	}
	inline const char* getSymbol (const uint32_t geomNum, const uint32_t atomNum) const {
		//since we return a constified pointer, we do not have to invalidate the object
		return &symbolVec[geomNum*natoms*3 + atomNum*3];
	}

	inline const atom getAtom(const uint32_t geomNum, const uint32_t atomNum) const {
		return atom(getSymbol(geomNum, atomNum), getX(geomNum, atomNum), getY(geomNum, atomNum), getZ(geomNum, atomNum));
	}
	xyzfile getXYZ (const uint32_t geomNum) const {
		xyzfile tmp;
		tmp.natoms = natoms;
		tmp.energy = energyVec[geomNum];
		tmp.extraline = extralineVec[geomNum];
		tmp.atomvec.reserve(natoms);
		for(uint32_t i=0; i<natoms; i++){
			tmp.atomvec.push_back(getAtom(geomNum, i));
		}
		return tmp;
	}

	//OTHER FUNCTIONS
	inline void reserve(uint32_t ngeom, uint32_t num_atoms){
		energyVec.reserve(ngeom);
		extralineVec.reserve(ngeom);
		symbolVec.reserve(ngeom*num_atoms*3);
		xcoordVec.reserve(ngeom*num_atoms);
		ycoordVec.reserve(ngeom*num_atoms);
		zcoordVec.reserve(ngeom*num_atoms);
	}
	inline void shrink_to_fit(){
		energyVec.shrink_to_fit();
		extralineVec.shrink_to_fit();
		symbolVec.shrink_to_fit();
		xcoordVec.shrink_to_fit();
		ycoordVec.shrink_to_fit();
		zcoordVec.shrink_to_fit();
	}
	inline void free(){
		isValidated = false;
		natoms = 0;
		energyVec = std::vector<double>();
		extralineVec = std::vector<std::string>();
		symbolVec = std::vector<char>();
		xcoordVec = std::vector<double>();
		ycoordVec = std::vector<double>();
		zcoordVec = std::vector<double>();
	}
};

class multiGEM : public multiXYZ{
    /// derived class of multiXYZ, also stores trajectory numbers and timesteps
	/// PLACEHOLDER
    public:
    std::vector<uint32_t> trajNumVec, trajStepVec;

	//CONSTRUCTORS
	//the default constructor just calls the default constructor of the base class
	multiGEM() : multiXYZ(){}
	
	//construct a multiGEM object from a single, multi-atom geometry
	// multiGEM(const uint32_t natoms_in, const double energy_in, const std::string& extraline_in, const char* const symbols_in,
			 // const double* const xcoord_in, const double* const ycoord_in, const double* const zcoord_in, const uint32_t trajNum_in,
			 // const uint32_t trajStep_in) : multiXYZ(natoms_in, energy_in, extraline_in, symbols_in, xcoord_in, ycoord_in, zcoord_in)
	// {
        // trajNumVec.push_back(trajNum_in);
        // trajStepVec.push_back(trajStep_in);
	// }
	
	//construct a multiGEM object from a single, multi-atom geometry (STL vector version)
	// multiGEM(const double energy_in, const std::string& extraline_in, const std::vector<char>& symbols_in,
			 // const std::vector<double>& xcoord_in, const std::vector<double>& ycoord_in, const std::vector<double>& zcoord_in,
			 // const uint32_t trajNum_in, const uint32_t trajStep_in) : multiXYZ(energy_in, extraline_in, symbols_in, xcoord_in, ycoord_in, zcoord_in)
	// {
        // trajNumVec.push_back(trajNum_in);
        // trajStepVec.push_back(trajStep_in);
	// }
	
	//construct a multiGEM object from a single, multi-atom geometry (STL array version)
	// template<std::size_t SIZE>
	// multiGEM(const double energy_in, const std::string& extraline_in, const std::array<char,SIZE*3>& symbols_in,
			 // const std::array<double,SIZE>& xcoord_in, const std::array<double,SIZE>& ycoord_in, const std::array<double,SIZE>& zcoord_in,
			 // const uint32_t trajNum_in, const uint32_t trajStep_in) : multiXYZ(energy_in, extraline_in, symbols_in, xcoord_in, ycoord_in, zcoord_in)
	// {
        // trajNumVec.push_back(trajNum_in);
        // trajStepVec.push_back(trajStep_in);
	// }
	
	//construct a multiGEM object from a single gemfile object
	// multiGEM(const gemfile& xyz_in, const uint32_t trajNum_in, const uint32_t trajStep_in): multiXYZ(static_cast<const xyzfile&>(xyz_in)){
        // trajNumVec.push_back(trajNum_in);
        // trajStepVec.push_back(trajStep_in);
	// }
	
	gemfile getGEM (const uint32_t geomNum) const {
		gemfile tmp;
		tmp.natoms = natoms;
		tmp.trajNum = trajNumVec[geomNum];
		tmp.trajStep = trajStepVec[geomNum];
		tmp.energy = energyVec[geomNum];
		tmp.extraline = extralineVec[geomNum];
		tmp.atomvec.reserve(natoms);
		for(uint32_t i=0; i<natoms; i++){
			tmp.atomvec.push_back(getAtom(geomNum, i));
		}
		return tmp;
	}

	//OTHER FUNCTIONS
	inline void reserve(uint32_t ngeom, uint32_t num_atoms){
		multiXYZ::reserve(ngeom, num_atoms);//call the same function in the base class
		trajNumVec.reserve(ngeom);
		trajStepVec.reserve(ngeom);
	}
	inline void shrink_to_fit(){
		multiXYZ::shrink_to_fit();
		trajNumVec.shrink_to_fit();
		trajStepVec.shrink_to_fit();
	}
	inline void free(){
		multiXYZ::free();
		trajNumVec = std::vector<uint32_t>();
		trajStepVec = std::vector<uint32_t>();
	}
};

class multiDOS : public multiGEM{
    /// derived class of multiGEM, also stores atom velocity vectors
	/// PLACEHOLDER
    public:
	std::vector<double> xvelVec, yvelVec, zvelVec;
	
};

class array2d{
    uint32_t tmp;
    uint32_t cols;
    //uint32_t rows;
    std::vector<double> vec;

public:
    double energy;
    //array2d(uint32_t col, uint32_t row) : tmp(col*row), cols(col), rows(row), vec(tmp, 0){}
    array2d(uint32_t col, uint32_t row) : tmp(col*row), cols(col), vec(tmp, 0){}
    const double& operator() (const uint32_t col, const uint32_t row) const{
        //assert((cols > col) && (rows > row)); //disabled check for speed
        return vec[cols*col + row];
    }
    double& operator() (const uint32_t col, const uint32_t row){
        //assert((cols > col) && (rows > row)); //disabled check for speed
        return vec[cols*col + row];
    }
    const std::vector<double>& get_vec() const {return vec;}
    std::vector<double>& get_vec() {return vec;}
    uint32_t size() const {return cols;}
};

class distmat{
public:
    std::vector<double> dists;
    double energy;
    distmat(uint32_t ndists) : dists(ndists, 0.0){}
};

class multiDistmat{
public:
    std::vector<double> energyVec;
	std::vector<double> distsVec;
	
	multiDistmat(const uint32_t numGeom, const uint32_t distsPerGeom, const uint32_t PERM_SYMM) :
		energyVec(numGeom), distsVec( static_cast<uint64_t>(numGeom) * static_cast<uint64_t>(PERM_SYMM) * static_cast<uint64_t>(distsPerGeom) ){}
	
	inline uint32_t size_noPerm() const {
		return energyVec.size();
	}
	
	template<uint32_t PERM_SYMM>
	inline uint64_t size_perm() const {
		return static_cast<uint64_t>(energyVec.size()) * static_cast<uint64_t>(PERM_SYMM);
	}
	
	template<uint32_t distsPerGeom, uint32_t PERM_SYMM>
	inline void check() const {
		const uint64_t geomCount = static_cast<uint64_t>(distsVec.size()) / ( static_cast<uint64_t>(distsPerGeom) * static_cast<uint64_t>(PERM_SYMM) );
		assert(energyVec.size() == geomCount);
	}
	
	template<uint32_t distsPerGeom>
	inline void reserve(const uint32_t res){
		energyVec.reserve(res);
		distsVec.reserve( static_cast<uint64_t>(res) * static_cast<uint64_t>(distsPerGeom) );
	}
	
	inline void free(){
		energyVec = std::vector<double>();
		distsVec = std::vector<double>();
	}
	// inline void resize(const uint64_t res, const uint64_t distsPerGeom){
		// energyVec.resize(res);
		// distsVec.resize(res*distsPerGeom);
	// }
};

struct gemDiff{
    uint64_t index;
    uint32_t trajNum, trajStep;
    double minEdiff, globalminEdiff;
    double smallestRMSD;
};

class multiXYZdiff{
public:
	std::vector<uint32_t> indexVec;
	std::vector<double> minEdiffVec, globalminEdiffVec, smallestRMSDvec;
	
	multiXYZdiff(const uint32_t numGeom) : indexVec(numGeom), minEdiffVec(numGeom), globalminEdiffVec(numGeom), smallestRMSDvec(numGeom){}
	
	inline uint32_t size() const {
		return indexVec.size();
	}
	inline void check() const {
		assert(indexVec.size() == minEdiffVec.size());
		assert(indexVec.size() == globalminEdiffVec.size());
		assert(indexVec.size() == smallestRMSDvec.size());
	}
	inline void reserve(const uint32_t res){
		indexVec.reserve(res);
		minEdiffVec.reserve(res);
		globalminEdiffVec.reserve(res);
		smallestRMSDvec.reserve(res);
	}
	//sorting related functions, they use selection sort
	inline void swapAll(const uint32_t i, const uint32_t j){
		std::swap(indexVec[i], indexVec[j]);
		std::swap(minEdiffVec[i], minEdiffVec[j]);
		std::swap(globalminEdiffVec[i], globalminEdiffVec[j]);
		std::swap(smallestRMSDvec[i], smallestRMSDvec[j]);
	}
	inline void sortByRMSD_desc(){
		check();
		const uint32_t numGeom = size();
		for(uint32_t i=0; i<(numGeom-1); i++){
			uint32_t largestRMSDindex = i;
			for(uint32_t j=i+1; j<numGeom; j++){
				if(smallestRMSDvec[j] > smallestRMSDvec[largestRMSDindex]) largestRMSDindex = j;
			}
			swapAll(i, largestRMSDindex);
		}
	}
	inline void sortByEdiff_desc(){
		check();
		const uint32_t numGeom = size();
		for(uint32_t i=0; i<(numGeom-1); i++){
			uint32_t largestEdiffIndex = i;
			for(uint32_t j=i+1; j<numGeom; j++){
				if(minEdiffVec[j] > minEdiffVec[largestEdiffIndex]) largestEdiffIndex = j;
			}
			swapAll(i, largestEdiffIndex);
		}
	}
	inline void sortByGlobalminEdiff_asc(){
		check();
		const uint32_t numGeom = size();
		for(uint32_t i=0; i<(numGeom-1); i++){
			uint32_t smallestGlobalminEdiffIndex = i;
			for(uint32_t j=i+1; j<numGeom; j++){
				if(globalminEdiffVec[j] < globalminEdiffVec[smallestGlobalminEdiffIndex]) smallestGlobalminEdiffIndex = j;
			}
			swapAll(i, smallestGlobalminEdiffIndex);
		}
	}
};

class multiGEMdiff : public multiXYZdiff{
public:
	std::vector<uint32_t> trajNumVec, trajStepVec;
	
	multiGEMdiff(const uint32_t numGeom) : multiXYZdiff(numGeom), trajNumVec(numGeom), trajStepVec(numGeom){}
	
	inline void check() const {
		multiXYZdiff::check();
		assert(indexVec.size() == trajNumVec.size());
		assert(indexVec.size() == trajStepVec.size());
	}
	inline void reserve(const uint32_t res){
		multiXYZdiff::reserve(res);
		trajNumVec.reserve(res);
		trajStepVec.reserve(res);
	}
};