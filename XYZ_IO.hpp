template <typename T>
inline void printXYZ(const T& geom, const uint32_t prec = 8){
    /// print the contents of an xyzfile-like struct to the screen
	/// the energy is printed by just printing the extraline in the object
	/// prec adjusts the precision of the coordinates
	std::cout << std::setprecision(prec);
    std::cout << geom.natoms;
    if constexpr (std::is_same<T,xyzfile>::value){
        std::cout << '\n';
    }else{
        std::cout << "   #" << geom.trajNum << '#' << geom.trajStep <<'\n';
    }
    std::cout << geom.extraline << '\n';
    for (uint32_t i=0; i<geom.natoms; i++){
		std::cout << geom.atomvec[i].symbol << "   ";
        if (geom.atomvec[i].symbol[1]==0) std::cout << ' '; //add extra space after 1 letter symbols, to keep fixed width
        std::cout << getSpacer(geom.atomvec[i].xcoord,3) << geom.atomvec[i].xcoord << ' ';
        std::cout << getSpacer(geom.atomvec[i].ycoord,3) << geom.atomvec[i].ycoord << ' ';
        std::cout << getSpacer(geom.atomvec[i].zcoord,3) << geom.atomvec[i].zcoord;
		if constexpr (std::is_same<T,dosfile>::value){ //print velocity vector in the same line, if present
			std::cout << "   ";
			std::cout << getSpacer(geom.velVec[i].xvel,3) << geom.velVec[i].xvel << ' ';
			std::cout << getSpacer(geom.velVec[i].yvel,3) << geom.velVec[i].yvel << ' ';
			std::cout << getSpacer(geom.velVec[i].zvel,3) << geom.velVec[i].zvel;
		}
		std::cout<<'\n';
    }
	std::cout<<std::flush; //not strictly necessary to do this, just to have a little bit smoother behavior
}

template <typename T>
int32_t CheckMultiXYZ_InternalConsistence (const std::vector<T>& inputVec){
	/// check the internal consistency of a vector of xyzfile-like objects
	/// this makes sure the number of atoms, the atom types and the atom order is the same in all objects in the vector
	/// also checks that the atom symbol strings are properly null-terminated
    for (uint32_t i = 0;  i < inputVec.size();  i++){
        if ((inputVec[i].natoms != inputVec.front().natoms)||(inputVec[i].natoms != inputVec[i].atomvec.size())){
            std::cout << "Error: geometry #" << i << "has an incorrect number of atoms" << std::endl;
            printXYZ(inputVec[i]);
            return -1;
        }
        bool atomerr=false;
        for (uint32_t j = 0;  j < inputVec.front().natoms;  j++){
            if (inputVec[i].atomvec[j].symbol[0] != inputVec.front().atomvec[j].symbol[0]) atomerr=true;
            if (inputVec[i].atomvec[j].symbol[1] != inputVec.front().atomvec[j].symbol[1]) atomerr=true;
            if (inputVec[i].atomvec[j].symbol[2] != 0) atomerr=true; //make sure the string is null terminated
        }
        if (atomerr == true){
            std::cout << "Error: geometry #" << i << "has incorrect atom order, atom type or is missing null termination" << std::endl;
            std::cout << "If this is the last point in the file, you probably have some empty lines at the end" << std::endl;
            printXYZ(inputVec[i]);
            return -2;
        }
    }
    return 0;
}

template <typename T>
int32_t CheckMultiXYZ_InternalConsistence (const T& input){
	/// Checks the internal consistency of a of multiXYZ-like object.
	/// This makes sure the atom types and the atom order is the same in all of the geometries.
	/// Also checks that the concatenated atom symbol strings are properly null-terminated.
	if (input.energyVec.size() != input.extralineVec.size()) return -1;
	if ((input.energyVec.size()*input.natoms*3) != input.symbolVec.size()) return -1;
	if ((input.energyVec.size()*input.natoms) != input.xcoordVec.size()) return -1;
	if ((input.energyVec.size()*input.natoms) != input.ycoordVec.size()) return -1;
	if ((input.energyVec.size()*input.natoms) != input.zcoordVec.size()) return -1;
	if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
		if (input.energyVec.size() != input.trajNumVec.size()) return -1;
		if (input.energyVec.size() != input.trajStepVec.size()) return -1;
	}
	if constexpr (std::is_same<T,multiDOS>::value){
		if ((input.energyVec.size()*input.natoms) != input.xvelVec.size()) return -1;
		if ((input.energyVec.size()*input.natoms) != input.yvelVec.size()) return -1;
		if ((input.energyVec.size()*input.natoms) != input.zvelVec.size()) return -1;
	}
    for (uint32_t i = 0;  i < input.energyVec.size();  i++){
        bool atomerr=false;
        for (uint32_t j = 0;  j < input.natoms;  j++){
            if (input.getSymbol(i,j)[0] != input.getSymbol(0,j)[0]) atomerr=true;
            if (input.getSymbol(i,j)[1] != input.getSymbol(0,j)[1]) atomerr=true;
            if (input.getSymbol(i,j)[2] != 0) atomerr=true; //make sure the string is null terminated
        }
        if (atomerr == true){
            std::cout << "Error: geometry #" << i << "has incorrect atom order, atom type or is missing null termination" << std::endl;
            std::cout << "If this is the last point in the file, you probably have some empty lines at the end" << std::endl;
            printXYZ(input.getXYZ(i));
            return -2;
        }
    }
    return 0;
}

template <typename T1, typename T2>
int32_t CheckMultiXYZ_Correspondance(const T1& A, const T2& B){
	/// Make sure the two multiXYZ-like objects are consistent with each other.
	/// Checks the number of atoms, the type of atoms and the atom order.
    if (A.natoms != B.natoms){
        std::cout << "Error: the objects have a different number of atoms" << std::endl;
        return -1;
    }
    for (uint32_t i = 0;  i < A.natoms;  i++){
        if ((A.getSymbol(0,i)[0] != B.getSymbol(0,i)[0]) || (A.getSymbol(0,i)[1] != B.getSymbol(0,i)[1])){
            std::cout << "Error: the objects have a different ordering or type of atoms" << std::endl;
            return -2;
        }
    }
    return 0;
}

template <typename T>
inline T parseXYZ(std::ifstream& infile, bool& reject, const double XYZPARSE_E_UPPER_LIMIT, const double XYZPARSE_E_LOWER_LIMIT, const double XYZPARSE_COORD_LIMIT){
    /// parses the contents of an already opened XYZ file into a struct
    reject = false;
    T pesPoint;
    atom atombuff; // temporary struct for reading in the atom coordinates
    infile >> pesPoint.natoms; // number of atoms in the xyz file
    assert(pesPoint.natoms > 0); // make sure we got the number of atoms successfully
	assert(pesPoint.natoms < 100); //hopefully this line turns into a bug one day :)
    if constexpr (std::is_same<T,xyzfile>::value){
		std::getline(infile, pesPoint.extraline); // get to the second line
    }else{
        std::getline(infile, pesPoint.extraline, '#'); //peel off the first hashtag
        std::getline(infile, pesPoint.extraline, '#'); //extract the trajectory number
        while((pesPoint.extraline[0]=='0')&&(pesPoint.extraline.size()>1)) pesPoint.extraline.erase(0, 1); //remove leading zeros
        pesPoint.trajNum = std::stoul(pesPoint.extraline);
        std::getline(infile, pesPoint.extraline); //extract the timestep number
        while((pesPoint.extraline[0]=='0')&&(pesPoint.extraline.size()>1)) pesPoint.extraline.erase(0, 1); //remove leading zeros
        pesPoint.trajStep = std::stoul(pesPoint.extraline);
    }
    std::getline(infile, pesPoint.extraline); // XYZ files contain a line with arbitrary data
    //XYZ files written by Molpro contain some text before the energy, this text needs to be removed
    if constexpr (std::is_same<T,xyzfile>::value){
        if (pesPoint.extraline.find(":") != std::string::npos) pesPoint.extraline.erase(0, 1+pesPoint.extraline.find(":"));
    }
    pesPoint.energy = std::stod(pesPoint.extraline);
    pesPoint.atomvec.reserve(pesPoint.natoms);
	//these constants need to be tuned for the particular PES and correlation level/basis set/ECP/etc.
    if ((pesPoint.energy < XYZPARSE_E_UPPER_LIMIT)&&(pesPoint.energy > XYZPARSE_E_LOWER_LIMIT)){
        for (uint32_t i=0; i<pesPoint.natoms; i++){
            infile >> atombuff.symbol;
            infile >> atombuff.xcoord;
            if (std::abs(atombuff.xcoord) > XYZPARSE_COORD_LIMIT) reject = true;
            infile >> atombuff.ycoord;
            if (std::abs(atombuff.ycoord) > XYZPARSE_COORD_LIMIT) reject = true;
            infile >> atombuff.zcoord;
            if (std::abs(atombuff.zcoord) > XYZPARSE_COORD_LIMIT) reject = true;
            pesPoint.atomvec.push_back(atombuff);
            if constexpr (std::is_same<T,dosfile>::value){
                atom_vel velbuff; // temporary struct for reading in the atom velocities
                infile >> velbuff.xvel;
                infile >> velbuff.yvel;
                infile >> velbuff.zvel;
                pesPoint.velVec.push_back(velbuff);
            }
        }
    }else{
        reject = true;
        for (uint32_t i=0; i<pesPoint.natoms; i++){
            infile >> atombuff.symbol;
            infile >> atombuff.xcoord;
            infile >> atombuff.ycoord;
            infile >> atombuff.zcoord;
            pesPoint.atomvec.push_back(atombuff);
            if constexpr (std::is_same<T,dosfile>::value){
                atom_vel velbuff; // temporary struct for reading in the atom velocities
                infile >> velbuff.xvel;
                infile >> velbuff.yvel;
                infile >> velbuff.zvel;
                pesPoint.velVec.push_back(velbuff);
            }
        }
    }
    return pesPoint;
}

int32_t InspectXYZfile(const std::string& inputFileName, bool* const trajMarks = nullptr,
						bool* const velvecs = nullptr, uint32_t* const estNumGeom = nullptr, uint32_t* const estExtralineLen = nullptr)
{
	///Inspects an XYZ-like file to figure out if it is a plain XYZ file or a GEM or DOS file.
	///Also reports the number of atoms, and an estimated number of geometries.
	///This function is not hardened against buffer overflow in any way. Do not use this on untrusted files.
	uint32_t natoms;
	std::string linebuf;
	double dummy;
	char dummy_c[3] = {0,0,0};
	std::ifstream infile;
	const int64_t infileSize = getFilesize(inputFileName);
    if (infileSize <= 10){ //a valid XYZ file below 11 bytes is impossible
		std::cout<<"Error: impossibly small XYZ file: " << inputFileName << std::endl;
		return -2;
	}
	infile.open(inputFileName);
    if (!infile.is_open()){
        std::cout << "Error: cannot open XYZ-like input file: " << inputFileName << std::endl;
        return -1;
    }
	infile >> natoms; //read the number of atoms in the first geometry
	
	//if the presence of trajectory numbers or velocity vectors is requested
	if ((trajMarks != nullptr)||(velvecs != nullptr)||(estExtralineLen != nullptr)){
		std::getline(infile, linebuf); //get to the second line and extract trajectory marks if present
		if (trajMarks != nullptr){
			if (linebuf.find('#') != std::string::npos){
				*trajMarks = true;
			}else{
				*trajMarks = false;
			}
		}
		//at this point we are at the beginning of the second line
		if ((velvecs != nullptr)||(estExtralineLen != nullptr)){
			std::getline(infile, linebuf); //skip the second line
			if (estExtralineLen != nullptr){
				*estExtralineLen = (11*linebuf.length())/10; //take 110% as our estimate
			}
			if (velvecs != nullptr){
				infile >> dummy_c; //skip the symbol of the first atom
				for (uint32_t i=0;i<3;i++) infile >> dummy; //skip the xyz coords of the first atom
				std::getline(infile, linebuf); //extract any remaining characters in the line
				//if there are no non-whitespace characters after the z coordinate, we have no velocities
				if (std::all_of(linebuf.cbegin(), linebuf.cend(), [](const unsigned char c){return std::isspace(c);})){
					*velvecs = false;
				}else{
					*velvecs = true;
				}
			}
		}
	}
	//estimate the number of geometries in the file
	if(estNumGeom != nullptr){
		uint32_t bound;
		if (velvecs != nullptr){ //if we scanned for velocities we have already gone past 1 atom
			bound = (natoms-1);
		}else{
			bound = natoms;
		}
		for(uint32_t j=0; j<bound; j++){
			infile >> dummy_c;
			for (uint32_t i=0;i<3;i++) infile >> dummy;
			if(*velvecs) for (uint32_t i=0;i<3;i++) infile >> dummy;
		}
		//try to estimate an upper bound, by taking 110% of the size of the fisrt geometry
		//the cast to int64 is necessary to handle files larger than ~2 GiB
		*estNumGeom = (11*(infileSize/(static_cast<int64_t>(infile.tellg())+1)))/10;
	}
	return natoms;
}

template <typename T>
int64_t LoadMultiXYZ(T& output, const std::string& inputFileName, const double XYZPARSE_E_UPPER_LIMIT = 1234.0,
					 const double XYZPARSE_E_LOWER_LIMIT = -1234.0, const double XYZPARSE_COORD_LIMIT = 1234.0, T* const rejectVec = nullptr){
	///Loads one or more xyz structures from a file into a multiXYZ object.
	///Since multiXYZ can only handle geometries with the same number of atoms, this function assumes that is satisfied by the file.
	///This can optionally filter the geometries, and save rejected geometries in a separate multiXYZ object.
	///Since this function works with "absolute" energies, the lower/upper limits need to be tuned for the particular PES and correlation level/basis set/ECP/etc.
	///This function is not hardened against buffer overflow in any way. Do not use this to load untrusted files.
	//TODO: investigate std::from_chars(...) performance instead of the implicit std::stod
	bool trajMarks, velvecs;
	uint32_t estNumGeom;
	uint32_t acceptCount = 0;
	uint32_t rejectCount = 0;
	std::string tmp_str;
	std::ifstream infile;
	
	const auto start = std::chrono::steady_clock::now(); //start timing the entire function
	const int64_t infileSize = getFilesize(inputFileName);
    assert(infileSize > 10); //a valid XYZ file below 11 bytes is impossible, this also catches if getFilesize fails
	const int32_t natoms_tmp = InspectXYZfile(inputFileName, &trajMarks, &velvecs, &estNumGeom);
	assert(natoms_tmp > 0);
	const uint32_t natoms = natoms_tmp;
	output.natoms = natoms;
	if (rejectVec != nullptr) rejectVec->natoms = natoms;
	output.reserve(estNumGeom, natoms); //preallocate some memory to avoid slow reallocations while reading
    infile.open(inputFileName);
    if (!infile.is_open()){
        std::cout << "Error: cannot open XYZ-like input file: " << inputFileName << std::endl;
        return -1;
    }
	std::cout << "Loading geometries from: " << inputFileName << '\n';
	std::cout << "Number of atoms: " << natoms << std::endl;
	std::cout << "Estimated number of geometries: " << estNumGeom << std::endl;
	if (trajMarks) std::cout << "Trajectory number and timestep data found"<<std::endl;
	if (velvecs) std::cout << "Velocity vectors found"<<std::endl;
	
	//temporary buffers for 1 complete geometry
	std::vector<char> symbols_tmp(natoms*3);
	std::vector<double> xcoords_tmp(natoms);
	std::vector<double> ycoords_tmp(natoms);
	std::vector<double> zcoords_tmp(natoms);
	
	std::vector<double> xvels_tmp(natoms);
	std::vector<double> yvels_tmp(natoms);
	std::vector<double> zvels_tmp(natoms);
	
	//main loop that reads in the XYZ-like file
	while(infile.tellg() < (infileSize-10)){
		uint32_t tmp_int;
		infile >> tmp_int; //read the number of atoms
		assert(tmp_int == natoms); //multiXYZ-like files can only store geometries with the same number of atoms
		
		//at this point we are still in the first line, right after the number of atoms
		if (trajMarks){ //if trajnum/step is present they need to be handled
			if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){ //if we are reading into an object that can store trajnum/step
				std::getline(infile, tmp_str, '#'); //peel off the first hashtag
				assert(tmp_str.size() < 1000); //sanity check
				std::getline(infile, tmp_str, '#'); //extract the trajectory number
				assert(tmp_str.size() < 1000); //sanity check
				while((tmp_str[0]=='0')&&(tmp_str.size()>1)) tmp_str.erase(0, 1); //remove leading zeros before the trajectory number
				output.trajNumVec.push_back(std::stoul(tmp_str)); //convert the string into integer and store it
				std::getline(infile, tmp_str); //extract the timestep number
				assert(tmp_str.size() < 1000); //sanity check
				while((tmp_str[0]=='0')&&(tmp_str.size()>1)) tmp_str.erase(0, 1); //remove leading zeros before the timestep number
				output.trajStepVec.push_back(std::stoul(tmp_str));
			}else{ //there is trajnum/step data but the object cannot store them, so get to the second line
				std::getline(infile, tmp_str);
				assert(tmp_str.size() < 1000); //sanity check
			}
		}else{ //no trajnum/step in the file
			if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){ //the object needs the data, so just give it a fake "magic" value
				output.trajNumVec.push_back(12345);
				output.trajStepVec.push_back(12345);
			}
			std::getline(infile, tmp_str); //get to the second line
			assert(tmp_str.size() < 1000); //sanity check
		}
		
		//at this point we are at the second line, there can be arbitrary data here, but usually XYZ-like files have an energy value here
		std::getline(infile, tmp_str); //read the second line
		assert(tmp_str.size() < 1000); //sanity check
		//at this point we are at the third line, this is the first line containing atom symbols, coords and potentially velocities
		
		//XYZ files written by Molpro contain some text before the energy, this text needs to be removed
		if (tmp_str.find(":") != std::string::npos) tmp_str.erase(0, 1+tmp_str.find(":"));
		
		//convert the string to double and check if it is inside the acceptable energy range
		if(const double energy_tmp = std::stod(tmp_str); (energy_tmp <= XYZPARSE_E_UPPER_LIMIT) && (energy_tmp >= XYZPARSE_E_LOWER_LIMIT)){
			//push the energy and the extraline into the output
			output.energyVec.push_back(energy_tmp);
			output.extralineVec.push_back(tmp_str);
			//read in the symbols and coords
			for(uint32_t i=0; i<natoms; i++){
				//yes, this is a buffer overflow hazard, do not open untrusted files
				infile>>&symbols_tmp[i*3];
				infile>>xcoords_tmp[i];
				infile>>ycoords_tmp[i];
				infile>>zcoords_tmp[i];
				//if the object can store velocities read that in as well
				if constexpr (std::is_same<T,multiDOS>::value){
					infile>>xvels_tmp[i];
					infile>>yvels_tmp[i];
					infile>>zvels_tmp[i];
				}
			}
			//check if any of the coordinates is too large
			//declare the comparison function via a lambda expression, this effectively declares an is_larger(double A){...} function
			auto is_larger = [XYZPARSE_COORD_LIMIT](const double A){return A > XYZPARSE_COORD_LIMIT;};
			if (std::any_of(xcoords_tmp.begin(), xcoords_tmp.end(), is_larger) || std::any_of(ycoords_tmp.begin(), ycoords_tmp.end(), is_larger) || std::any_of(zcoords_tmp.begin(), zcoords_tmp.end(), is_larger)){
				//at leat one coordinate is > the coord limit, reject the geometry
				rejectCount++;
				if (rejectVec != nullptr){ //if we are collecting rejected geometries
					//if the object can store trajnum/step they have already been pushed, and we need to copy them to the rejectVec and pop them off the main output
					if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
						rejectVec->trajNumVec.push_back(output.trajNumVec.back());
						rejectVec->trajStepVec.push_back(output.trajStepVec.back());
						output.trajNumVec.pop_back();
						output.trajStepVec.pop_back();
					}
					//copy, then pop off the energy and the extraline
					rejectVec->energyVec.push_back(output.energyVec.back());
					rejectVec->extralineVec.push_back(output.extralineVec.back());
					output.energyVec.pop_back();
					output.extralineVec.pop_back();
					//copy the data from the temporary array to the reject-output
					//std::vector::insert works a lot like a multi-element push_back in this arrangement
					rejectVec->symbolVec.insert(rejectVec->symbolVec.end(), symbols_tmp.begin(), symbols_tmp.end());
					rejectVec->xcoordVec.insert(rejectVec->xcoordVec.end(), xcoords_tmp.begin(), xcoords_tmp.end());
					rejectVec->ycoordVec.insert(rejectVec->ycoordVec.end(), ycoords_tmp.begin(), ycoords_tmp.end());
					rejectVec->zcoordVec.insert(rejectVec->zcoordVec.end(), zcoords_tmp.begin(), zcoords_tmp.end());
					if constexpr (std::is_same<T,multiDOS>::value){
						rejectVec->xvelVec.insert(rejectVec->xvelVec.end(), xvels_tmp.begin(), xvels_tmp.end());
						rejectVec->yvelVec.insert(rejectVec->yvelVec.end(), yvels_tmp.begin(), yvels_tmp.end());
						rejectVec->zvelVec.insert(rejectVec->zvelVec.end(), zvels_tmp.begin(), zvels_tmp.end());
					}
				}else{ //geometry rejected due to coord limit, we are not collecting rejects
					//if the object can store trajnum/step they have already been pushed, and we need to pop them off
					if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
						output.trajNumVec.pop_back();
						output.trajStepVec.pop_back();
					}
					output.energyVec.pop_back();
					output.extralineVec.pop_back();
				}
			}else{ //the coordinates are all <= the coord limit
				acceptCount++;
				//copy the data from the temporary array to the output
				//std::vector::insert works a lot like a multi-element push_back in this arrangement
				output.symbolVec.insert(output.symbolVec.end(), symbols_tmp.begin(), symbols_tmp.end());
				output.xcoordVec.insert(output.xcoordVec.end(), xcoords_tmp.begin(), xcoords_tmp.end());
				output.ycoordVec.insert(output.ycoordVec.end(), ycoords_tmp.begin(), ycoords_tmp.end());
				output.zcoordVec.insert(output.zcoordVec.end(), zcoords_tmp.begin(), zcoords_tmp.end());
				if constexpr (std::is_same<T,multiDOS>::value){
					output.xvelVec.insert(output.xvelVec.end(), xvels_tmp.begin(), xvels_tmp.end());
					output.yvelVec.insert(output.yvelVec.end(), yvels_tmp.begin(), yvels_tmp.end());
					output.zvelVec.insert(output.zvelVec.end(), zvels_tmp.begin(), zvels_tmp.end());
				}
			} // end of coordinate based rejection if
		}else{ //the energy is too low or too high, reject this geometry
			rejectCount++;
			if (rejectVec != nullptr){ //if we are collecting rejected geometries
				//if the object can store trajnum/step they have already been pushed, and we need to copy them to the rejectVec and pop them off the main output
				if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
					rejectVec->trajNumVec.push_back(output.trajNumVec.back());
					rejectVec->trajStepVec.push_back(output.trajStepVec.back());
					output.trajNumVec.pop_back();
					output.trajStepVec.pop_back();
				}
				//push the energy and the extraline into the rejectVec
				rejectVec->energyVec.push_back(energy_tmp);
				rejectVec->extralineVec.push_back(tmp_str);
				//read in the symbols and coords
				for(uint32_t i=0; i<natoms; i++){
					//yes, this is a buffer overflow hazard, do not open untrusted files
					infile>>&symbols_tmp[i*3];
					infile>>xcoords_tmp[i];
					infile>>ycoords_tmp[i];
					infile>>zcoords_tmp[i];
					//if the object can store velocities read that in as well
					if constexpr (std::is_same<T,multiDOS>::value){
						infile>>xvels_tmp[i];
						infile>>yvels_tmp[i];
						infile>>zvels_tmp[i];
					}
				}
				//copy the data from the temporary array to the reject-output
				//std::vector::insert works a lot like a multi-element push_back in this arrangement
				rejectVec->symbolVec.insert(rejectVec->symbolVec.end(), symbols_tmp.begin(), symbols_tmp.end());
				rejectVec->xcoordVec.insert(rejectVec->xcoordVec.end(), xcoords_tmp.begin(), xcoords_tmp.end());
				rejectVec->ycoordVec.insert(rejectVec->ycoordVec.end(), ycoords_tmp.begin(), ycoords_tmp.end());
				rejectVec->zcoordVec.insert(rejectVec->zcoordVec.end(), zcoords_tmp.begin(), zcoords_tmp.end());
				if constexpr (std::is_same<T,multiDOS>::value){
					rejectVec->xvelVec.insert(rejectVec->xvelVec.end(), xvels_tmp.begin(), xvels_tmp.end());
					rejectVec->yvelVec.insert(rejectVec->yvelVec.end(), yvels_tmp.begin(), yvels_tmp.end());
					rejectVec->zvelVec.insert(rejectVec->zvelVec.end(), zvels_tmp.begin(), zvels_tmp.end());
				}
			}else{ //the geometry is rejected due to energy, and we are not collecting it
				//if the object can store trajnum/step they have already been pushed, and we need to pop them off
				if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
					output.trajNumVec.pop_back();
					output.trajStepVec.pop_back();
				}
				//skip natom number of lines
				for(uint32_t i=0; i<natoms; i++) std::getline(infile, tmp_str);
			}
		} //end of energy-based rejection if
	} //end of main while loop
	
	//at this point all geometries have been read from the file, no need to close it, the destructor will
	std::cout<<acceptCount<<" geometries parsed and accepted in ";
	const auto end = std::chrono::steady_clock::now();
	std::cout << std::chrono::duration <double> (end - start).count() << " seconds" << std::endl << rejectCount;
    if (rejectVec != nullptr){
        std::cout<<" geometries rejected, but saved"<<std::endl;
    }else{
        std::cout<<" geometries rejected and discarded"<<std::endl;
    }
	
	// make sure all geometries have the same number of atoms and the same order of atoms, as the first geometry
	std::cout << "Checking accepted geometries for inconsistencies" << std::endl;
	if (const int32_t ret = CheckMultiXYZ_InternalConsistence(output); ret!=0) return ret;
	output.isValidated = true;
	std::cout<<"All accepted geometries checked"<<std::endl;
	if (rejectVec != nullptr){
		std::cout << "Checking rejected geometries for inconsistencies" << std::endl;
		if (const int32_t ret = CheckMultiXYZ_InternalConsistence(*rejectVec); ret!=0) return ret;
		rejectVec->isValidated = true;
		std::cout<<"All rejected geometries checked"<<std::endl;
	}
	return 0;
}

template <typename T>
int32_t LoadMultiXYZ(std::vector<T>& outputVec, const std::string& inputFileName, const double XYZPARSE_E_UPPER_LIMIT, const double XYZPARSE_E_LOWER_LIMIT,
                      const double XYZPARSE_COORD_LIMIT, double* const globalmin = nullptr, std::vector<T>* const rejectVec = nullptr){
    // LOAD POINTS USED IN FITTING
    // -------------------------------------------------------------
    // the input file must have the PES points in the XYZ format
    const auto start = std::chrono::steady_clock::now();
    if (globalmin != nullptr) *globalmin = 0.0;
    uint32_t rejectCount = 0; //needed to count rejects if we are discarding them
    const int64_t infileSize = getFilesize(inputFileName);
    assert(infileSize > 10); //a valid XYZ file below 11 bytes is impossible
    std::ifstream infile;
    infile.open(inputFileName);
    if (!infile.is_open()){
        std::cout << "Error: cannot open XYZ-like input file: " << inputFileName << std::endl;
        return -1;
    }
    // parse the input into a vector of XYZ file structs
    std::cout << "Loading geometries from: " << inputFileName << std::endl;
    while(infile.tellg() < (infileSize-10)){ //this probably does NOT work under Windows
        bool reject;
        const T tmpxyz = parseXYZ<T>(infile, reject, XYZPARSE_E_UPPER_LIMIT, XYZPARSE_E_LOWER_LIMIT, XYZPARSE_COORD_LIMIT);
        if (reject == true){
            if (rejectVec != nullptr){ //if the caller wants to keep rejects, in a separate object
                rejectVec->push_back(tmpxyz);
            }else{
                rejectCount++;
            }
        }else{
            outputVec.push_back(tmpxyz);
        }
        if (globalmin != nullptr) { //if the caller requests the energy of the lowest energy point from the dataset
			if (outputVec.size() > 0){ //otherwise this would crash if all geoms so far were rejected
				if (outputVec.back().energy < *globalmin) *globalmin = outputVec.back().energy;
			}
        }
        //if ((outputVec.size() % PROGRESS_GRANULARITY) == 0) std::cout<<outputVec.size()<<" PES points parsed"<<std::endl;
    }
    std::cout<<outputVec.size()<<" geometries parsed and accepted in ";
    const auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration <double> (end - start).count() << " seconds" <<std::endl;
    if (rejectVec != nullptr){
		if (rejectVec->size() != 0){
			std::cout<<rejectVec->size()<<" geometries rejected, but saved"<<std::endl;
		}else{
			std::cout<<"No geometries were rejected.\n";
		}
    }else{
		if (rejectCount != 0){
			std::cout<<rejectCount<<" geometries rejected and discarded"<<std::endl;
		}else{
			std::cout<<"No geometries were rejected.\n";
		}
    }
    // make sure all points have the same number of atoms and the same order of atoms, as the first point
    std::cout << "Checking accepted geometries for inconsistencies" << std::endl;
    if (int32_t ret = CheckMultiXYZ_InternalConsistence(outputVec); ret < 0) return ret;
    std::cout<<"All accepted points checked"<<std::endl;
    if ((rejectVec != nullptr)&&(rejectVec->size() != 0)){
        std::cout << "Checking rejected geometries for inconsistencies" << std::endl;
        if (int32_t ret = CheckMultiXYZ_InternalConsistence(*rejectVec); ret < 0) return ret;
        std::cout<<"All rejected points checked"<<std::endl;
    }
    return 0;
}

template <typename T>
inline T LoadSingleXYZ(const std::string& inputFileName, const double XYZPARSE_E_UPPER_LIMIT, const double XYZPARSE_E_LOWER_LIMIT, const double XYZPARSE_COORD_LIMIT){
    std::cout << "Loading geometries from: " << inputFileName << std::endl;
    std::ifstream infile;
    infile.open(inputFileName);
    if (!infile.is_open()){
        std::cout << "Error: cannot open XYZ-like input file: " << inputFileName << std::endl;
        assert(infile.is_open());
    }
    bool reject;
    const T tmp = parseXYZ<T>(infile, reject, XYZPARSE_E_UPPER_LIMIT, XYZPARSE_E_LOWER_LIMIT, XYZPARSE_COORD_LIMIT);
    assert(reject == false);
    return tmp;
}

template <typename T>
inline void WriteXYZatoms(std::ofstream& outfile, const T& pesPoint){
    ///write the atom symbols and coordinates of an xyzfile struct to an already opened output file stream
    assert(outfile.is_open());
    for (uint32_t i=0; i<pesPoint.natoms; i++){
        outfile<<pesPoint.atomvec[i].symbol;
        if (pesPoint.atomvec[i].symbol[1]==0) outfile<<' ';
        outfile<<' '<<pesPoint.atomvec[i].xcoord<<' '<<pesPoint.atomvec[i].ycoord<<"  "<<pesPoint.atomvec[i].zcoord<<'\n';
    }
}

template <typename T, bool writeTraj = true>
inline void WriteXYZ(const T& pesPoint, std::ofstream& outfile, const uint32_t prec = 8, const int32_t dummyTraj = -1, const int32_t dummyStep = -1){
    ///write the contents of an xyzfile struct to an already opened output file stream
    outfile << std::setprecision(prec) << std::fixed;
    outfile << pesPoint.natoms;
    if constexpr (std::is_same<T,xyzfile>::value){
		if ((dummyTraj >= 0) && (dummyStep >= 0)) {outfile << "   #" << dummyTraj << "#" << dummyStep;}
    }else{
        if constexpr (writeTraj) {outfile << "   #" << pesPoint.trajNum << "#" << pesPoint.trajStep;}
    }
	outfile << '\n';
    outfile << pesPoint.extraline << '\n';
    WriteXYZatoms(outfile, pesPoint);
	outfile << std::flush; //not strictly necessary to do this, just to have a little bit smoother behavior
}

template <typename T, bool writeTraj = true>
void WriteMultiXYZ(const std::string& filename, const std::vector<T>& xyzvec, const bool append = false, const uint32_t prec = 8, const int32_t dummyTraj = -1, const int32_t dummyStep = -1){
	if (xyzvec.size() != 0){ //avoid creating empty files
		std::ofstream outfile;
		std::cout << "Opening file " << filename << " to write " << xyzvec.size() << " geometries" << std::endl;
		if (append == true){
			outfile.open(filename, std::ofstream::out | std::ofstream::app);
			assert(outfile.is_open());
		}else{
			outfile.open(filename, std::ofstream::out | std::ofstream::trunc);
			assert(outfile.is_open());
		}
		for(uint32_t i=0; i < xyzvec.size(); i++){
			WriteXYZ<T,writeTraj>(xyzvec[i], outfile, prec, dummyTraj, dummyStep);
		}
		assert(outfile.good());
	}else{
		std::cout << "No geometries to be written to " << filename << std::endl;
	}
}

template <typename T>
void WriteSingleXYZ(const std::string& filename, const T& xyz){
    std::ofstream outfile;
    std::cout << "Opening file " << filename << " for writing" << std::endl;
    outfile.open(filename);
    assert(outfile.is_open());
    WriteXYZ(xyz, outfile);
    assert(outfile.good());
}
