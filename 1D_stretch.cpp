#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cassert>
#include <quadmath.h>

#include "misc_util.hpp"
#include "storageClasses.hpp"
#include "XYZ_IO.hpp"


class VECT{
	public:
	__float128 x,y,z;
	VECT(){}
	constexpr VECT(const double x_in, const double y_in, const double z_in) : x(x_in), y(y_in), z(z_in){}
	VECT(const atom atom_in) : x(atom_in.xcoord), y(atom_in.ycoord), z(atom_in.zcoord){}
};

void TranslateOriginToCoords(xyzfile& geom, const VECT& coords){
	///translates the origin of the Cartesian coordinate system to a given point
	for(uint32_t i=0; i<geom.natoms; i++){
		geom.atomvec[i].xcoord -= static_cast<double>(coords.x);
		geom.atomvec[i].ycoord -= static_cast<double>(coords.y);
		geom.atomvec[i].zcoord -= static_cast<double>(coords.z);
	}
}

VECT Substract(const VECT& a, const VECT& b){
	VECT c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}

VECT Scale(const VECT& a, const __float128 factor){
	VECT c;
	c.x = factor * a.x;
	c.y = factor * a.y;
	c.z = factor * a.z;
	return c;
}

auto CalcNorm(const VECT& a){
	//if constexpr (std::is_same<decltype(a.x), 
	return sqrtq(a.x*a.x + a.y*a.y + a.z*a.z);
}

VECT Normalize(const VECT& a){
	const auto norm = CalcNorm(a);
	return Scale(a, static_cast<decltype(a.x)>(1)/norm);
}

int main(int argc, char* argv[]){
	std::cout<<std::boolalpha; //printing bool variables should print true/false instead of 1/0
	const std::vector<std::string> cliArgs = getCLIargs(argc, argv);
	std::cout<< "*------------------------------------------------------------------------------*\n";
	std::cout<< "1D_scan_linear was started with the following command line arguments:\n";
	for(uint32_t i=0;i<cliArgs.size();i++) std::cout<< cliArgs.at(i)<< '\n';
	//assert((cliArgs.size() == 8)||(cliArgs.size() == 9)||(cliArgs.size() == 10)||(cliArgs.size() == 11));
	std::cout<< "*------------------------------------------------------------------------------*\n";
	
	//parse the arguments
	if(cliArgs.size() < 1) std::cout<<"Missing arguments: input geometry\n";
	const std::string inputGeomPath = cliArgs.at(0);
	std::cout<< "Input geometry file: "<< inputGeomPath<< '\n';
	const xyzfile input_geom = LoadSingleXYZ<xyzfile>(inputGeomPath, 123456.7, -123456.7, 1.0E+6);
	
	if(cliArgs.size() < 2) std::cout<<"Missing arguments: first atom of the scan axis\n";
	if(cliArgs.size() < 3) std::cout<<"Missing arguments: second atom of the scan axis\n";
	const auto scanAxisAtom1_Idx = std::stoul(cliArgs.at(1));
	const auto scanAxisAtom2_Idx = std::stoul(cliArgs.at(2));
	std::cout<< "First atom specifying the axis of the 1D scan : "<< scanAxisAtom1_Idx<< '\n';
	std::cout<< "Second atom specifying the axis of the 1D scan: "<< scanAxisAtom2_Idx<< '\n';
	
	//if(cliArgs.size() < 1) std::cout<<"Missing arguments: minimum distance on the scan \n";
	const double minDist = std::stod(cliArgs.at(3));
	const double maxDist = std::stod(cliArgs.at(4));
	//std::cout<< "Minimum distance between the second atom of the scan axis and the first atom of the moving fragment: "<< minDist<< '\n';
	std::cout<< "Minimum displacement of the moving fragment along the scan axis: "<< minDist<< '\n';
	std::cout<< "Maximum displacement of the moving fragment along the scan axis: "<< maxDist<< '\n';
	//std::cout<< "Maximum distance between the second atom of the scan axis and the first atom of the moving fragment: "<< maxDist<< '\n';
	const auto numPoints = std::stoul(cliArgs.at(5));
	std::cout<< "Number of points in the scan: "<< numPoints<< '\n';
	assert(numPoints >= 2);
	const double stepLen = (maxDist-minDist)/static_cast<double>(numPoints-1);
	std::cout<< "Calculated step length:                  "<< stepLen<< '\n';
	
	const bool monoatomicMovingFrag = [&]{
		if (cliArgs.at(6)=="Mono"){ return true; }else{ return false;}
	}();
	const bool diatomicMovingFrag = [&]{
		if (cliArgs.at(6)=="Di"){ return true; }else{ return false;}
	}();
	const bool triatomicMovingFrag = [&]{
		if (cliArgs.at(6)=="Tri"){ return true; }else{ return false;}
	}();
	std::cout<< "Is the moving fragment monoatomic? : "<< monoatomicMovingFrag<< '\n';
	std::cout<< "Is the moving fragment diatomic?   : "<< diatomicMovingFrag<< '\n';
	std::cout<< "Is the moving fragment triatomic?  : "<< triatomicMovingFrag<< '\n';
	assert( (cliArgs.at(6)=="Mono") || (cliArgs.at(6)=="Di") || (cliArgs.at(6)=="Tri") ); //type of moving fragment: monoatomic, diatomic, triatomic
	if(triatomicMovingFrag){ ///TODO: implement triatomic
		std::cout<<"Triatomic moving frgments are not implemented yet\n";
		abort();
	}
	
	const auto movingFragIdxArr = [&]{
		std::vector<uint64_t> arr;
		if (monoatomicMovingFrag){
			arr = { std::stoul(cliArgs.at(7)) };
			return arr;
		}
		if (diatomicMovingFrag){
			arr = { std::stoul(cliArgs.at(7)), std::stoul(cliArgs.at(8)) };
			return arr;
		}else{
			arr = { std::stoul(cliArgs.at(7)), std::stoul(cliArgs.at(8)), std::stoul(cliArgs.at(9)) };
			return arr;
		}
	}();
	std::cout<< "Atoms of the moving fragment:        ";
	for(const auto& i : movingFragIdxArr) std::cout<< i<< ", ";
	std::cout<< '\n';
	
	const bool alignFrag = [&]{
		uint32_t arg = 8;
		if(diatomicMovingFrag) arg += 1;
		if(triatomicMovingFrag) arg += 2;
		if(!( (cliArgs.at(arg)=="Align") || (cliArgs.at(arg)=="Noalign") ) ){
			std::cout<<"Error: fragment alignment must be either Align or Noalign\n";
			abort();
		}
		if (cliArgs.at(arg)=="Align"){ return true; }else{ return false;}
	}();
	std::cout<< "Should a fragment be aligned onto an alignment axis? : "<< alignFrag<< '\n';
	
	const auto alignAxisAtom1_Idx = [&]{
		if(!alignFrag){
			return static_cast<unsigned long int>(0);
		}else{
			uint32_t arg = 9;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			return std::stoul(cliArgs.at(arg));
		}
	}();
	const auto alignAxisAtom2_Idx = [&]{
		if(!alignFrag){
			return static_cast<unsigned long int>(0);
		}else{
			uint32_t arg = 10;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			return std::stoul(cliArgs.at(arg));
		}
	}();
	const auto alignDist = [&]{
		if(!alignFrag){
			return 0.0;
		}else{
			uint32_t arg = 11;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			return std::stod(cliArgs.at(arg));
		}
	}();
	//find out the number of atoms in the aligned fragment
	const bool monoatomicAlignedFrag = [&]{
		if(!alignFrag){
			return false;
		}else{
			uint32_t arg = 12;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			if (cliArgs.at(arg)=="MonoAlign"){ return true; }else{ return false;}
		}
	}();
	const bool diatomicAlignedFrag = [&]{
		if(!alignFrag){
			return false;
		}else{
			uint32_t arg = 12;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			if (cliArgs.at(arg)=="DiAlign"){ return true; }else{ return false;}
		}
	}();
	const bool triatomicAlignedFrag = [&]{
		if(!alignFrag){
			return false;
		}else{
			uint32_t arg = 12;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			if (cliArgs.at(arg)=="TriAlign"){ return true; }else{ return false;}
		}
	}();
	if(alignFrag){
		std::cout<< "Is the aligned fragment monoatomic? : "<< monoatomicAlignedFrag<< '\n';
		std::cout<< "Is the aligned fragment diatomic?   : "<< diatomicAlignedFrag<< '\n';
		std::cout<< "Is the aligned fragment triatomic?  : "<< triatomicAlignedFrag<< '\n';
		uint32_t arg = 12;
		if(diatomicMovingFrag) arg += 1;
		if(triatomicMovingFrag) arg += 2;
		assert( (cliArgs.at(arg)=="MonoAlign") || (cliArgs.at(arg)=="DiAlign") || (cliArgs.at(arg)=="TriAlign") ); //type of aligned fragment: monoatomic, diatomic, triatomic
		if(triatomicAlignedFrag){ ///TODO: implement triatomic
			std::cout<<"Triatomic aligned frgments are not implemented yet\n";
			abort();
		}
	}
	
	const auto alignedFragIdxArr = [&]{
		std::vector<uint64_t> arr;
		if(!alignFrag){
			return arr;
		}else{
			uint32_t arg = 13;
			if(diatomicMovingFrag) arg += 1;
			if(triatomicMovingFrag) arg += 2;
			if (monoatomicAlignedFrag){
				arr = { std::stoul(cliArgs.at(arg)) };
				return arr;
			}
			if (diatomicAlignedFrag){
				arr = { std::stoul(cliArgs.at(arg)), std::stoul(cliArgs.at(arg+1)) };
				return arr;
			}else{
				arr = { std::stoul(cliArgs.at(arg)), std::stoul(cliArgs.at(arg+1)), std::stoul(cliArgs.at(arg+2)) };
				return arr;
			}
		}
	}();
	std::cout<< "Atoms of the aligned fragment:        ";
	for(const auto& i : alignedFragIdxArr) std::cout<< i<< ", ";
	std::cout<< '\n';
	
	
	// const auto ref_geom = [&]{
		// ///TODO: implement aligning to an axis perpendicular to the alignment axis
		// ///TODO: implement aligning to an axis specified by three atoms (middle of the angle between them)
		// ///TODO: implement aligning to an axis specified by four atoms (axis specified by the fourth atom and the centroid of the triangle created by the first three)
		// ///TODO: implement aligning with a different distance
		// ///TODO: implement triatomic moving fragment alignment
        // auto tmp = input_geom;
		// if(alignMovingFrag){
			// const double alignDist = 0;
			// std::cout<< "First atom specifying the axis of the 1D scan:  "<< scanAxisAtom1_Idx<< '\n';
			// std::cout<< "Second atom specifying the axis of the 1D scan: "<< scanAxisAtom2_Idx<< '\n';
			// std::cout<< "Alignment distance between the second atom of the alignment axis and the first atom of the moving fragment: "<< alignDist<< '\n';
			
			// TranslateOriginToCoords(tmp, tmp.atomvec.at(alignAxisAtom2_Idx) );
			// const auto axis = Normalize( Substract(tmp.atomvec.at(alignAxisAtom2_Idx), tmp.atomvec.at(alignAxisAtom1_Idx)) );
			// tmp.atomvec.at(movingFragIdxArr.at(0)).xcoord = alignDist*static_cast<double>(axis.x);
			// tmp.atomvec.at(movingFragIdxArr.at(0)).ycoord = alignDist*static_cast<double>(axis.y);
			// tmp.atomvec.at(movingFragIdxArr.at(0)).zcoord = alignDist*static_cast<double>(axis.z);
			// if(!monoatomicMovingFrag){
				// const double bondLength_1 = static_cast<double>(CalcNorm( Substract(input_geom.atomvec.at(movingFragIdxArr.at(0)), input_geom.atomvec.at(movingFragIdxArr.at(1))) ));
				// tmp.atomvec.at(movingFragIdxArr.at(1)).xcoord = (alignDist + bondLength_1)*static_cast<double>(axis.x);
				// tmp.atomvec.at(movingFragIdxArr.at(1)).ycoord = (alignDist + bondLength_1)*static_cast<double>(axis.y);
				// tmp.atomvec.at(movingFragIdxArr.at(1)).zcoord = (alignDist + bondLength_1)*static_cast<double>(axis.z);
			// }
			// ///todo triatomic
		// }
        // TranslateOriginToCoords(tmp, tmp.atomvec.at(scanAxisAtom2_Idx) );
        // return tmp;
	// }();
	// std::cout<< "*------------------------------------------------------------------------------*\n";
	
	std::vector<std::remove_cv<decltype(input_geom)>::type> scanGeomVec;
	scanGeomVec.resize(numPoints);
	const auto scanAxis = Normalize( Substract(input_geom.atomvec.at(scanAxisAtom2_Idx), input_geom.atomvec.at(scanAxisAtom1_Idx)) );
	
	const auto ref_geom = [&]{
		auto tmp = input_geom;
		if(alignFrag){
			std::cout<< "First atom specifying the axis of the alignment : "<< alignAxisAtom1_Idx<< '\n';
			std::cout<< "Second atom specifying the axis of the alignment: "<< alignAxisAtom2_Idx<< '\n';
			std::cout<< "Alignment distance between the second atom of the alignment axis and the first atom of the aligned fragment: "<< alignDist<< '\n';
			
			TranslateOriginToCoords(tmp, tmp.atomvec.at(alignAxisAtom2_Idx) );
			const auto axis = Normalize( Substract(tmp.atomvec.at(alignAxisAtom2_Idx), tmp.atomvec.at(alignAxisAtom1_Idx)) );
			tmp.atomvec.at(alignedFragIdxArr.at(0)).xcoord = alignDist*static_cast<double>(axis.x);
			tmp.atomvec.at(alignedFragIdxArr.at(0)).ycoord = alignDist*static_cast<double>(axis.y);
			tmp.atomvec.at(alignedFragIdxArr.at(0)).zcoord = alignDist*static_cast<double>(axis.z);
			if(!monoatomicAlignedFrag){
				const double bondLength_1 = static_cast<double>(CalcNorm( Substract(input_geom.atomvec.at(alignedFragIdxArr.at(0)), input_geom.atomvec.at(alignedFragIdxArr.at(1))) ));
				tmp.atomvec.at(alignedFragIdxArr.at(1)).xcoord = (alignDist + bondLength_1)*static_cast<double>(axis.x);
				tmp.atomvec.at(alignedFragIdxArr.at(1)).ycoord = (alignDist + bondLength_1)*static_cast<double>(axis.y);
				tmp.atomvec.at(alignedFragIdxArr.at(1)).zcoord = (alignDist + bondLength_1)*static_cast<double>(axis.z);
			}
			///todo triatomic
		}
		TranslateOriginToCoords(tmp, tmp.atomvec.at(scanAxisAtom1_Idx) );
		if( monoatomicMovingFrag && (scanAxisAtom2_Idx == movingFragIdxArr.at(0)) ){
			tmp.atomvec.at(movingFragIdxArr.at(0)).xcoord = 0;
			tmp.atomvec.at(movingFragIdxArr.at(0)).ycoord = 0;
			tmp.atomvec.at(movingFragIdxArr.at(0)).zcoord = 0;
		}
		if( diatomicMovingFrag && (scanAxisAtom2_Idx == movingFragIdxArr.at(0)) ){
			tmp.atomvec.at(movingFragIdxArr.at(1)).xcoord -= tmp.atomvec.at(movingFragIdxArr.at(0)).xcoord;
			tmp.atomvec.at(movingFragIdxArr.at(1)).ycoord -= tmp.atomvec.at(movingFragIdxArr.at(0)).ycoord;
			tmp.atomvec.at(movingFragIdxArr.at(1)).zcoord -= tmp.atomvec.at(movingFragIdxArr.at(0)).zcoord;
			tmp.atomvec.at(movingFragIdxArr.at(0)).xcoord = 0;
			tmp.atomvec.at(movingFragIdxArr.at(0)).ycoord = 0;
			tmp.atomvec.at(movingFragIdxArr.at(0)).zcoord = 0;
		}
        return tmp;
	}();
	
	
	
	std::ofstream dist_out;
	dist_out.open("scan_dists.txt");
	assert(dist_out.is_open());
	for (uint32_t i=0; i<numPoints; i++){
		const auto disp = Scale(scanAxis, i*stepLen + minDist);
		auto tmp = ref_geom;
		tmp.atomvec.at(movingFragIdxArr.at(0)).xcoord += static_cast<double>(disp.x);
		tmp.atomvec.at(movingFragIdxArr.at(0)).ycoord += static_cast<double>(disp.y);
		tmp.atomvec.at(movingFragIdxArr.at(0)).zcoord += static_cast<double>(disp.z);
		if(!monoatomicMovingFrag){
			tmp.atomvec.at(movingFragIdxArr.at(1)).xcoord += static_cast<double>(disp.x);
			tmp.atomvec.at(movingFragIdxArr.at(1)).ycoord += static_cast<double>(disp.y);
			tmp.atomvec.at(movingFragIdxArr.at(1)).zcoord += static_cast<double>(disp.z);
		}
		if(triatomicMovingFrag){
			tmp.atomvec.at(movingFragIdxArr.at(2)).xcoord += static_cast<double>(disp.x);
			tmp.atomvec.at(movingFragIdxArr.at(2)).ycoord += static_cast<double>(disp.y);
			tmp.atomvec.at(movingFragIdxArr.at(2)).zcoord += static_cast<double>(disp.z);
		}
		scanGeomVec.at(i) = tmp;
		dist_out<< i*stepLen + minDist<< '\n';
		if( const auto actualDist = static_cast<double>( CalcNorm( Substract(tmp.atomvec.at(movingFragIdxArr.at(0)), ref_geom.atomvec.at(scanAxisAtom2_Idx)) ) ); std::abs( actualDist - (i*stepLen + minDist) ) > 1E-10){
			std::cout<<"Distance between atom #"<<movingFragIdxArr.at(0)<<" and atom #"<<scanAxisAtom2_Idx<<" is: "<< actualDist<< " , should be : "<< i*stepLen + minDist<< std::endl;
		}
	}
	assert(dist_out.good());
	WriteMultiXYZ("scan.xyz", scanGeomVec);
	return 0;
	/*
	
	
	//echo the arguments after parsing
	


	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	std::cout<< ": "<< inputGeomPath<< '\n';
	xyzfile ref_geom = LoadSingleXYZ<xyzfile>(REF_GEOM_PATH, 123456.7, -123456.7, 1.0E+6);
	TranslateOriginToCoords(ref_geom, CalcCOM(ref_geom));
	
	
	
	
	std::array<uint32_t,2> frag2_idx = {0, 5};
	constexpr uint32_t NATOMS = 8;
	constexpr __float128 APPROACH_STEP = 0.1;
	constexpr uint32_t ROTATIONS_PER_APPROACH_STEP = 10000;
	constexpr __float128 APPROACH_MIN_COMDIST = 1.0;
	constexpr double PRUNING_THRESH = 0.025;
	const std::string REF_GEOM_PATH = "ref_geom_MeOH_OHrad.xyz";
	
	xyzfile ref_geom = LoadSingleXYZ<xyzfile>(REF_GEOM_PATH, 123456.7, -123456.7, 1.0E+6);
	assert(ref_geom.natoms == NATOMS);
	{
		std::vector<xyzfile> rotatedGeomVec;
		auto PRNG = initPRNG();
		std::uniform_real_distribution<double> URD(0.0, 1.0);
		std::uniform_real_distribution<double> pertDist(-0.1, 0.1);
		while(true){
			for(uint32_t i=0; i<ROTATIONS_PER_APPROACH_STEP; i++){
				xyzfile rotated_geom = ref_geom;
				randPertGeom<NATOMS>(rotated_geom, PRNG, pertDist);
				{ //rotate the smaller reactant
					const ROTMAT R = GenRandRotmat(PRNG, URD);
					TranslateOriginToCoords(rotated_geom, CalcCOM<'S'>(rotated_geom, &frag2_idx));
					for(uint32_t j=0; j<ref_geom.natoms; j++){
						if (std::any_of(frag2_idx.cbegin(), frag2_idx.cend(), [&](uint32_t idx){ return idx == j;})){
							RotateAtom(rotated_geom.atomvec[j], R);
						}
					}
				}
				{ //rotate the larger reactant
					const ROTMAT R = GenRandRotmat(PRNG, URD);
					TranslateOriginToCoords(rotated_geom, CalcCOM<'F'>(rotated_geom, &frag2_idx));
					for(uint32_t j=0; j<ref_geom.natoms; j++){
						if (std::none_of(frag2_idx.cbegin(), frag2_idx.cend(), [&](uint32_t idx){ return idx == j;})){
							RotateAtom(rotated_geom.atomvec[j], R);
						}
					}
				}
				TranslateOriginToCoords(rotated_geom, CalcCOM(rotated_geom));
				rotatedGeomVec.push_back(rotated_geom);
			}
			
			TranslateOriginToCoords(ref_geom, CalcCOM<'F'>(ref_geom, &frag2_idx));
			const VECT frag2vec = CalcCOM<'S'>(ref_geom, &frag2_idx);
			if (const __float128 norm = CalcNorm(frag2vec); norm < APPROACH_MIN_COMDIST){
				break;
			}else{
				std::cout<<static_cast<double>(norm)<<std::endl;
				TranslateFrag<'S'>(ref_geom, Scale(frag2vec, -APPROACH_STEP/norm), &frag2_idx);
			}
		};
		WriteMultiXYZ("approach.xyz",rotatedGeomVec);
	}
	
	std::vector<xyzfile> pruned_geoms;
	{ //prune the generated geometries
		multiXYZ gen_geoms;
		assert(0==LoadMultiXYZ(gen_geoms, "approach.xyz", 123456.7, -123456.7, 1.0E+6));
		const auto PruningStart = std::chrono::steady_clock::now();
		PruneGems(gen_geoms, PRUNING_THRESH);
		const auto PruningFinish = std::chrono::steady_clock::now();
		std::cout << "Done in "<< std::chrono::duration <double> (PruningFinish - PruningStart).count() << " seconds" << std::endl;
		assert(CheckMultiXYZ_InternalConsistence(gen_geoms) == 0);
		gen_geoms.isValidated = true;
		std::cout<<gen_geoms.size()<<" geometries remain after pruning"<<std::endl<<std::endl;
		
		pruned_geoms.reserve(gen_geoms.size());
		for(uint32_t i=0; i<gen_geoms.size(); i++) pruned_geoms.push_back(gen_geoms.getXYZ(i));
	}
	WriteMultiXYZ("approach_pruned.xyz",pruned_geoms);
	
	return 0;*/
}