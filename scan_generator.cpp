#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

inline int CallBash_quiet(const std::string& cmd){
	const int ret = std::system(cmd.c_str());
	#pragma omp critical
	if (ret!=0) std::cout << "Return value of " << cmd << " : " << ret << std::endl;
	return ret;
}

std::string LoadString(const std::string& fname){
	std::ifstream input;
	input.open(fname);
	assert(input.is_open());
	std::string tmp;
	std::getline(input, tmp);
	return tmp;
}

auto CreateDir(const std::string& dirname){
	return CallBash_quiet("mkdir "+dirname);
}

int main(){
	const auto refgeomPath = LoadString("ref_geom_path.in");
	std::ifstream scanspecFile;
	scanspecFile.open("scan_specifications.in");
	assert(scanspecFile.is_open());
	do{
		std::string scanName, movingFragType, alignment, alignedFragType;
		uint32_t numScanParts = 0, scanAxis[2]={1234,1234}, movingAtoms[3]={1234,1234,1234}, alignAxis[2]={1234,1234}, alignedAtoms[3]={1234,1234,1234};
		double alignDist = 0.0;
		scanspecFile >> scanName;
		scanspecFile >> numScanParts;
		scanspecFile >> scanAxis[0];
		scanspecFile >> scanAxis[1];
		scanspecFile >> movingFragType;
		if( !((movingFragType=="Mono") || (movingFragType=="Di") || (movingFragType=="Tri")) ){ //type of moving fragment: monoatomic, diatomic, triatomic
			std::cout<<"Error: the type of the moving fragment must be Mono or Di or Tri"<<std::endl;
			std::cout<<"Please check the specification of "<<scanName<<std::endl;
			abort();
		}
		scanspecFile >> movingAtoms[0];
		if((movingFragType=="Di") || (movingFragType=="Tri")) scanspecFile >> movingAtoms[1];
		if(movingFragType=="Tri") scanspecFile >> movingAtoms[2];
		scanspecFile >> alignment;
		if(!( (alignment=="Align") || (alignment=="Noalign") ) ){
			std::cout<<"Error: fragment alignment must be either Align or Noalign"<<std::endl;
			abort();
		}
		if(alignment=="Align"){
			scanspecFile >> alignAxis[0];
			scanspecFile >> alignAxis[1];
			scanspecFile >> alignDist;
			scanspecFile >> alignedFragType;
			if( !((alignedFragType=="MonoAlign") || (alignedFragType=="DiAlign") || (alignedFragType=="TriAlign")) ){ //type of aligned fragment: monoatomic, diatomic, triatomic
				std::cout<<"Error: the type of the aligned fragment must be MonoAlign or DiAlign or TriAlign"<<std::endl;
				std::cout<<"Please check the specification of "<<scanName<<std::endl;
				abort();
			}
			scanspecFile >> alignedAtoms[0];
			if((alignedFragType=="DiAlign") || (alignedFragType=="TriAlign")) scanspecFile >> alignedAtoms[1];
			if(alignedFragType=="TriAlign") scanspecFile >> alignedAtoms[2];
		}
		assert(0==CreateDir(scanName));
		std::ofstream scanScriptFile;
		scanScriptFile.open(scanName+'/'+scanName+".sh");
		assert(scanScriptFile.is_open());
		//scanScriptFile<<"#!/bin/sh\n";
		for(uint32_t i=0; i<numScanParts; i++){
			double scanpartStart = 0.0, scanpartEnd = 0.0;
			uint32_t scanpartPoints = 0;
			scanspecFile >> scanpartStart;
			scanspecFile >> scanpartEnd;
			scanspecFile >> scanpartPoints;
			scanScriptFile<< "./a.out "<< refgeomPath<<' '<< scanAxis[0]<<' '<< scanAxis[1]<<' '<< scanpartStart<<' '<< scanpartEnd<<' '<< scanpartPoints<<' '<< movingFragType<<' '<< movingAtoms[0]<<' ';
			if((movingFragType=="Di") || (movingFragType=="Tri")) scanScriptFile<< movingAtoms[1]<<' ';
			if(movingFragType=="Tri") scanScriptFile<< movingAtoms[2]<<' ';
			scanScriptFile<< alignment;
			if(alignment=="Align"){
				scanScriptFile<<' '<< alignAxis[0]<<' '<< alignAxis[1]<<' '<< alignDist<<' '<< alignedFragType<<' '<< alignedAtoms[0];
				if((alignedFragType=="DiAlign") || (alignedFragType=="TriAlign")) scanScriptFile <<' '<< alignedAtoms[1];
				if(alignedFragType=="TriAlign") scanScriptFile <<' '<< alignedAtoms[2];
			}
			scanScriptFile<< std::endl;
			if(i==0){
				scanScriptFile<< "cat scan.xyz > "<< scanName<< ".xyz\n";
				scanScriptFile<< "cat scan_dists.txt > "<< scanName<< "_dists.txt\n";
			}else{
				scanScriptFile<< "cat scan.xyz >> "<< scanName<< ".xyz\n";
				scanScriptFile<< "cat scan_dists.txt >> "<< scanName<< "_dists.txt\n";
			}
		}
		scanScriptFile<< "rm scan.xyz scan_dists.txt"<<std::endl;
		assert(scanScriptFile.good());
		scanScriptFile.close();
		//the scan generation script is now done, run it to generate the actual scan points
		assert(0==CallBash_quiet("chmod +x "+scanName+'/'+scanName+".sh"));
		CallBash_quiet("cp -p ./* "+scanName+'/');
		assert(0==CallBash_quiet("cd "+scanName+"; ./"+scanName+".sh > "+scanName+"_gen.log"));
		assert(0==CallBash_quiet("truncate --size=-1 "+scanName+'/'+scanName+".xyz"));
		//copy the pesrecalc template into the scan directory
		assert(0==CallBash_quiet("cp -rp pesrecalc_template/* "+scanName+'/'));
	} while(!scanspecFile.eof());
	
	return 0;
}