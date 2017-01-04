#define USE_EIGEN

#include "config.h"
#include "suMesh.h"
#include "suVolume.h"
#include "suWriteSTL.h"
#include "suStructureEvolution.h"
#include <time.h>
#include <iostream>
#include "global.h"
//global parameters
suString inputFile, outputFile;
//int nOctreeDeep = 3;
//int nResolution = 3;         //resolution of MC reconstruct grid  (3-10)
//float  fThreshold = 2;          //threshold of metaball in MC reconstruct.
//int iteratorTimes = 1;
void helper()
{
	std::cout << "Evolve is a mesh structure optimizer. \n";
	std::cout << "Usage: \n    Evolve -i [input stl/obj/ply] -o [output stl/obj/ply]" << std::endl;
	std::cout << "Other parameter include -t -r -d" << std::endl;
	std::cout << std::endl;
	std::cout << "-t Specifythreshold of metaball in MC reconstruct [range: 0.5-3.0]." << std::endl;
	std::cout << "-r Specify resolution of MC reconstruct grid  [range: 3-10]." << std::endl;
	std::cout << "-d Specify octree deep [range:1-6]" << std::endl;
	std::cout << "-a Specify iterator times" << std::endl;
	std::cout << "-s Specify a config file. If this is specified, other parameters are ignored.";
}

void main(int argc, char* argv[])
{
	suCMDParser option(argv, argc);


	if (option.findOption("-h") || argc == 1)
	{
		helper();
		return;
	}

	////////////////////////////////////////// get parameters
	try {
		//get input file
		if (option.findOption("-i"))
		{
			inputFile = option.findOptionValue("-i");
			if (option.findOption("-o"))
			{
				outputFile = option.findOptionValue("-o");
			}
			else {
				throw std::exception("No output file!");
			}

		}
		else {
			throw std::exception("No input file!");
		}
	}
	catch (std::exception &e)
	{
		Console::WriteLine(suString(e.what()));
		helper();
		return;
	}

	int dValue = 0;
	int rValue = 0;
	float fValue = 0.0f;
	int aValue = 0;
	// find parameters
	try {
		if (option.findOption("-s"))
		{
			//from config file
			suString configFile = option.findOptionValue("-s");

			//add a class to parse and hold parameters
			suGlobal::gOnly().load(configFile);
			globalValue::globalValuePoint().loadValuesFromsuGlobal();
		}
		else
		{
			//from command line
			if (option.findOption("-d"))
			{
				dValue = option.findOptionValue("-d").ConvertToInteger();
			}
			if (option.findOption("-t"))
			{
				fValue = option.findOptionValue("-t").ConvertToFloat();
			}
			if (option.findOption("-r"))
			{
				rValue = option.findOptionValue("-r").ConvertToInteger();
			}
			if (option.findOption("-a"))
			{
				aValue = option.findOptionValue("-a").ConvertToInteger();
			}
			if (!dValue) globalValue::globalValuePoint().octreeDepth = dValue;
			if (!rValue) globalValue::globalValuePoint().mcResolution = rValue;
			if (fValue != 0) globalValue::globalValuePoint().metaballThreshold = fValue;
			if (!aValue) globalValue::globalValuePoint().iteratorTimes = aValue;
		}
		
	}
	catch (std::exception &e)
	{
		Console::WriteLine(e.what());
		Console::WriteLine("Using default parameters(d = 3, t = 2, r = 3)");
	}





	////////////////////////////////////////////////Optimization
	//read file
	Console::WriteLine("Reading " + inputFile);
	Console::AdvanceLine();

	SU::suVolume v;
	v.LoadMeshFromFile(inputFile.CString());

	//voxization
	time_t start, stop;
	start = time(NULL);


	v.PartitionSpace(globalValue::globalValuePoint().octreeDepth);
	stop = time(NULL);
	std::cout << "Use Time:" << (stop - start);
	std::vector<SU::OctNode*> newLeafNodes;   //store all leaf interior nodes
	std::vector<SU::OctNode*> &InterNodeArr = v.leafInternalNodes_;

	std::vector<SU::OctNode*>::iterator it = InterNodeArr.begin();
	for (; it != InterNodeArr.end(); it++)
	{
		if ((*it)->level_  < v.level_)
		{
			v.patitionToLevel(*it, v.level_, (*it)->label_, newLeafNodes);	 //add new node	pointer			
		}
		else {
			newLeafNodes.push_back(*it);  //copy original nodes pointer
		}
	}
	//Now we get a new leaf nodes array(boundary + interior)
	InterNodeArr = newLeafNodes;   //update v.leafInternalNodes_
								   //newLeafNodes.insert(newLeafNodes.end(), v.leafBoundaryNodes_.begin(), v.leafBoundaryNodes_.end());


	Console::WriteLine("Num of interior nodes: " + suString(v.leafInternalNodes_.size()));
	Console::AdvanceLine();


	//Run evolution
	suStructrueOptimizer optimizer;

	optimizer.set_iter_times(globalValue::globalValuePoint().iteratorTimes);
	optimizer.set_octree_level(globalValue::globalValuePoint().octreeDepth);

	optimizer.init(&v);

	optimizer.run();



	//Reconstruction, bool combination and save	
	std::vector<Eigen::Vector3f> bbox;

	Eigen::Vector3f minP;  minP << v.bbMin_[0], v.bbMin_[1], v.bbMin_[2]; bbox.push_back(minP);
	Eigen::Vector3f maxP;  maxP << v.bbMax_[0], v.bbMax_[1], v.bbMax_[2]; bbox.push_back(maxP);

	int nVolume_per_dim = pow(2, globalValue::globalValuePoint().octreeDepth);
	float volume_size = (maxP - minP).maxCoeff() / nVolume_per_dim;

	//generate metaballs on the model surface.
	std::vector<SU::METABALL> mballs;
	for (unsigned int i = 0; i <InterNodeArr.size(); i++)
	{
		SU::METABALL m;

		if (!(InterNodeArr[i]->out)
			&& InterNodeArr[i]->label_ == SU::INTERIOR_CELL)
		{
			SU::Point cP = InterNodeArr[i]->center();
			std::cout << "position = " << cP.x << ", " << cP.y << ", " << cP.z << std::endl;
			m.position << cP.x, cP.y, cP.z;
			m.squaredRadius = pow(volume_size / 2, 2);     //  metaball.r^2 = (volume box size / 2)^2													 

			mballs.push_back(m);
		}
	}

	Console::WriteLine("Output model in " + outputFile);
	int resolution = nVolume_per_dim * globalValue::globalValuePoint().mcResolution;

	//SU::write_metaball_to_stl(outputFile.CString(), /*samples*/mballs, fThreshold, resolution, bbox);

	v.saveVTK(outputFile.CString(), globalValue::globalValuePoint().octreeDepth);
	system("pause");
}