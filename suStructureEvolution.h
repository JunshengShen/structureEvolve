#pragma once
/*
    A volume structure evolution class.
	Morton code is used to index all voxels	
	*******
	\brief suStructrueOptimizer
	\todo
	\example   
	```
	suStructrueOptimizer  framework;
	framework.set_octree_level(5);
	framework.set_rule_lists(...);
	framework.init(internal_nodes);
	framework.run();
	```
   
**/
#include "suMesh.h"
#include "common.h"

//forecast defination
namespace SU
{
	struct OctNode;
	class suVolume;
};

struct point_strain
{
	double x_strain, y_strain, z_strain;
};
struct point_mises_strain
{
	float mises_strain;
	int point_morton;
};
//evolution framework
class suStructrueOptimizer
{
public:
	suStructrueOptimizer():nLevel_(0), nIterTimes_(999999999), pVolume_(0), fThresholdMC(1){};
	~suStructrueOptimizer() {};

	void set_octree_level(int level) { nLevel_ = level; }
	void set_iter_times(int nTimes) { nIterTimes_ = nTimes; }	
		
	int evolve();     //return num of iteration
	int init_morton_coding();  //generate a node array in nodeArr_, indexed with morton code
	int init_position_field(); 
	int init(SU::suVolume *pVolume);
	int run();
	void outOofemfile(std::string);
	void outForcedOofemFile(std::string);
	void useOofem(std::string);
	void readFeedback(std::string);
	static void outCoor(float, float, float, int, std::string coorAddress = "r:\\oofemOutFile.txt");
	void outForcedCoor(float, float, float, int, std::string coorAddress = "r:\\oofemOutFile.txt");
	//knowlege base
	void readRule() {};	
	void read_point_information(std::string address);
	void read_point_coor(std::string, float, float, float, int);
	bool export_stl_with_metaball(const char* fileName, std::vector<SU::OctNode*>& stlVector);
	void assignment(std::vector<SU::OctNode*>&);
	float return_max_strain(int, int);
	std::vector<int> forcedPoint;
private:
	SU::suVolume  *pVolume_;
	std::vector<SU::OctNode*> nodeArr_;  //morton coded node array 
	std::vector<SU::OctNode*> boundaryNodeArr_;            
	int nLevel_;     //octree level
	int nIterTimes_; //interation times      
	std::vector<float> all_point_mises_strain;
	std::vector<int> point_morton;
	void ass_point();
	int fThresholdMC ;
	void repair();
public:
	std::vector<point_mises_strain> sort_vector;
};


double trans(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8, char a9, char a10, char a11, char a12);


int asc2(char);