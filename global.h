#pragma once

/* a singleton class
 * to hold global states & parameters
**/
#include "config.h"

class suGlobal
{
public:
	suGlobal();
	~suGlobal();

	suString getParam(std::string para_name);
	bool load(suString conf_file);

	void clear();

	static suGlobal& gOnly();

private:
	static suGlobal *pOnly_;
	std::map<std::string, suString> paras_;
};

class globalValue{
public:
	int octreeDepth;
	int mcResolution;
	float metaballThreshold;
	std::string oofemPath;
	int iteratorTimes;
	float volumeFraction;
	float* loadRegion;
	static globalValue& globalValuePoint();
	int boundaryThickness;
	double maxStrain;
	void loadValuesFromsuGlobal();
	int currentIteratorTime;
	float density;
	float youngModules;
	float possionRatio;
	float force;
	float tAlpha;
	globalValue();
	~globalValue();
private:
	static globalValue *globalValues;
};