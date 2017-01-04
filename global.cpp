#include "global.h"

suGlobal* suGlobal::pOnly_ = 0;
globalValue* globalValue::globalValues = 0;
/*
* tools functions
**/
bool isNote(suString strLine)
{
	int idx = strLine.FindFirstIndex('#');
	if (idx != -1)
	{
		for (int i = 0; i < idx; i++)
		{
			if (strLine[i] != ' ')
				return false;
		}
		return true;
	}
	return false;

}

/* class suGlobal
* load parameters and manage them.
**/
suGlobal::suGlobal()
{

}

suGlobal::~suGlobal()
{
	clear();
}

suString suGlobal::getParam(std::string para_name)
{
	std::string strName = suString(para_name.c_str()).MakeLowercase().CString();
	const std::map<std::string, suString>::iterator it = paras_.find(strName);
	if (it != paras_.end())
	{
		return it->second;
	}

	return suString("");
}

bool suGlobal::load(suString conf_file)
{
	suVector<suString> _strParams;
	Utility::GetFileLines(conf_file, _strParams);

	paras_.clear();
	for (unsigned int i = 0; i < _strParams.Length(); i++)
	{
		if (!isNote(_strParams[i]))
		{
			suVector<suString> strKeys;
			_strParams[i].Partition('=', strKeys);
			if (strKeys.Length() == 2)
			{
				paras_[std::string(strKeys[0].MakeLowercase().CString())] = strKeys[1];
			}
		}
	}

	//dump params
	/*
	std::map<std::string, suString>::iterator it = paras_.begin();
	for (; it != paras_.end(); it++)
	std::cout << it->first.c_str() << " = " << it->second.CString() << std::endl;
	*/
	return true;
}

void suGlobal::clear()
{
	paras_.clear();
}

suGlobal & suGlobal::gOnly()
{
	if (!pOnly_)
	{
		pOnly_ = new suGlobal;
	}
	return *pOnly_;
}

globalValue & globalValue::globalValuePoint()
{
	if (!globalValues)
		globalValues = new globalValue;
	return *globalValues;
	// TODO: insert return statement here
}

void globalValue::loadValuesFromsuGlobal()
{
	globalValue::globalValuePoint().octreeDepth = suGlobal::gOnly().getParam("Octree Depth").ConvertToInteger();
	globalValue::globalValuePoint().mcResolution = suGlobal::gOnly().getParam("MC Resolution").ConvertToInteger();
	globalValue::globalValuePoint().metaballThreshold = suGlobal::gOnly().getParam("Metaball Threshold").ConvertToFloat();
	globalValue::globalValuePoint().oofemPath = suGlobal::gOnly().getParam("OOFEM Path").CString();
	globalValue::globalValuePoint().iteratorTimes = suGlobal::gOnly().getParam("Iteration Times").ConvertToInteger();
	globalValue::globalValuePoint().volumeFraction = suGlobal::gOnly().getParam("Volume Fraction").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[0] = suGlobal::gOnly().getParam("region x_min").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[1] = suGlobal::gOnly().getParam("region x_max").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[2] = suGlobal::gOnly().getParam("region y_min").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[3] = suGlobal::gOnly().getParam("region y_max").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[4] = suGlobal::gOnly().getParam("region z_min").ConvertToFloat();
	globalValue::globalValuePoint().loadRegion[5] = suGlobal::gOnly().getParam("region z_max").ConvertToFloat();
	globalValue::globalValuePoint().maxStrain = suGlobal::gOnly().getParam("max strain").ConvertToFloat();
	globalValue::globalValuePoint().boundaryThickness = suGlobal::gOnly().getParam("boundary thickness").ConvertToInteger();
	globalValue::globalValuePoint().density = suGlobal::gOnly().getParam("density").ConvertToFloat();
	globalValue::globalValuePoint().youngModules = suGlobal::gOnly().getParam("young modules").ConvertToFloat();
	globalValue::globalValuePoint().possionRatio = suGlobal::gOnly().getParam("possion ratio").ConvertToFloat();
	globalValue::globalValuePoint().force = suGlobal::gOnly().getParam("force").ConvertToFloat();
	globalValue::globalValuePoint().tAlpha = suGlobal::gOnly().getParam("tAlpha").ConvertToFloat();
}


globalValue::globalValue() :octreeDepth(3), mcResolution(5), metaballThreshold(1),
oofemPath("D:\\oofem\\build2.3\\debug\\"), iteratorTimes(1), volumeFraction(0.5),
currentIteratorTime(0)
{
	loadRegion = new float[6];
	for (int i = 0; i < 6; i++) {
		loadRegion[i] = 0.0;
	}
}

globalValue::~globalValue()
{
	delete[]loadRegion;
}
