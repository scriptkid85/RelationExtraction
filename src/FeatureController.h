#include <iostream>
#include <string>
#include <boost/unordered_map.hpp>
#include "fnv1.hpp"

using namespace std;

class FeatureController{
public:

	int relationNum;

	FeatureController();
	FeatureController(int n);
	double getFeatureValue(int relationLabel, string feature);
	void setFeatureValue(int relationLabel, string feature, double value);
	bool hasFeatureValue(int relationLabel, string feature);
	void addFeatureByOne(int relationLabel, string feature);
	void reduceFeatureByOne(int relationLabel, string feature);
	void clear();

private:
	static const int FEATURE_MAX_NUM = INT_MAX;
	static const int RELATION_MAX_NUM = 100;

	boost::unordered_map<string, double* , hash::fnv_1> featureTheta;
};
