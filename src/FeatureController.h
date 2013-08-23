#include <iostream>
#include <string>
#include <boost/unordered_map.hpp>
#include <fstream>
#include "fnv1.hpp"

using namespace std;

class FeatureController{
public:

	int relationNum;

	FeatureController();
	FeatureController(int n);

	void addFeatureCounter(string feature);
	void removeLessFrequentFeature(double firetime);

	double getFeatureValue(int relationLabel, string feature);
	int getFeatureCount(int relationLabel, string feature);
	void setFeatureValue(int relationLabel, string feature, double value);
	bool hasFeature(int relationLabel, string feature);
	void addFeatureByOne(int relationLabel, string feature);
	void reduceFeatureByOne(int relationLabel, string feature);

	double getBiFeatureValue(int relationLabel, string featureOne, string featureTwo);
	void setBiFeatureValue(int relationLabel, string featureOne, string featureTwo, double value);
	bool hasBiFeatureValue(int relationLabel, string featureOne, string featureTwo);
	void addBiFeatureByOne(int relationLabel, string featureOne, string featureTwo);
	void reduceBiFeatureByOne(int relationLabel, string featureOne, string featureTwo);

	void freezeSet();
	void unfreezeSet();

	void clear();
	void save(string file, int firetime);
	void load(string file);

private:
	static const int FEATURE_MAX_NUM = INT_MAX;
	static const int RELATION_MAX_NUM = 100;
	int setMutable_flag;

	boost::unordered_map<string, double* , hash::fnv_1> featureTheta;
};
