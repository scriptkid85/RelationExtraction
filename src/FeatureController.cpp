//============================================================================
// Name        : FeatureController.cpp
// Author      : Guanyu Wang
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include "FeatureController.h"
#include <boost/algorithm/string.hpp>
#include <string>

using namespace std;

//the first position in the array is the counter
FeatureController::FeatureController() {
	relationNum = RELATION_MAX_NUM;
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
	unfreezeSet();
}

FeatureController::FeatureController(int n) {
	relationNum = n;
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
	unfreezeSet();
}

void FeatureController::freezeSet() {
	cout << "Freeze the Feature Space" << endl;
	setMutable_flag = 0;
}

void FeatureController::unfreezeSet() {
	cout << "Unfreeze the Feature Space" << endl;
	setMutable_flag = 1;
}

void FeatureController::clear() {
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ++it)
		delete it->second;
	featureTheta.clear();
}

double FeatureController::getFeatureValue(int relationLabel, string feature) {

	if (!FeatureController::hasFeature(relationLabel, feature))
		return 0;
	double ret = (featureTheta.find(feature)->second)[relationLabel + 1];
	return ret;
}

int FeatureController::getFeatureCount(int relationLabel, string feature) {
	if (!FeatureController::hasFeature(relationLabel, feature))
		return 0;
	double ret = (featureTheta.find(feature)->second)[0];
	return ret;
}

void FeatureController::setFeatureValue(int relationLabel, string feature,
		double value) {
	if (!FeatureController::hasFeature(relationLabel, feature)){
		if(setMutable_flag == 1)
			featureTheta[feature] = new double[relationNum + 1]();
		else return;
	}
	featureTheta[feature][relationLabel + 1] = value;
}

void FeatureController::addFeatureByOne(int relationLabel, string feature) {
	if (!FeatureController::hasFeature(relationLabel, feature)) {
		if(setMutable_flag == 1){
			featureTheta[feature] = new double[relationNum + 1]();
			featureTheta[feature][relationLabel + 1] = 1;
		}
		else return;

	} else
		featureTheta[feature][relationLabel + 1] = this->getFeatureValue(
				relationLabel, feature) + 1;
}

void FeatureController::reduceFeatureByOne(int relationLabel, string feature) {
	if (!FeatureController::hasFeature(relationLabel, feature)) {
		if(setMutable_flag == 1){
			featureTheta[feature] = new double[relationNum + 1]();
			featureTheta[feature][relationLabel + 1] = -1;
		}
		else return;
	} else
		featureTheta[feature][relationLabel + 1] = this->getFeatureValue(
				relationLabel, feature) - 1;
}

bool FeatureController::hasFeature(int relationLabel, string feature) {

	if (featureTheta.find(feature) == featureTheta.end())
		return false;

	return true;
}



void FeatureController::addFeatureCounter(string feature){
	if (featureTheta.find(feature) == featureTheta.end()){
		if(setMutable_flag == 1){
			featureTheta[feature] = new double[relationNum + 1]();
			featureTheta[feature][0] = 1;
		}
		else return;
	}
	else featureTheta[feature][0] ++;
}



void FeatureController::removeLessFrequentFeature(double firetime) {
	int count = 0;
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
				featureTheta.begin(); it != featureTheta.end(); ++ it) {
		count ++;
	}
	cout << "removeLessFrequentFeature: total number of features before remove:" << count << endl;

	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ) {
		if((it -> second)[0] < firetime) {
			it = featureTheta.erase(it);
		}
		else{
			++ it;
		}
	}

	count = 0;
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
				featureTheta.begin(); it != featureTheta.end(); ++ it) {
		count ++;
	}
	cout << "removeLessFrequentFeature: total number after remove:" << count << endl;


}

void FeatureController::save(string filename, int firetime) {
	ofstream thetas;
	thetas.open(filename.c_str());
	int count = 0;
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ++it) {
//		if((it->second)[0] < firetime)continue;
		count ++;
		thetas << it->first;
		for (int i = 0; i < relationNum + 1; ++i)
			thetas << " " << (it->second)[i];
		thetas << endl;
	}
	cout << "total number of features saved:" << count << endl;
	thetas.close();
}

void FeatureController::load(string filename) {
	ifstream file(filename.c_str());
	string line;
	vector<string> featurevalue;
	getline(file, line);
	featurevalue.clear();
	boost::split(featurevalue, line, boost::is_any_of(" "));
	if (featurevalue.size() != relationNum + 2) {
		cout << "error loading featuretheta, mismatch of relation number, need: " << relationNum + 2 << ", but get: "<< featurevalue.size() << endl;
		exit(0);
	}
	featureTheta.clear();
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
	double *values = new double[relationNum + 1]();
	for (int i = 0; i < relationNum + 1; ++i) {
		values[i] = atof(featurevalue[i + 1].c_str());
	}
	featureTheta[featurevalue[0]] = values;
	while (getline(file, line)) {

		featurevalue.clear();
		boost::split(featurevalue, line, boost::is_any_of(" "));
		double *values = new double[relationNum + 1]();
		for (int i = 0; i < relationNum + 1; ++i) {
			values[i] = atof(featurevalue[i + 1].c_str());
		}
		featureTheta[featurevalue[0]] = values;
	}

}
