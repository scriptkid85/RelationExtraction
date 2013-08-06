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

FeatureController::FeatureController() {
	relationNum = RELATION_MAX_NUM;
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
}

FeatureController::FeatureController(int n) {
	relationNum = n;
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
}

void FeatureController::clear() {
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ++it)
		delete it->second;
	featureTheta.clear();
}

double FeatureController::getFeatureValue(int relationLabel, string feature) {

	if (!FeatureController::hasFeatureValue(relationLabel, feature))
		return 0;

	double ret = (featureTheta.find(feature)->second)[relationLabel];
	return ret;
}

void FeatureController::setFeatureValue(int relationLabel, string feature,
		double value) {
	if (!FeatureController::hasFeatureValue(relationLabel, feature))
		featureTheta[feature] = new double[relationNum]();
	featureTheta[feature][relationLabel] = value;
}

void FeatureController::addFeatureByOne(int relationLabel, string feature) {
	if (!FeatureController::hasFeatureValue(relationLabel, feature)) {
		featureTheta[feature] = new double[relationNum]();
		featureTheta[feature][relationLabel] = 1;
	} else
		featureTheta[feature][relationLabel] = this->getFeatureValue(
				relationLabel, feature) + 1;
}

void FeatureController::reduceFeatureByOne(int relationLabel, string feature) {
	if (!FeatureController::hasFeatureValue(relationLabel, feature)) {
		featureTheta[feature] = new double[relationNum]();
		featureTheta[feature][relationLabel] = -1;
	} else
		featureTheta[feature][relationLabel] = this->getFeatureValue(
				relationLabel, feature) - 1;
}

bool FeatureController::hasFeatureValue(int relationLabel, string feature) {

	if (featureTheta.find(feature) == featureTheta.end())
		return false;

	return true;
}

void FeatureController::save(string filename) {
	ofstream thetas;
	thetas.open(filename.c_str());
	for (boost::unordered_map<string, double*, hash::fnv_1>::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ++it) {
		thetas << it->first;
		for (int i = 0; i < relationNum; ++i)
			thetas << " " << (it->second)[i];
		thetas << endl;
	}
	thetas.close();
}

void FeatureController::load(string filename) {
	ifstream file(filename.c_str());
	string line;
	vector<string> featurevalue;
	getline(file, line);
	featurevalue.clear();
	boost::split(featurevalue, line, boost::is_any_of(" "));
	if (featurevalue.size() != relationNum + 1) {
		cout << "error loading featuretheta, mismatch of relation number"
				<< endl;
		exit(0);
	}
	featureTheta.clear();
	featureTheta.reserve(1000000);
	featureTheta.max_load_factor(0.1);
	double *values = new double[relationNum]();
	for (int i = 0; i < relationNum; ++i) {
		values[i] = atof(featurevalue[i + 1].c_str());
	}
	featureTheta[featurevalue[0]] = values;
	while (getline(file, line)) {

		featurevalue.clear();
		boost::split(featurevalue, line, boost::is_any_of(" "));
		double *values = new double[relationNum]();
		for (int i = 0; i < relationNum; ++i) {
			values[i] = atof(featurevalue[i + 1].c_str());
		}
		featureTheta[featurevalue[0]] = values;
	}

}
