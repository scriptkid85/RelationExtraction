//============================================================================
// Name        : FeatureController.cpp
// Author      : Guanyu Wang
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include "FeatureController.h"
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
	for(boost::unordered_map<string, double* , hash::fnv_1>::iterator it = featureTheta.begin(); it != featureTheta.end(); ++ it)
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
	if (!FeatureController::hasFeatureValue(relationLabel, feature)){
				featureTheta[feature] = new double[relationNum]();
				featureTheta[feature][relationLabel] = 1;
	}
	else featureTheta[feature][relationLabel] = this->getFeatureValue(
				relationLabel, feature) + 1;
}

void FeatureController::reduceFeatureByOne(int relationLabel, string feature) {
	if (!FeatureController::hasFeatureValue(relationLabel, feature)){
			featureTheta[feature] = new double[relationNum]();
			featureTheta[feature][relationLabel] = -1;
	}
	else featureTheta[feature][relationLabel] = this->getFeatureValue(
			relationLabel, feature) - 1;
}

bool FeatureController::hasFeatureValue(int relationLabel, string feature) {

	if (featureTheta.find(feature) == featureTheta.end())
		return false;

	return true;
}

