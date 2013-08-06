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
	for (int i = 0; i < relationNum; ++i) {
		boost::unordered_map<std::string, double, hash::fnv_1> features;
		features.max_load_factor(0.1);
		features.reserve(1000000);
		featureTheta[i] = features;
	}
}

FeatureController::FeatureController(int n) {
	relationNum = n;
	for (int i = 0; i < relationNum; ++i) {
		boost::unordered_map<std::string, double, hash::fnv_1> features;
		features.max_load_factor(0.1);
		features.reserve(1000000);
		featureTheta[i] = features;
	}
}

void FeatureController::clear() {
	for (boost::unordered_map<int, boost::unordered_map<std::string, double, hash::fnv_1> >::iterator it =
			featureTheta.begin(); it != featureTheta.end(); ++it) {
		it->second.clear();
	}
	featureTheta.clear();
}

double FeatureController::getFeatureValue(int relationLabel, string feature) {

	if (!FeatureController::hasFeatureValue(relationLabel, feature))
		return 0;

	double ret = featureTheta.find(relationLabel)->second.find(feature)->second;
	return ret;
}

void FeatureController::setFeatureValue(int relationLabel, string feature,
		double value) {
	(featureTheta[relationLabel])[feature] = value;
}

void FeatureController::addFeatureByOne(int relationLabel, string feature) {
	(featureTheta[relationLabel])[feature] = this->getFeatureValue(
			relationLabel, feature) + 1;
}

void FeatureController::reduceFeatureByOne(int relationLabel, string feature) {
	(featureTheta[relationLabel])[feature] = this->getFeatureValue(
			relationLabel, feature) - 1;
}

bool FeatureController::hasFeatureValue(int relationLabel, string feature) {

	boost::unordered_map<string, double, hash::fnv_1> featuresForOneRelation =
			featureTheta.find(relationLabel)->second;

	if (featuresForOneRelation.find(feature) == featuresForOneRelation.end())
		return false;

	return true;
}

