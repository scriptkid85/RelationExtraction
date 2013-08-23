//============================================================================
// Name        : RelationExtraction_UW.cpp
// Author      : Guanyu Wang
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include "FeatureController.h"
#include <vector>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <unistd.h>

#define MAX_RUN 50
#define RELATION_NUM 24

using namespace std;

typedef pair<string, double> FeatureValue;

FeatureController fc(RELATION_NUM + 1);
int run_num = MAX_RUN;
int fire_time_threshold = 0;

void printVector(vector<string> v) {
	for (size_t i = 0; i < v.size(); ++i) {
		cout << i << ":" << v[i] << ", ";
	}
	cout << endl;
}

void printVector(vector<int> v) {
	for (size_t i = 0; i < v.size(); ++i) {
		cout << i << ":" << v[i] << ", ";
	}
	cout << endl;
}

void printVector(vector<double> v) {
	for (size_t i = 0; i < v.size(); ++i) {
		cout << i << ":" << v[i] << ", ";
	}
	cout << endl;
}

double computeThetaExtract(int RelationLabel, vector<string> features) {
	double ret = 0;
	for (size_t feature = 0; feature < features.size(); ++feature) {
		ret += fc.getFeatureValue(RelationLabel, features[feature]);
	}
	return exp(ret);
}


double computeThetaExtractMatrix(int RelationLabel, vector<double> features) {
	double ret = 0;
	stringstream ss;
	for (size_t feature = 0; feature < features.size(); ++feature) {
		ss << feature;
		ret += features[feature] * fc.getFeatureValue(RelationLabel, ss.str());
	}
	return exp(ret);
}

string significantFeature(int RelationLabel, vector<string> features) {
	if (features.size() == 0)
		return "";
	double maxextract = fc.getFeatureValue(RelationLabel, features[0]);
	string ret = features[0];
	for (size_t feature = 1; feature < features.size(); ++feature) {
		double extract = fc.getFeatureValue(RelationLabel, features[feature]);
		if (extract > maxextract) {
			maxextract = extract;
			ret = features[feature];
		}
	}
	ret += "_";
	ostringstream strs;
	strs << maxextract;
	ret += strs.str();
	return ret;
}

vector<double> valueForRelations(vector<string> feature) {
	vector<double> ret;
	//always has none relation
	int numberofrelation = RELATION_NUM + 1;
	double thetaExtract = 0;
	for (int j = 0; j < numberofrelation; ++j) {
		thetaExtract = computeThetaExtract(j, feature);
		ret.push_back(thetaExtract);
	}
	return ret;
}

vector<string> sigFeatureForRelations(vector<string> feature) {
	vector<string> ret;
	//always has none relation
	int numberofrelation = RELATION_NUM + 1;
	for (int j = 0; j < numberofrelation; ++j) {
		ret.push_back(significantFeature(j, feature));
	}
	return ret;
}

vector<vector<int> > maximizeYZ(vector<vector<string> > featureOfOnePair) {
	vector<vector<int> > ret;
	vector<int> labelVector;
	size_t numberofsentence = featureOfOnePair.size();
	//always has none relation
	int numberofrelation = RELATION_NUM + 1;
	size_t sentence;
	vector<int> Y(numberofrelation, 0);
	vector<string> feature;
	for (sentence = 0; sentence < numberofsentence; ++sentence) {
		double max = -1, thetaExtract = 0;
		int label = 0;
		feature = featureOfOnePair[sentence];
		for (int j = 0; j < numberofrelation; ++j) {
			thetaExtract = computeThetaExtract(j, feature);
			if (thetaExtract > max) {
				max = thetaExtract;
				label = j;
			}
		}
		labelVector.push_back(label);
		Y[label] = 1;
	}

	ret.push_back(Y);
	ret.push_back(labelVector);

	return ret;
}

vector<int> maximizeZfast(vector<vector<string> > featureOfOnePair,
		vector<int> Y) {

	vector<int> ret(featureOfOnePair.size(), -1);
	vector<int> nonzeroLabel;
	vector<vector<double> > DP;
	size_t numberofrelation = Y.size();
	size_t numberofsentence = featureOfOnePair.size();

	//compute the dynamic programming table
	for (size_t i = 0; i < numberofrelation; ++i) {
		if (Y[i] == 0)
			continue;
		nonzeroLabel.push_back(i);
	}

	size_t numberofnonzeros = nonzeroLabel.size();
	if (numberofsentence < numberofnonzeros) {
		vector<int> empty;
		return empty;
	}
	//store all weight

	int id = 0;
	double max = -1;

	//assign a sentence to the nonzero labels first
	for (size_t i = 0; i < numberofnonzeros; ++i) {

		id = 0;
		int currentid = nonzeroLabel[i];

		max = computeThetaExtract(currentid, featureOfOnePair[0]);
		vector<double> dp;
		dp.push_back(max);

		for (size_t sentence = 1; sentence < numberofsentence; ++sentence) {
			double extract = computeThetaExtract(currentid,
					featureOfOnePair[sentence]);
			if (extract > max) {
				max = extract;
				id = sentence;
			}
			dp.push_back(extract);
		}
		ret[id] = currentid;
		DP.push_back(dp);
	}

	double extract;
	int best;
	for (size_t i = 0; i < numberofsentence; ++i) {
		int label = ret[i];
		if (label != -1)
			continue;

		best = 0;
		max = DP[0][i];
		extract = -1;
		for (size_t relation = 1; relation < numberofnonzeros; ++relation) {
			extract = DP[relation][i];
			if (extract > max) {
				max = extract;
				best = relation;
			}
		}
		ret[i] = nonzeroLabel[best];
	}

	for (size_t i = 0; i < DP.size(); i++) {
		DP[i].clear();
	}
	DP.clear();
	nonzeroLabel.clear();
	return ret;
}

vector<int> maximizeZ(vector<vector<string> > featureOfOnePair, vector<int> Y) {

	vector<int> ret(featureOfOnePair.size(), 0);
	vector<int> nonzeroLabel;
	vector<vector<double> > DP;

	//compute the dynamic programming table
	for (size_t i = 0; i < Y.size(); ++i) {
		if (Y[i] == 0)
			continue;
		nonzeroLabel.push_back(i);
	}

	if (featureOfOnePair.size() < nonzeroLabel.size()) {
//		cout
//				<< "Number of labels is larger than the number of mention sentences: "
//				<< endl;
		vector<int> empty;
		return empty;
	}

	vector<int> assignmentcount(Y.size(), 0);
	//store all weight
	vector<int> connected(nonzeroLabel.size(), 0);
	for (size_t i = 0; i < featureOfOnePair.size(); ++i) {
		vector<double> dp;
		double max = -1;
		for (size_t j = 0; j < nonzeroLabel.size(); ++j) {
			double thetaextract = computeThetaExtract(nonzeroLabel[j],
					featureOfOnePair[i]);
			dp.push_back(thetaextract);
			if (thetaextract > max) {
				max = thetaextract;
				ret[i] = nonzeroLabel[j];
			}
		}
		assignmentcount[ret[i]]++;
		DP.push_back(dp);
	}

	//approximation, find the relation has no assigned sentence
	for (size_t relation = 0; relation < nonzeroLabel.size(); relation++) {
		if (assignmentcount[nonzeroLabel[relation]] < 1) {
			int max = -1, id = 0;
			for (size_t i = 0; i < DP.size(); ++i) {
				if (assignmentcount[ret[i]] > 1 && DP[i][relation] > max) {
					max = DP[i][relation];
					id = i;
				}
			}
			assignmentcount[ret[id]]--;
			assignmentcount[nonzeroLabel[relation]] = 1;
			ret[id] = nonzeroLabel[relation];
		}
	}

	for (size_t i = 0; i < DP.size(); i++) {
		DP[i].clear();
	}
	DP.clear();
	assignmentcount.clear();
	nonzeroLabel.clear();

	return ret;
}

void UpdateTheta(vector<vector<string> > featureOfOnePair, vector<int> Zstar,
		vector<int> Z) {
	for (size_t i = 0; i < featureOfOnePair.size(); ++i) {
		for (size_t j = 0; j < featureOfOnePair[i].size(); ++j) {
			fc.addFeatureByOne(Zstar[i], featureOfOnePair[i][j]);
			fc.reduceFeatureByOne(Z[i], featureOfOnePair[i][j]);
		}
	}
}

//generate feature vector
vector<string> generateFV(string featurestring) {
	vector<string> ret, featurevalue;
	boost::split(ret, featurestring, boost::is_any_of(" "));
	for (size_t i = 0; i < ret.size(); ++i) {
		featurevalue.clear();
		boost::split(featurevalue, ret[i], boost::is_any_of(":"));
		if (featurevalue.size() != 2) {
			cout << "Invalid feature and value format: " << featurestring
					<< endl;
			continue;
		}
		ret[i] = featurevalue[0];
	}
	return ret;
}

vector<double> generateFVMatrix(string featurestring) {
	vector<string> features;
	vector<double> ret;
	boost::split(features, featurestring, boost::is_any_of(" "));
	for (size_t i = 0; i < ret.size(); ++i) {
		ret[i] = atof(features[i].c_str());
	}
	return ret;
}


//generate feature vector when do training
vector<string> generateTrainFV(string featurestring) {
	vector<string> ret = generateFV(featurestring);
	for(int i = 0; i < ret.size(); ++ i){
		fc.addFeatureCounter(ret[i]);
	}
	return ret;
}




vector<string> generateBiFV(string featurestring) {
	vector<string> ret, featurevalue, biret;
	boost::split(ret, featurestring, boost::is_any_of(" "));
	for (size_t i = 0; i < ret.size(); ++i) {
		featurevalue.clear();
		boost::split(featurevalue, ret[i], boost::is_any_of(":"));
		if (featurevalue.size() != 2) {
			cout << "Invalid feature and value format: " << featurestring
					<< endl;
			continue;
		}
		ret[i] = featurevalue[0];
	}
	for (size_t i = 0; i < ret.size(); ++i) {
		for (size_t j = i + 1; j < ret.size(); ++j) {
			string feature = "";
			if (ret[i].compare(ret[j]) > 0) {
				feature = ret[j];
				feature += "+";
				feature += ret[i];
			} else {
				feature = ret[i];
				feature += "+";
				feature += ret[j];
			}
			biret.push_back(feature);
		}
	}
	return biret;
}

//generate label vector
vector<int> generateLV(string labelstring) {
	vector<string> labels;
	vector<int> ret(RELATION_NUM + 1, 0);

	boost::split(labels, labelstring, boost::is_any_of(" "));
	for (size_t i = 0; i < labels.size(); ++i) {
		ret[atoi(labels[i].c_str())] = 1;
	}
	return ret;
}

void testSingle(string input) {
	int pos = input.find("--S--");
	vector<string> tokens;
	tokens.push_back(input.substr(0, pos));
	tokens.push_back(input.substr(pos + 5, input.length()));

	vector<vector<string> > featureOfOnePair;
	vector<string> fields;

	pos = tokens[1].find("| ");
	fields.push_back(tokens[1].substr(0, pos));
	fields.push_back(tokens[1].substr(pos + 2, tokens[1].length()));

	//get feature vector
	vector<string> features = generateFV(fields[1]);
	vector<vector<string> > fv;
	fv.push_back(features);
	cout << tokens[0] << ":";
	printVector(maximizeYZ(fv)[0]);
	fv[0].clear();
	fv.clear();
	fields.clear();
	tokens.clear();
}

void testMulti(string pair, vector<vector<string> > fv) {
	cout << pair << ":";
	printVector(maximizeYZ(fv)[0]);
}

bool vectorContain(vector<int> big, vector<int> small) {
	if (big == small)
		return true;
	bool flag = false;
	for (size_t i = 0; i < RELATION_NUM + 1; ++i) {
		if (big[i] == 1 && small[i] == 1)
			flag = true;
		else if (big[i] == 0 && small[i] == 1)
			return false;
	}
	return flag;
}

void testWithLabel(char *testfile) {

	int count = 1;
	int correct = 0;
	int partcorrect = 0;
	cout << "testFile: " << testfile << endl;
	ifstream file(testfile);

	string line;
	vector<string> labelAndInstance, tokens;
	tokens.clear();
	labelAndInstance.clear();
	getline(file, line);
	boost::split(labelAndInstance, line, boost::is_any_of("\t"));
	vector<int> Y = generateLV(labelAndInstance[0]);
	vector<int> Ytest;
	int pos = labelAndInstance[1].find("--S--");
	tokens.push_back(labelAndInstance[1].substr(0, pos));
	tokens.push_back(
			labelAndInstance[1].substr(pos + 5, labelAndInstance[1].length()));

	//get the pair
	string pair = tokens[0];
	vector<vector<string> > featureOfOnePair;

	pos = tokens[1].find("| ");

	//get feature vector
	vector<string> features = generateFV(
			tokens[1].substr(pos + 2, tokens[1].length()));
	featureOfOnePair.push_back(features);

	while (getline(file, line)) {
		tokens.clear();
		labelAndInstance.clear();
		boost::split(labelAndInstance, line, boost::is_any_of("\t"));
		int pos = labelAndInstance[1].find("--S--");
		tokens.push_back(labelAndInstance[1].substr(0, pos));
		tokens.push_back(
				labelAndInstance[1].substr(pos + 5,
						labelAndInstance[1].length()));

		//get the pair

		if (pair.compare(tokens[0]) == 0) {

			pos = tokens[1].find("| ");
			//get feature vector
			vector<string> features = generateFV(
					tokens[1].substr(pos + 2, tokens[1].length()));
			featureOfOnePair.push_back(features);
		} else {
			count++;
			Ytest = maximizeYZ(featureOfOnePair)[0];
			cout << pair << ": " << endl;
			cout << "test label: ";
			printVector(Ytest);
			cout << "true label: ";
			printVector(Y);
			if (Ytest == Y) {
				correct++;
				partcorrect++;
			} else if (vectorContain(Y, Ytest)) {
				partcorrect++;
			}
			Y.clear();
			Ytest.clear();
			pair = tokens[0];
			Y = generateLV(labelAndInstance[0]);
			featureOfOnePair.clear();
			pos = tokens[1].find("| ");
			//get feature vector
			vector<string> features = generateFV(
					tokens[1].substr(pos + 2, tokens[1].length()));
			featureOfOnePair.push_back(features);
		}
	}
	Y.clear();
	Ytest.clear();
	features.clear();
	tokens.clear();
	labelAndInstance.clear();

	cout << "exact-correctness: " << correct << "/" << count << "="
			<< correct / (double) count << endl;
	cout << "partial-correctness: " << partcorrect << "/" << count << "="
			<< partcorrect / (double) count << endl;
}

void testWithLabelwithBiFV(char *testfile) {

	int count = 1;
	int correct = 0;
	int partcorrect = 0;
	cout << "testFile: " << testfile << endl;
	ifstream file(testfile);

	string line;
	vector<string> labelAndInstance, tokens;
	tokens.clear();
	labelAndInstance.clear();
	getline(file, line);
	boost::split(labelAndInstance, line, boost::is_any_of("\t"));
	vector<int> Y = generateLV(labelAndInstance[0]);
	vector<int> Ytest;
	int pos = labelAndInstance[1].find("--S--");
	tokens.push_back(labelAndInstance[1].substr(0, pos));
	tokens.push_back(
			labelAndInstance[1].substr(pos + 5, labelAndInstance[1].length()));

	//get the pair
	string pair = tokens[0];
	vector<vector<string> > featureOfOnePair;

	pos = tokens[1].find("| ");

	//get feature vector
	vector<string> features = generateBiFV(
			tokens[1].substr(pos + 2, tokens[1].length()));
	featureOfOnePair.push_back(features);

	while (getline(file, line)) {
		tokens.clear();
		labelAndInstance.clear();
		boost::split(labelAndInstance, line, boost::is_any_of("\t"));
		int pos = labelAndInstance[1].find("--S--");
		tokens.push_back(labelAndInstance[1].substr(0, pos));
		tokens.push_back(
				labelAndInstance[1].substr(pos + 5,
						labelAndInstance[1].length()));

		//get the pair

		if (pair.compare(tokens[0]) == 0) {

			pos = tokens[1].find("| ");
			//get feature vector
			vector<string> features = generateBiFV(
					tokens[1].substr(pos + 2, tokens[1].length()));
			featureOfOnePair.push_back(features);
		} else {
			count++;
			Ytest = maximizeYZ(featureOfOnePair)[0];
			cout << pair << ": " << endl;
			cout << "test label: ";
			printVector(Ytest);
			cout << "true label: ";
			printVector(Y);
			if (Ytest == Y) {
				correct++;
				partcorrect++;
			} else if (vectorContain(Y, Ytest)) {
				partcorrect++;
			}
			Y.clear();
			Ytest.clear();
			pair = tokens[0];
			Y = generateLV(labelAndInstance[0]);
			featureOfOnePair.clear();
			pos = tokens[1].find("| ");
			//get feature vector
			vector<string> features = generateBiFV(
					tokens[1].substr(pos + 2, tokens[1].length()));
			featureOfOnePair.push_back(features);
		}
	}
	Y.clear();
	Ytest.clear();
	features.clear();
	tokens.clear();
	labelAndInstance.clear();

	cout << "exact-correctness: " << correct << "/" << count << "="
			<< correct / (double) count << endl;
	cout << "partial-correctness: " << partcorrect << "/" << count << "="
			<< partcorrect / (double) count << endl;
}

void testWithoutLabel(char *testfile) {

	cout << "testFile: " << testfile << endl;
	ifstream file(testfile);

	string line;
	vector<string> tokens;
	tokens.clear();
	getline(file, line);
	vector<int> Ytest;
	int pos = line.find("--S--");
	tokens.push_back(line.substr(0, pos));
	tokens.push_back(line.substr(pos + 5, line.length()));

	//get the pair
	string pair = tokens[0];
	vector<vector<string> > featureOfOnePair;

	pos = tokens[1].find("| ");

	//get feature vector
	vector<string> features = generateFV(
			tokens[1].substr(pos + 2, tokens[1].length()));
	featureOfOnePair.push_back(features);
	Ytest = maximizeYZ(featureOfOnePair)[0];
	cout << pair << ": " << endl;
	cout << "test label: ";
	printVector(Ytest);
	cout << "test label value: ";
	printVector(valueForRelations(featureOfOnePair[0]));
	printVector(sigFeatureForRelations(featureOfOnePair[0]));
	featureOfOnePair[0].clear();
	featureOfOnePair.clear();

	while (getline(file, line)) {
		tokens.clear();
		int pos = line.find("--S--");
		tokens.push_back(line.substr(0, pos));
		tokens.push_back(line.substr(pos + 5, line.length()));

		pair = tokens[0];
		pos = tokens[1].find("| ");
		//get feature vector
		vector<string> features = generateFV(
				tokens[1].substr(pos + 2, tokens[1].length()));
		featureOfOnePair.push_back(features);
		Ytest = maximizeYZ(featureOfOnePair)[0];
		cout << pair << ": " << endl;
		cout << "test label: ";
		printVector(Ytest);
		cout << "test label value: ";
		printVector(valueForRelations(featureOfOnePair[0]));
		printVector(sigFeatureForRelations(featureOfOnePair[0]));
		featureOfOnePair[0].clear();
		featureOfOnePair.clear();
		Ytest.clear();
	}
	Ytest.clear();
	features.clear();
	tokens.clear();

}

void testWithoutLabelwithBiFV(char *testfile) {

	cout << "testFile: " << testfile << endl;
	ifstream file(testfile);

	string line;
	vector<string> tokens;
	tokens.clear();
	getline(file, line);
	vector<int> Ytest;
	int pos = line.find("--S--");
	tokens.push_back(line.substr(0, pos));
	tokens.push_back(line.substr(pos + 5, line.length()));

	//get the pair
	string pair = tokens[0];
	vector<vector<string> > featureOfOnePair;

	pos = tokens[1].find("| ");

	//get feature vector
	vector<string> features = generateBiFV(
			tokens[1].substr(pos + 2, tokens[1].length()));
	featureOfOnePair.push_back(features);
	Ytest = maximizeYZ(featureOfOnePair)[0];
	cout << pair << ": " << endl;
	cout << "test label: ";
	printVector(Ytest);
	cout << "test label value: ";
	printVector(valueForRelations(featureOfOnePair[0]));
	printVector(sigFeatureForRelations(featureOfOnePair[0]));
	featureOfOnePair[0].clear();
	featureOfOnePair.clear();

	while (getline(file, line)) {
		tokens.clear();
		int pos = line.find("--S--");
		tokens.push_back(line.substr(0, pos));
		tokens.push_back(line.substr(pos + 5, line.length()));

		pair = tokens[0];
		pos = tokens[1].find("| ");
		//get feature vector
		vector<string> features = generateBiFV(
				tokens[1].substr(pos + 2, tokens[1].length()));
		featureOfOnePair.push_back(features);
		Ytest = maximizeYZ(featureOfOnePair)[0];
		cout << pair << ": " << endl;
		cout << "test label: ";
		printVector(Ytest);
		cout << "test label value: ";
		printVector(valueForRelations(featureOfOnePair[0]));
		printVector(sigFeatureForRelations(featureOfOnePair[0]));
		featureOfOnePair[0].clear();
		featureOfOnePair.clear();
		Ytest.clear();
	}
	Ytest.clear();
	features.clear();
	tokens.clear();

}

//input line format: relation_id relation_id...\tArg1--R--Arg2--S--SentenceID| f1:v1 f2:v2...
void train(char *trainfile) {
	cout << "trainfile: " << trainfile << endl;
	int run;
	for (run = 0; run < run_num; ++run) {
		cout << "run #: " << run << endl;
		if(run == 1)fc.removeLessFrequentFeature(fire_time_threshold);
		if(run == 1)fc.freezeSet();
		ifstream file(trainfile);

		string line, pair, labels;

		getline(file, line);
		vector<string> labelAndInstance, tokens;
		tokens.clear();
		labelAndInstance.clear();
		boost::split(labelAndInstance, line, boost::is_any_of("\t"));
		if (labelAndInstance.size() != 2) {
			cout << "Invalid label instance format(\\t split problem): " << line
					<< endl;
			exit(0);
		}

		//get the labels for this pair of entities
		labels = labelAndInstance[0];
		vector<int> Y = generateLV(labels);
		vector<int> Zstar;
		vector<vector<int> > YZ;

		int pos = labelAndInstance[1].find("--S--");
		tokens.push_back(labelAndInstance[1].substr(0, pos));
		tokens.push_back(
				labelAndInstance[1].substr(pos + 5,
						labelAndInstance[1].length()));

		//get the pair
		pair = tokens[0];
		vector<vector<string> > featureOfOnePair;
		vector<string> fields;

		pos = tokens[1].find("| ");
		fields.push_back(tokens[1].substr(0, pos));
		fields.push_back(tokens[1].substr(pos + 2, tokens[1].length()));

		//get feature vector
		vector<string> features = generateTrainFV(fields[1]);

		featureOfOnePair.push_back(features);

		int count = 1, errorcount = 0;
		while (getline(file, line)) {
			count++;
			//			cout << line << endl;
			//			cout << "\r";
			//			cout << countl;

			tokens.clear();
			labelAndInstance.clear();
			boost::split(labelAndInstance, line, boost::is_any_of("\t"));
			if (labelAndInstance.size() != 2) {
				cout << "Invalid label instance format: " << line << endl;
				continue;
			}

			pos = labelAndInstance[1].find("--S--");
			tokens.push_back(labelAndInstance[1].substr(0, pos));
			tokens.push_back(
					labelAndInstance[1].substr(pos + 5,
							labelAndInstance[1].length()));

			//if the pair does not changes
			if (pair.compare(tokens[0]) == 0) {

				pos = tokens[1].find("| ");
				fields.clear();
				fields.push_back(tokens[1].substr(0, pos));
				fields.push_back(tokens[1].substr(pos + 2, tokens[1].length()));

				//get feature vector
				vector<string> features = generateTrainFV(fields[1]);
				featureOfOnePair.push_back(features);

			}
			//if the pair changes, i.e. already parsed all mentions for this pair
			else {
				//process the featureOfOnePair for the previous pair
				YZ.clear();
				YZ = maximizeYZ(featureOfOnePair);

				//if Y prime is not equal to the distant labeled Y

				if (!(YZ[0] == Y)) {
					Zstar.clear();
					Zstar = maximizeZ(featureOfOnePair, Y);
					if (Zstar.size() == 0) {
						errorcount++;
//						cout << count << endl;
					}
					//					printVector(Zstar);
					//update theta
					else {

						UpdateTheta(featureOfOnePair, Zstar, YZ[1]);

					}
				}

				//update the pair and the label vector for the pair
				pair = tokens[0];
				labels = labelAndInstance[0];

				//clear label vector
				Y.clear();
				//clear the featureOfOnePair vector
				for (size_t i = 0; i < featureOfOnePair.size(); ++i) {
					featureOfOnePair[i].clear();
				}
				featureOfOnePair.clear();
				fields.clear();

				Y = generateLV(labels);
				pos = tokens[1].find("| ");
				fields.push_back(tokens[1].substr(0, pos));
				fields.push_back(tokens[1].substr(pos + 2, tokens[1].length()));

				//get feature vector
				vector<string> features = generateTrainFV(fields[1]);
				featureOfOnePair.push_back(features);
			}
		}
		Zstar.clear();
		YZ.clear();
		Y.clear();
	}

}


int main(int argc, char *args[]) {
	char c;
	int test_with_label_flag = 0,
		bifeature_flag = 0,
	opterr = 0;
	char *train_file = NULL, *test_file = NULL, *output_theta = NULL;
	while ((c = getopt(argc, args, "w:r:lbt:d:f:")) != -1){
		switch (c) {
		case 'l':
			test_with_label_flag ++;
			break;
		case 'b':
			bifeature_flag ++;
			break;
		case 't':
			train_file = optarg;
			break;
		case 'w':
			output_theta = optarg;
			break;
		case 'd':
			test_file = optarg;
			break;
		case 'f':
			fire_time_threshold = atoi(optarg);
			break;
		case 'r':
			run_num = atoi(optarg);
			break;
		case '?':
			if(optopt == 'b'){
				cout << "Using unifeatures." << endl;
				break;
			}
			else if (optopt == 'f')
				fprintf(stderr, "Option -%c requires an argument for the firetime threshold.\n", optopt);
		    else if (optopt == 't')
				fprintf(stderr, "Option -%c requires an argument for the training file.\n", optopt);
			else if (optopt == 'd')
				fprintf(stderr, "Option -%c requires an argument for the testing file.\n", optopt);
			else if (optopt == 'r')
				fprintf(stderr, "Option -%c requires an argument for the run times.\n", optopt);
			else if (optopt == 'w')
				fprintf(stderr, "Option -%c requires an argument to output the trained theta file.\n", optopt);
			else if (isprint(optopt))
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			break;

		}
	}



	if(!train_file){
		fprintf(stderr, "Need trainning file with -t.\n");
		return 1;
	}
	if(!test_file){
		fprintf(stderr, "Need testing file with -d.\n");
		return 1;
	}
	if(!output_theta){
		fprintf(stderr, "Need output_theta file with -w.\n");
		return 1;
	}


	cout << "Take training file: " << train_file << endl;
	cout << "Take testing file: " << test_file << endl;
	cout << "Write to theta file: " << output_theta << endl;
	cout << "Using firetime threshold: " << fire_time_threshold << endl;
	cout << "Using round of run: " << run_num << endl;


	if(test_with_label_flag == 1){
		cout << "Testing with given labels" << endl;
	}
	else cout << "Testing without given labels" << endl;

	if(bifeature_flag == 0){
		cout << "Training with Unifeatures" << endl;
	}
	else cout << "Training with Bifeatures" << endl;


	if(test_with_label_flag && !bifeature_flag){
		train(train_file);
		testWithLabel(test_file);
		fc.save(output_theta, 0);
		fc.clear();
	}
	else if(!test_with_label_flag && !bifeature_flag){
		train(train_file);
		testWithoutLabel(test_file);
		fc.save(output_theta, 0);
		fc.clear();
	}
}
