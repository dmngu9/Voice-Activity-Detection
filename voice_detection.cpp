/*
 * main.cpp
 *
 *  Created on: May 26, 2016
 *      Author: dmngu9
 */
#include "VAD.h"
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv){
	
	const int fs = 11025;
	const int WinSize = 256;
	const int order = 5;
	const double threshold = -6;
	const int uSize = 4;

	ifstream soundFile(argv[1]);
	vector<double> sound;

	if(soundFile.is_open()){
		double value;
		while(soundFile >> value){
			sound.push_back(value);
			if(sound.size() == 50176)
				break;
		}
		soundFile.close();
	}

	double* signal = &sound[0];

	VAD vad(WinSize,sound.size(),order,threshold);
	vector<double> outcome = vad.compute(signal);
	MatrixXd enFrame = vad.getNormalizedEnFrame();

	double maxLevel = 0;
	for(int i = 0; i < sound.size(); i++){
		if(maxLevel < abs(sound[i]))
			maxLevel = abs(sound[i]);
	}

	maxLevel += 0.01*maxLevel;
	vector<int> idx(outcome.size(),0);
	for(int i = 0; i < outcome.size(); i++){
		idx[i] = (outcome[i] > 2.5) ? 1 : 0;
	}

	VectorXd d = VectorXd::Zero(idx.size()-1);
	VectorXd vadStart, vadEnd;

	for(int i = 0; i < outcome.size()-1; i++){
		d(i) = idx[i+1] - idx[i];
		if(d(i) == 1){
			vadStart.conservativeResize(vadStart.size()+1);
			vadStart(vadStart.size()-1) = i;
		}
		else if (d(i) == -1){ 
			vadEnd.conservativeResize(vadEnd.size()+1);
			vadEnd(vadEnd.size()-1) = i;
		}
	}

	double q = (double) WinSize/fs;
	VectorXd temp = vadEnd - vadStart;
	VectorXd len = temp*q;
	vector<int> VAD_begin, VAD_end;
	for(int i = 0; i < len.size(); i++){
		if(len(i) >= (uSize*WinSize/fs)){
			VAD_begin.push_back(vadStart(i));
			VAD_end.push_back(vadEnd(i));
		}
	}

	//plot sound wave here 
	cout << enFrame.row(0) << endl;
	for(int i = 0; i < VAD_begin.size(); i++){
		double x_start = enFrame(0,VAD_begin[i]+1) + 0.5*order*WinSize;
		double x_end = enFrame(enFrame.rows()-1,VAD_end[i]+1) + 0.5*order*WinSize;
		// VectorXd x, y;
		// x << x_start, x_end, x_end, x_start, x_start;
		// y << maxLevel, maxLevel, -maxLevel, -maxLevel, maxLevel;

		stringstream ss;
    		ss << "example_" << i << ".txt";
    		ofstream myfile;
    		cout << ss.str();
		myfile.open (ss.str().c_str());
		myfile << x_start << "\t\t\t" << maxLevel << "\n";
		myfile << x_end << "\t\t\t" << maxLevel << "\n";
		myfile << x_end << "\t\t\t" << -maxLevel << "\n";
		myfile << x_start << "\t\t\t" << -maxLevel << "\n";
		myfile << x_start << "\t\t\t" << maxLevel << "\n";
		myfile.close();
	}
	return 0;
}



