/*
 * VAD.h
 *
 *  Created on: May 27, 2016
 *      Author: dmngu9
 */

#ifndef VAD_H_
#define VAD_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <fftw3.h>

using namespace std;
using namespace Eigen;

const double PI = 3.14;

class VAD{

private:
	int winSize;
	int signalSize;
	int NFFT2;
	int order;
	double threshold;
	MatrixXd enFrame;
	RowVectorXd hamming;
	VectorXd averageNoise;
	MatrixXd amplitude;

public:

	VAD(int winSize, int signalSize, int order,double threshold){
		this->winSize = winSize;
		this->signalSize = signalSize;
		this->order = order;
		this->threshold = threshold;
		this->enFrame = MatrixXd::Zero(this->winSize, this->signalSize/(this->winSize*0.5));
		this->hamming = createHammingWindow();
		this->NFFT2 = this->winSize/2;
	}

	~VAD(){}

	RowVectorXd createHammingWindow(){
		double alpha = 0.54;
		double beta = 0.46;
		RowVectorXd hamming = RowVectorXd::Zero(winSize);
		for(int i = 0; i < winSize; i++){
			hamming(i) = alpha - beta*cos((2*PI*i)/(winSize-1));
		}
		return hamming;
	}

	void buffer(double* signal){
		MatrixXd upper = MatrixXd::Zero(this->winSize/2,this->signalSize/(this->winSize*0.5));
		MatrixXd below = MatrixXd::Zero(this->winSize/2,this->signalSize/(this->winSize*0.5));

		int j = 0;
		int k = 0;
		for(int i =0; i < this->signalSize; i++){
			below(j,k) = signal[i];
			j++;
			if(j == this->winSize/2){
				j = 0;
				k++;
				if(k < this->signalSize/(this->winSize*0.5))
					upper.col(k) = below.col(k-1);
			}
		}

		this->enFrame.topRows(this->winSize/2)= upper;
		this->enFrame.bottomRows(this->winSize/2) = below;
	}

	VectorXd fft_calc(double* fft_input){
		VectorXd result = VectorXd::Zero(NFFT2);
		fftw_complex* fft_result;
		fftw_plan p;

		fft_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->winSize);
		p = fftw_plan_dft_r2c_1d(this->winSize, fft_input, fft_result,FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);

		for(int i = 0; i < this->NFFT2; i++){
			result(i) = sqrt(pow(fft_result[i][0],2)+ pow(fft_result[i][1],2));
		}
		fftw_free(fft_result);
		return result;
	}

	VectorXd computeNoiseAverageSpectrum(){
		VectorXd averageNoiseSpectrum;
		int wnum = this->enFrame.cols();
		VectorXd avgAmp = VectorXd::Zero(this->NFFT2);
		for(int i = 0; i < wnum; i++){
			VectorXd s = this->enFrame.col(i);//got 6 in each col
			double fft_input[this->winSize];//winsize 6
			for(int j = 0; j < this->winSize; j++){
				fft_input[j] = this->hamming(j) * s(j);//size of 6
			}

			VectorXd temp = fft_calc(fft_input);
			avgAmp += temp;
		}
		averageNoiseSpectrum = avgAmp/wnum;
		return averageNoiseSpectrum;
	}

	//signal is enFrame
	VectorXd getAmplitude(int index){
		VectorXd amp;
		if(amplitude.rows() > index){
			amp = amplitude.row(index);
		}
		else{
			VectorXd s = this->enFrame.col(index);
			double fft_input[this->winSize];
			for(int j = 0; j < this->winSize; j++){
				fft_input[j] = this->hamming(j) * s(j);//size of 6
			}
			amp = fft_calc(fft_input);
			amplitude.conservativeResize(amplitude.rows()+1, this->NFFT2);
			amplitude.row(index) = amp;
		}
		return amp;
	}

	VectorXd findMax(VectorXd& a, VectorXd& b){
		VectorXd result = VectorXd::Zero(this->NFFT2);
		for(int i = 0; i < this->NFFT2; i++){
			result(i) = (a(i) > b(i)) ? a(i) : b(i);
		}
		return result;
	}

	VectorXd ltse(int index){
		VectorXd maxmag = VectorXd::Zero(this->NFFT2);
		VectorXd maxamp;
		int i = index - order;
		while(i != index+order){
			VectorXd amp = getAmplitude(i);
			maxamp = findMax(amp,maxmag);
			i++;
		}
		return maxamp;
	}

	double ltsd(int index){
		if(index < (this->order) || (index+order >= this->enFrame.cols())){
			return 0.0;
		}
	
		VectorXd ltseOutput = ltse(index);
		ltseOutput = ltseOutput.array().square();
		VectorXd sp = ltseOutput.array()/this->averageNoise.array();
		double sum = 0;
		for(int i = 0; i < sp.size(); i++){
			sum += sp(i)/this->NFFT2;
		}		
		
		double result = 10 * log10(sum);

		if(result < this->threshold){
			this->averageNoise = 0.54 * this->averageNoise + (1-0.54)*sum*VectorXd::Ones(averageNoise.size());
		}
		return result;
	}

	vector<double> compute(double* signal){
		buffer(signal);
		int wnum = this->enFrame.cols();
		vector<double> ltsds;
		this->averageNoise = computeNoiseAverageSpectrum().array().square();
		for(int i = 0; i < wnum; i++){
			ltsds.push_back(ltsd(i));
		}
		return ltsds;
	}

	MatrixXd getNormalizedEnFrame(){
		MatrixXd upper = MatrixXd::Zero(this->winSize/2,this->signalSize/(this->winSize*0.5));
		MatrixXd below = MatrixXd::Zero(this->winSize/2,this->signalSize/(this->winSize*0.5));
		MatrixXd result = MatrixXd::Zero(this->winSize,this->signalSize/(this->winSize*0.5));

		int j = 0;
		int k = 0;
		for(int i = 1; i < this->signalSize+1; i++){
			below(j,k) = i;
			j++;
			if(j == this->winSize/2){
				j = 0;
				k++;
				if(k < this->signalSize/(this->winSize*0.5))
					upper.col(k) = below.col(k-1);
			}
		}

		result.topRows(this->winSize/2)= upper;
		result.bottomRows(this->winSize/2) = below;
		return result;
	}
};

#endif /* VAD_H_ */
