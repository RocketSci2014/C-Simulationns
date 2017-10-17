#include "stdafx.h"
#pragma once
#include "randomcup.h"
#include "stopwatch.h"
#include "global.h"
#include "QD.h"
using namespace std;

class smr :
	public randomcup, public stopwatch
{
protected:
	short int FlowDirection;	
	bool 	
		pass,			//if a photon burst has passed
		first,			//determine the first peak, the value of the first peak equals to the first wss peak
		keep;			//keep recording photon timings into pulse, for the case of MakeWss1

	int	
		CountDownTillAfterpulse,
		CountToNextFlowEntry,
		index,
		NumPassive,			//number of passive reverse, automatic recycling
		mindex,				//index of the molecule that passes the laser beam
		peakindex;			//the index of the peak
   
	double ttime[NUM_TIMES];
	std::vector<double> Pulse;	 //for maximum likehood, the timing of each photon
	unsigned short int Weights[Max_Binned_Photons][900];
	std::vector<int> Wss_Array;	 //weighted sliding values, length sig6/BinTime
	std::vector<QD*> QD_Array;	 //molecules in the nanochannel
	std::vector<double>::size_type tt;

//for recycling function	
	double
		NextWssTime, //time to make a new weighted sliding sum
		sum_of_elements,
		peak,		//peak for ML estimation
		hpeak,		//peak for maximum wss value
		hwss,		//peak value of the wss
		tE,			//the expeted crossing time
		PM1,		//lower limit for the delay time between two consecutive passages
		PM2,		//upper limit for the invader molecule control algorithm 
		PM3;
		
public:
	smr(void);
	~smr(void);
	void OpenFiles(void);  //write files
	void CloseFiles(void);  //close files
	void MonteCarlo(void);	//Monte Carlo simulation
	void InitializeWss(void);//initialize arrays for weighted sliding sum
	void GenerateWeights(void);//generate weights
	unsigned short int WeightedSlidingSum(unsigned short int PhotonCounts);//calculate the wss value for each wss period and increment wss by photon counts
	unsigned short int PhotonCount(void); //count the number of photons in a binning period
	bool MakeWss(void); //make wss and determine if a photon burst has passed 3 sigma
	bool MakeWss1(void);//make wss and dterminne if a photon burst has passed T0+T~T0+3T
	void WritePhoton(double time);//write the timing of a photon into the file
	void NextPhoton(double PhotonTiming,double PriorPhotonTiming);//apply conditions of SPAD dead time and after pulse to the photon timing
	void WriteComment(string comment); //write comment to the file
	void FindPeakMethod(double range); //set the range of WSS
};