#include "stdafx.h"
#pragma once
#include "randomcup.h"
#include "stopwatch.h"
#include "global.h"
#include "Molecule.h"
using namespace std;

class smr :
	public randomcup, public stopwatch
{
protected:
	short int FlowDirection;	
	bool 
	Pass,		      //if a photon burst has passed
	WSS_Increasing;

	MoleculeType Mtype;

	int	
	CountDownTillAfterpulse,
	CountToNextFlowEntryL,
	CountToNextFlowEntryR,
	index,
	MissedDetections,              //number of passive reverse, automatic recycling
	number_of_photons,//the number of photons in 1 ms
	total_number, //total number of photons
	RecycleCount,		//the index of the peak
	PeakPast; //if a peak has passed

	unsigned short int 
	PeakMax, //peak value of the wss
	WSS_LocalMax,
	WSS_LocalMin;
	std::vector<double> Pulse; //for maximum likehood, the timing of each photon
	//std::vector<unsigned short int> Weights;
	unsigned short int Weights[Max_Binned_Photons][900];
	std::vector<unsigned short int> Wss_Array; //weighted sliding values, length sig6/BinTime
	std::vector<Molecule*> Molecule_Array; //molecules in the nanochannel
	std::vector<double>::size_type tt;

//for recycling function	
	double
		NextBinTime,            //time to make a new weighted sliding sum
		Peak,          //peak for ML estimation
		MaxTime,        //peak for maximum wss value		
		tE,				//the expected time of detection	
		ttime[NUM_TIMES];
		

#ifdef Blend
	int NType1,
		NType2;
#endif
		
public:
	smr(void);
	~smr(void);
	void OpenFiles(void);  //write files
	void CloseFiles(void);  //close files
	void MonteCarlo(void);
	void InitializeWss(void);//initialize arrays for weighted sliding sum
	void GenerateWeights(void);//generate weights
	unsigned short int WeightedSlidingSum(unsigned short int PhotonCounts);//calculate the wss value for each wss period and increment wss by photon counts
	unsigned short int PhotonCount(void); //count the number of photons in a binning period
	void MakeWss(void); //make wss and determine if a photon burst has passed 3 sigma
	void MakeWss1(void);//use method 5 to make WSS
	void WritePhoton(double time);//write the timing of a photon into the file
	void NextPhoton(double PhotonTiming,double PriorPhotonTiming);//apply conditions of SPAD dead time and after pulse to the photon timing
	void WriteComment(string comment);
	void FindPeakMethod(double range);
	void FindPeakMethod(void);
	void Binning(void);
};

