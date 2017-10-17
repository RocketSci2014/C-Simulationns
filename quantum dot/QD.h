#pragma once
#include "randomcup.h"
#include "global.h"
class QD :
	public randomcup
{
public:
    double timing;			//timing of the next photophysics event
	double decaypath[5];	//time of the upcoming photophysics events 0.decay1 1. decay2 2. decay3 3. decayA 4. excitation2
	QDState state;			//state of the molecule
	int position;			//position of the molecule in a nanochannel
	void FindNextExcitation(double ttime1);//find the timing of the next excitation after one photophysics event
	double PhotoPhysics(void);				//To generate the photophysics, ttime is the end time of the current loop 
	QD(double ttime2,int position1);		//generate QD and find the next excitation time of the molecule
	~QD(void);
};

