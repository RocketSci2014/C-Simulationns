#include "StdAfx.h"
#include "QD.h"
#include "global.h"
using namespace std;

QD::QD(double ttime2,int position1)
{
	InitializeTimeToNextExcitation1();	
	InitializeTimeToNextExcitation2();
	InitializeTimeToDecay3();
	InitializeTimeToDecayA();
	position=position1;
	state=S0;
	FindNextExcitation(ttime2);
	for(int ii=0;ii<5;ii++)
		decaypath[ii]=INFTY;
}


void QD::FindNextExcitation(double ttime1)
{
	if( position > -BEAM_SIZE && position< BEAM_SIZE)
	{		
		timing=ttime1+(Period-fmod( ttime1,Period ))+TimeToNextExcitation1( position); 	// Note: the random time is added to the present time tt[index] so the routine can be called no matter what process is being executed
		state=S1;
	}
	else
   		timing=INFTY;	
	return;
}


double QD::PhotoPhysics(void)
{
	double PhotonTiming=INFTY;
	switch (state){	 
				case S0: //groundstate of QD
					FindNextExcitation(timing);
				    break;
				case S1: //S1
					decaypath[0]=TimeToDecay1();
					if( position > -BEAM_SIZE && position< BEAM_SIZE)
					decaypath[4]=(Period-fmod( timing,Period ))+TimeToNextExcitation2( position);
					else
						decaypath[4]=INFTY;
					if(decaypath[0]<decaypath[4])
					{
						state=S0;
						timing+=decaypath[0];
						if (DetectPhoton()==true)
							PhotonTiming=timing+PhotonTimingError();//time to emit a photon
					}
					else
					{
						state=S2;
						timing+=decaypath[4];
					}
					break;
				case S2: //S2
					decaypath[1]=TimeToDecay2();
					decaypath[3]=TimeToDecayA(TauFactor());
					if (decaypath[1]<decaypath[3])
					{
						state=S1;
						timing+=decaypath[1];
						if (DetectPhoton()==true)
							PhotonTiming=timing+PhotonTimingError();//time to emit a photon
					}
					else
					{
						state=S3;
						timing+=decaypath[3];
					}
					break;
				case S3: //S3
					decaypath[2]=TimeToDecay3(TauFactor());
					state=S1;
					timing+=decaypath[2];
					break;
					}
	return PhotonTiming;
}