/*====Generate random numbers for the Monte-Carlo simulation, for blending of two kinds of molecules, the generation functions for diffusion are all different===*/	
#include "stdafx.h"
#include "global.h"
#pragma once

class randomcup
{
protected:
VSLStreamStatePtr stream1;      //the stream for generating random numbers

/*arrays for generating random steps*/
int* ir;						//time to next diffusion
double* r;						//the case of continous laser
float* rf;						//random step of a molecule
double* Random_TimeToNextDiffuseEntry;
int* Random_TimeToNextFlowEntry;
short int* Random_Diffuse;		//steps of diffusing
short int*	Random_PositionOfNextDiffuseEntry;
short int J_Diffuse;
short int J_TimeToNextDiffuseEntry;
short int J_PositionOfNextDiffuseEntry;
short int J_TimeToNextFlowEntry;
double* Random_TimeToNextBackground;
double* Random_TimeToAfterpulse;
double* Random_PhotonTimingError;
int* Random_PhotonsUntilAnAfterpulse;
short int J_TimeToNextBackground;
short int J_TimeToAfterpulse;
short int J_PhotonsUntilAnAfterpulse;
short int J_PhotonTimingError;

/*QD*/
double** Random_TimeToNextExcitation1;
double** Random_TimeToNextExcitation2;
short int J_TimeToNextExcitation1[BEAM_SIZE];
short int J_TimeToNextExcitation2[BEAM_SIZE];
double PE1[BEAM_SIZE];	//excitation from S0 to S1
double PE2[BEAM_SIZE];	//excitation from S1 to S2
double* Random_TimeToDecay1;
double* Random_TimeToDecay2;
double** Random_TimeToDecay3;
double** Random_TimeToDecayA;
short int* Random_TauFactor;
bool* Random_DetectPhoton;
short int J_TimeToDecay1;
short int J_TimeToDecay2;
short int J_TimeToDecay3[400];//400 is for the equally distributed tau and tauA
short int J_TimeToDecayA[400];
short int J_TauFactor;
short int J_DetectPhoton;
int e1[L_grid2];			//used to map grid coordinates to array positions for molecule excitation for beam 1

public:
	randomcup(void);
	~randomcup(void);
	int TimeToNextFlowEntry( void );
	double TimeToNextDiffuseEntry( void );
	short int Diffuse(void);
	short int PositionOfNextDiffuseEntry( void );
	double TimeToNextBackground(void);
	double PhotonTimingError(void);
	double TimeToAfterpulse(void);
	int PhotonsUntilAnAfterpulse( void );

	/*photophysics of QD*/
	double TimeToDecay1(void);		//S1
	double TimeToDecay2(void);		//S2
	double TimeToDecay3(short int);	//S3, input is the power factor 
	double TimeToDecayA(short int);	//angus ionization
	void InitializeTimeToDecay3(void);
	void InitializeTimeToDecayA(void);
	short int TauFactor(void);		//to generate a uniformly distributed factor
	double tau[400];				//decay from ionized state
	double tauA[400];				//Auger ionization
	double TimeToNextExcitation1(int);
	double TimeToNextExcitation2(int);
	void InitializeTimeToNextExcitation1(void);
	void InitializeTimeToNextExcitation2(void);
	bool DetectPhoton(void);		//detetor's photon detection efficiency	
};