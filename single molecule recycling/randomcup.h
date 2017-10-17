/*====Generate random numbers for the Monte-Carlo simulation, for blending of two kinds of molecules, the generation functions for diffusion are all different===*/	
#include "stdafx.h"
#pragma once
#ifndef RANDOMCUP_H
#define RANDOMCUP_H

class randomcup
{
protected:
VSLStreamStatePtr stream1;      //the stream for generating random numbers
/*arrays for generating random steps*/
int* ir;//time to next diffusion
double* r;//the case of continous laser
float* rf;//random step of a molecule
double* Random_TimeToNextDiffuseEntryL;
double* Random_TimeToNextDiffuseEntryR;
int* Random_TimeToNextFlowEntryL;
int* Random_TimeToNextFlowEntryR;
short int* Random_Diffuse;//steps of diffusing
short int*	Random_PositionOfNextDiffuseEntryL;
short int*	Random_PositionOfNextDiffuseEntryR;
short int J_Diffuse;
short int J_TimeToNextDiffuseEntryL;
short int J_TimeToNextDiffuseEntryR;
short int J_PositionOfNextDiffuseEntryL;
short int J_PositionOfNextDiffuseEntryR;
short int J_TimeToNextFlowEntryL;
short int J_TimeToNextFlowEntryR;
double* Random_TimeToSingletDecay;
double* Random_TimeToTripletDecay;
double* Random_TimeToNextBackground;
double* Random_TimeToAfterpulse;
double* Random_PhotonTimingError;
double** Random_TimeToNextExcitation1;
int* Random_PhotonsUntilAnAfterpulse;
short int J_TimeToNextExcitation1[BEAM_SIZE];
MoleculeState* Random_DecayPath;
short int J_TimeToSingletDecay;
short int J_TimeToTripletDecay;
short int J_TimeToNextBackground;
short int J_TimeToAfterpulse;
short int J_PhotonsUntilAnAfterpulse;
short int J_PhotonTimingError;
short int J_DecayPath;
double PE1[BEAM_SIZE];			/*  dimensionless	Probability of excitation per laser pulse for either laser at each grid point offset from center*/
#ifdef Blend
short int* Random_Diffuse1;
short int*	Random_PositionOfNextDiffuseEntryL1;
short int*	Random_PositionOfNextDiffuseEntryR1;
double* Random_TimeToNextDiffuseEntryL1;
double* Random_TimeToNextDiffuseEntryR1;
MoleculeType* Random_TypeOfNextFlowEntryL;
MoleculeType* Random_TypeOfNextFlowEntryR;
short int J_Diffuse1;
short int J_PositionOfNextDiffuseEntryL1;
short int J_PositionOfNextDiffuseEntryR1;
short int J_TimeToNextDiffuseEntryL1;
short int J_TimeToNextDiffuseEntryR1;
short int J_MoleculeTypeL;
short int J_MoleculeTypeR;

#endif
int e1[L_grid2];					/* used to map grid coordinates to array positions for molecule excitation for beam 1 */ 

public:
	randomcup(void);
	~randomcup(void);
	int TimeToNextFlowEntryL( void );
	int TimeToNextFlowEntryR( void );
	double TimeToNextDiffuseEntryL( void );
	double TimeToNextDiffuseEntryR( void );
	short int Diffuse(void);
	short int PositionOfNextDiffuseEntryL( void );
	short int PositionOfNextDiffuseEntryR( void );
	MoleculeState DecayPath(void);
	double TimeToSingletDecay(void);
	double TimeToTripletDecay(void);
	double TimeToNextBackground(void);
	double PhotonTimingError(void);
	double TimeToNextExcitation1(int);
	void InitializeTimeToNextExcitation(void);
	double TimeToAfterpulse(void);
	int PhotonsUntilAnAfterpulse( void );

#if defined Blend
	double TimeToNextDiffuseEntryL1( void );
	double TimeToNextDiffuseEntryR1( void );
	short int Diffuse1(void);
	short int PositionOfNextDiffuseEntryL1( void );
	short int PositionOfNextDiffuseEntryR1( void );
	MoleculeType TypeOfNextFlowEntryL( void );
	MoleculeType TypeOfNextFlowEntryR( void );
#endif
	
};

#endif