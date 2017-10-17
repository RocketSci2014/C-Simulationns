
//For prredefine, if use Wss1, the blend molecule can be applied; 

#pragma once
#ifndef GLOBAL_H
#define GLOBAL_H
#include <stdio.h>
#include <tchar.h>
#include<iostream>
#include<fstream>
#include <iomanip>
#include<vector>
#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include <stdio.h>
#include "mkl.h"
#include <time.h>
#include <ctime>
#include <sys/types.h>
#include <sys/timeb.h>
#include <string.h>
#include <string>
#include <direct.h>
#include <sstream>

/*define constants*/
/*define functionalities of the simulation*/
//#define WssRecord		1				//record all the wss values
#define Blend			1				//blend two different kinds of molecules into the same solution, use the self adaptive fitting algorithm and wss to discriminate different values of diffusion coefficients 
//#define Location		1				//record the locatioins of a single molecule when it comes into the laser beam
#define Mono			1				//only one molecule is in the recycling, this function is imcompatible with the funcction blend!!
#define Method1		1				//binning
//#define Method4		1				//+/- 3sig around the WSS Peak
//#define Method5			1				//average photon times at the situation of muliple crossings


/*define length of the simulation, parameters for smr*/
#define SimulationTime	    10         //time of simulation
#ifdef Method1
#define BinTime				1.0e-3			//the period for binning is 1 ms
#else
#define	BinTime				10.0e-6		      //the period to calculate weighted sliding sum
#endif
#define	ReversalDelay		30.0e-3          //the reversal delay in smr
#define MaxMissedDetects    4             //number of automatic recycles
#define Max_Binned_Photons	16				//the upper limit of the lookup table

/*define states in the Markov chain*/
#define NextDiffuseEntryL	0
#define NextDiffuseEntryR	1
#define NextDiffuse			2
#define NextFlow			3
#define NextBackground		4
#define SetFlow				5
#if defined Blend
#define NextDiffuseEntryL1	6
#define NextDiffuseEntryR1	7
#define NextDiffuse1		8
#define NUM_TIMES 9				//=MAX_MOLECULES + 8=dimension of array ttime[]
#else
#define NUM_TIMES 6					//=MAX_MOLECULES +6=dimension of array ttime[]
#endif

/*state of the molecule*/
enum MoleculeState {GroundState,SingletWithPhotonDue,SingletWithNoPhoton,TripletState,Bleached};  //state of a molecule

/*types of the molecules*/
enum MoleculeType {Type1,Type2};      //type of a molecule, in the simulation we can use a mixture of two species of molecules

/*parameters of nanochannel*/
#define L_grid  5000                  //Rochester paper-increase the length of the nanochannel
#define L_grid2 10001                  //Rochester paper
#define Delx 1.0e-8						//Grid spacing 
#define pi 3.14159265
#define Avagadro 6.02e23                   //Avagadro constant
#define ChannelWidth1 6.5e-7                //Rochester paper section of the channel
#define ChannelWidth2 4.5e-7                //width and height of the channel
#define P_Entry 0.763540130047191       //probability of entering the nanochannel
#define INFTY			9.9e30          //infinity

/* parameters of photophysics*/ 
#define BEAM_SIZE 150					//Number of grid points from center out to 3 sigma for either beam=1.5 * Waist / Delx; 144 or 150
#ifdef Mono
#define PROB_Bleached				0.0
#define PROB_TripletState			0.0
#else
#define PROB_Bleached				1.0e-6 //the times of passes can be 500; if 1.0e-6, passes gonna be 50
//#define PROB_Bleached				0.0
#define PROB_TripletState			1.0e-3 //the tipical value of triplet quenchinng
//#define PROB_TripletState			0.0
#endif
#define PROB_SingletWithPhotonDue	0.008 //the photon detectionn efficiency is 0.8%



/*define derived variables*/
const double
#ifdef Mono
	ConcentrationL=0.0,
	ConcentrationR=0.0,
#else
    ConcentrationL=10.0e-12,				// M		sample concentration 100 pM the left reservoir
	ConcentrationR=10.0e-12,
#endif
	C0L=1.0e3*Avagadro*ConcentrationL*Delx*ChannelWidth1*ChannelWidth2,		// dimensionless	sample concentration per grid point left reservoir
	C0R=1.0e3*Avagadro*ConcentrationR*Delx*ChannelWidth1*ChannelWidth2,		//right reservoir
	C1L=P_Entry * C0L,			// dimensionless	Probability of new molecule entry per diffusion time step left reservoir
	C1R=P_Entry * C0R,
    VFmax=3.4e-4,						// Rochester paper m/s		max electrokinetic flow speed 100v 3.4/4.0
	D=0.98e-10,
    Dt_Diffuse=Delx * Delx /(2.0 * D),	// s	    Time step for diffusion 
	Dt_Flow=Delx / VFmax,				// s		Time step for flow 
    Waist=1e-6,						// m		Beam waist for each beam 0.96/1.0
	ClockFrequency=1.0e8,               // for continuous laser, it is the trigering frequency of the piezo clock
	Power=98.0e-6,                     //Rochester paper  
	sigma_a=4.40e-20,					//Absorption cross section of Atto 532 m^2
	E_lambda=3.73e-19,				//J		Photon energy for  Atto 532, emision at 553 nm 
	SingletLifetime=3.8e-9,             //Rochester paper dsDNA13 singlet lifetime
	TripletLifetime=1.0e-6,				// s		Phosphorescence Lifetime 
	SPAD_Sigma=1.27e-10,				// s		Standard deviation of SPAD timing error 
	SPAD_Offset=3.0 * SPAD_Sigma,		// s		Offset of SPAD timing 
	MeanTimeBetweenBackground=0.00167,    //the background is justified into 1000/s
	ProbAfterpulse=0.005,		// dimensionless	Probability of an Afterpulse 
	SPADDeadTime=4.0e-8,				// s		SPAD dead time 
	MeanTimeToAfterpulse=1.0e-7,		 // s		Mean time until an afterpulse is 100 ns
	sig=1.5e-3,							//s		Non-ideal width for the WSS weights
	factor_sig=2.0/pow(sig,2),
	sig3=3*sig,
	sig6=6*sig,                         //for weighted sliding sum
	PM1=-1.0*ReversalDelay,				//lower limit for the delay time between two consecutive passages
	PM2=ReversalDelay;					//upper limit for the invader molecule control algorithm 

const unsigned short int
	Wss_Th=1750,
	Fluctuation=64,			//half the height of WSS peak
	UpperLimit=(unsigned short int)(14500.0*Power/98.0e-6);				//Peaks above UpperLimit indicate an Invader

const int
	 wss_sig=int(sig/BinTime),
	 wss_sig3=3*wss_sig,
	 wss_sig6=6*wss_sig;	

/*define parameters used in the blend of species*/
#ifdef Blend
const double
	ConcentrationL1=10.0e-12,
	ConcentrationR1=10.0e-12,
	C0L1=1.0e3*Avagadro*ConcentrationL1*Delx*ChannelWidth1*ChannelWidth2,		// dimensionless	sample concentration per grid point left reservoir
	C0R1=1.0e3*Avagadro*ConcentrationR1*Delx*ChannelWidth1*ChannelWidth2,		//right reservoir
	C2L=P_Entry*C0L1,
	C2R=P_Entry*C0R1,
	D1=1.3e-10,							//diffuson coefficient of dsDNA15mers molecules, used in blend molecule simulation
	Dt_Diffuse1=Delx * Delx /(2.0 * D1), // s	    Time step for diffusion, for the second kind of molecule 
	RatioL=ConcentrationL/(ConcentrationL+ConcentrationL1),// ratio between concentrations of specie one and specie two
	RatioR=ConcentrationR/(ConcentrationR+ConcentrationR1);
#endif

#endif