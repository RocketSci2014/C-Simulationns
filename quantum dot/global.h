
/*The simulation for quantum dots excited by 2PE and recycled in a nanochannel, parameters are defined according to the values from previous papers*/

#pragma once
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


/*define parameters of the simulation, parameters for wss and smr*/
#define SimulationTime	    10000         //time of simulation
#define	BinTime				10.0e-6		      //the period to calculate weighted sliding sum
#define	ReversalDelay		30.0e-3          //the reversal delay in smr
#define Wss_Th				1750			//600 for the width of WSS, by taking the maximum value of the background
#define Fluctuation			64			//half the height of WSS peak
#define Auto_Recycle        10             //number of automatic recycles
#define Max_Binned_Photons	16				//the upper limit of the lookup table



/*define states in the Markov chain*/
#define NextDiffuseEntry	0
#define NextDiffuse			1
#define NextFlow			2
#define NextBackground		3
#define SetFlow				4
#define NUM_TIMES			5				//=MAX_MOLECULES +6=dimension of array ttime[]

/*state of the molecule*/
enum QDState{S0,S1,S2,S3};//S0 is ground state, S1 is 1 e-h, S2 is 2 e-h, S3 is ionized

/*parameters of nanochannel*/
#define L_grid  10000                  //Rochester paper-increase the length of the nanochannel
#define L_grid2 20001                  //Rochester paper
#define Delx 1.0e-8						//Grid spacing 
#define pi 3.14159265
#define Avagadro 6.02e23                   //Avagadro constant
#define ChannelWidth1 6.5e-7                //Rochester paper section of the channel
#define ChannelWidth2 4.5e-7                //width and height of the channel
#define P_Entry 0.763540130047191       //probability of entering the nanochannel
#define INFTY			9.9e30          //infinity
#define BEAM_SIZE 26		//Number of grid points from center out to 3 sigma for either beam=1.5 * Waist / Delx; 

/*define derived variables and parameters for QD*/
const double
	tau1=4.0e-9,					//4 ns
	tau2=1.0e-9,					//1ns
	Efficiency=0.05,				//photon detection efficieny 5%, have to relist the parameters
	sigma1=7.9e-51,					//absorption cross section 2pe
	sigma2=2.0e-22,
	power=0.5e-3,
	waist=0.17e-6,					//obtained from fitting the ACF of Rhodamine B
	E_lambda=2.49e-19,				//photon energy for 800 nm wavelength photon
	Concentration=10e-12,			// M		sample concentration 10 pM
	C0=1.0e3*Avagadro*Concentration*Delx*ChannelWidth1*ChannelWidth2,		// dimensionless	sample concentration per grid point 
	C1=P_Entry * C0,				// dimensionless	Probability of new molecule entry per diffusion time step
    VFmax=11e-4,					// Rochester paper m/s		max electrokinetic flow speed 100v make tD=2tF
	D=2.4e-10,	                    //stokes-eintein relation, D=kT/6Pigr, k=1.38e-23 m^2.kg.s_2.K^-1 g=8.9e-4 Pa.s r=2 nm
    Dt_Diffuse=Delx * Delx /(2.0 * D),	// s	    Time step for diffusion 
	Dt_Flow=Delx / VFmax,				// s		Time step for flow 
	ClockFrequency=1.0e8,               // for continuous laser, it is the trigering frequency of the piezo clock
	SPAD_Sigma=1.27e-10,				// s		Standard deviation of SPAD timing error 
	SPAD_Offset=3.0 * SPAD_Sigma,		// s		Offset of SPAD timing 
	MeanTimeBetweenBackground=0.00167,    //the background is justified into 1000/s
	ProbAfterpulse=0.005,		// dimensionless	Probability of an Afterpulse 
	SPADDeadTime=4.0e-8,				// s		SPAD dead time 
	MeanTimeToAfterpulse=1.0e-7,		// s		Mean time until an afterpulse is 100 ns
	sig=waist/VFmax,					//s 1 sigma in the unit of time
	factor_sig=2.0/pow(sig,2),
	sig3=3*sig,
	sig6=6*sig,                         //for weighted sliding sum
	Period=13.2e-9;						//period of the pulse laser, frequency 76 MHz

const int
	 wss_sig=int(sig/BinTime),
	 wss_sig3=3*wss_sig,
	 wss_sig6=6*wss_sig;	