#include "stdafx.h"
#include "global.h"
#include "randomcup.h"
using namespace std;


const int seed1=10;						// the seed for generating random numbers 		
#define NUM 1024						/*			Number of Random Numbers to generate at each time */
#define NUM1 NUM-1

randomcup::randomcup(void)
{
	vslNewStream( &stream1, VSL_BRNG_MT19937, seed1 );//initialize stream
/*initialize arrays used in random number generation*/
  this->Random_TimeToNextDiffuseEntry=new double[NUM];
  this->Random_TimeToNextFlowEntry=new int[NUM];
  this->Random_Diffuse=new short int[NUM];
  this->Random_PositionOfNextDiffuseEntry=new short int [NUM];
  this->J_Diffuse=NUM1;
  this->J_TimeToNextDiffuseEntry=NUM1;
  this->J_PositionOfNextDiffuseEntry=NUM1;
  this->J_TimeToNextFlowEntry=NUM1;
  this->J_TimeToNextBackground=NUM1;
  this->J_TimeToAfterpulse=NUM1;
  this->J_PhotonsUntilAnAfterpulse=NUM1;
  this->J_PhotonTimingError=NUM1;

  /*photophysics parameters of the QD*/
  this->J_TimeToDecay1=NUM1;
  this->J_TimeToDecay2=NUM1;
  for (int i=0;i<400;i++)
  {
	  this->J_TimeToDecay3[i]=short int(NUM1);	
	  this->J_TimeToDecayA[i]=short int(NUM1);	
  }
  this->J_TauFactor=NUM1;
  this->J_DetectPhoton=NUM1;
  this->Random_TimeToDecay1=new double [NUM];
  this->Random_TimeToDecay2=new double [NUM];
  this->Random_TimeToDecay3=new double* [400];
  this->Random_TimeToDecayA=new double* [400];
  for (int i=0;i<400;i++)
  {
	  this->Random_TimeToDecay3[i]=new double [NUM];
	  this->Random_TimeToDecayA[i]=new double [NUM];
  }
  this->Random_TauFactor=new short int [NUM];
  this->Random_DetectPhoton=new bool [NUM];
    this->Random_TimeToNextExcitation1=new double* [BEAM_SIZE];
	this->Random_TimeToNextExcitation2=new double* [BEAM_SIZE];
    for (int i=0;i<BEAM_SIZE;i++)
	{
		this->Random_TimeToNextExcitation1[i]=new double [NUM];
		this->Random_TimeToNextExcitation2[i]=new double [NUM];
		this->J_TimeToNextExcitation1[i]=short int(NUM1);	
		this->J_TimeToNextExcitation2[i]=short int(NUM1);	
	}

  this->ir=new int [NUM];
  this->r=new double [NUM];
  this->rf=new float[NUM];
  this->Random_TimeToNextBackground=new double [NUM];
  this->Random_TimeToAfterpulse=new double [NUM];
  this->Random_PhotonTimingError=new double [NUM];
  this->Random_PhotonsUntilAnAfterpulse=new int [NUM];	      
  for( int i=0; i<L_grid2; i++ ) {e1[i]=BEAM_SIZE-1;}//define e1
      e1[L_grid]=0;
  for(int i=1; i<BEAM_SIZE;i++)
	  e1[L_grid+i]=e1[L_grid-i]=i;
}

randomcup::~randomcup(void)
{
	/*delete arrays*/
	delete [] Random_TimeToNextDiffuseEntry;
	delete [] Random_TimeToNextFlowEntry;
	delete [] Random_Diffuse;
	delete [] Random_PositionOfNextDiffuseEntry;
	delete [] ir;
	delete [] r;
	delete [] rf;
	delete [] Random_TimeToNextBackground;
	delete [] Random_TimeToAfterpulse;
	delete [] Random_PhotonTimingError;
	delete [] Random_PhotonsUntilAnAfterpulse;
	delete [] Random_TimeToDecay1;
	delete [] Random_TimeToDecay2;
	for (int i=0;i<400;i++)
	 {
		delete[] Random_TimeToDecay3[i];
		delete[] Random_TimeToDecayA[i];
	 }
	delete [] Random_TimeToDecay3;
	delete [] Random_TimeToDecayA;
	delete [] Random_TauFactor;
	delete [] Random_DetectPhoton;
	for (int i=0;i<BEAM_SIZE;i++)
	 {
		delete[] Random_TimeToNextExcitation1[i];
		delete[] Random_TimeToNextExcitation2[i];
	 }
	delete [] Random_TimeToNextExcitation1;
	delete [] Random_TimeToNextExcitation2;

    vslDeleteStream(&stream1);
}


int randomcup::TimeToNextFlowEntry( void )
{
   if(J_TimeToNextFlowEntry==NUM1) 
   {
	   vslNewStream( &stream1, VSL_BRNG_MT19937, seed1 );
	   J_TimeToNextFlowEntry=-1;
	if( C0<1.0e-20 ) {
		for( int iz=0;iz<NUM;iz++ ) {
			Random_TimeToNextFlowEntry[iz]=2147483647; }}//infinity
	else
	   viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF,stream1, NUM, Random_TimeToNextFlowEntry, C0 );	/* Generating */
   }
   J_TimeToNextFlowEntry++;
   return( Random_TimeToNextFlowEntry[J_TimeToNextFlowEntry]);
}


double randomcup::TimeToNextDiffuseEntry( void )
{
	int
		iz;

   if(J_TimeToNextDiffuseEntry==NUM1) 
   {
	   J_TimeToNextDiffuseEntry=-1;
		if( C0<1.0e-20 ) {
			for( int iz=0;iz<NUM;iz++ ) {
				Random_TimeToNextDiffuseEntry[iz]=INFTY; }}//infinity
		else
		{
			viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, ir, C1 );	/* Generating */
			for(iz=0;iz<NUM;iz++) { Random_TimeToNextDiffuseEntry[iz]=Dt_Diffuse*(double)ir[iz]; }
		}
   }
   J_TimeToNextDiffuseEntry++;
   return( Random_TimeToNextDiffuseEntry[J_TimeToNextDiffuseEntry] );
}


short int randomcup::Diffuse(void)
{
	const float
		PB[13]={
0.382924922548026,
0.624655260005155,
0.866385597462284,
0.926983133405366,
0.987580669348448,
0.993557705595188,
0.999534741841929,
0.999763973247840,
0.999993204653751,
0.999996583337313,
0.999999962020875,
0.999999980970278,
1.0 /*0.999999999919680 */ };
	const short int
		IB[13]={0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6 };
	short int
		i,
		j;

   if(J_Diffuse==NUM1) 
   {
	   J_Diffuse=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_Diffuse[i]=IB[j];
	   }
   }
   J_Diffuse++;
   return( Random_Diffuse[J_Diffuse] );
}
short int randomcup::PositionOfNextDiffuseEntry( void )
{
	const float
		PB[12]={
0.404066599711270,
0.808133199422541,
0.895625174513805,
0.983117149605069,
0.991249444832068,
0.999381740059067,
0.999686395503419,
0.999991050947771,
0.999995500604795,
0.999999950261819,
0.999999975130910,
1.000000000000000
 };
	const  int
		IB[12]={-L_grid, L_grid, 1-L_grid, L_grid-1, 2-L_grid, L_grid-2, 3-L_grid, L_grid-3, 4-L_grid, L_grid-4, 5-L_grid, L_grid-5 };
	short int
		i,
		j;

   if(J_PositionOfNextDiffuseEntry==NUM1) 
   {
	   J_PositionOfNextDiffuseEntry=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); // Generating 
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_PositionOfNextDiffuseEntry[i]=IB[j];
	   }
   }
   J_PositionOfNextDiffuseEntry++;
   return( Random_PositionOfNextDiffuseEntry[J_PositionOfNextDiffuseEntry] );
}


double randomcup::TimeToDecay1( void )
{
   if(J_TimeToDecay1==NUM1) 
   {
	   J_TimeToDecay1=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToDecay1, 0.0, tau1 );	// Generating 
   }
   J_TimeToDecay1++;
   return( Random_TimeToDecay1[J_TimeToDecay1] );
}

double randomcup::TimeToDecay2( void )
{
   if(J_TimeToDecay2==NUM1) 
   {
	   J_TimeToDecay2=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToDecay2, 0.0, tau2 );	// Generating 
   }
   J_TimeToDecay2++;
   return( Random_TimeToDecay2[J_TimeToDecay2] );
}

short int randomcup::TauFactor(void)
{
	if(J_TauFactor==NUM1) 
   {
	   J_TauFactor=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 4.0 ); //Generating the random numbers for the blend type of molecules
	   for (int ii=0;ii<NUM;ii++)
	   {
		  Random_TauFactor[ii]=short int (floor(rf[ii]/0.01));
	   }
	}
   J_TauFactor++;
   return Random_TauFactor[J_TauFactor];
}

void randomcup::InitializeTimeToDecay3(void)
{	
	for(int ii=0;ii<400;ii++)
	{
		tau[ii]=exp((ii*0.01-2.0)*2.3);//-2 -> 2 ln10=2.3
		vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, tau[ii] );	// Generating 
	    for(int iz=0;iz<NUM;iz++) { Random_TimeToDecay3[ii][iz]=r[iz]; }
	    J_TimeToDecay3[ii]=-1;
	}
	return;
}

void randomcup::InitializeTimeToDecayA(void)
{	
	for(int ii=0;ii<400;ii++)
	{
		tauA[ii]=2*exp((ii*0.01-8)*2.3);//-5 -> -1
		vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, tauA[ii] );	// Generating 
	    for(int iz=0;iz<NUM;iz++) { Random_TimeToDecayA[ii][iz]=r[iz]; }
	    J_TimeToDecayA[ii]=-1;
	}
	return;
}

double randomcup::TimeToDecay3( short int factor)
{
   if(J_TimeToDecay3[factor]==NUM1) 
   {
	   J_TimeToDecay3[factor]=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, tau[factor] );	// Generating 
	    for( int iz=0;iz<NUM;iz++ ) { 
		   Random_TimeToDecay3[factor][iz]=r[iz]; }
   }
   J_TimeToDecay3[factor]++;
   return( Random_TimeToDecay3[factor][J_TimeToDecay3[factor]] );
}

double randomcup::TimeToDecayA( short int factor)
{
   if(J_TimeToDecayA[factor]==NUM1) 
   {
	   J_TimeToDecayA[factor]=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, tauA[factor] );	/* Generating */
	    for( int iz=0;iz<NUM;iz++ ) { 
		   Random_TimeToDecayA[factor][iz]=r[iz]; }
   }
   J_TimeToDecayA[factor]++;
   return( Random_TimeToDecayA[factor][J_TimeToDecayA[factor]] );
}



double randomcup::TimeToNextBackground( void )
{
   if(J_TimeToNextBackground==NUM1) 
   {
	   J_TimeToNextBackground=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToNextBackground, 0.0, MeanTimeBetweenBackground );	// Generating 
   }
   J_TimeToNextBackground++;
   return( Random_TimeToNextBackground[J_TimeToNextBackground] );
}

double randomcup::PhotonTimingError( void )
{
   if(J_PhotonTimingError==NUM1) 
   {
	   J_PhotonTimingError=-1;
	   vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream1, NUM, Random_PhotonTimingError, SPAD_Offset, SPAD_Sigma );			// Generating 
   }
   J_PhotonTimingError++;
   return( Random_PhotonTimingError[J_PhotonTimingError] );
}

double randomcup::TimeToNextExcitation1(int position ) 	
{

    int grid_position=e1[position+L_grid];
	if(J_TimeToNextExcitation1[grid_position]==NUM1) 
   {
	   J_TimeToNextExcitation1[grid_position]=-1;
	   viRngGeometric( VSL_METHOD_IGEOMETRIC_ICDF,stream1,NUM,ir,PE1[grid_position] );
	   for( int iz=0;iz<NUM;iz++ ) { 
		   Random_TimeToNextExcitation1[grid_position][iz]=Period*(double)ir[iz]; } 
   }
   J_TimeToNextExcitation1[grid_position]++;
   return( Random_TimeToNextExcitation1[grid_position][J_TimeToNextExcitation1[grid_position]] );
}

double randomcup::TimeToNextExcitation2(int position ) 	
{

    int grid_position=e1[position+L_grid];
	if(J_TimeToNextExcitation2[grid_position]==NUM1) 
   {
	   J_TimeToNextExcitation2[grid_position]=-1;
	   viRngGeometric( VSL_METHOD_IGEOMETRIC_ICDF,stream1,NUM,ir,PE2[grid_position] );
	   for( int iz=0;iz<NUM;iz++ ) { 
		   Random_TimeToNextExcitation2[grid_position][iz]=Period*(double)ir[iz]; } 
   }
   J_TimeToNextExcitation2[grid_position]++;
   return( Random_TimeToNextExcitation2[grid_position][J_TimeToNextExcitation2[grid_position]] );
}

void randomcup::InitializeTimeToNextExcitation1( void )
{
	double
		x,
		PE0;
	PE0=power/(E_lambda * pi * waist * waist);
	PE0=PE0*PE0*sigma1*Period;
	for(int ii=0;ii<BEAM_SIZE;ii++ ) 
	{
		x=(Delx*(double)ii)/waist;
		x*=x;
		x*=-2.0;
		PE1[ii]=PE0*exp( x );
		viRngGeometric( VSL_METHOD_IGEOMETRIC_ICDF,stream1,NUM,ir,PE1[ii] );	
		for(int iz=0;iz<NUM;iz++ ) { 
			Random_TimeToNextExcitation1[ii][iz]=Period*(double)ir[iz]; }
		J_TimeToNextExcitation1[ii]=-1;
	}
   return;
}


void randomcup::InitializeTimeToNextExcitation2( void )
{
	double
		x,
		PE0;
	PE0=(sigma2*power*Period)/(E_lambda*pi*waist*waist) ; //transitions per laser pulse
	for(int ii=0;ii<BEAM_SIZE;ii++ ) 
	{
		x=(Delx*(double)ii)/waist;
		x*=x;
		x*=-2.0;
		PE2[ii]=PE0*exp( x );
		viRngGeometric( VSL_METHOD_IGEOMETRIC_ICDF,stream1,NUM,ir,PE2[ii] );	
		for(int iz=0;iz<NUM;iz++ ) { 
			Random_TimeToNextExcitation2[ii][iz]=Period*(double)ir[iz]; }
		J_TimeToNextExcitation2[ii]=-1;
	}
   return;
}

 bool randomcup::DetectPhoton( void )
{
   if(J_DetectPhoton==NUM1) 
   {
	   J_DetectPhoton=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); //Generating the random numbers for the blend type of molecules
	   for (int ii=0;ii<NUM;ii++)
	   {
		   if (rf[ii]<Efficiency)
			   Random_DetectPhoton[ii]=true;
		   else
			   Random_DetectPhoton[ii]=false;
	   }
   }
   J_DetectPhoton++;
   return Random_DetectPhoton[J_DetectPhoton];
};

double randomcup::TimeToAfterpulse( void )
{
   if(J_TimeToAfterpulse==NUM1) 
   {
	   J_TimeToAfterpulse=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToAfterpulse, 0.0, MeanTimeToAfterpulse );	// Generating 
   }
   J_TimeToAfterpulse++;
   return( Random_TimeToAfterpulse[J_TimeToAfterpulse] );
}

int randomcup::PhotonsUntilAnAfterpulse( void )
{
	if(J_PhotonsUntilAnAfterpulse==NUM1) 
	{
		J_PhotonsUntilAnAfterpulse=-1;
		viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, Random_PhotonsUntilAnAfterpulse, ProbAfterpulse );	// Generating 
	}
	J_PhotonsUntilAnAfterpulse++;
	return( Random_PhotonsUntilAnAfterpulse[J_PhotonsUntilAnAfterpulse] );
}