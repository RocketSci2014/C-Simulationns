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
  this->Random_TimeToNextDiffuseEntryL=new double[NUM];
  this->Random_TimeToNextDiffuseEntryR=new double[NUM];
  /*for the blend of species*/
#if defined Blend
  this->Random_TimeToNextDiffuseEntryL1=new double[NUM];
  this->Random_TimeToNextDiffuseEntryR1=new double[NUM];
  this->J_TimeToNextDiffuseEntryL1=NUM1;
  this->J_TimeToNextDiffuseEntryR1=NUM1;
  this->Random_PositionOfNextDiffuseEntryL1=new short int [NUM];
  this->Random_PositionOfNextDiffuseEntryR1=new short int [NUM];
  this->J_PositionOfNextDiffuseEntryL1=NUM1;
  this->J_PositionOfNextDiffuseEntryR1=NUM1;
  this->Random_Diffuse1=new short int[NUM];
  this->J_Diffuse1=NUM1;
  this->Random_TypeOfNextFlowEntryL=new MoleculeType[NUM];
  this->Random_TypeOfNextFlowEntryR=new MoleculeType[NUM];
#endif
  this->Random_TimeToNextFlowEntryL=new int[NUM];
  this->Random_TimeToNextFlowEntryR=new int[NUM];
  this->Random_Diffuse=new short int[NUM];
  this->Random_PositionOfNextDiffuseEntryL=new short int [NUM];
  this->Random_PositionOfNextDiffuseEntryR=new short int [NUM];
  this->J_Diffuse=NUM1;
  this->J_TimeToNextDiffuseEntryL=NUM1;
  this->J_TimeToNextDiffuseEntryR=NUM1;
  this->J_PositionOfNextDiffuseEntryL=NUM1;
  this->J_PositionOfNextDiffuseEntryR=NUM1;
  this->J_TimeToNextFlowEntryL=NUM1;
  this->J_TimeToNextFlowEntryR=NUM1;
  this->ir=new int [NUM];
  this->r=new double [NUM];
  this->rf=new float[NUM];
  this->Random_TimeToSingletDecay=new double [NUM];
  this->Random_TimeToTripletDecay=new double [NUM];
  this->Random_TimeToNextBackground=new double [NUM];
  this->Random_TimeToAfterpulse=new double [NUM];
  this->Random_PhotonTimingError=new double [NUM];
  this->Random_TimeToNextExcitation1=new double* [BEAM_SIZE];
  for (int i=0;i<BEAM_SIZE;i++)
		this->Random_TimeToNextExcitation1[i]=new double [NUM];
  this->Random_PhotonsUntilAnAfterpulse=new int [NUM];
  for(int jj=0;jj<BEAM_SIZE;jj++)		
		this->J_TimeToNextExcitation1[jj]=short int(NUM1);		      
		this->Random_DecayPath=new MoleculeState [NUM];
		this->J_TimeToSingletDecay=NUM1;
		this->J_TimeToTripletDecay=NUM1;
		this->J_TimeToNextBackground=NUM1;
		this->J_TimeToAfterpulse=NUM1;
		this->J_PhotonsUntilAnAfterpulse=NUM1;
		this->J_PhotonTimingError=NUM1;
		this->J_DecayPath=NUM1;
  for( int i=0; i<L_grid2; i++ ) {e1[i]=BEAM_SIZE-1;}//define e1
      e1[L_grid]=0;
  for(int i=1; i<BEAM_SIZE;i++)
	  e1[L_grid+i]=e1[L_grid-i]=i;
}

randomcup::~randomcup(void)
{
	delete [] Random_TimeToNextDiffuseEntryL;
	delete [] Random_TimeToNextDiffuseEntryR;
#if defined Blend
	delete [] Random_TimeToNextDiffuseEntryL1;
	delete [] Random_TimeToNextDiffuseEntryR1;
	delete [] Random_PositionOfNextDiffuseEntryL1;
	delete [] Random_PositionOfNextDiffuseEntryR1;
	delete [] Random_Diffuse1;
	delete [] Random_TypeOfNextFlowEntryL;
	delete [] Random_TypeOfNextFlowEntryR;
#endif
	delete [] Random_TimeToNextFlowEntryL;
	delete [] Random_TimeToNextFlowEntryR;
	delete [] Random_Diffuse;
	delete [] Random_PositionOfNextDiffuseEntryL;
	delete [] Random_PositionOfNextDiffuseEntryR;
	delete [] ir;
	delete [] r;
	delete [] rf;
	delete [] Random_TimeToSingletDecay;
	delete [] Random_TimeToTripletDecay;
	delete [] Random_TimeToNextBackground;
	delete [] Random_TimeToAfterpulse;
	delete [] Random_PhotonTimingError;
	 for (int i=0;i<BEAM_SIZE;i++)
		delete[] Random_TimeToNextExcitation1[i];
	delete [] Random_TimeToNextExcitation1;
	delete [] Random_PhotonsUntilAnAfterpulse;
	delete [] Random_DecayPath;
    vslDeleteStream(&stream1);
}


int randomcup::TimeToNextFlowEntryL( void )
{
   if(J_TimeToNextFlowEntryL==NUM1) 
   {
	   vslNewStream( &stream1, VSL_BRNG_MT19937, seed1 );
	   J_TimeToNextFlowEntryL=-1;
	if( C0L<1.0e-20 ) {
		for( int iz=0;iz<NUM;iz++ ) {
			Random_TimeToNextFlowEntryL[iz]=2147483647; }}//infinity
	else
	   viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF,stream1, NUM, Random_TimeToNextFlowEntryL, C0L );	/* Generating */
   }
   J_TimeToNextFlowEntryL++;
   return( Random_TimeToNextFlowEntryL[J_TimeToNextFlowEntryL]);
}

int randomcup::TimeToNextFlowEntryR( void )
{
   if(J_TimeToNextFlowEntryR==NUM1) 
   {
	   vslNewStream( &stream1, VSL_BRNG_MT19937, seed1 );
	   J_TimeToNextFlowEntryR=-1;
	if( C0R<1.0e-20 ) {
		for( int iz=0;iz<NUM;iz++ ) {
			Random_TimeToNextFlowEntryR[iz]=2147483647; }}//infinity
	else
	   viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF,stream1, NUM, Random_TimeToNextFlowEntryR, C0R );	/* Generating */
   }
   J_TimeToNextFlowEntryR++;
   return( Random_TimeToNextFlowEntryR[J_TimeToNextFlowEntryR]);
}

double randomcup::TimeToNextDiffuseEntryL( void )
{
	int
		iz;

   if(J_TimeToNextDiffuseEntryL==NUM1) 
   {
	   J_TimeToNextDiffuseEntryL=-1;
		if( C0L<1.0e-20 ) {
			for( int iz=0;iz<NUM;iz++ ) {
				Random_TimeToNextDiffuseEntryL[iz]=INFTY; }}//infinity
		else
		{
			viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, ir, C1L );	/* Generating */
			for(iz=0;iz<NUM;iz++) { Random_TimeToNextDiffuseEntryL[iz]=Dt_Diffuse*(double)ir[iz]; }
		}
   }
   J_TimeToNextDiffuseEntryL++;
   return( Random_TimeToNextDiffuseEntryL[J_TimeToNextDiffuseEntryL] );
}

double randomcup::TimeToNextDiffuseEntryR( void )
{
	int
		iz;

   if(J_TimeToNextDiffuseEntryR==NUM1) 
   {
	   J_TimeToNextDiffuseEntryR=-1;
		if( C0R<1.0e-20 ) {
			for( int iz=0;iz<NUM;iz++ ) {
				Random_TimeToNextDiffuseEntryR[iz]=INFTY; }}//infinity
		else
		{
			viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, ir, C1R );	/* Generating */
			for(iz=0;iz<NUM;iz++) { Random_TimeToNextDiffuseEntryR[iz]=Dt_Diffuse*(double)ir[iz]; }
		}
   }
   J_TimeToNextDiffuseEntryR++;
   return( Random_TimeToNextDiffuseEntryR[J_TimeToNextDiffuseEntryR] );
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
short int randomcup::PositionOfNextDiffuseEntryL( void )
{
	const float
		PB[6]={
0.808133199827390,
0.983117150015168,
0.999381740374242,
0.999991051158142,
0.999999950367008,
1.000000000000000
 };
	const  int
		IB[6]={-L_grid, 1-L_grid, 2-L_grid, 3-L_grid, 4-L_grid, 5-L_grid};
	short int
		i,
		j;

   if(J_PositionOfNextDiffuseEntryL==NUM1) 
   {
	   J_PositionOfNextDiffuseEntryL=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_PositionOfNextDiffuseEntryL[i]=IB[j];
	   }
   }
   J_PositionOfNextDiffuseEntryL++;
   return( Random_PositionOfNextDiffuseEntryL[J_PositionOfNextDiffuseEntryL] );
}

short int randomcup::PositionOfNextDiffuseEntryR( void )
{
	const float
		PB[6]={
0.808133199827390,
0.983117150015168,
0.999381740374242,
0.999991051158142,
0.999999950367008,
1.000000000000000
 };
	const  int
		IB[12]={L_grid, L_grid-1, L_grid-2, L_grid-3, L_grid-4, L_grid-5 };
	short int
		i,
		j;

   if(J_PositionOfNextDiffuseEntryR==NUM1) 
   {
	   J_PositionOfNextDiffuseEntryR=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_PositionOfNextDiffuseEntryR[i]=IB[j];
	   }
   }
   J_PositionOfNextDiffuseEntryR++;
   return( Random_PositionOfNextDiffuseEntryR[J_PositionOfNextDiffuseEntryR] );
}

MoleculeState randomcup::DecayPath( void )
{
	const float
		PB[4]={
               1.0-(PROB_SingletWithPhotonDue+PROB_TripletState+PROB_Bleached),
               1.0-(PROB_TripletState+PROB_Bleached),
               1.0-(PROB_Bleached),
               1.0
              };

	const MoleculeState
		IB[4]={SingletWithNoPhoton, SingletWithPhotonDue, TripletState, Bleached }; /* 0,1,2,3 */ 
	short int
		i,
		j;

   if(J_DecayPath==NUM1) 
   {
	   J_DecayPath=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_DecayPath[i]=IB[j];
	   }
   }
   J_DecayPath++;
   return( Random_DecayPath[J_DecayPath] );
}

double randomcup::TimeToSingletDecay( void )
{
   if(J_TimeToSingletDecay==NUM1) 
   {
	   J_TimeToSingletDecay=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToSingletDecay, 0.0, SingletLifetime );	/* Generating */
   }
   J_TimeToSingletDecay++;
   return( Random_TimeToSingletDecay[J_TimeToSingletDecay] );
}

double randomcup::TimeToTripletDecay( void )
{
   if(J_TimeToTripletDecay==NUM1) 
   {
	   J_TimeToTripletDecay=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToTripletDecay, 0.0, TripletLifetime );	/* Generating */
   }
   J_TimeToTripletDecay++;
   return( Random_TimeToTripletDecay[J_TimeToTripletDecay] );
}

double randomcup::TimeToNextBackground( void )
{
   if(J_TimeToNextBackground==NUM1) 
   {
	   J_TimeToNextBackground=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToNextBackground, 0.0, MeanTimeBetweenBackground );	/* Generating */
   }
   J_TimeToNextBackground++;
   return( Random_TimeToNextBackground[J_TimeToNextBackground] );
}

double randomcup::PhotonTimingError( void )
{
   if(J_PhotonTimingError==NUM1) 
   {
	   J_PhotonTimingError=-1;
	   vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream1, NUM, Random_PhotonTimingError, SPAD_Offset, SPAD_Sigma );				/* Generating */
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
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, PE1[grid_position] );	/* Generating */
	   for(int iz=0;iz<NUM;iz++) { Random_TimeToNextExcitation1[grid_position][iz]=r[iz]; }
   }
   J_TimeToNextExcitation1[grid_position]++;
   return( Random_TimeToNextExcitation1[grid_position][J_TimeToNextExcitation1[grid_position]] );
}


void randomcup::InitializeTimeToNextExcitation( void )
{
	double
		x,
		PE0;
	PE0=(E_lambda * pi * Waist * Waist)/(sigma_a  * Power ) ; //transitions per second continous pulse
	for(int ii=0;ii<BEAM_SIZE;ii++)
	{
		x=(Delx * (double)ii)/Waist;
		x *=x;
		x *=-2.0;
		PE1[ ii ]=double(PE0 / exp( x ));
		vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, r, 0.0, PE1[ii] );	/* Generating */
	    for(int iz=0;iz<NUM;iz++) { Random_TimeToNextExcitation1[ii][iz]=r[iz]; }
	    J_TimeToNextExcitation1[ii]=-1;
	}
   return;
}

double randomcup::TimeToAfterpulse( void )
{
   if(J_TimeToAfterpulse==NUM1) 
   {
	   J_TimeToAfterpulse=-1;
	   vdRngExponential( VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream1, NUM, Random_TimeToAfterpulse, 0.0, MeanTimeToAfterpulse );	/* Generating */
   }
   J_TimeToAfterpulse++;
   return( Random_TimeToAfterpulse[J_TimeToAfterpulse] );
}

int randomcup::PhotonsUntilAnAfterpulse( void )
{
	if(J_PhotonsUntilAnAfterpulse==NUM1) 
	{
		J_PhotonsUntilAnAfterpulse=-1;
		viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, Random_PhotonsUntilAnAfterpulse, ProbAfterpulse );	/* Generating */
	}
	J_PhotonsUntilAnAfterpulse++;
	return( Random_PhotonsUntilAnAfterpulse[J_PhotonsUntilAnAfterpulse] );
}

#if defined Blend
 MoleculeType randomcup::TypeOfNextFlowEntryL( void )
{
   if(J_MoleculeTypeL==NUM1) 
   {
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); //Generating the random numbers for the blend type of molecules
	   for (int ii=0;ii<NUM;ii++)
	   {
		   if (rf[ii]<RatioL)
			   Random_TypeOfNextFlowEntryL[ii]=Type1;
		   else
			   Random_TypeOfNextFlowEntryL[ii]=Type2;
	   }
   }
   J_MoleculeTypeL++;
   return Random_TypeOfNextFlowEntryL[J_MoleculeTypeL];
};

 MoleculeType randomcup::TypeOfNextFlowEntryR( void )
{
   if(J_MoleculeTypeR==NUM1) 
   {
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); //Generating the random numbers for the blend type of molecules
	   for (int ii=0;ii<NUM;ii++)
	   {
		   if (rf[ii]<RatioR)
			   Random_TypeOfNextFlowEntryR[ii]=Type1;
		   else
			   Random_TypeOfNextFlowEntryR[ii]=Type2;
	   }
   }
   J_MoleculeTypeR++;
   return Random_TypeOfNextFlowEntryR[J_MoleculeTypeR];
};

double randomcup::TimeToNextDiffuseEntryL1( void )
{
	int
		iz;

   if(J_TimeToNextDiffuseEntryL1==NUM1) 
   {
	   J_TimeToNextDiffuseEntryL1=-1;
	   viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, ir, C2L );	/* Generating */
	   for(iz=0;iz<NUM;iz++) { Random_TimeToNextDiffuseEntryL1[iz]=Dt_Diffuse1*(double)ir[iz]; }
   }
   J_TimeToNextDiffuseEntryL1++;
   return( Random_TimeToNextDiffuseEntryL1[J_TimeToNextDiffuseEntryL1] );
}

double randomcup::TimeToNextDiffuseEntryR1( void )
{
	int
		iz;

   if(J_TimeToNextDiffuseEntryR1==NUM1) 
   {
	   J_TimeToNextDiffuseEntryR1=-1;
	   viRngGeometric( VSL_RNG_METHOD_GEOMETRIC_ICDF, stream1, NUM, ir, C2R );	/* Generating */
	   for(iz=0;iz<NUM;iz++) { Random_TimeToNextDiffuseEntryR1[iz]=Dt_Diffuse1*(double)ir[iz]; }
   }
   J_TimeToNextDiffuseEntryR1++;
   return( Random_TimeToNextDiffuseEntryR1[J_TimeToNextDiffuseEntryR1] );
}

short int randomcup::Diffuse1(void)
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

   if(J_Diffuse1==NUM1) 
   {
	   J_Diffuse1=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_Diffuse1[i]=IB[j];
	   }
   }
   J_Diffuse1++;
   return( Random_Diffuse1[J_Diffuse1] );
}

short int randomcup::PositionOfNextDiffuseEntryL1( void )
{
	const float
		PB[6]={
0.808133199827390,
0.983117150015168,
0.999381740374242,
0.999991051158142,
0.999999950367008,
1.000000000000000
 };
	const  int
		IB[6]={-L_grid, 1-L_grid, 2-L_grid, 3-L_grid, 4-L_grid, 5-L_grid};
	short int
		i,
		j;

   if(J_PositionOfNextDiffuseEntryL1==NUM1) 
   {
	   J_PositionOfNextDiffuseEntryL1=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_PositionOfNextDiffuseEntryL1[i]=IB[j];
	   }
   }
   J_PositionOfNextDiffuseEntryL1++;
   return( Random_PositionOfNextDiffuseEntryL1[J_PositionOfNextDiffuseEntryL1] );
}

short int randomcup::PositionOfNextDiffuseEntryR1( void )
{
	const float
		PB[6]={
0.808133199827390,
0.983117150015168,
0.999381740374242,
0.999991051158142,
0.999999950367008,
1.000000000000000
 };
	const  int
		IB[12]={L_grid, L_grid-1, L_grid-2, L_grid-3, L_grid-4, L_grid-5 };
	short int
		i,
		j;

   if(J_PositionOfNextDiffuseEntryR1==NUM1) 
   {
	   J_PositionOfNextDiffuseEntryR1=-1;
	   vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream1, NUM, rf, 0.0, 1.0 ); /* Generating */
	   for( i=0;i<NUM;i++)
	   {
		   j=0;
		   while( rf[i]>PB[j] ) j++;
		   Random_PositionOfNextDiffuseEntryR1[i]=IB[j];
	   }
   }
   J_PositionOfNextDiffuseEntryR1++;
   return( Random_PositionOfNextDiffuseEntryR1[J_PositionOfNextDiffuseEntryR1] );
}

#endif