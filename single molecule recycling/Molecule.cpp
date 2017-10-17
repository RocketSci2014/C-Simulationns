#include "StdAfx.h"
#include "Molecule.h"
#include "global.h"

Molecule::Molecule(double ttime2,int position1,MoleculeType type1)
{
	InitializeTimeToNextExcitation();	
	position=position1;
	FindNextExcitation(ttime2);
	type=type1;
	state=GroundState;
}


Molecule::~Molecule(void)
{

}


void Molecule::FindNextExcitation(double ttime1)
{
	if( position > -BEAM_SIZE && position< BEAM_SIZE)
	{
		timing = ttime1 + ( TimeToNextExcitation1(position) );
		}
	else
   		timing=INFTY;	
	return;
}


double Molecule::PhotoPhysics(void)
{
	double PhotonTiming=INFTY;
	switch (state){	 
				case GroundState: // After GroundState, Molecule becomes excited; decide decay path that molecule should take
					switch( DecayPath() )
					{
					case SingletWithPhotonDue: 
						state=SingletWithPhotonDue;
						timing +=TimeToSingletDecay();
						break;
					case SingletWithNoPhoton: 
						state=SingletWithNoPhoton;
						timing +=TimeToSingletDecay();
						break;
					case TripletState: 
						state=TripletState;
						timing +=TimeToTripletDecay();
						break;
					case Bleached: 
						state=Bleached;
						timing=INFTY;
						break;
					}
					break;
				case SingletWithPhotonDue:
					PhotonTiming=timing+PhotonTimingError();//time to emit a photon
				case SingletWithNoPhoton:
				case TripletState:
					state=GroundState;
					FindNextExcitation(timing);//find next excitation based on the urrent timing
					break;
				case Bleached:
					std::cout<<"\n Error: Simulation time can't reach the bleaching time!"<<std::endl;
					break;
					}
	return PhotonTiming;
}