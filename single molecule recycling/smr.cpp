#include "stdafx.h"
#include "global.h"
#include "smr.h"
using namespace std;

/* used by idamin() */ 
const int n=NUM_TIMES;
const int incx=1;
ofstream PhotonData;//timings of photons
ofstream Comment;
ofstream MoleculeData;//type, wss Peak for each molecule

#ifdef Location
ofstream Locations;
#endif

#ifdef WssRecord
ofstream WssRecords;
#endif


smr::smr(void)
{   
	smr::OpenFiles();
	InitializeWss();
	FlowDirection=1;   //flow direction is positive
   for(int i=0; i<NUM_TIMES; i++ ) {ttime[i]=INFTY;}
   ttime[NextFlow]=Dt_Flow;
   CountDownTillAfterpulse=PhotonsUntilAnAfterpulse();
   ttime[NextBackground]=TimeToNextBackground(); 
   ttime[NextDiffuseEntryL]=TimeToNextDiffuseEntryL();
   ttime[NextDiffuseEntryR]=TimeToNextDiffuseEntryR();
   CountToNextFlowEntryL=TimeToNextFlowEntryL();
   CountToNextFlowEntryR=TimeToNextFlowEntryR();
#ifdef Mono
   	ttime[NextDiffuse]=Dt_Diffuse;
#endif
#if defined Blend
   ttime[NextDiffuseEntryL1]=TimeToNextDiffuseEntryL1();
   ttime[NextDiffuseEntryR1]=TimeToNextDiffuseEntryR1();
#endif
#ifdef Method1
   	NextBinTime=0.0;
#else
   NextBinTime=BinTime;//the starting point of wss is 1xBinTime
#endif
   index=idamin(&n, ttime, &incx);
   index--;
	number_of_photons=0;
	total_number=0;
	Peak=0.0;
	MaxTime=0.0;  
	tE=0.0;
	PeakMax=0;         
	Pass=false;
#ifdef Method5
	MissedDetections=1;
#else
	MissedDetections=0;
#endif
	PeakPast=0;
	RecycleCount=0;	
	WSS_LocalMax=0;
	WSS_LocalMin=65535;
	WSS_Increasing=true;
#ifdef Blend
	NType1=0;
	NType2=0;
#endif
}


smr::~smr(void)
{
	WriteComment;//write comment
	CloseFiles();
	Wss_Array.clear();
}

void smr::WriteComment(string comment)
{
	Comment<<comment;
}

void smr::OpenFiles(void)
{   	
	std::string directory=stopwatch::currentDateTime();
    std::string root="C:/Users/bwang/Documents/Bo Wang/simulation codes/single molecule recycling-4/data/";
	directory=root+directory;
	mkdir(directory.c_str());
	PhotonData.open(directory+"/PhotonData.bin",std::ios::binary);
	Comment.open(directory+"/Comment.txt");
    MoleculeData.open(directory+"/MoleculeData.txt");

#ifdef Location
   Locations.open(directory+"/Locations.txt");
#endif

#ifdef WssRecord
   WssRecords.open(directory+"/WssRecords.txt");
#endif

}

void smr::CloseFiles(void)
{
	PhotonData.close();
	Comment.close();
	MoleculeData.close();

#ifdef Location
	Locations.close();
#endif

#ifdef WssRecord
	WssRecords.close();
#endif

}

void smr::InitializeWss(void)
{
	Wss_Array.reserve(wss_sig6);
	GenerateWeights();
	for (int ii=0;ii<wss_sig6;ii++)
		Wss_Array.push_back(0.0);
}

void smr::GenerateWeights(void)
{
	double temp,temp1;
	for(int ii=0;ii<wss_sig6;ii++)
	{
		temp=ii-wss_sig3;
		temp=temp*temp;
		temp1=wss_sig*wss_sig;
		Weights[0][ii]=unsigned short int(128*exp(temp*(-0.5/temp1))+0.5);//the peak of the weights is 128
			for(int jj=1;jj<Max_Binned_Photons;jj++)
				{
					Weights[jj][ii]=(jj+1)*Weights[0][ii];
				}
	}

}

unsigned short int smr::WeightedSlidingSum(unsigned short int PhotonCounts)
{
	switch(PhotonCounts)
	{
	    case 0:
		break;
		default:
			if (PhotonCounts<=Max_Binned_Photons)
				for (int ii=0;ii<wss_sig6;ii++)
				    Wss_Array[ii]+=Weights[PhotonCounts-1][ii];
			else
				for (int ii=0;ii<wss_sig6;ii++)
				    Wss_Array[ii]+=PhotonCounts*Weights[0][ii];
		break;
	}
	//calculate the wss value by taking the last element in the array
	unsigned short int WssValue=Wss_Array[0];
	Wss_Array.erase(Wss_Array.begin());
	Wss_Array.push_back(0);
    return WssValue;
}

unsigned short int smr::PhotonCount(void)//assign pulse to the beginning of pulse when starting the Monte-Carlo program
{
	unsigned short int PhotonCounts=0;	
	if (!Pulse.empty())
		tt++;//make the array counter at the right place
		while (tt<Pulse.size()&&Pulse[tt]<NextBinTime)
		{
			PhotonCounts++;
			tt++;
		}
		tt--;
	return PhotonCounts;
}


void smr::FindPeakMethod(double range)
{
	int ii=0;
	int jj=0;
	double sum_of_elements=0.0;
	while(ii<Pulse.size()&&Pulse[ii]<MaxTime-sig3+range)//photons should be within the range of hpeak-range tp hpeak+range
	{
		if(Pulse[ii]>MaxTime-sig3-range)
		{
			sum_of_elements+=Pulse[ii];
			jj++;
		}
		ii++;
	}
	Peak=sum_of_elements/jj;//calculate the ML estimated peak from the range, which is 3 sigma in this simulation. 
	//The algorithm in this simulation doesn't look for a local minimum as the end of the burst, instead it directly averages the photon times within the plus and minus 3 sigma range around the Peak of WSS.
}

void smr::FindPeakMethod(void)
{
	int ii=0;
	double sum_of_elements=0.0;
	while(ii<Pulse.size()&&Pulse[ii]<MaxTime)//photons should be within the range of hpeak-range tp hpeak+range
	{
		sum_of_elements+=Pulse[ii];
		ii++;
	}
	Peak=sum_of_elements/ii;//this algorithm pairs with the MakeWss1
}

void smr::Binning(void)
{
	double delay;
	if (number_of_photons>7)
	{
		if(Pass==false)
		{
			total_number+=number_of_photons;
			if (total_number>40)
			{
				Peak=NextBinTime;
				Pass=true;
				total_number=0;
			}	
		}
	}
	else
	{
		if (Pass==true)
			{
				Pass=false;
				delay=Peak-tE;
				if (delay>PM1&&delay<PM2||RecycleCount==0)
				{
					MissedDetections=0;
					ttime[SetFlow]=Peak+ReversalDelay; //also set values of Peak hpeak and hwss to 0 in the monte-carlo function
					tE=Peak+2*ReversalDelay;	
				}
				else
				{	
					MissedDetections=MaxMissedDetects+1;
					ttime[SetFlow]=Peak; //terminate the recyle when finding a invader molecule	

				}
			}
	}
		number_of_photons=0;
		NextBinTime+=BinTime;
}

void smr::MakeWss(void)
{
	unsigned short int wss=WeightedSlidingSum(PhotonCount());
	double delay;
#ifdef WssRecord
	WssRecords<<wss<<"\n";
#endif
	if(!Pulse.empty())
	{
		if (wss<=Wss_Th-Fluctuation)
		{
			if (Pass==true)
			{
				FindPeakMethod(3.0*sig);//try the 3.0 sigma range
				Pass=false;
				delay=Peak-tE;
				if (delay>PM1&&delay<PM2||RecycleCount==0)
				{	
					MissedDetections=0;
					ttime[SetFlow]=Peak+ReversalDelay; //also set values of Peak hpeak and hwss to 0 in the monte-carlo function
					tE=Peak+2*ReversalDelay;					
				}
				else
				{	
					MissedDetections=MaxMissedDetects+1;
					ttime[SetFlow]=Peak; //terminate the recyle when finding a invader molecule							
				}
				PeakMax=0;
			}
		}
		else if (wss>Wss_Th)
		{
			if (Pass==false)
			{
				while (!Pulse.empty()&&Pulse[0]<NextBinTime-sig6)
				{
				   Pulse.erase(Pulse.begin());
				   tt--;
				}
				  Pass=true;
			}

			if (wss>PeakMax)
			{
				PeakMax=wss;
				MaxTime=NextBinTime;
			}
		}
	}
	NextBinTime+=BinTime;
}

void smr::MakeWss1(void)
{
	unsigned short int WSS=WeightedSlidingSum(PhotonCount());
	if( WSS_Increasing ) {									// Alternately find peaks and valleys in WSS
		PeakPast=0;											// We are not yet waiting for a peak to become passed because the WSS is increasing
		if( WSS>WSS_LocalMax ) {
			MaxTime=NextBinTime;
			WSS_LocalMax=WSS; }
		else {
			if( WSS<WSS_LocalMax-Fluctuation ) {			// If the WSS is now decreasing...
				WSS_LocalMin=WSS;
				WSS_Increasing=false; } } }
	else { //WSS_Decreasing
		if( WSS<WSS_LocalMin ) {
			WSS_LocalMin=WSS; 
			if( NextBinTime>MaxTime+sig3 ) {	// If the time has stepped past the time of the WSS maximum plus 3 sigma_t
				PeakMax=WSS_LocalMax;
				PeakPast++; } }									// We are no longer waitig for the peak to become passed
		else {
			if( WSS>WSS_LocalMin+Fluctuation ) {				// If the WSS is now increasing... we can now utilize the last local maximum peak ...
				WSS_LocalMax=WSS;
				WSS_Increasing=true; 
				PeakMax=WSS_LocalMax;
				PeakPast++; }									// We are no longer waitig for the peak to become passed
			else {
				if( NextBinTime>MaxTime+sig3 ) {
					PeakMax=WSS_LocalMax;
					PeakPast++; } } } }
			if( PeakPast==1&&PeakMax>Wss_Th ) {						// Was there a new WSS peak >threshold ?
				if( PeakMax>UpperLimit ) {									//	If the peak was too large, it was from two molecules passing together, so there is an invader
					MissedDetections=MaxMissedDetects+1; 					// ...set MissedDetections to Max, and ... 
					ttime[SetFlow]=ttime[index];							// ...schedule to reload a new molecule
					 }  
				else {
					if( MissedDetections>0 ) {									// Was this the first peak since the flow was reversed ?
//						First Peak		
						while (!Pulse.empty()&&Pulse[0]<MaxTime-sig6)
						{
							Pulse.erase(Pulse.begin());
							tt--;
						}
						FindPeakMethod();
						MissedDetections=0; 
						ttime[SetFlow]=Peak+ReversalDelay; }
					else {
						//	Repeated peak
						FindPeakMethod();
						ttime[SetFlow]=Peak+ReversalDelay; } }
				}										// ...set the time for the next flow switch to now, so that the next flow switch will reload a new molecule

			NextBinTime+=BinTime;								// Finally, we advance the PhotonBinTime
}

void smr::WritePhoton(double time)
{
    unsigned long int ClockCount;
	ClockCount=(unsigned long int)fmod( time*ClockFrequency, 4294967296.0 );
	PhotonData.write((char*)&ClockCount,sizeof(ClockCount));
}

void smr::NextPhoton(double PhotonTiming,double PriorPhotonTiming)
{
	if(PhotonTiming>PriorPhotonTiming+SPADDeadTime)
	{
#ifdef Method1
		number_of_photons++;
#else
		Pulse.push_back(PhotonTiming);
#endif
		if (Pulse.size()==1)
			tt=0;
		if( --CountDownTillAfterpulse<=0 )  /*There is a possibility that this photon gives an afterpulse */ 
			{
				CountDownTillAfterpulse=PhotonsUntilAnAfterpulse();
#ifdef Method1				
				WritePhoton(PhotonTiming+TimeToAfterpulse());
#else
				Pulse.push_back(PhotonTiming+TimeToAfterpulse());
				WritePhoton(Pulse.back());
#endif

			}
		else
			WritePhoton(PhotonTiming);
    }
}

void smr::MonteCarlo(void)
{
	double CurrentTime;//current system time
	double PhotonTiming;//the timing of the current photon
	int itemp;//to store the temperary value of diffuse steps
	double PriorPhotonTiming=0;//to store the timing of the last photon, initial value is set to 0
	Molecule* NewMolecule; //temperary storage for the new molecule
#ifdef Mono
	NewMolecule=new Molecule(ttime[index],0,Type1);//place a molecule at the center of the nanochannel at the beginning
	Molecule_Array.push_back(NewMolecule);
#endif
	while (ttime[index]<SimulationTime)
	{
			/*do photophysics for each molecule*/
		if (!Molecule_Array.empty())
		{
		   for (int ii=0;ii<Molecule_Array.size();ii++)
		   {
			   CurrentTime=Molecule_Array[ii]->timing;
			   while (CurrentTime<ttime[index])	
			   {
			       PhotonTiming=Molecule_Array[ii]->PhotoPhysics();//update current timing of photophysics
				   if(PhotonTiming<INFTY)//for each molecule, do photophysics until exceed the system time
				   {
						NextPhoton(PhotonTiming,PriorPhotonTiming);
						PriorPhotonTiming=PhotonTiming;
				   }
					CurrentTime=Molecule_Array[ii]->timing;//update current timing of photophysics
			   }
		   }
		}
		/*make WSS for the upcoming photons*/

		while (NextBinTime<=ttime[index])//make wss until the end of the iteration
#ifdef Method1
			Binning();
#elif defined Method4
			MakeWss();
#elif defined Method5
			MakeWss1();
#endif
		  /*do Monte-Carlo simulation, switch between Next Diffuse entry,next diffuse, next flow, set flow, and background*/
		   switch(index)
		   {
			      case NextDiffuseEntryL:
					  NewMolecule=new Molecule(ttime[index],PositionOfNextDiffuseEntryL( ),Type1);
					  Molecule_Array.push_back(NewMolecule);//push back a new molecule into the array	
#if defined Blend
					NType1++;
					if( NType1==1 ) 
					   {// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
						   ttime[NextDiffuse]=ttime[NextDiffuseEntryL] + Dt_Diffuse;	
					   }
#else 
					  if( Molecule_Array.size()==1 ) 
					   {
						   ttime[NextDiffuse]=ttime[NextDiffuseEntryL] + Dt_Diffuse;	// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
					   }
#endif
					   ttime[NextDiffuseEntryL]	+=TimeToNextDiffuseEntryL();
					   break;
#if defined Blend
		          case NextDiffuseEntryL1:
			   	      NewMolecule=new Molecule(ttime[index],PositionOfNextDiffuseEntryL1( ),Type2);
					  Molecule_Array.push_back(NewMolecule);//push back a new molecule into the array	
				      NType2++;
					if( NType2==1 ) 
					{// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
						ttime[NextDiffuse1]=ttime[NextDiffuseEntryL1] + Dt_Diffuse1;	
					}
				   ttime[NextDiffuseEntryL1]+=TimeToNextDiffuseEntryL1();
				   break;
#endif
				   case NextDiffuseEntryR:
					  NewMolecule=new Molecule(ttime[index],PositionOfNextDiffuseEntryR( ),Type1);
					  Molecule_Array.push_back(NewMolecule);//push back a new molecule into the array	
#if defined Blend
					NType1++;
					if( NType1==1 ) 
					   {// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
						   ttime[NextDiffuse]=ttime[NextDiffuseEntryR] + Dt_Diffuse;	
					   }
#else 
					  if( Molecule_Array.size()==1 ) 
					   {
						   ttime[NextDiffuse]=ttime[NextDiffuseEntryR] + Dt_Diffuse;	// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
					   }
#endif
					   ttime[NextDiffuseEntryR]	+=TimeToNextDiffuseEntryR();
					   break;
#if defined Blend
		          case NextDiffuseEntryR1:
			   	      NewMolecule=new Molecule(ttime[index],PositionOfNextDiffuseEntryR1( ),Type2);
					  Molecule_Array.push_back(NewMolecule);//push back a new molecule into the array	
				      NType2++;
					if( NType2==1 ) 
					{// If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
						ttime[NextDiffuse1]=ttime[NextDiffuseEntryR1] + Dt_Diffuse1;	
					}
				   ttime[NextDiffuseEntryR1]+=TimeToNextDiffuseEntryR1();
				   break;
#endif
				   case NextDiffuse:	
					   if (!Molecule_Array.empty())
					   {
						   int jj=0;
						   while (jj<Molecule_Array.size())
						   {
#ifdef Blend
							   if(Molecule_Array[jj]->type==Type1) 
								   itemp=Diffuse();
							   else
								   itemp=0;
#else 
							   itemp=Diffuse();
#endif
							   if( itemp !=0 )
							   {
								   if(( Molecule_Array[jj]->position<0 &&  Molecule_Array[jj]->position+itemp >=0)||( Molecule_Array[jj]->position>0 &&  Molecule_Array[jj]->position+itemp <=0))
									{ 	//determine the type of the passed by molecule	
										Mtype=Type1; 
									}
								   Molecule_Array[jj]->position+=itemp;
								   if( Molecule_Array[jj]->state==GroundState )				   
								      Molecule_Array[jj]->FindNextExcitation(ttime[index]);
								   if( Molecule_Array[jj]->position>L_grid || Molecule_Array[jj]->position<-L_grid )
								   {
#ifdef Mono
									   Molecule_Array[jj]->position=0;
									   Molecule_Array[jj]->state=GroundState;
									   Molecule_Array[jj]->FindNextExcitation(ttime[index]);
									   jj++;
#else
									   delete Molecule_Array[jj];
									   Molecule_Array.erase(Molecule_Array.begin()+jj);//remove the molecule, and the length of the array is reduced by one
#ifdef Blend
									   NType1--;
#endif
#endif
								   }
								   else
									   jj++;
							   }
							   else
								   jj++;
						   }
#ifdef Blend
						   if(NType1==0)
							   ttime[NextDiffuse]=INFTY;
#else
						   if(Molecule_Array.empty())
							   ttime[NextDiffuse]=INFTY;
#endif
					   }
					   ttime[NextDiffuse]+=Dt_Diffuse;
					   break;
#ifdef Blend					   
	                case NextDiffuse1:	
					   if (!Molecule_Array.empty())
					   {
						   int jj=0;
						   while (jj<Molecule_Array.size())
						   {
							   if(Molecule_Array[jj]->type==Type2) 
								   itemp=Diffuse1();
							   else
								   itemp=0;
							   if( itemp !=0 )
							   {
								   if(( Molecule_Array[jj]->position<0 &&  Molecule_Array[jj]->position+itemp >=0)||( Molecule_Array[jj]->position>0 &&  Molecule_Array[jj]->position+itemp <=0))
							        { 	
										Mtype=Type2; //index of the molecule in the Molecule_Array
								    }
								   Molecule_Array[jj]->position+=itemp;
								   if( Molecule_Array[jj]->state==GroundState )				   
								      Molecule_Array[jj]->FindNextExcitation(ttime[index]);
								   if( Molecule_Array[jj]->position>L_grid || Molecule_Array[jj]->position<-L_grid )
								   {
									   delete Molecule_Array[jj];
									   Molecule_Array.erase(Molecule_Array.begin()+jj);//remove the molecule, and the length of the array is reduced by one
									   NType2--;
								   }
								   else
									   jj++;
							   }
							   else
								   jj++;
						   }

						   if(NType2==0)
							   ttime[NextDiffuse1]=INFTY;
					   }
					   ttime[NextDiffuse1]+=Dt_Diffuse1;
					   break;

#endif
					case NextFlow:			// 3 
							   switch( FlowDirection )
							   { // Note: case 0: should never occur as ttime[NextFlow] would equal INFTY
							   case 1:
								   if(!Molecule_Array.empty())
								   {
  									   int kk=0;
									   while (kk<Molecule_Array.size())
									   {
											Molecule_Array[kk]->position++;
#ifdef Location
											if(Molecule_Array[kk]->position>-5000&&Molecule_Array[kk]->position<5000)  //3 sigma range of laser beam, the molecule staying in recycle
												Locations<<Molecule_Array[kk]->position<<"\t"<<ttime[NextFlow]<<"\t"<<FlowDirection<<"\n";
#endif
											if( Molecule_Array[kk]->position==0 )
											{ 	
												Mtype=Molecule_Array[kk]->type; //index of the molecule in the Molecule_Array
											}
											if( Molecule_Array[kk]->state==GroundState )				   
												Molecule_Array[kk]->FindNextExcitation(ttime[index]);
											if( Molecule_Array[kk]->position>L_grid )
											{
#ifdef Mono
											   Molecule_Array[kk]->position=0;
											   Molecule_Array[kk]->state=GroundState;
											   Molecule_Array[kk]->FindNextExcitation(ttime[index]);
											   kk++;
#else
#ifdef Blend
											   if (Molecule_Array[kk]->type==Type1)
												   NType1--;
											   else
												   NType2--;
#endif
											   delete Molecule_Array[kk];
												Molecule_Array.erase(Molecule_Array.begin()+kk);//remove the molecule, and the length of the array is reduced by one
#endif
											}
											else
												kk++;
									   }	
#ifdef Blend
										if (NType1==0)			 
											ttime[NextDiffuse]=INFTY;   
										if (NType2==0)
											ttime[NextDiffuse1]=INFTY;
#else
										  if (Molecule_Array.empty())
											ttime[NextDiffuse]=INFTY;
#endif
								   }
								   if( CountToNextFlowEntryL==0 )
								   {									   									   
#ifdef Blend
										NewMolecule=new Molecule(ttime[index],-L_grid,TypeOfNextFlowEntryL());
										if(NewMolecule->type==Type1)
										{
											NType1++;
											if (NType1==1)
											ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;
										}
										else
										{
											NType2++;
											if(NType2==1)
											ttime[NextDiffuse1]=ttime[index] + Dt_Diffuse1;
										}
#else
									   NewMolecule=new Molecule(ttime[index],-L_grid,Type1);
									   if( Molecule_Array.empty())
									   {
										   ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;	
									   }
#endif
									   Molecule_Array.push_back(NewMolecule);
									   CountToNextFlowEntryL=TimeToNextFlowEntryL();
								   }
								   CountToNextFlowEntryL--;					/* If there is flow, count down the Flow time steps until the next flow entry */ 
								   break;
							   case -1:
   								   if(!Molecule_Array.empty())
								   {
  									   int dd=0;
									   while (dd<Molecule_Array.size())
									   {
											Molecule_Array[dd]->position--;
#ifdef Location
											if(Molecule_Array[dd]->position>-5000&&Molecule_Array[dd]->position<5000) //3 sigma range of laser beam
												Locations<<Molecule_Array[dd]->position<<"\t"<<ttime[NextFlow]<<"\t"<<FlowDirection<<"\n";
#endif
											if( Molecule_Array[dd]->position==0 )
											{ 	
												Mtype=Molecule_Array[dd]->type;; //index of the molecule in the Molecule_Array
											}
											if( Molecule_Array[dd]->state==GroundState )				   
												Molecule_Array[dd]->FindNextExcitation(ttime[index]);
											if( Molecule_Array[dd]->position<-L_grid )
											{
#ifdef Mono
											   Molecule_Array[dd]->position=0;
											   Molecule_Array[dd]->state=GroundState;
											   Molecule_Array[dd]->FindNextExcitation(ttime[index]);
											   dd++;
#else
#ifdef Blend
											   if (Molecule_Array[dd]->type==Type1)
												   NType1--;
											   else
												   NType2--;
#endif
											   delete Molecule_Array[dd];
											   Molecule_Array.erase(Molecule_Array.begin()+dd);//remove the molecule, and the length of the array is reduced by one
#endif
											}
											else
												dd++;
									   }	
#ifdef Blend
										if (NType1==0)			 
											ttime[NextDiffuse]=INFTY;   
										if (NType2==0)
											ttime[NextDiffuse1]=INFTY;
#else
										  if (Molecule_Array.empty())
											ttime[NextDiffuse]=INFTY;
#endif
								   }
								   if( CountToNextFlowEntryR==0 )
								   {

#ifdef Blend
										NewMolecule=new Molecule(ttime[index],L_grid,TypeOfNextFlowEntryR());
										if(NewMolecule->type==Type1)
										{
											NType1++;
											if (NType1==1)
											ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;
										}
										else
										{
											NType2++;
											if(NType2==1)
											ttime[NextDiffuse1]=ttime[index] + Dt_Diffuse1;
										}
#else
									   NewMolecule=new Molecule(ttime[index],L_grid,Type1);
									   if( Molecule_Array.empty())
									   {
										   ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;	
									   }
#endif
									   Molecule_Array.push_back(NewMolecule);
									   CountToNextFlowEntryR=TimeToNextFlowEntryR();
								   }
								   CountToNextFlowEntryR--;					/* If there is flow, count down the Flow time steps until the next flow entry */ 
								   break;
								}
							   ttime[NextFlow]+=Dt_Flow;
							   break;
				   case SetFlow:		/* 5 */
								if(MissedDetections==0)
								{
									if(RecycleCount==0)
										MoleculeData<<"Type\t"<<Mtype<<"\n";
									MoleculeData<<std::setprecision(12)<<Peak<<"\n";
								}
							    if (MissedDetections<=MaxMissedDetects)
								{									
									if(FlowDirection==1)
									{
										FlowDirection=-1;
									}
									else if (FlowDirection==-1)
									{
										FlowDirection=1;
									}																	
									tE=ttime[SetFlow]+ReversalDelay;
									ttime[NextFlow]=ttime[SetFlow]+0.5*(ttime[SetFlow]-(ttime[NextFlow]-Dt_Flow));
									RecycleCount++;
									MissedDetections++;
									ttime[SetFlow]+=ReversalDelay*4*double(MissedDetections);//extenderd auto recycling delay
									//ttime[SetFlow]=ttime[SetFlow]+2.0*ReversalDelay;
								}
								else
								{
									FlowDirection=1;
									RecycleCount=0;
#ifdef Method5
									MissedDetections=1;
#endif
									ttime[SetFlow]=INFTY;
								}								
							   break;
					case NextBackground:		// 4 
						       NextPhoton(ttime[index],PriorPhotonTiming);
						       PriorPhotonTiming=ttime[index];
							   ttime[NextBackground]+=TimeToNextBackground();
							   break;
	}
	index=idamin(&n, ttime, &incx);
	index--;
  }
}
