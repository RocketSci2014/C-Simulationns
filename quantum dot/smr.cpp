#include "stdafx.h"
#include "global.h"
#include "smr.h"
using namespace std;

/* used by idamin() */ 
const int n=NUM_TIMES;
const int incx=1;

ofstream PhotonData;	//timings of photons
ofstream Comment;

smr::smr(void)
{   
	smr::OpenFiles();
	InitializeWss();
	FlowDirection=1;			//flow direction is positive
   for(int i=0; i<NUM_TIMES; i++ ) {ttime[i]=INFTY;}
   ttime[NextFlow]=Dt_Flow;
   ttime[NextDiffuse]=Dt_Diffuse;
   CountDownTillAfterpulse=PhotonsUntilAnAfterpulse();
   ttime[NextBackground]=TimeToNextBackground(); 
   ttime[NextDiffuseEntry]=TimeToNextDiffuseEntry();
   CountToNextFlowEntry=TimeToNextFlowEntry();
   NextWssTime=BinTime;			//the starting point of wss is 1xBinTime
   index=idamin(&n, ttime, &incx);
   index--;
	PM1=-1.0*ReversalDelay;		//lower limit for the delay time between two consecutive passages
	PM2=ReversalDelay;			//upper limit for the invader molecule control algorithm 
	sum_of_elements=0;
	peak=0.0;
	hpeak=0.0;        
	hwss=0.0;         
	tE=0.0;
	peakindex=0;
	pass=false;
	first=true;
	keep=false;
	NumPassive=0;
	mindex=-1;					//initialize the index
}


smr::~smr(void)
{
	WriteComment;				//write comment
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
    std::string root="C:/Users/bwang/Documents/Bo Wang/simulation codes/QD 1D/data/";
	directory=root+directory;
	mkdir(directory.c_str());
	PhotonData.open(directory+"/PhotonData.bin",std::ios::binary);
   Comment.open(directory+"/Comment.txt");
}

void smr::CloseFiles(void)
{
	PhotonData.close();
	Comment.close();

}


void smr::InitializeWss(void)
{
	this->Wss_Array.reserve(wss_sig6);
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
		Weights[0][ii]=unsigned short int(128*exp(temp*(-0.5/temp1))+0.5); //the peak of the weights is 128
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
	}
	unsigned short int WssValue=Wss_Array[0]; //calculate the wss value by taking the last element in the array
	Wss_Array.erase(Wss_Array.begin());
	Wss_Array.push_back(0);
    return WssValue;
}

unsigned short int smr::PhotonCount(void)
	/*assign pulse to the beginning of pulse when starting the Monte-Carlo program*/
{
	unsigned short int PhotonCounts=0;	
	if (!Pulse.empty())
		tt++;		//make the array counter at the right place
		while (tt<Pulse.size()&&Pulse[tt]<NextWssTime)
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
	while(ii<Pulse.size()&&Pulse[ii]<hpeak-sig3+range)//photons should be within the range of hpeak-range tp hpeak+range
	{
		if(Pulse[ii]>hpeak-sig3-range)
		{
			sum_of_elements+=Pulse[ii];
			jj++;
		}
		ii++;
	}
	peak=sum_of_elements/jj;						//calculate the ML estimated peak from the range
	sum_of_elements=0;
}


bool smr::MakeWss(void)
{
	int wss=WeightedSlidingSum(PhotonCount());
	double delay;
	bool PassIndicator=false;
	if(!Pulse.empty())
	{
		if (wss<=Wss_Th-Fluctuation)
		{
			if (pass==true)
			{
				FindPeakMethod(3.0*sig);				//try the 3.0 sigma range
				pass=false;
				PassIndicator=true;
				delay=peak-tE;
				if (delay>PM1&&delay<PM2||peakindex==0)//conditions for setting ttime[SetFlow]
				{
					ttime[SetFlow]=peak+ReversalDelay; //also set values of peak hpeak and hwss to 0 in the monte-carlo function
					tE=peak+2.0*ReversalDelay;
					NumPassive=0;					
				}
				else
				{
					ttime[SetFlow]=INFTY;				//terminate the recyle when finding a invader molecule
					peakindex=0;
				}
				hwss=0;
			}
		}
		else if (wss>Wss_Th)
		{
			if (pass==false)
			{
				while (!Pulse.empty()&&Pulse[0]<NextWssTime-sig6)
				{
				   Pulse.erase(Pulse.begin());
				   tt--;
				}
				  pass=true;
			}

			if (wss>hwss)
			{
				hwss=wss;
				hpeak=NextWssTime;
			}
		}
	}
	NextWssTime+=BinTime;
	return PassIndicator; 
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
		Pulse.push_back(PhotonTiming);
		if (Pulse.size()==1)
			tt=0;
		if( --CountDownTillAfterpulse<=0 )  /*There is a possibility that this photon gives an afterpulse */ 
			{
				CountDownTillAfterpulse=PhotonsUntilAnAfterpulse();
				Pulse.push_back(PhotonTiming+TimeToAfterpulse());
				WritePhoton(Pulse.back());
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
	QD* NewQD; //temperary storage for the new molecule
	while (ttime[index]<SimulationTime)
	{
			/*do photophysics for each QD*/
		if (!QD_Array.empty())
		{
		   for (int ii=0;ii<QD_Array.size();ii++)
		   {
			   CurrentTime=QD_Array[ii]->timing;
			   while (CurrentTime<ttime[index])	
			   {
			       PhotonTiming=QD_Array[ii]->PhotoPhysics();//update current timing of photophysics
				   if(PhotonTiming<INFTY)					//for each molecule, do photophysics until exceed the system time
				   {
						NextPhoton(PhotonTiming,PriorPhotonTiming);
						PriorPhotonTiming=PhotonTiming;
				   }
					CurrentTime=QD_Array[ii]->timing;		//update current timing of photophysics
			   }
		   }
		}
		 /*do Monte-Carlo simulation, switch between Next Diffuse entry,next diffuse, next flow, set flow, and background*/
		   switch(index)
		   {
			      case NextDiffuseEntry: //0
					  NewQD=new QD(ttime[index],PositionOfNextDiffuseEntry( ));
					  QD_Array.push_back(NewQD);//push back a new molecule into the array	
					  if( QD_Array.size()==1 ) 
					   {
						   ttime[NextDiffuse]=ttime[NextDiffuseEntry] + Dt_Diffuse;	//If there were  previously no molecules and one has just entered, we need to set the times for NextDiffuse 
					   }
					   ttime[NextDiffuseEntry]	+=TimeToNextDiffuseEntry();
					   break;
				   case NextDiffuse:	//1
					   if (!QD_Array.empty())
					   {
						   int jj=0;
						   while (jj<QD_Array.size())
						   {
							   itemp=Diffuse();
							   if( itemp !=0 )
							   {
								   QD_Array[jj]->position+=itemp;
								   if( QD_Array[jj]->state==S0 )				   
								      QD_Array[jj]->FindNextExcitation(ttime[index]);
								   if( QD_Array[jj]->position>L_grid || QD_Array[jj]->position<-L_grid )
								   {
									   delete QD_Array[jj];
									   QD_Array.erase(QD_Array.begin()+jj);//remove the molecule, and the length of the array is reduced by one
								   }
								   else
									   jj++;
							   }
							   else
								   jj++;
						   }
						   if(QD_Array.empty())
							   ttime[NextDiffuse]=INFTY;
					   }
					   ttime[NextDiffuse]+=Dt_Diffuse;
					   break;
					case NextFlow:			// 2
							   switch( FlowDirection )
							   { /* Note: case 0: should never occur as ttime[NextFlow] would equal INFTY*/
							   case 1:
								   if(!QD_Array.empty())
								   {
  									   int kk=0;
									   while (kk<QD_Array.size())
									   {
											QD_Array[kk]->position++;
											if( QD_Array[kk]->state==S0 )				   
												QD_Array[kk]->FindNextExcitation(ttime[index]);
											if( QD_Array[kk]->position>L_grid )
											{

											   delete QD_Array[kk];
												QD_Array.erase(QD_Array.begin()+kk);//remove the molecule, and the length of the array is reduced by one
											}
											else
												kk++;
									   }	

										  if (QD_Array.empty())
											ttime[NextDiffuse]=INFTY;
								   }
								   if( CountToNextFlowEntry==0 )
								   {									   									   
									   NewQD=new QD(ttime[index],-L_grid);
									   if( QD_Array.empty())
									   {
										   ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;	
									   }
									   QD_Array.push_back(NewQD);
									   CountToNextFlowEntry=TimeToNextFlowEntry();
								   }
								   CountToNextFlowEntry--;					//If there is flow, count down the Flow time steps until the next flow entry 
								   break;
							   case -1:
   								   if(!QD_Array.empty())
								   {
  									   int dd=0;
									   while (dd<QD_Array.size())
									   {
											QD_Array[dd]->position--;

											if( QD_Array[dd]->state==S0 )				   
												QD_Array[dd]->FindNextExcitation(ttime[index]);
											if( QD_Array[dd]->position<-L_grid )
											{
											   delete QD_Array[dd];
											   QD_Array.erase(QD_Array.begin()+dd);//remove the molecule, and the length of the array is reduced by one
											}
											else
												dd++;
									   }	

										  if (QD_Array.empty())
											ttime[NextDiffuse]=INFTY;
								   }
								   if( CountToNextFlowEntry==0 )
								   {

									   NewQD=new QD(ttime[index],L_grid);
									   if( QD_Array.empty())
									   {
										   ttime[NextDiffuse]=ttime[index] + Dt_Diffuse;	
									   }
									   QD_Array.push_back(NewQD);
									   CountToNextFlowEntry=TimeToNextFlowEntry();
								   }
								   CountToNextFlowEntry--;					//If there is flow, count down the Flow time steps until the next flow entry 
								   break;
								}
							   ttime[NextFlow]+=Dt_Flow;
							   break;
				   case SetFlow:		//3
							   if(FlowDirection==1)
							   {
								   FlowDirection=-1;
							   }
							   else if (FlowDirection==-1)
							   {
								   FlowDirection=1;
							   }
							   peakindex++;
							   tE=ttime[SetFlow]+ReversalDelay;
							   if (NumPassive<Auto_Recycle)
							   {
								 ttime[SetFlow]=ttime[SetFlow]+2.0*ReversalDelay;
								 NumPassive++;
							   }
							   else
							   {
								  ttime[SetFlow]=INFTY;
								  peakindex=0;
							   }
							   break;
					case NextBackground:		//5 
						       NextPhoton(ttime[index],PriorPhotonTiming);
						       PriorPhotonTiming=ttime[index];
							   ttime[NextBackground]+=TimeToNextBackground();
							   break;
	}
	index=idamin(&n, ttime, &incx);
	index--;
  }
}
