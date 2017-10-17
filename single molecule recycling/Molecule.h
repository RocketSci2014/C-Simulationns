#pragma once
#ifndef MOLECULE_H
#define MOLECULE_H
#include "global.h"
#include "randomcup.h"
using namespace std;
class Molecule:
	public randomcup
{
public:
	double timing;//timing of the next photophysics event
	MoleculeState state;//state of the molecule
	MoleculeType type;//type of the molecule
	int position;//position of the molecule in a nanochannel

	Molecule(double ttime2,int position1,MoleculeType type1);//generate molecule and find the next excitation time of the molecule
	void FindNextExcitation(double ttime1);//find the timing of the next excitation after one photophysics event
	double PhotoPhysics(void);//To generate the photophysics, ttime is the end time of the current loop 
	~Molecule(void);
};

#endif