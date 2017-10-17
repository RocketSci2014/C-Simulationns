// single molecule recycling-1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "global.h"
#include "randomcup.h"
#include "stopwatch.h"
#include "smr.h"
#include"Molecule.h"
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	string comment;
	stopwatch a;
	smr b;
	b.MonteCarlo();
	comment="DNA 15mers, simulation time 1000s, concentration 0 pM, no photobleach, new parameters";
	b.WriteComment(comment);
	std::cout<<a.get_time()<<std::endl;
	return 0;
}

