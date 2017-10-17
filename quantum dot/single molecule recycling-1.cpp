// single molecule recycling-1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "global.h"
#include "randomcup.h"
#include "stopwatch.h"
#include "smr.h"
#include"QD.h"
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	string comment;
	stopwatch a;
	smr b;
	b.MonteCarlo();
	comment="2PE 10000s power law";
	b.WriteComment(comment);
	std::cout<<a.get_time()<<std::endl;
	return 0;
}

