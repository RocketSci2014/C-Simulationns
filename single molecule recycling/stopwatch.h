#include "stdafx.h"
#pragma once
#ifndef STOPWATCH_H
#define STOPWATCH_H
class stopwatch
{
private:
    clock_t start;
public:
	stopwatch(void): start(clock()){} //start counting time;
	double get_time();
	std::string currentDateTime();
	~stopwatch(void);
};

#endif
