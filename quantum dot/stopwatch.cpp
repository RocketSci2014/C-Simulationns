#include "stdafx.h"
#include "global.h"
#include "stopwatch.h"


double stopwatch::get_time()
{
	return double(clock()-start)/CLOCKS_PER_SEC;
}
std::string stopwatch::currentDateTime()
{
	time_t     now = time(0);
    struct tm  *tstruct=localtime(&now);
	std::ostringstream time_stream;
    //strftime(buf, sizeof(buf), "%Y-%m-%d-%X", &tstruct);
    time_stream<<(tstruct->tm_year+1900)<<"-"<<(tstruct->tm_mon+1)<<"-"<<tstruct->tm_mday<<"@"<<tstruct->tm_hour<<"-"<<tstruct->tm_min<<"-"<<tstruct->tm_sec;
    return time_stream.str();
}
stopwatch::~stopwatch(void)
{
}
