#include<string.h>
using namespace std;

/**
 * Contains information about current processed segment (chromosome, start and end of the region)
 */
 
 class CurrentSegment {
    public:
		string chr;
    	int start;
    	int end;

    	CurrentSegment(string chr, int start, int end) {
        	this->chr = chr;
        	this->start = start;
        	this->end = end;
    };
};
