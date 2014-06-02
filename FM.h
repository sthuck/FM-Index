/*  FM-Index - Text Index
 *  Copyright (C) 2011  Matthias Petri
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. */
 
#ifndef FM_H
#define	FM_H
#define size_uchar 256

#include <stdarg.h>

#include "util.h"
#include <list>

#include "wt.hpp"
#include "bit_vectors.hpp"

using namespace std;
using namespace sdsl;

#define DEFAULT_SAMPLERATE      64
#define RRR_SAMPLERATE			20

class FM {
	typedef  wt<unsigned char*,rrr_vector<>,rrr_vector<>::rank_1_type,rrr_vector<>::select_1_type,rrr_vector<>::select_0_type> myWt;
//	typedef wt_blcd<> myWt;   //uncompressed Wavelet tree

public:
    FM(uint8_t* T,uint32_t n,uint8_t debug,uint32_t samplerate);
    void build(uint8_t* T,uint32_t n,uint32_t samplerate);
    static FM* load(char* filename);
    int32_t save(char* filename);
    uint8_t* remap0(uint8_t* T,uint32_t n);
    uint32_t count(uint8_t* pattern,uint32_t m);
    uint32_t* locate(uint8_t* pattern,uint32_t m,uint32_t* matches);
    uint32_t getSize();
    uint8_t* extract(uint32_t start,uint32_t stop);
    uint8_t* reconstructText(uint32_t* n);
    float getSizeN();
    virtual ~FM();
    uint32_t backSearchHelper(uint8_t myChar, myWt* bwt, uint32_t& sp, uint32_t& ep);
    uint32_t backSearchHelperNoUpdate(uint8_t myChar, myWt* bwt, uint32_t sp, uint32_t ep);
    int searchHelper(uint8_t myChar, uint32_t& sp, uint32_t& ep, myWt* bwt_r, uint32_t& sp_r, uint32_t& ep_r) ;
    uint32_t* locateAfterSearch(uint32_t sp,uint32_t ep,uint32_t* matches);
    list<uint32_t*>* Search1Error(uint8_t* pattern,uint32_t m);
    list<uint32_t*>* Search2Error(uint8_t* pattern,uint32_t m);
    int doNTimesForwardSearch(uint8_t* pattern, uint32_t pos, uint32_t endPostion, uint32_t& sp, uint32_t& ep, uint32_t& sp_r, uint32_t& ep_r );
    int doNTimesBackSearch(uint8_t* pattern, uint32_t pos, uint32_t endPostion, uint32_t& sp, uint32_t& ep, uint32_t& sp_r, uint32_t& ep_r );
    bool SanityCheck(uint8_t* pattern,uint32_t m,uint8_t error_number);


    void deb();
    void debr();

public:
	static void info(const char *format,...)
	{
		if(FM::verbose == 1) {
			va_list vargs;
			va_start (vargs, format);
			vfprintf (stderr, format, vargs);
			fprintf (stderr, "\n");
		}
	}
	static int verbose;

private:
    FM();
private:

  uint32_t sigma;
  uint32_t samplerate;
  int32_t I;
  uint32_t n;
  uint32_t C[size_uchar+1];
  uint8_t remap[size_uchar];
  uint8_t* remap_reverse;
  uint32_t* suffixes;    //This is ssa
  uint32_t* positions;
  rrr_vector<>* sampled;  //This is Mark[0]
  myWt *T_bwt;
  myWt *T_bwt_reverse;
  int32_t* SA;
  int32_t* SA_reverse;
  uint8_t* X;
  uint8_t* X_reverse;
  uint8_t debug;
};


#endif	/* FM_H */

