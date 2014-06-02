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

#include "FM.h"
#include "util.h"
#include "divsufsort.h"
#include "wavelet_trees.hpp"
#include "int_vector.hpp"

#include <stack>
#include <algorithm>
#define pushRange() backup->push(sp); backup->push(ep); backup->push(sp_r); backup->push(ep_r);
#define popRange() ep_r=backup->top();backup->pop();sp_r=backup->top();backup->pop();ep=backup->top();backup->pop();sp=backup->top();backup->pop();

void* malloc(size_t num);

int FM::verbose = 0;

FM::FM(uint8_t* T,uint32_t N,uint8_t _debug=0,uint32_t _samplerate = DEFAULT_SAMPLERATE)
	: samplerate(_samplerate),n(N),debug(_debug) {

    /* 0 terminate */
    if(T[N-1] != 0) {
        T = (uint8_t*) safe_realloc(T,(N+1) * sizeof(uint8_t));
        T[N] = 0;
        this->n++;
    }

    build(T,n,samplerate);
}

FM::FM()  :
		sigma(0),
        samplerate(0),
		I(0),
        n(0),
		remap_reverse(NULL),
		suffixes(NULL),
		positions(NULL),
		sampled(new rrr_vector<>),
		T_bwt(new myWt()),
		T_bwt_reverse(new myWt()),
		SA(NULL),
		SA_reverse(NULL),
		X(NULL),
		X_reverse(NULL),
        debug(0)
		{}

FM::~FM() {
	free(remap_reverse);
	free(suffixes);
	free(positions);
	delete T_bwt;
	delete T_bwt_reverse;
	delete sampled;
	if (this->debug) {
		free(this->X);
		free(this->X_reverse);
		free(this->SA);
		free(this->SA_reverse);
	}
}

uint32_t
FM::getSize() {
    uint32_t bytes = 0;

    bytes += sizeof(this->n);
    bytes += sizeof(this->samplerate);
    bytes += sizeof(this->sigma);
    bytes += sizeof(this->I);
    bytes += sizeof(this->remap);
    bytes += sizeof(this->C);
    bytes += this->sigma * sizeof(uint8_t); /* remap_reverse */
    bytes += ((n/samplerate)+1) * sizeof(uint32_t); /* suffixes */
	bytes += ((n/samplerate)+2) * sizeof(uint32_t); /* positions */
    bytes += sdsl::util::get_size_in_bytes(*(this->sampled));
    bytes += sdsl::util::get_size_in_bytes(*(this->T_bwt));
    bytes += sdsl::util::get_size_in_bytes(*(this->T_bwt_reverse));

    return bytes;
}

float
FM::getSizeN() {
    uint32_t bytes = getSize();
    return (float)(bytes)/(float)(n);
}

uint8_t*
FM::remap0(uint8_t* T,uint32_t n) {
    uint8_t* X;
    uint32_t i,j,size=0;
    uint32_t freqs[size_uchar];

    for(i=0;i<size_uchar;i++) freqs[i]=0;
    for(i=0;i<n;i++) if(freqs[T[i]]++==0) size++;

    this->sigma=size;

    // remap alphabet
    if (freqs[0]>1) {i=1;sigma++;} //test if some character of T is zero, we already know that text[n-1]='\0'
    else i=0;

    remap_reverse = (uint8_t*) malloc(size*sizeof(uint8_t));
    for(j=0;j<size_uchar;j++) {
      if(freqs[j]!=0) {
        remap[j]=i;
        remap_reverse[i++]=j;
      }
    }
    // remap text
    X = (uint8_t*)malloc(n * sizeof(uint8_t));
    for(i=0;i<n-1;i++) // the last character must be zero
      X[i]=remap[T[i]];
    X[n-1]=0;		   // Aviad: I think it's a bug, malloc doesn't initialize
    return X;
}

void
FM::build(uint8_t* T,uint32_t n,uint32_t samplerate) {
    uint8_t* X;
    uint8_t* X_bwt;
    int32_t* SA;
    uint32_t i,prev,tmp,start,stop;
    float elapsed;
    start = gettime();
	info("building index.");

    /* remap if 0 in text */
    info("- remapping alphabet.");
    X = remap0(T,n);
    free(T);

    /* create cumulative counts */
    info("- creating cumulative counts C[].");
    for (i=0;i<size_uchar+1;i++) C[i]=0;
    for (i=0;i<n;++i) C[X[i]]++;
    prev=C[0];C[0]=0;
    for (i=1;i<size_uchar+1;i++) {
      tmp = C[i];
      C[i]=C[i-1]+prev;
      prev = tmp;
    }

    /* perform k-BWT */
    info("- performing bwt.");
    SA = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
    if( divsufsort(X,SA,n) != 0 ) {
        fatal("error divsufsort");
    }

    /* sample SA for locate() */
    info("- sample SA locations.");
    suffixes = (uint32_t*) safe_malloc( ((n/samplerate)+1) * sizeof(uint32_t));

    bit_vector *_sampled = new bit_vector();
    _sampled->resize(n);

    tmp = 0;
    for(i=0;i<n;i++) {
        if( SA[i] % samplerate == 0) {
            suffixes[tmp] = SA[i];
            (*_sampled)[i]=1;
            tmp++;
        } else (*_sampled)[i]=0;
    }

    this->sampled = new rrr_vector<>(*_sampled,RRR_SAMPLERATE);
    delete(_sampled);

	/* sample SA for display() */
	positions = (uint32_t*) safe_malloc( ((n/samplerate)+2) * sizeof(uint32_t));
    for (i=0;i<this->n;i++) {
        if (SA[i] % samplerate == 0) positions[SA[i]/samplerate] = i;
	}
    positions[(this->n-1)/samplerate+1] = positions[0];

    info("- creating bwt output.");
    X_bwt = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
    for(i=0;i<n;i++) {
        if(SA[i]==0) {
            X_bwt[i] = X[n-1];
            this->I = i;
        } else X_bwt[i] = X[SA[i]-1];
    }

    info("- create RRR wavelet tree over bwt.");
    T_bwt = new myWt(X_bwt,n);
    //free(X_bwt);
    if (debug)
    	this->SA=SA;
    else
    	free(SA);

    uint8_t* X_reverse = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
	for (i=0;i<=n-2;i++){
		X_reverse[n-2-i]=X[i];
	}
	X_reverse[n-1]=X[n-1];

	if (debug)
	    	this->X=X;
	    else
	    	free(X);

	int32_t* SA_reverse = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
	if( divsufsort(X_reverse,SA_reverse,n) != 0 ) {
		 fatal("error divsufsort");
	}

   uint8_t* X_bwt_reverse = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
   for(i=0;i<n;i++) {
	   if(SA_reverse[i]==0) {
		   X_bwt_reverse[i] = X_reverse[n-1];
	   } else X_bwt_reverse[i] = X_reverse[SA_reverse[i]-1];
   }

    if (!debug) {
    	free(SA_reverse);
    	free(X_reverse);
    }
    else {
    	this->SA_reverse=SA_reverse;
    	this->X_reverse=X_reverse;
    }

    info("- create RRR wavelet tree over bwt_r.");
    T_bwt_reverse = new myWt(X_bwt_reverse,n);
   //free(X_bwt_reverse);

    stop = gettime();
    elapsed = (float)(stop-start)/1000000;

    /* build aux data */
    info("build FM-Index done. (%.3f sec)",elapsed);

    uint32_t bytes;
	info("space usage:");
	bytes = sigma * sizeof(uint8_t);
	info("- remap_reverse: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = sizeof(this->C);
	info("- C: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = ((n / samplerate) + 1) * sizeof(uint32_t);
	info("- Suffixes: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = ((n / samplerate) + 2) * sizeof(uint32_t);
	info("- Positions: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = sdsl::util::get_size_in_bytes(*this->sampled);
	info("- Sampled: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = sdsl::util::get_size_in_bytes(*T_bwt);
	info("- T_bwt: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	bytes = sdsl::util::get_size_in_bytes(*T_bwt_reverse);
	info("- T_bwt_reverse: %d bytes (%.2f\%)", bytes, (float) bytes / getSize() * 100);
	info("input Size n = %lu bytes\n", this->n);
	info("index Size = %lu bytes (%.2f n)", getSize(), getSizeN());
}


int32_t
FM::save(char* filename) {
    std::ofstream f;
    f.open(filename,std::ios::out | std::ios::binary);

	info("writing FM Index to file '%s'",filename);
    if(f.is_open()) {
        f.write(reinterpret_cast<char*>(&samplerate),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&sigma),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&I),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&n),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(C),sizeof(uint32_t)*(size_uchar+1));
        f.write(reinterpret_cast<char*>(remap),sizeof(uint8_t)*size_uchar);
        f.write(reinterpret_cast<char*>(remap_reverse),sizeof(uint8_t)*sigma);
        f.write(reinterpret_cast<char*>(suffixes),sizeof(uint32_t)*((n/samplerate)+1));
		f.write(reinterpret_cast<char*>(positions),sizeof(uint32_t)*((n/samplerate)+2));
        T_bwt->serialize(f);
        sampled->serialize(f);
        T_bwt_reverse->serialize(f);
        f.write(reinterpret_cast<char*>(&debug),sizeof(uint8_t));
        if (debug) {
        	f.write(reinterpret_cast<char*>(SA),sizeof(uint32_t)*n);
        	f.write(reinterpret_cast<char*>(SA_reverse),sizeof(uint32_t)*n);
        	f.write(reinterpret_cast<char*>(X),sizeof(uint8_t)*n);
        	f.write(reinterpret_cast<char*>(X_reverse),sizeof(uint8_t)*n);
        }
        f.close();
    } else return 1;

    return 0;
}

FM*
FM::load(char* filename) {
	FM* newIdx = new FM();
	std::ifstream f;
	f.open(filename, std::ios::in | std::ios::binary);

	if (f.is_open()) {
		info("loading FM Index from file '%s'", filename);
		f.read(reinterpret_cast<char*>(&newIdx->samplerate), sizeof(uint32_t));
		f.read(reinterpret_cast<char*>(&newIdx->sigma), sizeof(uint32_t));
		f.read(reinterpret_cast<char*>(&newIdx->I), sizeof(uint32_t));
		f.read(reinterpret_cast<char*>(&newIdx->n), sizeof(uint32_t));
		f.read(reinterpret_cast<char*>(newIdx->C), sizeof(uint32_t) * (size_uchar + 1));
		f.read(reinterpret_cast<char*>(newIdx->remap), sizeof(uint8_t) * size_uchar);
		newIdx->remap_reverse = (uint8_t*) safe_malloc(sizeof(uint8_t) * (newIdx->sigma));
		f.read(reinterpret_cast<char*>(newIdx->remap_reverse), sizeof(uint8_t) * newIdx->sigma);
		newIdx->suffixes = (uint32_t*) safe_malloc(sizeof(uint32_t) * ((newIdx->n / newIdx->samplerate) + 1));
		f.read(reinterpret_cast<char*>(newIdx->suffixes), sizeof(uint32_t) * ((newIdx->n / newIdx->samplerate) + 1));
		newIdx->positions = (uint32_t*) safe_malloc(sizeof(uint32_t) * ((newIdx->n / newIdx->samplerate) + 2));
		f.read(reinterpret_cast<char*>(newIdx->positions), sizeof(uint32_t) * ((newIdx->n / newIdx->samplerate) + 2));
		newIdx->T_bwt->load(f);
		newIdx->sampled->load(f);
		newIdx->T_bwt_reverse->load(f);

		f.read(reinterpret_cast<char*>(&newIdx->debug), sizeof(uint8_t));
		if (newIdx->debug) {
			newIdx->SA = (int32_t*) safe_malloc(sizeof(int32_t) * newIdx->n);
			f.read(reinterpret_cast<char*>(newIdx->SA), sizeof(int32_t) * newIdx->n);
			newIdx->SA_reverse = (int32_t*) safe_malloc(sizeof(int32_t) * newIdx->n);
			f.read(reinterpret_cast<char*>(newIdx->SA_reverse), sizeof(int32_t) * newIdx->n);
			newIdx->X = (uint8_t*) safe_malloc(sizeof(uint8_t) * newIdx->n);
			f.read(reinterpret_cast<char*>(newIdx->X), sizeof(uint8_t) * newIdx->n);
			newIdx->X_reverse = (uint8_t*) safe_malloc(sizeof(uint8_t) * newIdx->n);
			f.read(reinterpret_cast<char*>(newIdx->X_reverse), sizeof(uint8_t) * newIdx->n);
		}

		f.close();
		info("samplerate = %d", newIdx->samplerate);
		info("sigma = %d", newIdx->sigma);
		info("I = %d", newIdx->I);
		info("n = %d", newIdx->n);
	} else {
		delete newIdx;
		return NULL;
	}

	return newIdx;
}

uint32_t*
FM::locateAfterSearch(uint32_t sp,uint32_t ep,uint32_t* matches) {
    uint32_t* locations;
    uint8_t c;
    uint32_t i;
    sdsl::rrr_rank_support<>* rs = new rrr_rank_support<>();
    sdsl::util::init_support<rrr_rank_support<>,rrr_vector<> >(*rs,sampled);
    if (sp<=ep) {
        /* determine positions */
        *matches = ep-sp+1;
        uint32_t locate=0;
        locations= (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
        i=sp;
        int32_t j,dist,rank;
        while (i<=ep) {
            j=i,dist=0;
            while (!((*sampled)[j])) {
                c = (*T_bwt)[j];
                rank = T_bwt->rank(j+1,c)-1;
                j = C[c]+rank; // LF-mapping
                ++dist;
            }

            locations[locate]=suffixes[rs->rank(j+1)-1]+dist;
            locate++;
            ++i;
        }
        /* locations are in SA order */
        std::sort(locations,locations+(*matches));
        delete (rs);
        return locations;
    } else {
      /* no matches */
      *matches = 0;
      delete(rs);
      return NULL;
    }
    delete(rs);
    return locations;
}

int FM::doNTimesBackSearch(uint8_t* pattern, uint32_t pos, uint32_t endPostion,
		uint32_t& sp, uint32_t& ep, uint32_t& sp_r, uint32_t& ep_r) {

	if (sp > ep || sp_r > ep_r)
		return 0;
	int res = 1;
	for (uint32_t l = pos; (l >= endPostion && res); l--) {
		res = searchHelper(remap[pattern[l]], sp_r, ep_r, this->T_bwt, sp, ep);
		if (l==0)
			break;  //l is uint
	}
	return res;
}

int FM::doNTimesForwardSearch(uint8_t* pattern, uint32_t pos,
		uint32_t endPostion, uint32_t& sp, uint32_t& ep, uint32_t& sp_r,
		uint32_t& ep_r) {

	if (sp > ep || sp_r > ep_r)
		return 0;
	int res = 1;
	for (uint32_t l = pos; (l <= endPostion && res); l++) {
		res = searchHelper(remap[pattern[l]], sp, ep, this->T_bwt_reverse, sp_r,
				ep_r);
	}
	return res;
}

/**
 * This function implements a Bi-directional BWT search.
 * forward or backward: it depends on what params you pass to it.
 *
 * searchHelper(c,sp,ep,T_bwt_reverse,sp_r,ep_r)) ==== THIS IS FORWARD SEARCH
 * searchHelper(c,sp_r,ep_r,T_bwt,sp,ep)) ==== THIS IS BACKWARD SEARCH
 *
 *
 * Receives: myChar - Should be already mapped!
 * bwt, start_postion, and end_position
 * and bwt_reverse, start_postion_reverse,end_postion_reverse
 *
 * Returns: 0 if failed, 1 if success
 * Changes: sp,ep,sp_r,ep_r
 */
int FM::searchHelper(uint8_t myChar, uint32_t& sp, uint32_t& ep, myWt* bwt_r, uint32_t& sp_r, uint32_t& ep_r) {

	uint32_t old_ep_r = ep_r, old_sp_r = sp_r, x;
	if (!this->backSearchHelper(myChar, bwt_r, sp_r, ep_r)) return 0;

	/**
	 * now we update sp and ep, according to sp_r and ep_r.
	 * THIS IS THE CORE OF 2-WAY SEARCH
	 */
	x = 0;
	for (int j = 0; j < myChar; j++) {
		x += backSearchHelperNoUpdate(j, bwt_r, old_sp_r, old_ep_r);
	}
	sp = sp + x;
	ep = sp + (ep_r - sp_r);
	return 1;
}


/**
 * This function does a backward search on a specific bwt array, and updates parameters: sp, ep.
 * The bwt array can be the "normal" bwt or the reversed bwt.
 */
uint32_t
FM::backSearchHelper(uint8_t myChar, myWt* bwt, uint32_t& sp, uint32_t& ep) {

	 sp = C[myChar] + bwt->rank(sp-1+1,myChar); /* LF Mapping */
	 ep = C[myChar] + bwt->rank(ep+1,myChar)-1; /* LF Mapping */



	 if (sp<=ep) {
	       return ep-sp+1;
	     } else {
	       return 0;
	     }
}
/**
 * This function does a backward search, and DOES NOT updates parameters: sp, ep.
 */
uint32_t
FM::backSearchHelperNoUpdate(uint8_t myChar, myWt* bwt, uint32_t sp, uint32_t ep) {

	sp = C[myChar] + bwt->rank(sp-1+1,myChar); /* LF Mapping */
	ep = C[myChar] + bwt->rank(ep+1,myChar)-1; /* LF Mapping */


	 if (sp<=ep) {
	       return ep-sp+1;
	     } else {
	       return 0;
	     }
}


/**
 * Length of Pattern is at least 2! sigma>2; (more than two chars at least + '0')
 */
list<uint32_t*>*
FM::Search1Error(uint8_t* pattern,uint32_t m) {
	list<uint32_t*>* resultsList = new list<uint32_t*>();
	uint32_t* singleResult;
	stack<uint32_t>* backup = new stack<uint32_t>();
	uint32_t e,x;
	x = (m%2==0) ? m/2-1 : m/2;
	m=m-1;

	//First case, error is in [0..x].

	uint32_t sp = C[remap[pattern[m]]];
	uint32_t ep = C[remap[pattern[m]]+1]-1;
	uint32_t sp_r=sp;
	uint32_t ep_r=ep;
	//ugly code :(
	if (doNTimesBackSearch(pattern,m-1,x+1,sp,ep,sp_r,ep_r)) 				//got [x+1..m]
		for (uint32_t i=x;;i--) {											//i=location of e
			pushRange();
			if (doNTimesBackSearch(pattern,x,i+1,sp,ep,sp_r,ep_r))			//got [i+1..x..m])
				for (e=1;e<sigma;e++) {										//running through all the letters
					pushRange();
					if (e!=remap[pattern[i]]) 								//searching one char, e: for all e in Sigma, where pattern[i]!=e. So in total - e[i+1..x..m]
						if (searchHelper(e,sp_r,ep_r,T_bwt,sp,ep))			// e[i+1..x..m]
							if ((i==0) || doNTimesBackSearch(pattern,i-1,0,sp,ep,sp_r,ep_r)) {     //[0..i-1]e[i+1...m])
								singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
								singleResult[0]=sp;
								singleResult[1]=ep;
								resultsList->push_back(singleResult);
							 }
					popRange();
				}
			popRange()
			if (i==0) break;
		}

	//Now, for second option: Error is in  [x..m]. resetting everything.
	backup->empty();
	sp = C[remap[pattern[0]]];
	ep = C[remap[pattern[0]]+1]-1;
	sp_r=sp;
	ep_r=ep;

	if (doNTimesForwardSearch(pattern,1,x,sp,ep,sp_r,ep_r))		  				//got [0..x]
		for (uint32_t i=x+1;i<=m;i++) {
			pushRange();
			if (doNTimesForwardSearch(pattern,x+1,i-1,sp,ep,sp_r,ep_r))			//[0..x..i-1]
				for (e=1;e<sigma;e++) {											//running through all the letters
					pushRange();
					if (e!=remap[pattern[i]]) {
						if(searchHelper(e,sp,ep,T_bwt_reverse,sp_r,ep_r))   	//[0..i-1]e
							if(doNTimesForwardSearch(pattern,i+1,m,sp,ep,sp_r,ep_r)) {	//[0..i-1]e[i+1...m]
								singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
								singleResult[0]=sp;
								singleResult[1]=ep;
								resultsList->push_back(singleResult);
							 }
					}
					popRange();
				}
			popRange()
		}
	delete(backup);
	return resultsList;
}


list<uint32_t*>*
FM::Search2Error(uint8_t* pattern,uint32_t m) {
	list<uint32_t*>* resultsList = new list<uint32_t*>();
	uint32_t* singleResult;
	uint32_t i;
	uint32_t e, e2, s1, s2;
	s1 = m/3; s2 = m - s1--;  s2--;  m=m-1;
	stack<uint32_t>* backup = new stack<uint32_t>();

	//First case, mismatches in the first two parts.
	uint32_t sp = C[remap[pattern[m]]];
	uint32_t ep = C[remap[pattern[m]]+1]-1;
	uint32_t sp_r = sp; uint32_t ep_r = ep;

	if (doNTimesBackSearch(pattern, m-1, s2+1, sp, ep, sp_r, ep_r))	//got [s2+1...m]
		for (uint32_t j = s2; j >= 1; j--) {						//j = location of e1
			pushRange();
			if (doNTimesBackSearch(pattern, s2, j+1, sp, ep, sp_r, ep_r)) //[j+1..s2..m]
				for (e = 1; e < sigma; e++) {
					pushRange();
					if (e != remap[pattern[j]])
						if (searchHelper(e, sp_r, ep_r, T_bwt, sp, ep)) 	//got e1[j+1...m]
							for (uint32_t l=j-1;;l--) {
								pushRange();
								if (doNTimesBackSearch(pattern,j-1,l+1,sp,ep,sp_r,ep_r))
										for (e2=1;e2<sigma;e2++) {
											pushRange();
											if (e2!=remap[pattern[l]])
												if (searchHelper(e2, sp_r, ep_r, T_bwt, sp, ep))
													if (l==0 || doNTimesBackSearch(pattern,l-1,0,sp,ep,sp_r,ep_r)) {
														singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
														singleResult[0]=sp;
														singleResult[1]=ep;
														resultsList->push_back(singleResult);
													}
											popRange();
										}
								popRange();
								if (l==0) break;
							}
					popRange();
				}
			popRange();
		}

	backup->empty();
	//Case B: Both mismatches occuer in last part
	//This is a lot like how I picture hell.
	sp=C[remap[pattern[0]]];
	ep=C[remap[pattern[0]]+1]-1;
	sp_r=sp; ep_r=ep;

	if (doNTimesForwardSearch(pattern,1,s2,sp,ep,sp_r,ep_r)) 			 								//[0..s2]
		for (i=s2+1;i<=m-1;i++) {											 							//i=position of e1
			pushRange();
			if (doNTimesForwardSearch(pattern,s2+1,i-1,sp,ep,sp_r,ep_r))		  						//[0..s2..i-1]
				for (e=1;e<sigma;e++) {											  						//running through all the letters
					pushRange();
					if (e!=remap[pattern[i]])
						if (searchHelper(e,sp,ep,T_bwt_reverse,sp_r,ep_r))   							//[0..i-1]e
							for (uint32_t j=i+1;j<=m;j++) {												//j=position of e2
								pushRange();
								if (doNTimesForwardSearch(pattern,i+1,j-1,sp,ep,sp_r,ep_r))  			//[0..i-1]e1[i+1..j-1]
									for (e2=1;e2<sigma;e2++) {
										pushRange();
										if (e2!=remap[pattern[j]]) {
											if (searchHelper(e2,sp,ep,T_bwt_reverse,sp_r,ep_r))			//[0..i-1]e1[i+1..j-1]e2
												if (doNTimesForwardSearch(pattern,j+1,m,sp,ep,sp_r,ep_r)) {  	//0..i-1]e1[i+1..j-1]e2[j+1...m]
													singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
													singleResult[0]=sp;
													singleResult[1]=ep;
													resultsList->push_back(singleResult);
												 }
										}
										popRange();
									}
								popRange();
							}
					popRange();
				}
			popRange();
		}
	//Case C: The mismatches in the second and third parts.
	backup->empty();
	sp=C[remap[pattern[0]]];
	ep=C[remap[pattern[0]]+1]-1;
	sp_r=sp;ep_r=ep;

	if (doNTimesForwardSearch(pattern,1,s1,sp,ep,sp_r,ep_r))  //[0..s1]
		for (i=s1+1;i<=s2;i++) {
			pushRange();
			if (doNTimesForwardSearch(pattern,s1+1,i-1,sp,ep,sp_r,ep_r))  //[0..s1..i-1]
				for (e=1;e<sigma;e++) {
					pushRange();
					if (e!=remap[pattern[i]])
						if (searchHelper(e,sp,ep,T_bwt_reverse,sp_r,ep_r))
							if(doNTimesForwardSearch(pattern,i+1,s2,sp,ep,sp_r,ep_r))
								for (uint32_t j=s2+1;j<=m;j++) {
									pushRange();
									if (doNTimesForwardSearch(pattern,s2+1,j-1,sp,ep,sp_r,ep_r))  // [0..s1..i-1]e[i+1..s2...j-1]
										for (e2=1;e2<sigma;e2++) {
											pushRange();
											if (e2!=remap[pattern[j]])
												if (searchHelper(e2,sp,ep,T_bwt_reverse,sp_r,ep_r))
													if (doNTimesForwardSearch(pattern,j+1,m,sp,ep,sp_r,ep_r)) {
														singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
														singleResult[0]=sp;
														singleResult[1]=ep;
														resultsList->push_back(singleResult);
													}
											popRange();
										}
									popRange();
								}
					popRange();
				}
			popRange();
		}

	//Case D: mimatches occur in the first or last parts
	backup->empty();
	sp=C[remap[pattern[s1+1]]];
	ep=C[remap[pattern[s1+1]]+1]-1;
	sp_r=sp; ep_r=ep;

	if (doNTimesForwardSearch(pattern,s1+2,s2,sp,ep,sp_r,ep_r))  //[s1+1..s2]
		for (i=s1;;i--) {
			pushRange();
			if (doNTimesBackSearch(pattern,s1,i+1,sp,ep,sp_r,ep_r))  //[i+1...s1...s2]
				for (e=1;e<sigma;e++) {
					pushRange();
					if (e!=remap[pattern[i]])
						if (searchHelper(e,sp_r,ep_r,T_bwt,sp,ep))
							if ((i==0) || doNTimesBackSearch(pattern,i-1,0,sp,ep,sp_r,ep_r)) // [0..]e[s1..s2]
								for (uint32_t j=s2+1;j<=m;j++) {
									pushRange();
									if (doNTimesForwardSearch(pattern,s2+1,j-1,sp,ep,sp_r,ep_r))
										for (e2=1;e2<sigma;e2++) {
											pushRange();
											if (e2!=remap[pattern[j]])
												if (searchHelper(e2,sp,ep,T_bwt_reverse,sp_r,ep_r))
													if (doNTimesForwardSearch(pattern,j+1,m,sp,ep,sp_r,ep_r)) {
														singleResult = (uint32_t*)malloc(2*sizeof(uint32_t));
														singleResult[0]=sp;
														singleResult[1]=ep;
														resultsList->push_back(singleResult);
													}
											popRange();
										}
									popRange();
								}
					popRange();
				}
			popRange();
			if (i==0) break; //i is unsigned
		}


	delete(backup);
	return resultsList;
}

bool
FM::SanityCheck(uint8_t* pattern,uint32_t m,uint8_t error_number) {
	if (m<=error_number)
		return false;
	int result=true;
	for (uint32_t i=0;i<m && result;i++) {
		uint8_t ch =remap[pattern[i]];
		if (!ch)
			result=false;
	}
	return result;
}

uint32_t
FM::count(uint8_t* pattern,uint32_t m) {
    uint8_t c = remap[pattern[m-1]]; /* map pattern to our alphabet */
    uint32_t i=m-1;
    uint32_t j = 1;

    uint32_t sp = C[c]; /* starting range in M from p[m-1] */
    uint32_t ep = C[c+1]-1;

	/* while there are possible occs and pattern not done */
    while (sp<=ep && i>=1) {
      c = remap[pattern[--i]]; /* map pattern to our alphabet */
      sp = C[c] + T_bwt->rank(sp-1+1,c); /* LF Mapping */
      ep = C[c] + T_bwt->rank(ep+1,c)-1; /* LF Mapping */
      j++;
    }

    if (sp<=ep) {
      return ep-sp+1;
    } else {
      return 0;
    }
}

uint32_t*
FM::locate(uint8_t* pattern,uint32_t m,uint32_t* matches) {
    uint32_t* locations;
    uint8_t c =  remap[pattern[m-1]];
    uint32_t i=m-1;
    sdsl::rrr_rank_support<>* rs = new rrr_rank_support<>();
    sdsl::util::init_support<rrr_rank_support<>,rrr_vector<> >(*rs,sampled);

    /* count occs */
    uint32_t sp = C[c];
    uint32_t ep = C[c+1]-1;
    while (sp<=ep && i>=1) {
      c =  remap[pattern[--i]];
      sp = C[c] + T_bwt->rank(sp-1+1,c);
      ep = C[c] + T_bwt->rank(ep+1,c)-1;
    }

    if (sp<=ep) {
        /* determine positions */
        *matches = ep-sp+1;
        uint32_t locate=0;
        locations= (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
        i=sp;
        int32_t j,dist,rank;
        while (i<=ep) {
            j=i,dist=0;
            while (!((*sampled)[j])) {
				c = (*T_bwt)[j];
				rank = T_bwt->rank(j+1,c)-1;
                j = C[c]+rank; // LF-mapping
                ++dist;
            }
            uint32_t temp = suffixes[rs->rank(j+1)-1];
            locations[locate]= temp + dist;
            locate++;
            ++i;
        }
        /* locations are in SA order */
        std::sort(locations,locations+(*matches));
        delete(rs);
        return locations;
    } else {
      /* no matches */
      *matches = 0;
      delete(rs);
      return NULL;
    }
    delete(rs);
    return locations;
}

uint8_t*
FM::extract(uint32_t start,uint32_t stop)
{
    uint8_t* T;
	uint32_t m,j,skip,todo,dist;
	uint8_t c;
	sdsl::rrr_rank_support<> *rs = new rrr_rank_support<>();
	sdsl::util::init_support<rrr_rank_support<>,rrr_vector<> >(*rs,sampled);

	/* last text pos is n-2 */
	if(stop > (this->n-1) ) stop = n-2;
    if(start > stop) {
		return NULL;
	}

	m = stop-start+1; /* snippet len */
	T = (uint8_t*) safe_malloc( (m+1) * sizeof(uint8_t)  );

	/* determine start pos of backwards search */
	j = positions[(stop/samplerate)+1];

	/* determine distance from start pos to the text snippet we want */
	if ((stop/samplerate+1) == ((n-1)/samplerate+1))
	   skip = n-1 - stop;
	else
		if (((samplerate-stop)%samplerate)==0)
			skip= samplerate-1;
		else
			//skip = (samplerate-stop)%samplerate-1;
			skip = suffixes[rs->rank(j+1)-1]-stop-1;

	/* start the backwards search */
	todo = m;
	dist = 0;
	while(todo>0) {
		c = (*T_bwt)[j];
		j = C[c] + T_bwt->rank(j+1,c)-1;

		/* check if we are at the snippet */
		if(dist>=skip) {
			c = remap_reverse[c];
			T[todo-1] = c;
			todo--;
		}
		dist++;
	}

	/* terminate */
	T[m] = 0;
	delete(rs);
    return T;
}

uint8_t*
FM::reconstructText(uint32_t* size)
{
    uint8_t* T;
    uint8_t c;
    uint32_t j,i;

    T = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );

    j = I; /* I is sa[I] = 0 -> last sym in T */
    for(i=0;i<n;i++) {
        c = (*T_bwt)[j]; /* L[j] = c */
        T[n-i-1] = remap_reverse[c]; /* undo sym mapping */
        j = C[c]+T_bwt->rank(j+1,c)-1; /* LF-mapping: j = LF[j] */
    }

	if(T[n-1] == 0) *size = n-1;
    else *size = n;

    return T;
}



/*
 * Debug Functions
 */

void FM::deb() {
	printf("\n");
	for (uint32_t i=0;i<this->n;i++) {
		uint32_t sa = this->SA[i];
		printf("%d\t%d\t",i,sa);
		uint8_t bwt = (*(this->T_bwt_reverse))[i];
		if (bwt==0)
			printf("-NULL-\t");
		else
			printf("%c\t",remap_reverse[bwt]);
		for (uint32_t j=sa; j<this->n;j++) {
			uint8_t ch = X[j];
			if (ch==0)
				printf("-NULL-");
			else
				printf("%c",remap_reverse[ch]);
		}
		printf("\n");
	}
}

void FM::debr() {
	printf("\n");
	for (uint32_t i=0;i<this->n;i++) {
		uint32_t sa = this->SA_reverse[i];
		printf("%d\t%d\t",i,sa);
		uint8_t bwt = (*(this->T_bwt_reverse))[i];
		if (bwt==0)
			printf("-NULL-\t");
		else
			printf("%c\t",remap_reverse[bwt]);

		for (uint32_t j=sa; j<this->n;j++) {
			uint8_t ch = X_reverse[j];
			if (ch==0)
				printf("-NULL-");
			else
				printf("%c",remap_reverse[ch]);
		}
		printf("\n");
	}
}

