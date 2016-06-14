
#include <vector>
#include "./helpers.hpp"


#include <algorithm>
#include <numeric>
#include<boost/foreach.hpp>

#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>

using namespace  boost::accumulators;



template<class Innercode, class Outercode, class DNAmap>
class EnDecode {
	public:
	Innercode innercode;
	Outercode outercode;
	DNAmap dnamap;
	vector< vector< typename Innercode::Symbol > > groundtruth;
	unsigned l;
	unsigned m;
    unsigned numblocks; 

	EnDecode(){};
	EnDecode(const Innercode& innercode_, const Outercode& outercode_, const DNAmap& dnamap_,unsigned l_){
		innercode = innercode_;
		outercode = outercode_;
		dnamap = dnamap_;
		l = l_;
		assert(innercode.n>l);
		m = innercode.k - l;
	
	};

	void encode(string& str, vector<string>& urn);

	void decode(string& str, const vector<string>& drawnseg);
};




////////////////////////////
template<class Innercode, class Outercode, class DNAmap>
void EnDecode<Innercode,Outercode,DNAmap>::encode(string& str, vector<string>& urn){

    //const int prime = pfe::prime; 

    // convert the string to a representation in GF(47) 
    vector<typename Innercode::Symbol> data_b47; 
    char2pfe<typename Innercode::Symbol>(str,data_b47); // data in pfe vector 
 
    // add the filelength at the beginning: 
    vector<typename Innercode::Symbol> data_b47_lv = tobase<typename Innercode::Symbol>(str.size(),47,6); 
    data_b47.insert(data_b47.begin(),data_b47_lv.begin(),data_b47_lv.end()); // not very efficient.. 

	// determine the number of blocks
    if(data_b47.size() % (m*outercode.k) != 0) { 
        numblocks = data_b47.size()/(m*outercode.k)+1; 
		unsigned oldsize = data_b47.size();
		data_b47.resize(numblocks*m*outercode.k); 
    	// fill the rest with zeros.. 
		for(unsigned i=oldsize;i<data_b47.size();++i) data_b47[i] = typename Innercode::Symbol(0);
	} else numblocks = data_b47.size()/(m*outercode.k); 
     
    urn.resize(numblocks*outercode.n); // the number of DNA segments 
 
    unsigned b = 0; 
 
    for(unsigned b=0;b<numblocks;++b){ // run over blocks 
     
        vector<typename Outercode::Symbol> infvec(outercode.k); // outer information vector 
    	
		
		for(unsigned i=0;i<outercode.k;++i) {
			vector<typename Innercode::Symbol> tmpefevec(data_b47.begin()+b*m*outercode.k+i*m,data_b47.begin()+b*m*outercode.k+i*m+m);  
			assert(tmpefevec.size() == m);
			infvec[i] = typename Outercode::Symbol(tmpefevec);
			//cout << infvec[i] << endl;
		}
        
		vector<typename Outercode::Symbol> c; // outer codeword 
		outercode.RSencode(infvec,c); // encode the information, obtain outer codeword
 

        for(unsigned i =0;i<outercode.n;++i){ // run over elements in a block 
             
            // compute and append index of each  
            vector<typename Innercode::Symbol> indexpfe = tobase<typename Innercode::Symbol>(i+outercode.n*b,47,l); 
			vector<typename Innercode::Symbol> infveci = c[i].tovec();
			//append index at back: 
			//infveci.insert(infveci.end(),indexpfe.begin(),indexpfe.end());
			//append index at front:
			infveci.insert(infveci.begin(),indexpfe.begin(),indexpfe.end());
			assert(infveci.size()== innercode.k );
			//cout << "infvec+index: " << infveci << endl;

            // encode the inner code, write it to urn[i] 
            vector<typename Innercode::Symbol> icw;  
			innercode.RS_shortened_encode(infveci,icw);
			dnamap.cw2frag(icw,urn[i+outercode.n*b]); 
        }
		cout << "encoded block " << b << " of " << numblocks << endl;
    } 
	
	
	// return numblocks;
	
} 


template<class T>
bool pairCompare(const T & x, const T & y) {
  return x.second < y.second; 
}


////////////////////////////
template<class Innercode, class Outercode, class DNAmap>
void EnDecode<Innercode,Outercode,DNAmap>::decode(string& str, const vector<string>& drawnseg){ 


	//// inner decoding

	// to keep track    
	pair<unsigned,unsigned> erctric(0,0);
	vector<int> frag_freqs(numblocks*outercode.n, 0); // to count the frequencies of each DNA sequence
	unsigned ic_decerr = 0; // counts the inner decoding errors
	unsigned errctr = 0;
	unsigned numfrag = 0;

	typedef typename Innercode::Symbol pfe;


	vector<vector<pfe> > segm_ordered(numblocks*outercode.n, vector<pfe>());

	typedef map<vector<pfe>, unsigned > seqctrmap;
	vector< seqctrmap > segm_ordered_map(numblocks*outercode.n, seqctrmap() );
	vector<unsigned> curerr(numblocks*outercode.n, innercode.n ); // error of the element in segm_ordered_map;

	// run over all sequences
	BOOST_FOREACH(const string& sLine, drawnseg){

		vector<pfe> rececw; // the received block 
		pair<unsigned,unsigned> erctrictmp;
		vector<pfe> irec(innercode.k); // recovered information from inner codeword

		dnamap.frag2cw(rececw,sLine); // map DNA segment to pfe vector
		erctrictmp = innercode.RS_shortened_decode(irec,rececw);
		

		if(erctrictmp.second != innercode.n_u){ // else we have a decoding error
			erctric.first += erctrictmp.first;
			erctric.second += erctrictmp.second;
			unsigned ind = frombase<pfe>(vector<pfe>(irec.begin(),irec.begin()+l),47);
			if(ind < frag_freqs.size()){
				if(erctrictmp.second < curerr[ind]){ // discard list, 
					curerr[ind] = erctrictmp.second;
					segm_ordered_map[ind].empty();
					segm_ordered_map[ind].insert(make_pair(irec,1));
				} else if (erctrictmp.second == curerr[ind]){
					typename seqctrmap::iterator cur = segm_ordered_map[ind].find(irec);
					if(cur != segm_ordered_map[ind].end()){ // segment already appeared  
						cur->second++;  
					} else {
						segm_ordered_map[ind].insert(make_pair(irec,1));
					}
				
				}

				// compare to groundtruth whether there is a decoding error
				
				if(! groundtruth.empty() ){ 
					if(groundtruth[ind] != irec) {
						errctr++;
					} else {
						// this fragment has been decoded without error
						frag_freqs[ind]++;
					}
				}
			}else{ // if the index is not in the appropriate range, there must be an error
				errctr ++;
			}
		} else {ic_decerr++;}
		numfrag++;
	}




for(unsigned ind=0;ind<segm_ordered_map.size();++ind){
	if(! segm_ordered_map[ind].empty() ){
	
		// take the one (out of all segments with smallest error) with largest frequency as the estimate
		typename seqctrmap::iterator max = max_element(segm_ordered_map[ind].begin(), segm_ordered_map[ind].end(), pairCompare<typename seqctrmap::value_type> );
		segm_ordered[ind] = max->first;

/*
		if( segm_ordered[ind] != groundtruth[ind]){ 
			cout << ind << "	" << curerr[ind] <<  ": ";
			for (seqctrmap::iterator it = segm_ordered_map[ind].begin();it != segm_ordered_map[ind].end();++it){
				
				cout << "(" << it->second << ", " << (groundtruth[ind] == it->first) << ") ";
			}
			cout << endl;
		}
*/	
	}
}
	



cout << "fragm. error prob: " << (float) errctr / (float) numfrag << endl;
cout << "numer of sequences: " << numfrag << endl;
cout << "inner code: " << (float)erctric.first / (float) numfrag << " erasures on average corrected" << endl;
cout << "inner code: " << (float)erctric.second/ (float) numfrag << " errors on average corrected" << endl;
cout << "inner code, dec. errors: " << ic_decerr << endl;




	unsigned total_erasure_ctr = 0;
	unsigned total_error_ctr = 0;
	unsigned total_blockerror_ctr = 0;
	//// determine the number of errors and erasures in each block
	for(unsigned b=0;b<numblocks;++b){
		unsigned block_erasure_ctr = 0;
		unsigned block_error_ctr = 0;
		for(unsigned i =0;i<outercode.n;++i){ // run over elements in a block 
			if( ! segm_ordered[b*outercode.n+i].empty()){
            	//compare to groundtruth
				if(! groundtruth.empty() ){ 
					if(groundtruth[b*outercode.n+i] != segm_ordered[b*outercode.n+i]) block_error_ctr++;
				}
			} // else we have an erasure
            else{
                block_erasure_ctr++;
            }	
		}
		total_erasure_ctr += block_erasure_ctr;
		total_error_ctr += block_error_ctr;

		cout << b <<"-> errors: " << block_error_ctr << "	erasures: " << block_erasure_ctr <<"	total:" << block_error_ctr*2+block_erasure_ctr << endl;
		
		if(block_error_ctr*2+block_erasure_ctr > outercode.n-outercode.k){
			total_blockerror_ctr++;
		}
		
	}




	cout << "num block errors: " << total_blockerror_ctr << endl;
	cout << "OC: symbol err. prob.: " << (float) total_error_ctr / (float) (outercode.n*numblocks) << endl;
	cout << "OC: symbol era. prob.: " << (float) total_erasure_ctr / (float) (outercode.n*numblocks) << endl;

	cout << "OC: symbol err. prob.: " << total_error_ctr  << endl;
	cout << "OC: symbol era. prob.: " << total_erasure_ctr << endl;
	//cout << "(" << frac << "," << (double) total_error_ctr/ (double) numit << ")" << endl;

	// decode the blocks

	unsigned filelengthstr;
	for(unsigned b=0;b<numblocks;++b){ // run over blocks 
        vector<typename Outercode::Symbol >reccw(outercode.n);
		unsigned erasure_ctr = 0;
		for(unsigned i =0;i<outercode.n;++i){ // run over elements in a block 
			if( ! segm_ordered[b*outercode.n+i].empty()){
				vector<pfe> tmppfevec = vector<pfe>(segm_ordered[b*outercode.n + i].begin()+l, segm_ordered[b*outercode.n + i].end() );
				reccw[i] = typename Outercode::Symbol( tmppfevec );
			} // else we have an erasure
			else{
				reccw[i] = typename Outercode::Symbol();
				erasure_ctr++;
			}
		}

		cout << erasure_ctr << " many erasures in block " << b << " of length " << outercode.n <<  endl;
		
		
		// decode outer codword
    
		
		vector<typename Outercode::Symbol> infvecrec;
		//pair<unsigned,unsigned> erctroc = RS_decode_spec<efe,dftefe>(infvecrec,reccw,a,k,n,dftefe_d);
		pair<unsigned,unsigned> erctroc = outercode.RS_decode_spec(infvecrec,reccw);
		cout << "decoded block " << b << endl;
		cout << "outer code: " << erctroc.first << " errasures, " << erctroc.second << " errors corrected"<< endl;
		

		for(unsigned i=0;i<outercode.k;++i){
			string tmp;
			vector<pfe> tmpinfvec;
			if( !((b==0) && (i==0)) ){
				tmpinfvec = infvecrec[i].tovec(); 
			}else{
				//determine filelength first
				tmpinfvec = infvecrec[i].tovec(); 
        		filelengthstr = frombase<pfe>(vector<pfe>(tmpinfvec.begin(),tmpinfvec.begin()+6  ),47); 
				tmpinfvec = vector<pfe>(tmpinfvec.begin()+6,tmpinfvec.end());
			}
			pfe2char<pfe>(tmp,tmpinfvec);
			str.append(tmp);
		}
	}


	str.resize(filelengthstr);


}




