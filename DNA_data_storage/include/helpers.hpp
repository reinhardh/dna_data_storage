


#ifndef HELPERSDNA 
#define HELPERSDNA

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <boost/bimap.hpp>
#include <string>
#include <map>
#include<vector>

using namespace std;


template<class vectt>
unsigned hamming_distance(vectt& v1,vectt& v2){
	assert(v1.size()=v2.size());
	unsigned hd = 0;
	for(unsigned i=0;i<v1.size();++i)
		if(v1[i] != v2[i]) hd++;
	return hd;
}

template<class vectt>
void flipvecdir(vectt& v){
	vectt tmp = v;
	for(unsigned i=0;i<v.size();++i){
		v[i] = tmp[v.size()-1-i];
	}
}


// convert decimal number to a number with another base
template<class base_t>
vector<base_t> tobase(unsigned nu, unsigned base, unsigned size){
	
	vector<base_t> nn(size,(base_t) 0 ); 
	if(nu == 0) return nn;
	
	unsigned div = nu;
	unsigned i = 0;
	while(div != 0){
		nn[i] = base_t(div % base); // rest
		div = div/base;	
		i++;
	}
	return nn;
};


// convert vector of numbers in some basis to 
template<class base_t>
unsigned frombase(const vector<base_t>& nn,unsigned base){
	int nu = 0;
	//base_t adf =  nn[0];
	//int adfi =  (int) adf;
	unsigned basepow = 1;
	for(unsigned i=0;i<nn.size();++i) {
		
		base_t adf = nn[i];
		nu += ((unsigned)adf) * basepow;//((unsigned) (adf)) * ((int) pow((double) base,(int) i));
		basepow *= base;
	}
	return nu;
};

// converts each 2 characters in data_char to 3 elements of GF(47), which are appended at 
template<class pfe>
void char2pfe(string& data_char, vector<pfe>& data_b47){
    if((data_char.size() % 2) != 0 )data_char.push_back('\n');
   
    data_b47.resize(0);
	// convert vector of char to vector with base 47;
    // every 2 char's are coded to three pfe's 
    for(unsigned i=0;i<data_char.size()/2; ++i){
        // the unsigned char type below is important, otherwise errors occur when converting..
		unsigned indec = frombase(vector<unsigned char>(data_char.begin()+i*2, data_char.begin()+i*2+2),256);
        vector<char> adf(data_char.begin()+i*2, data_char.begin()+i*2+2);
        //for(unsigned j=0;j<adf.size();++j) cout << adf[j];
        //cout << endl;
        
        vector<pfe> tmpvec = tobase<pfe>(indec,47,3);
        data_b47.insert(data_b47.begin()+i*3,tmpvec.begin(),tmpvec.end());
    	//cout << i*3 << "  " << data_b47.size() << endl;
	}
};


template<class pfe>
void pfe2char(string& data_char, vector<pfe>& data_b47){
    
    assert(data_b47.size() % 3  == 0);

    data_char.resize(2*(data_char.size()/3) );
    // convert vector of char to vector with base 47;
    // every 2 char's are coded to three pfe's 
    for(unsigned i=0;i<data_b47.size()/3; ++i){
        unsigned indec = frombase(vector<pfe>(data_b47.begin()+i*3, data_b47.begin()+i*3+3),47);
        //vector<char> adf(data_b47.begin()+i*3, data_b47.begin()+i*3+3);
        //for(unsigned j=0;j<adf.size();++j) cout << adf[j];
        //cout << endl;
        
        vector<char> tmpvec = tobase<char>(indec,256,2);
        data_char.insert(data_char.begin()+i*2,tmpvec.begin(),tmpvec.end());
    }

};



/*
maps each element of GF(47) to a string with letters {A,C,G,T} where 
*/

template<class PFE>
class DNAmap{
private: 
	//typedef boost::bimap< PFE , std::string > mapdna_type;
	//mapdna_type mapdna;
	
	map<PFE,string> pfetostr;
	map<string,PFE> strtopfe;

public: 
	DNAmap(){
		//initialize map
		unsigned prime = 47;
		
		char nucl[] = "ACGT"; // nucleotides	

		vector<string> allpos(4*4*4);
		for(unsigned i=0;i<allpos.size();++i){
			char cur[] = "AAA";
			cur[0] = nucl[i % 4];
			cur[1] = nucl[(i/4) % 4];
			cur[2] = nucl[(i/16) % 4];
			allpos[i] = string(cur);
		}
		
		unsigned j = 0;
		for(unsigned i=0;i<prime;++i){
			while( allpos[j][1] == allpos[j][2]) j++;
			//mapdna.insert( typename mapdna_type::value_type(PFE(i), allpos[j] ));
			pfetostr[PFE(i)] = allpos[j];
			strtopfe[allpos[j]] = PFE(i);
			
			j++;
		}
		//print_map(mapdna.right, " ", cout);
	}

	/*
	void printmap(){
		for(auto it = pfetostr.cbegin(); it != pfetostr.cend(); ++it){
    		std::cout << it->first << " " << it->second << std::endl;
		}
	}
	*/

	// codeword to DNA fragment
	void cw2frag(const vector<PFE>& cw, string& frag ){
		frag.resize(0);
		for(unsigned i=0;i<cw.size();++i){
			frag.append(pfetostr[cw[i]]); //mapdna.left.at(cw[i]); 
		}
	}
	
	// DNA fragment to codeword
	void frag2cw(vector<PFE>& cw, const string& frag ){
		assert( (frag.size() % 3) == 0   );
		cw.resize(frag.size()/3);
		for(unsigned i=0;i<cw.size();++i){ 
		
			//char charartmp[] = "AAA"; // nucleotides	
			//charartmp[0] = frag[i*3];
			//charartmp[1] = frag[i*3+1];
			//charartmp[2] = frag[i*3+2];
			//string tmp(charartmp);
			string tmp = string(frag.begin()+i*3, frag.begin()+i*3+3  );
				
			//print_map(mapdna.right, " ", cout);
			
			//cout <<  mapdna.right.at(tmp )  << endl; 
			//cw[i] = PFE( mapdna.right.at(tmp ) ); 
			cw[i] = strtopfe[tmp];
			//cout << cw[i] << endl;
		}
	}







};

#endif
