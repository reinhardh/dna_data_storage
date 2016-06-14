/*
Program to translate text to DNA and vice versa
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/program_options.hpp>
#include "../include/PFE.hpp"
#include "../include/EFE.hpp"
#include "../include/DFT.hpp"
#include "../include/helpers.hpp"
#include "../include/encodedecode.hpp"
#include "../include/ReedSolomon.hpp"
#include <string>
#include <fstream>
#include <streambuf>


using namespace std;



// define static variables 

typedef int el_t;
template<> int PFE<el_t>::prime = 47;
template<> boost::bimap<el_t,el_t> PFE<el_t>::exp_el = boost::bimap<el_t,el_t>();

typedef PFE<el_t> pfe;

template<> polynomial<pfe> EFE<pfe>::prim_poly = polynomial<pfe>();
template<> unsigned EFE<pfe>::m = 30;

 template<class cf_t> 
    ostream &operator<<(ostream &stream, const vector<cf_t>& x) 
    { 
	for(unsigned i = 0; i<x.size(); ++i) stream << x[i] << "  ";
	return stream; // must return stream 
    }; 



namespace po = boost::program_options;

int main(int ac, char* av[])
{

// command line options

int numblocks;
string infile;
string outfile;

int opt;
po::options_description desc("Allowed options");
desc.add_options()
    ("help", "produce help message")
    ("encode", "encode")
    ("decode", "decode")
    ("disturb", "draw uniformly at random from the input lines, add errors to each line")
	("input",po::value<string>(&infile)->default_value(""),"inputfile")	
	("output",po::value<string>(&outfile)->default_value(""),"outputfile")	
	("numblocks",po::value<int>(&numblocks)->default_value(0),"numblocks")	
;

po::variables_map vm;
po::store(po::parse_command_line(ac, av, desc), vm);
po::notify(vm);    

if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
}




const unsigned l = 3;
const unsigned m = 30; // EFE<pfe>::m; // GF(p^m)
DNAmap<pfe> dnamap;

// initialize the inner code	
pfe::initialize_exp_el(33);
typedef EFE<pfe> efe;

const unsigned M = 50000;
const int prime = PFE<el_t>::prime;

const unsigned N = 39;
const unsigned N_u = 46; // the underlying length of the shortened inner code
const unsigned K = 33;


// primitive element of the inner code 
const pfe A = pfe(33);
typedef DFT_PRIM<pfe> dftpfe;
dftpfe dftpfe_d(N,A);

RScode<pfe,dftpfe> innercode(N,K,A,dftpfe_d,N_u);


// initialize the outer code
const unsigned n = 713; // DFT length, divides 47^33-1
const unsigned P = 23;
const unsigned Q = 31; // P*Q = n
const unsigned k = 594;//	

EFE<pfe>::m = m;

// initialize the primitive polynomial
pfe ppvv[m+1] = {1,13,15,22,45,  19,34,10,17,5,41,21,  12,3,23,17,27,8,  34,32,7,40,41,1,  32,26,24,32,37,5,1}; // for m=30	
vector<pfe> ppv(ppvv, ppvv+m+1);
EFE<pfe>::prim_poly = polynomial<pfe>(ppv);
// initialize the element of order n, the Fourier kernel  
pfe aa[m] = {26,35,41,0,24,3,31,8,12,10,9,24,44,11,24,9,43,14,26,32,7,13,15,29,25,38,8,24,40,32}; // element of order 713
vector<pfe> avv(aa, aa+m); // 
efe a = efe(avv); // this is an element of order n

typedef DFT_FFT<efe> dftefe;
dftefe dftefe_d(n,a,P,Q); // Fourier transform for the outer code 

RScode<efe,dftefe> outercode(n,k,a,dftefe_d);

// encoder/decoder 
EnDecode<RScode<pfe,dftpfe>, RScode<efe,dftefe>, DNAmap<pfe> > endecode(innercode,outercode,dnamap,l);



/////////////////////// encode 
if (vm.count("encode")) {
	
	cout << "start encoding.." << endl;

	if(infile == "" || outfile ==""){
		cout << "in/outfile not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;


	// read data
	std::ifstream t(infile.c_str());
	std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());

	// encode
	vector<string> urn(n*numblocks);
	endecode.encode(str, urn);
	numblocks = endecode.numblocks; 
	
	cout << "encoded " << str.size() << " Bytes to " << numblocks << " blocks, resulting in "
	<< n*numblocks << " DNA segments of length " << N << " each." << endl;
    
	ofstream out;
	out.open(outfile.c_str());
	for(unsigned i=0;i<urn.size();++i) out << urn[i] << endl;
	out.close();
	return 0;
}

/////////////////////// decode 
if (vm.count("decode")) {
	//cout << "generate DFT lookup table"<< endl; 
	typedef DFT_FFT<efe> dftefe;
	dftefe dftefe_d(n,a,P,Q); // Fourier transform for the outer code 
	
	if(infile == "" || outfile =="" || numblocks==0 ) {
		cout << "in/outfile/numblocks not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;
	cout << "numblocks: "<<numblocks << endl; 

	vector<string> drawnseg;

	
	string sLine = "";
	ifstream in;
	in.open(infile.c_str());
	while (!in.eof()){
		getline(in, sLine);
		drawnseg.push_back(sLine);
	}
	drawnseg.resize(drawnseg.size()-1); // erase the last, empty line

	string recstr;
	cout << "start decode.." << endl;	
	endecode.numblocks = numblocks;
	endecode.decode(recstr, drawnseg);
	
	ofstream out;
	out.open(outfile.c_str());
	out << recstr;
	out.close();
	
}
////////////////////////// disturb
if (vm.count("disturb")) {
	
	if(infile == "" || outfile =="") {
		cout << "in/outfile not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;
	
	vector<string> urn;
	
	
	string sLine = "";
	ifstream in;
	in.open(infile.c_str());
	while (!in.eof()){
		getline(in, sLine);
		urn.push_back(sLine);
	}
	urn.resize(urn.size()-1); // erase the last, empty line
	
	boost::mt19937 rng; 	
	boost::uniform_int<> unif(0,urn.size()-1); // distribution that maps to 0,..,urn.size()-1
	boost::uniform_int<> unif_N(0,N-1); // distribution that maps to 0,..,N-1
	boost::uniform_int<> unif_4(0,4-1); // distribution that maps to 0,..,N-1
	boost::bernoulli_distribution<> bern(0.013);
	boost::bernoulli_distribution<> faircoin(0.5);
	
	char nucl[] = "ACGT";

	vector<string> drawnseg(M);
	
	//write data
	ofstream out;
	out.open(outfile.c_str());
	
	// draw M times
	string tmpstr;
	for(unsigned i=0;i<M;++i){
		unsigned randind = unif(rng);
		tmpstr = urn[randind];
		
		// introduce errors
		
		for(unsigned j=0;j<2;++j){
			tmpstr[unif_N(rng)] = nucl[unif_4(rng)];
		}
	
		unsigned ctr = 0;
		for(unsigned j=0;j<tmpstr.size();++j)
			if(bern(rng)) {
				tmpstr[j] = nucl[unif_4(rng)];
				ctr++;
			}
			//cout << ctr << endl;
		//if(faircoin(rng)) flipvecdir(tmpstr); // flip every second 

		// write to file
		out << tmpstr;
		if(i!= M-1) out << endl;
	}
	out.close();

}

}
