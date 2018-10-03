/* annotation_mutations.cpp
 * 
 * Takes a reference FASTA file and a SAM alignment file, and outputs for each read the 
 * number of total mutations, substitutions, deletions, and insertions per read
 *
 * ./annotate_mutations REFERENCE ALIGNMENTS
 * 
 * Seungsoo Kim
 * September 27, 2018
 */

#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;

vector<pair<int,string>> parsecigar(string cigar) {
	string s = cigar;
	smatch sm;
	regex e(R"((\d+)([MDI]))");
	vector<pair<int,string>> parsedcigar;
	while (regex_search(s,sm,e)) {
		parsedcigar.push_back(make_pair(stoi(sm.str(1)),sm.str(2)));
		s = sm.suffix().str();
	}
	return parsedcigar;
}

vector<string> parsemd(string md) {
	string s = md;
	smatch sm;
	regex e(R"(\d+|\D+)");
	vector<string> parsedmd;
	while (regex_search(s,sm,e)) {
		parsedmd.push_back(sm.str(0));
		s = sm.suffix().str();
	}
	return parsedmd;
}

int main(int argc, char** argv) {

	string line;
	string refseq;
	ifstream refseqf(argv[1]);
	getline(refseqf,line);
	getline(refseqf,refseq);
	refseqf.close();

	ifstream input(argv[2]);
	// variables to store entries in SAM file
	string readname;
	int flag;
	string alnchr;
	int alnstart;
	int mapq;
	string cigar;
	string matechr;
	int matealnstart;
	int insertlen;
	string seq;
	string quals;
	string tags;
	string md;
	while (getline(input,line)) {
		if (line.compare(0,1,"@")==0) {
			continue;
		}
		stringstream ss(line);
		// read in entries of SAM file
		ss >> readname;
		ss >> flag;
		ss >> alnchr;
		ss >> alnstart;
		ss >> mapq;
		ss >> cigar;
		ss >> matechr;
		ss >> matealnstart;
		ss >> insertlen;
		ss >> seq;
		ss >> quals;
		getline(ss,tags);
		// minimum MAPQ
		if (mapq < 1) {
			continue;
		}

		// read through tags until MD string is found
		stringstream tagss(tags);
		while (!tagss.eof()) {
			tagss >> md;
			if (md.compare(0,5,"MD:Z:") == 0) {
				md = md.substr(5);
				break;
			}
		}
		// parse CIGAR string and store as vector
		vector<pair<int,string>> cigarvector = parsecigar(cigar);
		// parse MD string and store as vector
		vector<string> mdvector = parsemd(md);
		// counter variables for position in CIGAR vector and MD vector
		int cigariter = 0;
		int mditer = 0;
		// vector to store mutations found
		vector<string> muts;
		// counter variables for position in read and reference sequence
		int readpos = 0;
		int refpos = alnstart - 1;
		// counter variables for mutation counts
		int mutcount = 0;
		int subcount = 0;
		int delcount = 0;
		int inscount = 0;
	
		// special case where alignment starts after beginning of reference - infer deletion
		if (alnstart > 1) {
			muts.push_back("1D" + to_string(alnstart-1));
			delcount++;
			mutcount++;
		}

		// loop until end of reference sequence is reached
		while (readpos < seq.length()) {
			// if substitution(s)
			if (!isdigit(mdvector[mditer][0]) && (mdvector[mditer].compare(0,1,"^") != 0) && (cigarvector[cigariter].second.compare("M") == 0)) {
				for (int i = 0; i < mdvector[mditer].length(); i++) {
					muts.push_back(to_string(refpos + 1) + seq[readpos]);
					refpos++;
					readpos++;
					subcount++;
					mutcount++;
					cigarvector[cigariter].first--;
				}
				mditer++;
			}
			// if deletion
			else if (cigarvector[cigariter].second.compare("D") == 0) {
				muts.push_back(to_string(refpos + 1) + "D" + to_string(cigarvector[cigariter].first));
				refpos += cigarvector[cigariter].first;
				cigariter++;
				mditer++;
				delcount++;
				mutcount++;
			}
			// if insertion
			else if (cigarvector[cigariter].second.compare("I") == 0) {
				muts.push_back(to_string(refpos + 1) + "+" + seq.substr(readpos,cigarvector[cigariter].first));
				readpos += cigarvector[cigariter].first;
				cigariter++;
				inscount++;
				mutcount++;
			}
			// if match
			else {
				int mdint = stoi(mdvector[mditer]);
				if (cigarvector[cigariter].first < mdint) {
					refpos += cigarvector[cigariter].first;
					readpos += cigarvector[cigariter].first;
					mdvector[mditer] = to_string(mdint-cigarvector[cigariter].first);
					cigariter++;
				}
				else if (cigarvector[cigariter].first > mdint) {
					refpos += mdint;
					readpos += mdint;
					cigarvector[cigariter].first -= mdint;
					mditer++;
				}
				else {
					refpos += mdint;
					readpos += mdint;
					cigariter++;
					mditer++;
				}
			}
		}
		if (refpos < refseq.length()) {
			muts.push_back(to_string(refpos + 1) + "D" + to_string(refseq.length() - refpos));
			mutcount++;
			delcount++;
		}
		cout << readname << '\t' << mutcount << '\t' << subcount << '\t' << delcount << '\t' << inscount << '\t';
		if (muts.size() > 0) {
			for (int i = 0; i < muts.size() - 1; i++) {
				cout << muts[i] << ',';
			}
			cout << muts.back() << '\n';
		}
		else {
			cout << "WT" << '\n';
		}
	}
	input.close();
}
