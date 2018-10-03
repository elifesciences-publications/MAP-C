#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>

using namespace std;

int main(int argc, char** argv) {
//  ifstream input(argv[1]);
  string line;

  while (getline(cin,line)) {

    stringstream s(line);
    string chr;
    int alnpos;
    string cigar;
    s >> chr;
    s >> alnpos;
    s >> cigar;

    int nts;
    char type;
    int readend = alnpos;
    stringstream cig(cigar);
    while (cig >> nts, cig >> type, !cig.eof()) {
      if (type=='M' | type=='D') {
        readend += nts;
      }
    }
    readend--;
    cout << chr << '\t' << alnpos << '\t' << readend << '\n'; 
//    cout << strain << '\t' << umi << '\t' << alnpos << '\t' << cigar << '\t' << startpos << '\t' << endpos << '\t' << seq.substr(0,startpos) << '[' << seq.substr(startpos,endpos-startpos) << ']' << seq.substr(endpos) << '\n'; 

  }
//  input.close();


}
