#include "bio.h"

#include <iostream>
using std::cout; using std::cin; using std::endl;
#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <algorithm>
using std::swap;
using std::reverse;
using std::distance;
using std::find;



/*
This function should return true if and only if
every character in the input is one of ATCG.
*/
bool is_valid_DNA_sequence(const string & input) {
  long check=0;
  for (auto chr : input){
    if (chr =='A' or chr =='T' or chr =='C' or chr =='G'){
      check++;
    }
  }
  if (check==input.length())
    return true;
  else{
    return false;
  }
}

//This function should reverse all nucleotide within sequence.
string reverseStr(string &str){
    int n = str.length();
    for (int i=0; i<n/2; i++)
       swap(str[i], str[n-i-1]);
    return str;
}

//This function should replace all required nucleotide within sequence.
string rplcStr(string str){
  for(int i =0;i<str.length();i++){
    if (str[i]=='A')
      str[i] = 'T';
    else if (str[i]=='T')
      str[i] = 'A';
    else if (str[i]=='C')
      str[i] = 'G';
    else if (str[i]=='G')
      str[i] = 'C';
  }
  return str;
}

/*
This function should replace all required nucleotide within RNA sequence to
get orinail RNA sequence or antiparallel RNA sequence.
*/
string tranTU(string str){
  for(int i =0;i<str.length();i++){
    if (str[i]=='T')
      str[i] = 'U';
  }
  return str;
}

vector<vector<string>> get_offset_sequence_in_vector(const string & input){
    string original=input,off0,off1,off2;
    vector<string> pff0,pff1,pff2;
    vector<vector<string>> onetype;
    off0 = original;
    off1 = off0.substr(1, off0.length());
    off2 = off1.substr(1, off1.length());
    for (int i=0;i<off0.length();i+=3){
      if (off0.substr(i,3).length()==3)
        pff0.push_back(off0.substr(i,3));
    }
    for (int i=0;i<off1.length();i+=3){
      if (off1.substr(i,3).length()==3)
        pff1.push_back(off1.substr(i,3));
    }
    for (int i=0;i<off2.length();i+=3){
      if (off2.substr(i,3).length()==3)
        pff2.push_back(off2.substr(i,3));
    }
    onetype.push_back(pff0);
    onetype.push_back(pff1);
    onetype.push_back(pff2);
    return onetype;
}

string translation(string input){
  string result;
  vector<string> codons {"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
  "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
  "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
  "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
  "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
  "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
  "GUG", "UAG", "UGA", "UAA"};
  vector<string> acid {"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
  "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
  "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
  "P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
  "Y", "V", "V", "V", "V", "*", "*", "*"};
  vector<string>::iterator it;
  it=find(codons.begin(),codons.end(),input);
  int pos = distance(codons.begin(), it);
  result = acid[pos];
  return result;
}

vector<int> positions_in_string(string amino,char wanted){
  vector<int> position;
  for(int i=0;i<amino.length(); i++){
    if(amino[i]==wanted)
      position.push_back(i);
  }
  return position;
}


/*
This function should calculate the reverse complement DNA sequence.

The first argument is the sequence, the second argument is a pointer to
an empty string, which you should modify to store the result.

This is obtained by reversing the input sequence and swaping each
nucleotide/letter with it's complement:
A <-> T
C <-> G

Example:
input = AAATTCGGGG
reverse = GGGGCTTAAA
reverse complement = CCCCGAATTT
*/
void get_reverse_complement_sequence(const string & input,  string * const output) {
  string original = input;
  string reverse;
  reverse = reverseStr(original);
  *output = rplcStr(reverse);
}


/*
This function should return the RNA transcript from a DNA sequence.

A RNA transcript is the reverse complement of the DNA sequence, but RNA
has U (uracil) instead of T (thiamine).

Make sure you don't have redundant code with the
get_reverse_complement_sequence function.
*/
string get_RNA_transcript(const string &input) {
    string result;
    get_reverse_complement_sequence(input, &result);
    result = tranTU(result);
    return result;
}

/*
This function should return a vector of vector of strings with each possible RNA
reading frame from the given DNA sequence.

There are three possible reading frames (because the genetic code has three
nucleotides per amino acid) in each direction (you can also transcribe DNA in
the reverse complement direction, called the antiparallel strand).

Order the sequences like so:
1: Original (0 offset)
2: Original (1 offset)
3: Original (2 offset)
4: Antiparallel (0 offset)
5: Antiparallel (1 offset)
6: Antiparallel (2 offset)

With in the input sequence of: AATTCCCGAAA
Original RNA transcript = UUUCGGGAAUU
Antiparallel RNA transcript = AAUUCCCGAAA

The offsets (starting at pos 0, 1, and 2) of the two RNA transcripts
UUUCGGAAUU
UUCGGGAAUU
UCGGGAAUU
AAUUCCCGAAA
AUUCCCGAAA
UUCCCGAAA

Instead of returning a vector of 6 strings, break each string into a vector
of length 3 strings (called codons) These codons will be useful
for the next translation step.

UUUGCCCAAUU -> {"UUU", "CGG", "GAA"}
// drop any remaining letters that don't fill a codon
UUGCCCAAUU -> {"UUC", "GGG", "AAU"}
UGCCCAAUU -> {"UCG", "GGA", "AUU"}
AAUUCCCGAAA -> {"AAU", "UCC", "CGA"}
AUUCCCGAAA -> {"AUU", "CCC", "GAA"}
UUCCCGAAA -> {"UUC", "CCG", "AAA"}

*/
vector<vector<string>> get_reading_frames_as_codons(const string &input) {
    vector<vector<string>> result,voriginal,vanti;
    string original,antiparallel;
    original = get_RNA_transcript(input);
    antiparallel = tranTU(input);
    voriginal = get_offset_sequence_in_vector(original);
    vanti = get_offset_sequence_in_vector(antiparallel);
    for(int i=0;i<3;i++){
      result.push_back(voriginal[i]);
    }
    for(int i=0;i<3;i++){
      result.push_back(vanti[i]);
    }
    return result;
}

/*
This function translates/converts a vector<string> (vector of codons) into a
string of amino acids using the genetic code
(see https://en.wikipedia.org/wiki/Genetic_code).

For example, the codons:
{"UUU", "GCC", "CAA"}
translates to:
F (Phenylalanine), A (Alanine), Q (Glutamine)
abreviated:
FAQ

To make your lives easier, here's a list of the possible codons:
"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA"

And there corresponding amino acids ("*" represents STOP codons,
more on them later):

"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*"
*/
string translate(const vector<string> &codon_sequence) {
    string result;
    for (auto codons:codon_sequence){
      result += translation(codons);
    }
    return result;
}

/*
This function takes a DNA sequence and returns the longest possible
amino acid sequence / protein that is encoded by that sequence
(open reading frame). A valid open reading frame begins with the
codon AUG (the amino acid, Methionine (M)) and runs until a stop codon (*)
is encountered. There may be multiple open reading frames in a sequence, and
you need to check all six reading frames in order given by
get_reading_frames_as_codons. If there are ties for longest, favor the first
one found.

Return the longest open reading frame as an amino acid sequence. It must start
with an 'M' and end with a '*' with no other '*''s within.
*/
string get_longest_open_reading_frame(const string &DNA_sequence) {
    string longest, before;
    string sub;
    vector<string> result;
    for (vector<string> frame : get_reading_frames_as_codons(DNA_sequence)) {
      before =  translate(frame);
      if ((before.find('M')!=string::npos) and (before.find('*')!=string::npos)){
        vector<int> positionM, positionS;
        positionM=positions_in_string(before,'M');
        positionS=positions_in_string(before,'*');
        for (auto M :positionM){
          for (auto star: positionS){
            if (star>M){
              sub = before.substr(M,star-M+1);
              result.push_back(sub);
              break;
            }
          }
        }
      }
      else{
        longest = "";
      }
    }
    sort(result.begin(), result.end(), [](const string &s1, const string &s2) {return s1.size() < s2.size(); });
    if (result.empty()==false){
      longest = result[result.size()-1];
    }
    return longest;
}
