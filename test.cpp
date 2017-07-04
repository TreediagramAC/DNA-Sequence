#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

void translation(string input){
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
    cout<<acid[pos]<<endl;
}

int main(){
  string dna;
  cin>>dna;
  translation(dna);
}
