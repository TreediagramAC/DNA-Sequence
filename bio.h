#ifndef BIO_H
#define	BIO_H
#include <string>
#include <vector>
using std::string;
using std::vector;

bool is_valid_DNA_sequence(const string & input);
string reverseStr(string &str);
string rplcStr(string str);
void get_reverse_complement_sequence(const string & input,  string * const output);
string get_RNA_transcript(const string &input);
vector<vector<string>> get_reading_frames_as_codons(const string &input);
string translate(const vector<string> &codon_sequence);
string get_longest_open_reading_frame(const string &DNA_sequence);
string tranTU(string str);
vector<vector<string>> get_offset_sequence_in_vector(const string & input);
string translation(string input);

#endif
