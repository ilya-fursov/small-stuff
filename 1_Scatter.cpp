/*
 * Scatter.cpp
 *
 *  Created on: 30.01.2016
 *      Author: FursovIV
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

vector<string> STATES(10);
const int UNQ_OFFSET = 64;
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

template <typename T>
int TIndex(T val, vector<T> ARR);					// index of 'val' in 'ARR', -1 if not found
template <typename T>
int VecIndex(const vector<T> &vec, const vector<vector<T> > &DB);		// index of row 'vec' in matrix 'DB', -1 if not found
int StrIndex(const string &S, const vector<string> &VECS);	// index of S in VECS[], -1 if not found
															// more precisely, VECS[i] can be only a substring in S
															// search goes from beginning to end, the first match is returned
int StrIndexINV(const string &S, const vector<string> &VECS);	// as StrIndex, but search goes from end to beginning, the first match is returned
bool ReadTokenComm(FILE *F, char **str, bool &new_line, char *str0, const int str0_len);
															// reads a token from the file (delimited by ' ', '\t', '\r', '\n'), dropping "--..." comments
															// returns true on success, false on failure/EOF
															// the token is saved to str
															// set new_line = true in the first call, then the function will manage it
															// str0 is a working array, it should have been allocated
vector<string> TokenizeStr(char *str);				// tokenizes str according to delimiters ' ', '\t', '\r', '\n'
														// the input str will be corrupted!
vector<string> TokenizeStrComm(const char *str, int SZ);	// removes "--" comment and tokenizes as TokenizeStr; SZ = size of 'str'
vector<string> TokenizeStrConstr(char *str, int N);	// tokenizes as TokenizeStr, but not more than 'N' times;
														// the remainder of 'str' is put to (N+1)'th element: result[N]
														// the input str will be corrupted!
bool IsComm(const char *str);		// true if 'str' = spaces/tabs [--[something]]
void StripChars(char *str, const char *rm);			// removes everything from 'str' starting from the first occurrence of any symbols of 'rm'
char *StripChars2(char *str, const char *rm);		// works as StripChars, and also removes leading/trailing spaces/tabs,
													// the result is returned (it resides in the memory block of 'str'),
													// 'str' itself gets its trailing symbols removed
int GetState(const char *s);	// 0 - 's' contains nothing of below
								// 1 - #ENDHEADER found, 2 - #STARTFOOTER
								// 3 - BULLEYE, 4 - ENDBULLEYE
								// 5 - REDUCE, 6 - ENDREDUCE
								// 7 - RESTORE, 8 - ENDRESTORE
								// this function is case insensitive
								// if #3 - #8 is commented, it is treated as #0

string GatherFile(vector<string> fn_in, vector<FILE*> in, FILE *out, vector<int> di, vector<int> dj, vector<int> dk, const char *egrid, vector<vector<int> > &RESTORE_DB);
													// converts sector coords of the main section of 'in' to global coords of 'out'
													// egrid is supplied to 'out'
													// MASTER section is taken from in[0]
													// RESTORE_DB is a 'database' for all loaded restore lines, to handle duplication
													// returns a message

string ScatterFile(FILE *in, vector<FILE*> out, vector<int> regnum, vector<string> egrid, vector<int> di, vector<int> dj, vector<int> dk, vector<int> NI, vector<int> NJ, vector<int> NK);
											// converts global coords of the main section of 'in' to sector coords of 'out'
											// sector coords are checked for MIN and MAX (NI, NJ, NK)
											// regnum is the global number of each sector (for #UNIQUE)
											// egrid is supplied to 'out'
											// returns a message

void ParseCoordsBulleye(vector<string> VEC, int &i, int &j, int &k1, int &k2, string name);
void ParseCoordsRestore(vector<string> VEC, int &i1, int &j1, int &i2, int &j2, int &k1, int &k2);
void ParseCoordsRestore(vector<string> VEC, int &i1, int &j1, int &i2, int &j2, int &k1, int &k2, int &wid);
void ToUpper(string &s);
void PinToBoundary(int i1, int j1, int &i2, int &j2, int NI, int NJ);	// (i1, j1) is should be inside the boundary, (i2, j2) - outside
																		// (i2, j2) is changed to the boundary intersection preserving the line direction
int GetStrNum(string str0, const char *needle);	// if str0 == "*needle X", then X is returned
string RptRegs(vector<int> regnum, vector<bool> nonunq_reg);	// string of comma-separated 'regnum' for which 'nonunq_reg' == true
//---------------------------------------------------------------------------------------
void Gather(vector<string> args)
{
	size_t L = args.size();
	if (L < 2)
		throw string("need at least 2 lines with parameters");

	vector<vector<int> > RESTORE_DB;

	const int Buff = 1024;
	char msg[Buff];
	char fn_out[Buff];
	char egrid_out[Buff];
	char fn_in[Buff];

	if (sscanf(args[L-1].c_str(), "%s %s", fn_out, egrid_out) != 2)						// works fine with tabs
		throw string("output specification (file_name, egrid_name) is a mess");

	FILE *out = fopen(fn_out, "w");

	vector<FILE*> in(L-1);
	vector<string> fnames(L-1);
	vector<int> di(L-1);
	vector<int> dj(L-1);
	vector<int> dk(L-1);

	if (out == 0)
		throw string("cannot open output file ") + fn_out;

	try
	{
		for (size_t i = 0; i < L-1; i++)		// read each input file
		{
			if (sscanf(args[i].c_str(), "%s %d %d %d", fn_in, &di[i], &dj[i], &dk[i]) != 4)		// works fine with tabs
			{
				sprintf(msg, "mess in input specification # %d (file_name, delta_i, delta_j, delta_k)", int(i+1));
				throw string(msg);
			}

			fnames[i] = fn_in;
			in[i] = fopen(fn_in, "r");
			if (in[i] == 0)
				throw string("cannot open input file ") + fn_in;

			cout << "Reading " << fn_in << "\n";
		}

		cout << GatherFile(fnames, in, out, di, dj, dk, egrid_out, RESTORE_DB) << "\n";

		for (size_t i = 0; i < L-1; i++)
		{
			fclose(in[i]);
			in[i] = 0;
		}
	}
	catch (...)
	{
		for (size_t i = 0; i < L-1; i++)
			fclose(in[i]);
		fclose(out);
		throw;
	}

	fclose(out);
}
//---------------------------------------------------------------------------------------
void Scatter(vector<string> args)
{
	size_t L = args.size();
	if (L < 2)
		throw string("need at least 2 lines with parameters");

	const int Buff = 1024;
	char msg[Buff];
	char fn_out[Buff];
	char egrid_out[Buff];
	char fn_in[Buff];

	if (sscanf(args[0].c_str(), "%s", fn_in) != 1)
		throw string("input specification (file_name) is a mess");

	FILE *in = fopen(fn_in, "r");
	if (in == 0)
		throw string("cannot open input file ") + fn_in;

	vector<FILE*> out(L-1);
	vector<int> regnum(L-1);
	vector<string> egrid(L-1);
	vector<int> di(L-1);
	vector<int> dj(L-1);
	vector<int> dk(L-1);
	vector<int> NI(L-1);
	vector<int> NJ(L-1);
	vector<int> NK(L-1);
	try
	{
		for (size_t i = 0; i < L-1; i++)		// open each output file
		{
			if (sscanf(args[i+1].c_str(), "%s %d %s %d %d %d %d %d %d", fn_out, &regnum[i], egrid_out, &di[i], &dj[i], &dk[i], &NI[i], &NJ[i], &NK[i]) != 9)
			{
				sprintf(msg, "mess in output specification # %d (file_name, reg_num, egrid_name, delta_i, delta_j, delta_k, NI, NJ, NK)", int(i+1));
				throw string(msg);
			}

		    out[i] = fopen(fn_out, "w");
		    egrid[i] = egrid_out;

			if (out[i] == 0)
				throw string("cannot open output file ") + fn_out;
		}

		cout << "Reading " << fn_in << "\n";
		cout << ScatterFile(in, out, regnum, egrid, di, dj, dk, NI, NJ, NK);

		for (size_t i = 0; i < L-1; i++)
		{
			fclose(out[i]);
			out[i] = 0;
		}
	}
	catch (...)
	{
		fclose(in);
		for (size_t i = 0; i < L-1; i++)
			fclose(out[i]);
		throw;
	}

	fclose(in);
}
//---------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	STATES[0] = "#MASTER";
	STATES[1] = "#ENDMASTER";
	STATES[2] = "BULLEYE";
	STATES[3] = "ENDBULLEYE";
	STATES[4] = "REDUCE";
	STATES[5] = "ENDREDUCE";
	STATES[6] = "RESTORE";
	STATES[7] = "ENDRESTORE";
	STATES[8] = "#UNIQUE";
	STATES[9] = "#NONUNIQUE";

	if (argc != 2)
	{
		cerr << "Expected 1 argument: control file name\n";
		return 1;
	}

	FILE *FileContr = fopen(argv[1], "r");	// control file
	const int Buff = 4096;
	char str0[Buff];

	vector<string> CMDS(2);
	CMDS[0] = "GATHER";
	CMDS[1] = "SCATTER";

	vector<string> ARGS;	// arguments read from the control file

	int cmd_ind = -1;		// -1 - search mode, 0 - gather, 1 - scatter

	if (FileContr != 0)
	{
		try
		{
			// read the control file
			while (!feof(FileContr))
			{
				if (fgets(str0, Buff, FileContr) != 0)
				{
					char *str1 = StripChars2(str0, "\r\n");

					if (str1[0] != 0)	// non-empty string
					{
						if (cmd_ind == -1)	// looking for the command
						{
							string str2 = str1;
							ToUpper(str2);
							cmd_ind = StrIndex(str2, CMDS);
							if (cmd_ind == -1)
								throw string("a command should be specified: SCATTER or GATHER");
						}
						else				// main mode: reading settings
							ARGS.push_back(str1);
					}
				}
				else
					break;
			}

			if (cmd_ind == -1)
				throw string("a command should be specified: SCATTER or GATHER");

			// perform operations
			if (cmd_ind == 0)
				Gather(ARGS);
			if (cmd_ind == 1)
				Scatter(ARGS);

			fclose(FileContr);
		}
		catch (const string &e)
		{
			fclose(FileContr);
			cerr << "ERROR: " << e.c_str() << "\n";
		}
	}
	else
	{
		cerr << "Cannot open control file " << argv[1] << "\n";
		return 1;
	}

	cout << "Finished\n";
	return 0;
}
//---------------------------------------------------------------------------------------
template <typename T>
int TIndex(T val, vector<T> ARR)
{
	for (size_t i = 0; i < ARR.size(); i++)
		if (ARR[i] == val)
			return i;

	return -1;
}
//---------------------------------------------------------------------------------------
template <typename T>
int VecIndex(const vector<T> &vec, const vector<vector<T> > &DB)
{
	size_t SZ = vec.size();
	for (size_t i = 0; i < DB.size(); i++)
	{
		bool equal = true;
		if (SZ != DB[i].size())
			throw string("vector size mismatch in VecIndex");
		for (size_t j = 0; j < SZ; j++)
			if (vec[j] != DB[i][j])
			{
				equal = false;
				break;
			}
		if (equal)
			return i;
	}

	return -1;
}
//---------------------------------------------------------------------------------------
int StrIndex(const string &S, const vector<string> &VECS)
{
	for (size_t i = 0; i < VECS.size(); i++)
		if (S.find(VECS[i]) != string::npos)
			return i;

	return -1;
}
//---------------------------------------------------------------------------------------
int StrIndexINV(const string &S, const vector<string> &VECS)
{
	for (int i = VECS.size()-1; i >= 0 ; i--)
		if (S.find(VECS[i]) != string::npos)
			return i;

	return -1;
}
//---------------------------------------------------------------------------------------
bool ReadTokenComm(FILE *F, char **str, bool &new_line, char *str0, const int str0_len)
{
	static const char COMM[] = "--";			// comment beginning
	static const char DELIM[] = " \t\r\n";		// delimiters

	*str = 0;

	while (*str == 0)
	{
		if (new_line)
		{
			if (fgets(str0, str0_len, F) != 0)	// read the line
			{
				// remove the comment
				char *comm_ind = strstr(str0, COMM);
				if (comm_ind != 0)		// comment found
					comm_ind[0] = 0;	// set end-of-line at the comment start

				new_line = false;

				// get the first token
				*str = strtok(str0, DELIM);
			}
			else
				return false;
		}
		else
			*str = strtok(0, DELIM);

		if (*str == 0)
			new_line = true;
	}

	return true;
}
//---------------------------------------------------------------------------------------
vector<string> TokenizeStr(char *str)
{
	static const char DELIM[] = " \t\r\n";		// delimiters
	vector<string> res;

	char *str0 = strtok(str, DELIM);
	while (str0 != 0)
	{
		res.push_back(str0);
		str0 = strtok(0, DELIM);
	}

	return res;
}
//---------------------------------------------------------------------------------------
vector<string> TokenizeStrComm(const char *str, int SZ)
{
	static const char COMM[] = "--";			// comment beginning

	char *str1 = new char[SZ];
	memcpy(str1, str, SZ);

	// remove the comment
	char *comm_ind = strstr(str1, COMM);
	if (comm_ind != 0)		// comment found
		comm_ind[0] = 0;	// set end-of-line at the comment start

	vector<string> res = TokenizeStr(str1);
	delete [] str1;

	return res;
}
//---------------------------------------------------------------------------------------
vector<string> TokenizeStrConstr(char *str, int N)
{
	static const char DELIM[] = " \t\r\n";		// delimiters
	vector<string> res;
	char *str_end = str + strlen(str);			// points to the end of 'str'

	if (N > 0)
	{
		char *str0 = strtok(str, DELIM);
		if (str0 != 0)
			res.push_back(str0);
		int count = 1;

		while (str0 != 0 && count < N)
		{
			str0 = strtok(0, DELIM);
			if (str0 != 0)
				res.push_back(str0);
			count++;
		}

		if (str0 != 0 && (str0 + strlen(str0) + 1 < str_end))
			res.push_back(str0 + strlen(str0) + 1);		// remainder
	}
	else
		res.push_back(str);

	return res;
}
//---------------------------------------------------------------------------------------
bool IsComm(const char *str)
{
	static const char COMM[] = "--";

	char *str1 = const_cast<char*>(str) + strspn(str, " \t");	// leading space/tab removal
	if (*str1 == 0)
		return true;

	char *comm_ind = strstr(str1, COMM);	// find "--"

	return comm_ind == str1;
}
//---------------------------------------------------------------------------------------
void StripChars(char *str, const char *rm)
{
	char *aux = strpbrk(str, rm);
	if (aux != 0)
		*aux = 0;
}
//---------------------------------------------------------------------------------------
char *StripChars2(char *str, const char *rm)
{
	StripChars(str, rm);

	char *res = str + strspn(str, " \t");		// leading space/tab removal

	for (int i = strlen(res)-1; i >= 0; i--)	// trailing space/tab removal
		if (res[i] == ' ' || res[i] == '\t')
			res[i] = 0;
		else
			break;

	return res;
}
//---------------------------------------------------------------------------------------
int GetState(const char *s)
{
	string S = s;
	ToUpper(S);

	int res = StrIndexINV(S, STATES) + 1;
	if (res >= 3 && res <= 8 && IsComm(s))
		res = 0;

	return res;
}
//---------------------------------------------------------------------------------------
string GatherFile(vector<string> fn_in, vector<FILE*> in, FILE *out, vector<int> di, vector<int> dj, vector<int> dk, const char *egrid, vector<vector<int> > &RESTORE_DB)
{
	const int Buff = 2048;
	char str0[Buff];
	char str_print[Buff];

	string MSG = "";
	int count_found_restores = 0;

// коды new_st
//	[1] = "#MASTER";
//	[2] = "#ENDMASTER";
//	[3] = "BULLEYE";
//	[4] = "ENDBULLEYE";
//	[5] = "REDUCE";
//	[6] = "ENDREDUCE";
//	[7] = "RESTORE";
//	[8] = "ENDRESTORE";
//	[9] = "#UNIQUE";
//	[10] = "#NONUNIQUE";

	string kwd_buff = "";			// RESTORE first writes to this buffer
	int kwd_buff_size = 0;			// number of lines in kwd_buff

	string main_buff = "";			// main buffer, which is flushed on #MASTER, #ENDMASTER, eof

	size_t L = in.size();
	size_t contr = 0;		// file currently being read
	bool FEOF = false;		// 'true' when all files reach their eof

	// current state we are in
	vector<int> state(L);	// 0 - master, 1 - main part (neutral), 2 - bulleye, 3 - reduce, 4 - restore [5 - in the footer : obsolete]
	vector<string> make_title(L);
	for (size_t i = 0; i < L; i++)
	{
		state[i] = 1;
		make_title[i] = "";
	}

	while (!FEOF)
	{
		bool pass_contr = false;	// when 'true', control is passed to the next file
		while (!feof(in[contr]) && !pass_contr)
		{
			if (fgets(str0, Buff, in[contr]) == 0)			// read a line from input
				break;

			StripChars(str0, "\r\n");
			int new_st = GetState(str0);	// info on the new state from the new line (0..10)
			string str1;					// will be written to output

			if (new_st == 10)	// #NONUNIQUE
				throw fn_in[contr] + string(": found a record which is potentially non-unique across the sectors, GATHER cannot continue.\n"
							 "Leave the record only in one segment [seg-num] with mark #UNIQUE [seg-num]\n"
							 "Details:\n") + str0;
			switch (state[contr])
			{
				//-------------------
				case 0:		// master
					if (new_st != 2 && new_st != 0)
						throw string("unexpected keyword ") + STATES[new_st-1] + " in the MASTER section (missing #ENDMASTER ?) in " + fn_in[contr] ;

					if (contr == 0)
					{
						str1 = string(str0) + "\r\n";
						string STR1 = str1;
						ToUpper(STR1);
						if (STR1.find("INPUTGRIDBIN") != string::npos)
							str1 = string("INPUTGRIDBIN  ") + egrid + ".EGRID\r\n";

						str1 = make_title[contr] + str1;
						make_title[contr] = "";
					}
					else
						str1 = "";

					if (new_st == 2)
					{
						state[contr] = 1;
						pass_contr = true;

						sprintf(str_print, "\r\n------------------------\r\n-- %s\r\n------------------------\r\n\r\n", fn_in[contr].c_str());
						make_title[contr] = str_print;
					}

					break;
				//-------------------
				case 1:		// main, neutral
					str1 = make_title[contr];
					make_title[contr] = "";
					if (new_st == 0)
					{
						vector<string> tok = TokenizeStrComm(str0, Buff);
						if (tok.size() > 0)
							MSG += string("warning (") + fn_in[contr] + ("): processing not possible (keyword will be copied unchanged): ") + tok[0] + "\n";
					}
					if (new_st == 1)
					{
						pass_contr = true;
						state[contr] = 0;
					}
					if (new_st == 3)
						state[contr] = 2;
					if (new_st == 5)
						state[contr] = 3;
					if (new_st == 7)
						state[contr] = 4;
					if (new_st == 2 || new_st == 4 || new_st == 6 || new_st == 8)
						throw string("unpaired keyword ") + STATES[new_st-1] + " (" + fn_in[contr] + ")";
					if (new_st >= 9)
						throw string("unexpected keyword ") + STATES[new_st-1] + " in the main section (" + fn_in[contr] + ")";

					if (new_st == 7)						// start filling the buffer 'kwd_buff'
					{
						kwd_buff = string(str0) + "\r\n";
						kwd_buff_size = 1;
					}
					else if (new_st != 1)
						str1 += string(str0) + "\r\n";
					else if (contr == 0)
						make_title[contr] = string(str0) + "\r\n";

					break;
				//-------------------
				case 2:		// bulleye
					if (new_st == 4)
					{
						state[contr] = 1;
						str1 = string(str0) + "\r\n";
					}
					else if (new_st > 0 && new_st <= 8)
						throw string("unexpected keyword ") + STATES[new_st-1] + " in BULLEYE section (" + fn_in[contr] + ")";
					else
					{		// converting coordinates
						if (!IsComm(str0))
						{
							vector<string> VEC = TokenizeStrConstr(str0, 5);
							int I, J, K1, K2;
							ParseCoordsBulleye(VEC, I, J, K1, K2, "BULLEYE");
							I += di[contr];
							J += dj[contr];
							K1 += dk[contr];
							K2 += dk[contr];
							sprintf(str_print, "%-6s\t%d\t%d\t%d\t%d\t%s\r\n", VEC[0].c_str(), I, J, K1, K2, VEC[5].c_str());
							str1 = str_print;
						}
						else
							str1 = string(str0) + "\r\n";
					}
					break;
				//-------------------
				case 3:		// reduce
					if (new_st == 6)
					{
						state[contr] = 1;
						str1 = string(str0) + "\r\n";
					}
					else if (new_st > 0 && new_st <= 8)
						throw string("unexpected keyword ") + STATES[new_st-1] + " in REDUCE section (" + fn_in[contr] + ")";
					else
					{		// converting coordinates
						if (!IsComm(str0))
						{
							vector<string> VEC = TokenizeStrConstr(str0, 5);
							int I, J, K1, K2;
							ParseCoordsBulleye(VEC, I, J, K1, K2, "REDUCE");	// same pattern for coordinates as in BULLEYE
							I += di[contr];
							J += dj[contr];
							K1 += dk[contr];
							K2 += dk[contr];
							sprintf(str_print, "%-6s\t%d\t%d\t%d\t%d\t%s\r\n", VEC[0].c_str(), I, J, K1, K2, VEC[5].c_str());
							str1 = str_print;
						}
						else
							str1 = string(str0) + "\r\n";
					}
					break;
				//-------------------
				case 4:		// restore
					if (new_st == 8)
					{
						state[contr] = 1;
						// flush the buffer
						if (kwd_buff_size > 1)
							str1 = kwd_buff + string(str0) + "\r\n";

						kwd_buff = "";
						kwd_buff_size = 0;
					}
					else if (new_st > 0)
						throw string("unexpected keyword ") + STATES[new_st-1] + " in RESTORE section (" + fn_in[contr] + ")";
					else
					{		// converting coordinates
						if (!IsComm(str0))
						{
							vector<string> VEC = TokenizeStrConstr(str0, 7);
							int I1, J1, I2, J2, K1, K2, WID;
							ParseCoordsRestore(VEC, I1, J1, I2, J2, K1, K2, WID);
							I1 += di[contr];
							J1 += dj[contr];
							I2 += di[contr];
							J2 += dj[contr];
							K1 += dk[contr];
							K2 += dk[contr];

							vector<int> VEC_RESTORE(7);
							VEC_RESTORE[0] = I1;
							VEC_RESTORE[1] = J1;
							VEC_RESTORE[2] = I2;
							VEC_RESTORE[3] = J2;
							VEC_RESTORE[4] = K1;
							VEC_RESTORE[5] = K2;
							VEC_RESTORE[6] = WID;
							int ind = VecIndex(VEC_RESTORE, RESTORE_DB);
							if (ind == -1)	// not found
							{
								RESTORE_DB.push_back(VEC_RESTORE);
								if (VEC.size() > 7)
									sprintf(str_print, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\r\n", I1, J1, I2, J2, K1, K2, WID, VEC[7].c_str());
								else
									sprintf(str_print, "%d\t%d\t%d\t%d\t%d\t%d\t%d\r\n", I1, J1, I2, J2, K1, K2, WID);

								kwd_buff += str_print;
								kwd_buff_size += 1;
							}
							else
								count_found_restores++;
						}
						else
							str1 = string(str0) + "\r\n";
					}
					break;
				//-------------------
				case 5:		// former footer
					throw string("found illegal state 5");
			}

			main_buff += str1;		// write to buffer
		}

		contr++;	// pass control to the next file
		if (contr >= L)
		{
			contr = 0;

			// check sync & dump
			FEOF = true;
			for (size_t i = 0; i < L; i++)
			{
				if (state[i] != state[0])
					throw string("#MASTER and #ENDMASTER are not synchronized between ") + fn_in[0] + " and " + fn_in[i];
				if (!feof(in[i]))
					FEOF = false;
			}

			fprintf(out, "%s", main_buff.c_str());
			main_buff = "";
		}
	}

	if (count_found_restores > 0)
	{
		sprintf(str_print, "Information: %d duplicate RESTORE lines were suppressed\n", count_found_restores);
		MSG += str_print;
	}

	return MSG;
}
//---------------------------------------------------------------------------------------
string ScatterFile(FILE *in, vector<FILE*> out, vector<int> regnum, vector<string> egrid, vector<int> di, vector<int> dj, vector<int> dk, vector<int> NI, vector<int> NJ, vector<int> NK)
{
	const int Buff = 2048;
	char str0[Buff];
	char str_print[Buff];

	string MSG = "";

	// коды new_st
	//	[1] = "#MASTER";
	//	[2] = "#ENDMASTER";
	//	[3] = "BULLEYE";
	//	[4] = "ENDBULLEYE";
	//	[5] = "REDUCE";
	//	[6] = "ENDREDUCE";
	//	[7] = "RESTORE";
	//	[8] = "ENDRESTORE";
	//	[9] = "#UNIQUE";
	//	[10] = "#NONUNIQUE";

	size_t L = out.size();
	vector<string> str1(L);					// str1[i] will be written to out[i]
	vector<string> kwd_buff(L);				// BULLEYE, REDUCE, RESTORE first write to this buffer
	vector<int> kwd_buff_size(L);			// number of lines in kwd_buff

	vector<bool> nonunq_regnum(L);			// regions for which non-unique records exist
	bool non_unique_found = false;

	for (size_t i = 0; i < L; i++)
	{
		str1[i] = "";
		kwd_buff[i] = "";
		kwd_buff_size[i] = 0;
		nonunq_regnum[i] = false;
	}

	// current state we are in
	int state = 1;	// 0 - master, 1 - main part (neutral), 2 - bulleye, 3 - reduce, 4 - restore [5 - in the footer: obsolete]
	while (!feof(in))
	{
		if (fgets(str0, Buff, in) == 0)			// read a line from input
			break;

		StripChars(str0, "\r\n");
		string str0_copy = str0;
		int new_st = GetState(str0);	// info on the new state from the new line (0..10)

		switch (state)
		{
			//-------------------
			case 0:		// master
				if (new_st != 2 && new_st != 0)
					throw string("unexpected keyword ") + STATES[new_st-1] + " in the MASTER section (missing #ENDMASTER ?)";
				if (new_st == 2)
					state = 1;

				for (size_t i = 0; i < L; i++)
				{
					string str_aux = string(str0);
					str1[i] = str_aux + "\r\n";
					ToUpper(str_aux);
					if (str_aux.find("INPUTGRIDBIN") != string::npos)
						str1[i] = string("INPUTGRIDBIN  ") + egrid[i] + ".EGRID\r\n";
					if (new_st == 2)
						str1[i] += "\r\n";
				}
				break;
			//-------------------
			case 1:		// main, neutral
				if (new_st == 0)
				{
					vector<string> tok = TokenizeStrComm(str0, Buff);
					if (tok.size() > 0)
						MSG += string("warning: processing not possible (keyword will be copied unchanged): ") + tok[0] + "\n";
				}
				if (new_st == 1)
					state = 0;
				if (new_st == 3)
					state = 2;
				if (new_st == 5)
					state = 3;
				if (new_st == 7)
					state = 4;
				if (new_st == 2 || new_st == 4 || new_st == 6 || new_st == 8)
					throw string("unpaired keyword ") + STATES[new_st-1];
				if (new_st >= 9)
					throw string("unexpected keyword ") + STATES[new_st-1] + " in the main section";

				for (size_t i = 0; i < L; i++)
					if (new_st == 0 && !IsComm(str0))
						str1[i] = string(str0) + "\r\n\r\n";		// only unrecognized keywords which are not comments, are taken
					else if (new_st == 1)
						str1[i] = string(str0) + "\r\n";			// master
					else
						str1[i] = "";

				if (new_st == 3 || new_st == 5 || new_st == 7)	// for main keywords - start filling the buffer
					for (size_t i = 0; i < L; i++)
					{
						kwd_buff[i] = string(str0) + "\r\n";
						kwd_buff_size[i] = 1;
					}

				break;
			//-------------------
			case 2:		// bulleye
				if (new_st == 4)
				{
					state = 1;
					// flush the buffer
					for (size_t i = 0; i < L; i++)
					{
						if (kwd_buff_size[i] > 1)
							str1[i] = kwd_buff[i] + string(str0) + "\r\n\r\n";

						kwd_buff[i] = "";
						kwd_buff_size[i] = 0;
					}
				}
				else if (new_st > 0 && new_st <= 8)
					throw string("unexpected keyword ") + STATES[new_st-1] + " in BULLEYE section";
				else
				{		// converting coordinates
					if (!IsComm(str0))
					{
						int reg = -1;		// index in 'regnum' for selective writing only to 'reg'-segment
						int glob_reg;
						if (new_st == 9)
						{
							glob_reg = GetStrNum(str0, "#UNIQUE");
							reg = TIndex(glob_reg, regnum);
							if (reg == -1)
							{
								sprintf(str_print, "warning (BULLEYE): sector specified in #UNIQUE (%d) was not found, the line will be skipped\n", glob_reg);
								MSG += (string)str_print + "Details:\n" + str0_copy + "\n";
							}
						}

						vector<string> VEC = TokenizeStrConstr(str0, 5);
						int I, J, K1, K2;
						ParseCoordsBulleye(VEC, I, J, K1, K2, "BULLEYE");

						int count_regs = 0;		// count segments which take the record
						for (size_t i = 0; i < L; i++)	// 1 - count segments
						{
							int Isec = I - di[i];
							int Jsec = J - dj[i];
							int K1sec = K1 - dk[i];
							int K2sec = K2 - dk[i];

							if (Isec >= 1     && Jsec >= 1     && K1sec >= 1     && K2sec >= 1 &&
								Isec <= NI[i] && Jsec <= NJ[i] && K1sec <= NK[i] && K2sec <= NK[i] &&
								(new_st != 9 || (new_st == 9 && (int)i == reg)))
							{
								count_regs++;
							}
						}
						if (new_st == 9 && count_regs == 0 && reg != -1)
						{
							sprintf(str_print, "warning (BULLEYE): sector specified in #UNIQUE (%d) is inconsistent with the record's coordinates; "
											   "make corrections or the record will be skipped during SCATTER\n", glob_reg);
							MSG += (string)str_print + "Details:\n" + str0_copy + "\n";
						}

						for (size_t i = 0; i < L; i++)	// 2 - take record
						{
							int Isec = I - di[i];
							int Jsec = J - dj[i];
							int K1sec = K1 - dk[i];
							int K2sec = K2 - dk[i];

							if (Isec >= 1     && Jsec >= 1     && K1sec >= 1     && K2sec >= 1 &&
								Isec <= NI[i] && Jsec <= NJ[i] && K1sec <= NK[i] && K2sec <= NK[i] &&
								(new_st != 9 || (new_st == 9 && (int)i == reg)))
							{
								sprintf(str_print, "%-6s\t%d\t%d\t%d\t%d\t%s", VEC[0].c_str(), Isec, Jsec, K1sec, K2sec, VEC[5].c_str());
								string str_pad = str_print;
								if (count_regs > 1)
								{
									int pad_len = UNQ_OFFSET - 1 - strlen(str_print);
									if (pad_len > 0)
										str_pad += string(pad_len, ' ');
									str_pad += "  -- #NONUNIQUE";
									non_unique_found = true;
									nonunq_regnum[i] = true;
								}

								kwd_buff[i] += str_pad + "\r\n";
								kwd_buff_size[i] += 1;
							}
						}
					}
				}
				break;
			//-------------------
			case 3:		// reduce
				if (new_st == 6)
				{
					state = 1;
					// flush the buffer
					for (size_t i = 0; i < L; i++)
					{
						if (kwd_buff_size[i] > 1)
							str1[i] = kwd_buff[i] + string(str0) + "\r\n\r\n";

						kwd_buff[i] = "";
						kwd_buff_size[i] = 0;
					}
				}
				else if (new_st > 0 && new_st <= 8)
					throw string("unexpected keyword ") + STATES[new_st-1] + " in REDUCE section";
				else
				{		// converting coordinates
					if (!IsComm(str0))
					{
						int reg = -1;		// index in 'regnum' for selective writing only to 'reg'-segment
						int glob_reg;
						if (new_st == 9)
						{
							glob_reg = GetStrNum(str0, "#UNIQUE");
							reg = TIndex(glob_reg, regnum);
							if (reg == -1)
							{
								sprintf(str_print, "warning (REDUCE): sector specified in #UNIQUE (%d) was not found, the line will be skipped\n", glob_reg);
								MSG += (string)str_print + "Details:\n" + str0_copy + "\n";
							}
						}

						vector<string> VEC = TokenizeStrConstr(str0, 5);
						int I, J, K1, K2;
						ParseCoordsBulleye(VEC, I, J, K1, K2, "REDUCE");	// same pattern for coordinates as in BULLEYE

						int count_regs = 0;		// count segments which take the record
						for (size_t i = 0; i < L; i++)	// 1 - count segments
						{
							int Isec = I - di[i];
							int Jsec = J - dj[i];
							int K1sec = K1 - dk[i];
							int K2sec = K2 - dk[i];

							if (Isec >= 1     && Jsec >= 1     && K1sec >= 1     && K2sec >= 1 &&
								Isec <= NI[i] && Jsec <= NJ[i] && K1sec <= NK[i] && K2sec <= NK[i] &&
								(new_st != 9 || (new_st == 9 && (int)i == reg)))
							{
								count_regs++;
							}
						}
						if (new_st == 9 && count_regs == 0 && reg != -1)
						{
							sprintf(str_print, "warning (REDUCE): sector specified in #UNIQUE (%d) is inconsistent with the record's coordinates; "
											   "make corrections or the record will be skipped during SCATTER\n", glob_reg);
							MSG += (string)str_print + "Details:\n" + str0_copy + "\n";
						}

						for (size_t i = 0; i < L; i++)	// 2 - take record
						{
							int Isec = I - di[i];
							int Jsec = J - dj[i];
							int K1sec = K1 - dk[i];
							int K2sec = K2 - dk[i];

							if (Isec >= 1     && Jsec >= 1     && K1sec >= 1     && K2sec >= 1 &&
								Isec <= NI[i] && Jsec <= NJ[i] && K1sec <= NK[i] && K2sec <= NK[i] &&
								(new_st != 9 || (new_st == 9 && (int)i == reg)))
							{
								sprintf(str_print, "%-6s\t%d\t%d\t%d\t%d\t%s", VEC[0].c_str(), Isec, Jsec, K1sec, K2sec, VEC[5].c_str());
								string str_pad = str_print;
								if (count_regs > 1)
								{
									int pad_len = UNQ_OFFSET - 1 - strlen(str_print);
									if (pad_len > 0)
										str_pad += string(pad_len, ' ');
									str_pad += "  -- #NONUNIQUE";
									non_unique_found = true;
									nonunq_regnum[i] = true;
								}

								kwd_buff[i] += str_pad + "\r\n";
								kwd_buff_size[i] += 1;
							}
						}
					}
				}
				break;
			//-------------------
			case 4:		// restore
				if (new_st == 8)
				{
					state = 1;
					// flush the buffer
					for (size_t i = 0; i < L; i++)
					{
						if (kwd_buff_size[i] > 1)
							str1[i] = kwd_buff[i] + string(str0) + "\r\n\r\n";

						kwd_buff[i] = "";
						kwd_buff_size[i] = 0;
					}
				}
				else if (new_st > 0)
					throw string("unexpected keyword ") + STATES[new_st-1] + " in RESTORE section";
				else
				{		// converting coordinates
					if (!IsComm(str0))
					{
						vector<string> VEC = TokenizeStrConstr(str0, 6);
						int I1, J1, I2, J2, K1, K2;
						ParseCoordsRestore(VEC, I1, J1, I2, J2, K1, K2);
						for (size_t i = 0; i < L; i++)
						{
							int I1sec = I1 - di[i];
							int J1sec = J1 - dj[i];
							int I2sec = I2 - di[i];
							int J2sec = J2 - dj[i];
							int K1sec = K1 - dk[i];
							int K2sec = K2 - dk[i];
							bool ins1 = I1sec >= 1 && I1sec <= NI[i] && J1sec >= 1 && J1sec <= NJ[i];
							bool ins2 = I2sec >= 1 && I2sec <= NI[i] && J2sec >= 1 && J2sec <= NJ[i];
							if (K1sec >= 1 && K1sec <= NK[i] && K2sec >= 1 && K2sec <= NK[i] && (ins1 || ins2))
							{
								if (!ins1 || !ins2)		// one point is outside: find boundary point
								{
									if (ins1)
										PinToBoundary(I1sec, J1sec, I2sec, J2sec, NI[i], NJ[i]);
									else
										PinToBoundary(I2sec, J2sec, I1sec, J1sec, NI[i], NJ[i]);
								}

								sprintf(str_print, "%d\t%d\t%d\t%d\t%d\t%d\t%s\r\n", I1sec, J1sec, I2sec, J2sec, K1sec, K2sec, VEC[6].c_str());
								kwd_buff[i] += str_print;
								kwd_buff_size[i] += 1;
							}
						}
					}
				}
				break;
			//-------------------
			case 5:		// former footer
				throw string("found illegal state 5");
		}

		for (size_t i = 0; i < L; i++)
		{
			fprintf(out[i], "%s", str1[i].c_str());
			str1[i] = "";
		}
	}

	if (non_unique_found)
		MSG += string("warning: some BULLEYE/REDUCE records were placed in more than one segment.\n"
					 "These records were marked with #NONUNIQUE and should be corrected before any future GATHER can work.\n"
					 "To make correction, manually replace #NONUNIQUE with #UNIQUE [seg-num] to attribute the record\n"
					 "to segment [seg-num]. Currently segments with #NONUNIQUE are ") + RptRegs(regnum, nonunq_regnum) + "\n";
	return MSG;
}
//---------------------------------------------------------------------------------------
void ParseCoordsBulleye(vector<string> VEC, int &i, int &j, int &k1, int &k2, string name)
{
	if (VEC.size() != 6)
		throw string("invalid number of items in ") + name + " string";
	int count = 0;
	count += sscanf(VEC[1].c_str(), "%d", &i);
	count += sscanf(VEC[2].c_str(), "%d", &j);
	count += sscanf(VEC[3].c_str(), "%d", &k1);
	count += sscanf(VEC[4].c_str(), "%d", &k2);
	if (count < 4)
		throw string("cannot parse coordinates in ") + name + " string";
}
//---------------------------------------------------------------------------------------
void ParseCoordsRestore(vector<string> VEC, int &i1, int &j1, int &i2, int &j2, int &k1, int &k2)
{
	if (VEC.size() != 7)
		throw string("invalid number of items in RESTORE string");
	int count = 0;
	count += sscanf(VEC[0].c_str(), "%d", &i1);
	count += sscanf(VEC[1].c_str(), "%d", &j1);
	count += sscanf(VEC[2].c_str(), "%d", &i2);
	count += sscanf(VEC[3].c_str(), "%d", &j2);
	count += sscanf(VEC[4].c_str(), "%d", &k1);
	count += sscanf(VEC[5].c_str(), "%d", &k2);
	if (count < 6)
		throw string("cannot parse coordinates in RESTORE string");
}
//---------------------------------------------------------------------------------------
void ParseCoordsRestore(vector<string> VEC, int &i1, int &j1, int &i2, int &j2, int &k1, int &k2, int &wid)
{
	if (VEC.size() != 7 && VEC.size() != 8)
		throw string("invalid number of items in RESTORE string");
	int count = 0;
	count += sscanf(VEC[0].c_str(), "%d", &i1);
	count += sscanf(VEC[1].c_str(), "%d", &j1);
	count += sscanf(VEC[2].c_str(), "%d", &i2);
	count += sscanf(VEC[3].c_str(), "%d", &j2);
	count += sscanf(VEC[4].c_str(), "%d", &k1);
	count += sscanf(VEC[5].c_str(), "%d", &k2);
	count += sscanf(VEC[6].c_str(), "%d", &wid);
	if (count < 7)
		throw string("cannot parse coordinates in RESTORE string");
}
//---------------------------------------------------------------------------------------
void ToUpper(string &s)
{
	transform(s.begin(), s.end(), s.begin(), ::toupper);
}
//---------------------------------------------------------------------------------------
void PinToBoundary(int i1, int j1, int &i2, int &j2, int NI, int NJ)
{
	double t = numeric_limits<double>::max();
	double r;
	const double eps = 1e-5;

	r = double(NI - i1)/(i2 - i1);
	if (r > 0 && r < t)
		t = r;

	r = double(1 - i1)/(i2 - i1);
	if (r > 0 && r < t)
		t = r;

	r = double(NJ - j1)/(j2 - j1);
	if (r > 0 && r < t)
		t = r;

	r = double(1 - j1)/(j2 - j1);
	if (r > 0 && r < t)
		t = r;

	if (!(1-eps <= double(i1) + double(i2 - i1)*t && double(i1) + double(i2 - i1)*t <= NI+eps &&
		  1-eps <= double(j1) + double(j2 - j1)*t && double(j1) + double(j2 - j1)*t <= NJ+eps))
	{
		t = 0;
	}

	i2 = floor(double(i1) + double(i2 - i1)*t + 0.5);
	j2 = floor(double(j1) + double(j2 - j1)*t + 0.5);

	if (i2 < 1)
		i2 = 1;
	if (i2 > NI)
		i2 = NI;
	if (j2 < 1)
		j2 = 1;
	if (j2 > NJ)
		j2 = NJ;
}
//---------------------------------------------------------------------------------------
int GetStrNum(string str0, const char *needle)
{
	ToUpper(str0);
	const char *s0 = str0.c_str();
	s0 = strstr(s0, needle);

	if (s0 == 0)
		throw string("needle not found in GetStrNum");

	int res;
	string templ = string(needle) + " %d";
	if (sscanf(s0, templ.c_str(), &res) != 1)
		throw string("number not recognized after ") + needle;

	return res;
}
//---------------------------------------------------------------------------------------
string RptRegs(vector<int> regnum, vector<bool> nonunq_reg)
{
	const int Buff = 512;
	char str[Buff];

	if (regnum.size() != nonunq_reg.size())
		throw string("vector size mismatch in RptRegs");

	string res = "";
	for (size_t i = 0; i < regnum.size(); i++)
		if (nonunq_reg[i])
		{
			if (res.size() > 0)
				sprintf(str, ", %d", regnum[i]);
			else
				sprintf(str, "%d", regnum[i]);

			res += str;
		}

	return res;
}
//---------------------------------------------------------------------------------------
