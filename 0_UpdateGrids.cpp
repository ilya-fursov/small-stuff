#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>

#define IJK(ARR) ARR[(i-1) + (j-1)*NI + (k-1)*Nij]		// macro for grid cell accessing

using namespace std;
//---------------------------------------------------------------------------------------
// the main calculation function
//---------------------------------------------------------------------------------------
string Calculation(vector<double*> geo1, vector<double*> GEO2,
				   vector<double*> gdm1, vector<double*> GDM2, size_t NI, size_t NJ, size_t NK,
				   vector<string> geo_names, vector<string> gdm_names)
{	// geo1 -  old geomodel: NTG  PORO PERMEABILITY SWAT IRR.WATERSATURATION SOIL
	// GEO2 -  new geomodel: NTG  PORO PERMEABILITY SWAT IRR.WATERSATURATION SOIL
	// gdm1 -	old GDM: ACTNUM NTG PORO PERMX SWATINIT SWL SOWCR
	// GDM2 -	new GDM: ACTNUM NTG PORO PERMX SWATINIT SWL SOWCR (output)
	// NI, NJ, NK - model dimensions
	// geo_names - names of geo grids, gdm_names - names of gdm grids
	// return: "" on success, error message on failure

	const int Buff = 4096;
	char msg[Buff];
	size_t Nij = NI*NJ;

	size_t count_SWL_patch = 0;		// ad hoc patching
	double SWL_patch = 0;
	size_t count_SOWCR_patch = 0;
	double SOWCR_patch = 0;

	FILE *err_rpt = fopen("Zero_cell_report.txt", "w");		// file for reporting problem cells
	bool see_err_rpt = false;

	// GRID cell [i, j, k] is accessed by GRID[(i-1) + (j-1)*NI + (k-1)*Nij], or IJK(GRID)
	// where i = 1..NI, j = 1..NJ, k = 1..NK
	try
	{
		for (size_t k = 1; k <= NK; k++)	// loop order: k, j, i - for better locality
			for (size_t j = 1; j <= NJ; j++)
				for (size_t i = 1; i <= NI; i++)
				{
					// ACTNUM
					int act_geo1 = IJK(geo1[0]);
					int act_GEO2 = IJK(GEO2[0]);
					int act_gdm_1 = IJK(gdm1[0]);
					int ntg_gdm_1 = IJK(gdm1[1]);
					act_gdm_1 = act_gdm_1 && ntg_gdm_1;		// note: gdm1 cell is considered active if both ACTNUM and NTG are 1
					int ACT = act_GEO2 - act_geo1 + act_gdm_1;
					if (ACT < 0) ACT = 0;
					if (ACT > 1) ACT = 1;
					IJK(GDM2[0]) = ACT;

					// NTG
					IJK(GDM2[1]) = ACT;


					//*****************************************
					// some very ad-hoc patching: remove zero & negative values in gdm1
					// 0	  1	  2	   3	 4		  5	  6
					// ACTNUM NTG PORO PERMX SWATINIT SWL SOWCR
					if (act_gdm_1)
					{
						if (IJK(gdm1[5]) <= 0 && j > 196)	// SWL
						{
							IJK(gdm1[5]) = -6.887 * IJK(gdm1[2]) + 1.786;
							SWL_patch += IJK(gdm1[5]);
							count_SWL_patch++;
						}
						if (IJK(gdm1[6]) <= 0 && k <= 60)	// SOWCR BV5
						{
							IJK(gdm1[6]) = 2.8043 * IJK(gdm1[2]) - 0.3196;
							SOWCR_patch += IJK(gdm1[6]);
							count_SOWCR_patch++;
						}
						if (IJK(gdm1[6]) <= 0 && k > 60)	// SOWCR BV6
						{
							IJK(gdm1[6]) = 4.866 * IJK(gdm1[2]) - 0.73;
							SOWCR_patch += IJK(gdm1[6]);
							count_SOWCR_patch++;
						}
					}
					// end of the very ad-hoc section
					//*****************************************

					// continuous props: consider different combinations of active cells
					if (act_geo1 && act_GEO2 && act_gdm_1)			// all three are active
						for (int c = 1; c <= 5; c++)
						{
							if (IJK(geo1[c]) == 0)
							{
								sprintf(msg, "ERROR: Old-geomodel grid %s has 0 value in active cell %zd %zd %zd\n", geo_names[c].c_str(), i, j, k);
								throw msg;
							}
							if (IJK(gdm1[c+1]) == 0)
							{
								fprintf(err_rpt, "Old-GDM grid	%s	has 0 value in active cell	%zd	%zd	%zd\r\n", gdm_names[c+1].c_str(), i, j, k);
								see_err_rpt = true;
							}
							if (IJK(GEO2[c]) == 0)
							{
								fprintf(err_rpt, "New-geo grid	%s	has 0 value in active cell	%zd	%zd	%zd\r\n", geo_names[c].c_str(), i, j, k);
								see_err_rpt = true;
							}

							IJK(GDM2[c+1]) = IJK(gdm1[c+1]) * IJK(GEO2[c]) / IJK(geo1[c]);
						}
					else if (!act_geo1 && act_GEO2)					// geo1 (old) - inactive, geo2 (new) - active, gdm1 - any
						for (int c = 1; c <= 5; c++)
						{
							if (IJK(GEO2[c]) == 0)
							{
								fprintf(err_rpt, "New-geo grid	%s	has 0 value in active cell	%zd	%zd	%zd\r\n", geo_names[c].c_str(), i, j, k);
								see_err_rpt = true;
							}

							IJK(GDM2[c+1]) = IJK(GEO2[c]);
						}
					else if (!act_geo1 && !act_GEO2 && act_gdm_1)	// only gdm1 (old) is active
						for (int c = 1; c <= 5; c++)
						{
							if (IJK(gdm1[c+1]) == 0)
							{
								fprintf(err_rpt, "Old-GDM grid	%s	has 0 value in active cell	%zd	%zd	%zd\r\n", gdm_names[c+1].c_str(), i, j, k);
								see_err_rpt = true;
							}

							IJK(GDM2[c+1]) = IJK(gdm1[c+1]);
						}
					else											// all that remains
						for (int c = 1; c <= 5; c++)
							IJK(GDM2[c+1]) = 0;

					// some corrections
					if (ACT)	//ACTNUM NTG PORO PERMX SWATINIT SWL SOWCR
					{
						if (IJK(GDM2[4]) > 1)		// SWATINIT patch
							IJK(GDM2[4]) = 1;
						if (IJK(GDM2[6]) < 0.001)	// SOWCR patch
							IJK(GDM2[6]) = 0.001;
						if (IJK(GDM2[5]) + IJK(GDM2[6]) >= 0.999)		// SWL + SOWCR >= 0.999 --> patch SWL
							IJK(GDM2[5]) = 0.999 - IJK(GDM2[6]);
					}
				}
	}
	catch (const char *e)
	{
		fclose(err_rpt);
		return string(e);
	}

	fclose(err_rpt);
	if (see_err_rpt)
		cout << "WARNING! Problem cells found. See details in Zero_cell_report.txt\n";

	// ad hoc patching
	if (count_SWL_patch > 0)
		SWL_patch /= count_SWL_patch;
	if (count_SOWCR_patch > 0)
		SOWCR_patch /= count_SOWCR_patch;
	if (count_SWL_patch > 0)
		cout << "MESSAGE: in old simulation model corrected " << count_SWL_patch << " zero/negative SWL values (average corrected value " << SWL_patch << ")\n";
	if (count_SOWCR_patch > 0)
		cout << "MESSAGE: in old simulation model corrected " << count_SOWCR_patch << " zero/negative SOWCR values (average corrected value " << SOWCR_patch << ")\n";

	return "";
}
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

int StrIndex(const string &S, const vector<string> &VECS);	// index of S in VECS[], -1 if not found
															// more precisely, VECS[i] can be only a substring in S
string ReadGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2);	// reads a number of grids from file
															// data = vector<double[len]>, all memory should be allocated prior to calling,
															// S1[i], S2 - start and end of grid[i] which is loaded to data[i]
															// returns "" on success, error message on failure
void WriteGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, vector<string> pat, int slen);	// writes a number of grids to file
															// data = vector<double[len]>
															// S1[i], S2 - start and end of grid[i] which is taken from data[i]
															// pat[i] - pattern string for writing 'doubles', slen - target length of each line
void WriteGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, int slen);
															// same as above, but automatically sets 'pat'
void WriteGridsTest(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, int slen);	// version for performance benchmark
bool ReadTokenComm(FILE *F, char **str, bool &new_line, char *str0, const int str0_len);
															// reads a token from the file (delimited by ' ', '\t', '\r', '\n'), dropping "--..." comments
															// returns true on success, false on failure/EOF
															// the token is saved to str
															// set new_line = true in the first call, then the function will manage it
															// str0 is a working array, it should have been allocated
vector<string> TokenizeStr(char *str);				// tokenizes str according to delimiters ' ', '\t', '\r', '\n'
													// the input str will be corrupted!
template <typename T>
size_t DiffValInd(T *arr, size_t i0, size_t len);		// returns first index i1 after i0: arr[i1] != arr[i0]; if not found, returns 'len'
														// arr should have size 'len'
void InitGrid(vector<double*> &G, size_t LEN);			// each double* is set to new double[LEN]
void FreeGrid(vector<double*> &G);						// delete each double*
void StripChars(char *str, const char *rm);			// removes everything from 'str' starting from the first occurrence of any symbols of 'rm'
string RptGridRead(const char *fname, vector<string> grids);	// report message about the grids read
//---------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	if (argc != 2)
	{
		cerr << "Expected 1 argument: control file name\n";
		return 1;
	}

	FILE *FileContr = fopen(argv[1], "r");	// control file
	const int Buff = 1024;
	char str0[Buff];
	char file_geo1[Buff], file_geo2[Buff], file_GDM_gr1[Buff], file_GDM_rp1[Buff], file_GDM_GR2[Buff], file_GDM_RP2[Buff];

	if (FileContr != 0)
	{
		const int NGEO = 6, NIGRID = 4, NIRP = 3, NOUTGRID = 4, NOUTRP = 3;		// numbers of input and output grids
		int ni(0), nj(0), nk(0);	// model dimensions

		fgets(str0, Buff, FileContr);
		sscanf(str0, "DIMS %d %d %d", &ni, &nj, &nk);

		size_t Nijk = ni*nj*nk;

		fgets(str0, Buff, FileContr);
		fgets(str0, Buff, FileContr);		// reading INPUT
		fgets(file_geo1, Buff, FileContr);	// geomodel files: old and new
		fgets(file_geo2, Buff, FileContr);
		StripChars(file_geo1, "\r\n");
		StripChars(file_geo2, "\r\n");

		fgets(str0, Buff, FileContr);		// grids read from geomodels: NTG PORO PERMEABILITY SWAT IRR.WATERSATURATION SOIL
		vector<string> GridsGeo = TokenizeStr(str0);

		fgets(str0, Buff, FileContr);
		fgets(file_GDM_gr1, Buff, FileContr);
		StripChars(file_GDM_gr1, "\r\n");
		fgets(str0, Buff, FileContr);
		vector<string> GridsGDM_gr = TokenizeStr(str0);	// "grid" grids: ACTNUM NTG PORO PERMX

		fgets(str0, Buff, FileContr);
		fgets(file_GDM_rp1, Buff, FileContr);
		StripChars(file_GDM_rp1, "\r\n");
		fgets(str0, Buff, FileContr);
		vector<string> GridsGDM_rp = TokenizeStr(str0);	// "rp" grids: SWATINIT SWL SOWCR

		fgets(str0, Buff, FileContr);
		fgets(str0, Buff, FileContr);		// reading OUTPUT
		fgets(file_GDM_GR2, Buff, FileContr);
		StripChars(file_GDM_GR2, "\r\n");
		fgets(str0, Buff, FileContr);
		vector<string> GridsGDM_gr_out = TokenizeStr(str0);	// "grid" grids: ACTNUM NTG PORO PERMX

		fgets(str0, Buff, FileContr);
		fgets(file_GDM_RP2, Buff, FileContr);
		StripChars(file_GDM_RP2, "\r\n");
		fgets(str0, Buff, FileContr);
		vector<string> GridsGDM_rp_out = TokenizeStr(str0);	// "rp" grids: SWATINIT SWL SOWCR

		if (GridsGeo.size() != (size_t)NGEO || GridsGDM_gr.size() != (size_t)NIGRID || GridsGDM_rp.size() != (size_t)NIRP ||
			GridsGDM_gr_out.size() != (size_t)NOUTGRID || GridsGDM_rp_out.size() != (size_t)NOUTRP)
		{
			fclose(FileContr);
			cerr << "Incorrect number of grids defined in the control file\n";
			return 1;
		}

		vector<double*> Geo1(NGEO);		// declare arrays of grids
		vector<double*> Geo2(NGEO);
		vector<double*> GDMgrid(NIGRID);
		vector<double*> GDMrp(NIRP);
		vector<double*> GDMgridOut(NOUTGRID);
		vector<double*> GDMrpOut(NOUTRP);

		InitGrid(Geo1, Nijk);			// allocate memory
		InitGrid(Geo2, Nijk);
		InitGrid(GDMgrid, Nijk);
		InitGrid(GDMrp, Nijk);
		InitGrid(GDMgridOut, Nijk);
		InitGrid(GDMrpOut, Nijk);

		string msg = "";
		msg += ReadGrids(file_geo1, Nijk, Geo1, GridsGeo, "/");		// read grids from files
		if (msg == "")
			cout << RptGridRead(file_geo1, GridsGeo) << "\n";

		msg += ReadGrids(file_geo2, Nijk, Geo2, GridsGeo, "/");
		if (msg == "")
			cout << RptGridRead(file_geo2, GridsGeo) << "\n";

		msg += ReadGrids(file_GDM_gr1, Nijk, GDMgrid, GridsGDM_gr, "/");
		if (msg == "")
			cout << RptGridRead(file_GDM_gr1, GridsGDM_gr) << "\n";

		msg += ReadGrids(file_GDM_rp1, Nijk, GDMrp, GridsGDM_rp, "/");
		if (msg == "")
			cout << RptGridRead(file_GDM_rp1, GridsGDM_rp) << "\n";

		if (msg != "")
		{
			fclose(FileContr);
			cerr << msg.c_str();
			return 1;
		}

		cout << "Calculating...\n";

		// concat the grids from GDMs for convenience
		GDMgrid.insert(GDMgrid.end(), GDMrp.begin(), GDMrp.end());
		GDMgridOut.insert(GDMgridOut.end(), GDMrpOut.begin(), GDMrpOut.end());
		GridsGDM_gr.insert(GridsGDM_gr.end(), GridsGDM_rp.begin(), GridsGDM_rp.end());

		msg = Calculation(Geo1, Geo2, GDMgrid, GDMgridOut, ni, nj, nk, GridsGeo, GridsGDM_gr);
		cout << msg;

		// write to files
		GDMgridOut.erase(GDMgridOut.begin() + NOUTGRID, GDMgridOut.begin() + NOUTGRID + NOUTRP);	// return the vector to its original length
		cout << "Writing " << file_GDM_GR2 << "\n";
		WriteGrids(file_GDM_GR2, Nijk, GDMgridOut, GridsGDM_gr_out, "/", 80);
		cout << "Writing " << file_GDM_RP2 << "\n";
		WriteGrids(file_GDM_RP2, Nijk, GDMrpOut, GridsGDM_rp_out, "/", 80);

		FreeGrid(Geo1);			// free memory
		FreeGrid(Geo2);
		FreeGrid(GDMgrid);
		//FreeGrid(GDMrp);		// effectively, freed above
		FreeGrid(GDMgridOut);
		FreeGrid(GDMrpOut);

		fclose(FileContr);
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
int StrIndex(const string &S, const vector<string> &VECS)
{
	for (size_t i = 0; i < VECS.size(); i++)
		if (S.find(VECS[i]) != string::npos)
			return i;

	return -1;
}
//---------------------------------------------------------------------------------------
string ReadGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2)
{
	const int Buff = 4096;
	char strmsg[Buff], str0[Buff];
	char *str = 0;
	FILE *File = fopen(file, "r");

	size_t GridCount(0);	// counts total grids already read
	size_t ValCount(0);		// counts values read for the current grid
	size_t c;				// index within currently read grid
	bool seek_beg = true;
	bool new_line = true;
	int ind;				// index within the grid names array

	if (File == 0)
		return (string)"Cannot open " + file + "\n";

	while (ReadTokenComm(File, &str, new_line, str0, Buff))	// reads a token to str, ignoring comments
	{
		string S = str;
		if (seek_beg)
		{
			ind = StrIndex(S, S1);
			if (ind != -1)								// found the starting string for grid-i
			{
				seek_beg = false;
				ValCount = 0;
				c = 0;
				continue;
			}
		}
		if (!seek_beg && S.find(S2) != string::npos)	// found the ending string
		{
			seek_beg = true;
			if (ValCount < len)
			{
				fclose(File);
				sprintf(strmsg, " grid contains less values (%zd) than NI*NJ*NK (%zd)\n", ValCount, len);
				return string(file) + ": " + S1[ind] + string(strmsg);
			}

			GridCount++;
			if (GridCount >= S1.size())		// all grids have been read
			{
				fclose(File);
				return "";
			}

			continue;
		}
		if (!seek_beg)	// reading the main data
		{
			size_t cnt;
			double d;
			int read;

			read = sscanf(str, "%zd*%lg", &cnt, &d);
			if (read == 2)	// RPT*VAL successfully read
			{
				ValCount += cnt;
				if (ValCount > len)
				{
					fclose(File);
					sprintf(strmsg, " grid contains more values than NI*NJ*NK = %zd\n", len);
					return string(file) + ": " + S1[ind] + string(strmsg);
				}
				for (size_t i = 0; i < cnt; i++)
				{
					data[ind][c] = d;
					c++;
				}
			}
			else
			{
				read = sscanf(str, "%lg", &d);
				if (read == 1) // VAL successfully read
				{
					ValCount += 1;
					if (ValCount > len)
					{
						fclose(File);
						sprintf(strmsg, " grid contains more values than NI*NJ*NK = %zd\n", len);
						return string(file) + ": " + S1[ind] + string(strmsg);
					}
					data[ind][c] = d;
					c++;
				}
				else
				{
					fclose(File);
					sprintf(strmsg, " grid contains non-numeric symbol %s\n", str);
					return string(file) + ": " + S1[ind] + string(strmsg);
				}
			}
		}
	}
	fclose(File);

	if (!seek_beg)
	{
		if (ValCount < len)
		{
			if (ind == -1)
				return "Index error!\n";

			sprintf(strmsg, " grid contains less values [%zd] than NI*NJ*NK [%zd]\n", ValCount, len);
			return string(file) + ": " + S1[ind] + string(strmsg);
		}
		GridCount++;
	}

	if (GridCount < S1.size())
	{
		sprintf(strmsg, "Only %zd grid(s) found out of %zd\n", GridCount, S1.size());
		return string(file) + ": " + string(strmsg);
	}


	return "";
}
//---------------------------------------------------------------------------------------
void WriteGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, vector<string> pat, int slen)
{	// ver.3	1.39

	if (data.size() != S1.size())
	{
		cerr << "WRONG ARRAY SIZES IN WriteGrids\n";
		return;
	}

	ofstream F(file, ios_base::out);

	const int Buff = 128;
	char str[Buff], str0[Buff];
	int count;

	for (size_t i = 0; i < S1.size(); i++)	// different grids
	{
		F << S1[i];

		F << "\r\n ";
		count = 1;
		for (size_t c = 0; c < len; c++)	// different grid values
		{
			size_t c1 = DiffValInd(data[i], c, len);
			if (c1 > c+1)	// more than one value to write
			{
				count += sprintf(str0, " %zd*", c1-c);
				F << str0;
			}
			else
			{
				count += 1;
				F << " ";
			}

			count += sprintf(str, pat[i].c_str(), data[i][c]);
			F << str;

			if (count > slen)
			{
				count = 1;
				F << "\r\n ";
			}

			c = c1 - 1;
		}

		F << " " << S2 << "\r\n\r\n";
	}

	F.close();
}
//---------------------------------------------------------------------------------------
void WriteGrids(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, int slen)
{
	vector<string> pat(S1.size());
	for (size_t i = 0; i< S1.size(); i++)
		if (S1[i] == "ACTNUM")
			pat[i] = "%.0f";
		else
			pat[i] = "%.4f";

	WriteGrids(file, len, data, S1, S2, pat, slen);
}
//---------------------------------------------------------------------------------------
void WriteGridsTest(const char *file, size_t len, vector<double*> data, vector<string> S1, string S2, int slen)
{
	ofstream F(file, ios_base::out);

	for (size_t i = 0; i < S1.size(); i++)	// different grids
	{
		F << S1[i] << "\r\n";

		int cnt = 0;
		for (size_t c = 0; c < len; c++)	// different grid values
		{
			F << data[i][c] << " ";
			cnt++;
			if (cnt == slen)
			{
				F << "\r\n";
				cnt = 0;
			}
		}

		F << " " << S2 << "\r\n\r\n";
	}

	F.close();
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
template <typename T>
size_t DiffValInd(T *arr, size_t i0, size_t len)
{
	for (size_t i = i0; i < len; i++)
		if (arr[i] != arr[i0])
			return i;

	return len;
}
//---------------------------------------------------------------------------------------
void InitGrid(vector<double*> &G, size_t LEN)
{
	for (size_t i = 0; i < G.size(); i++)
		G[i] = new double[LEN];
}
//---------------------------------------------------------------------------------------
void FreeGrid(vector<double*> &G)
{
	for (size_t i = 0; i < G.size(); i++)
		delete [] G[i];
}
//---------------------------------------------------------------------------------------
void StripChars(char *str, const char *rm)
{
	char *aux = strpbrk(str, rm);
	if (aux != 0)
		*aux = 0;
}
//---------------------------------------------------------------------------------------
string RptGridRead(const char *fname, vector<string> grids)
{
	string res = fname;
	res += ": loaded grids ";
	for (size_t i = 0; i < grids.size()-1; i++)
		res += grids[i] + ", ";

	return res + grids[grids.size()-1];
}
//---------------------------------------------------------------------------------------
