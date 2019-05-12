#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 6)
	{
		cerr << "Expected 5 arguments: file0, file1, file2, grid_start, grid_end.\n"
				"Lines in file0 between grid_start and grid_end (first occurence only) will be written to file2.\n"
				"Remaining lines in file0 will be written to file1.\n";
		return 1;
	}

	ifstream F0(argv[1]);
	if (!F0.good())
	{
		cerr << "Error opening input file " << argv[1] << "\n";
		return 1;
	}
	
	ofstream F1(argv[2]);
	ofstream F2(argv[3]);
	
	const int BUFF = 4096;
	char buff[BUFF];
	string S;
	bool write1 = true;
	bool search1 = true;
	
	while (!F0.eof())
	{
		F0.getline(buff, BUFF);
		if (search1 || !write1)
			S = buff;
		
		if (search1)	// looking for the starting line
		{
			size_t ind = S.find(argv[4]);
			if (ind != string::npos)	// starting line found
			{
				write1 = false;
				search1 = false;
			}
		}
		
		if (write1)
			F1 << buff << "\n";
		else
			F2 << buff << "\n";
			
		if (!write1)	// looking for the ending line
		{
			size_t ind = S.find(argv[5]);
			if (ind != string::npos)	// ending line found
				write1 = true;
		}
	}
	
	return 0;	
}