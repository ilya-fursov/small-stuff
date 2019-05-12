#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char **argv)
{
	const int Buff = 1024;
	char str[Buff];
	FILE *File;
	
	if (argc != 2 && argc != 5 && argc != 8)
	{
		cerr << "Expected 1 argument: file -> Program counts number of values;\n";
		cerr << "Or 4 arguments: file, NI, NJ, NK -> Program reads all values;\n";
		cerr << "Or 7 arguments: file, NI, NJ, NK, i, j, k -> Program reports value at (i, j, k);\n";
		return 1;
	}	
	
	double *Vals = 0;
	int ni(0), nj(0), nk(0);
	size_t Nij, Nijk;
	bool read_val = false;
	if (argc == 5 || argc == 8)
	{
		int count = 0;
		count += sscanf(argv[2], "%d", &ni);
		count += sscanf(argv[3], "%d", &nj);
		count += sscanf(argv[4], "%d", &nk);
		if (count != 3 || ni == 0 || nj == 0 || nk == 0)
		{
			cerr << "Three positive integer values expected: NI, NJ, NK.\n";
			return 1;
		}
		read_val = true;
		Nij = ni*nj;
		Nijk = ni*nj*nk;
		Vals = new double[Nijk];
	}
	
	int i1(0), j1(0), k1(0);
	if (argc == 8)
	{
		int count = 0;
		count += sscanf(argv[5], "%d", &i1);
		count += sscanf(argv[6], "%d", &j1);
		count += sscanf(argv[7], "%d", &k1);
		if (count != 3 || i1 <= 0 || j1 <= 0 || k1 <= 0 || i1 > ni || j1 > nj || k1 > nk)
		{
			cerr << "Three integer indices expected: i=1..NI, j=1..NJ, k=1..NK.\n";
			return 1;
		}
	}
	
	File = fopen(argv[1], "r");
	size_t TotCount(0);
	size_t c = 0;
	while (!feof(File))
	{
		if (fscanf(File, " %s ", str))
		{
			int cnt;
			double d;
			int read;
			
			read = sscanf(str, "%d*%lg", &cnt, &d);
			if (read == 2)	// RPT*VAL successfully read
			{
				TotCount += cnt;
				if (read_val)
				{
					if (TotCount > Nijk)
					{
						cerr << "File contains more values than NI*NJ*NK = " << Nijk << "\n";
						fclose(File);	
						delete [] Vals;
						return 1;
					}
					for (size_t i = 0; i < cnt; i++)
					{
						Vals[c] = d;
						c++;
					}
				}
				
				//cout << cnt << "*" << d << endl;	// DEBUG
			}
			else
			{
				read = sscanf(str, "%lg", &d);
				if (read == 1) // VAL successfully read
				{
					TotCount += 1;
					if (read_val)
					{
						if (TotCount > Nijk)
						{
							cerr << "File contains more values than NI*NJ*NK = " << Nijk << "\n";
							fclose(File);	
							delete [] Vals;
							return 1;
						}
						Vals[c] = d;
						c++;
					}				
					
					//cout << d << endl;	// DEBUG
				}
			}
		}
	}
	fclose(File);	
	
	cout << "Total values: " << TotCount << endl;
	
	if (read_val && TotCount < Nijk)
		cerr << "Total values read is less than NI*NJ*NK = " << Nijk << "\n";
	else if (argc == 8)
		cout << "Value at " << i1 << ", " << j1 << ", " << k1 << " -> " << Vals[(i1-1) + (j1-1)*ni + (k1-1)*Nij];
		
	delete [] Vals;	
	return 0;
}