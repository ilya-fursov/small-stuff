#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>

//------------------------------------------------------------------------------------------
class Exception : public std::exception
{
protected:
	std::string msg;
public:
	Exception() : msg(""){};
	Exception(std::string s) : msg(s){};						
	~Exception() noexcept {};
	const char *what() const noexcept {return msg.c_str();};
};
//------------------------------------------------------------------------------------------
std::string ToUpper(const std::string &s)
{
	std::string res = s;

	for (auto &i : res)
		i = toupper(i);

	return res;
}
//------------------------------------------------------------------------------------------
template <class ContainerT>
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters, const bool trimEmpty)
{
   std::string::size_type pos, lastPos = 0;
   tokens = ContainerT();
   while(true)
   {
      pos = str.find_first_of(delimiters, lastPos);
      if(pos == std::string::npos)
      {
         pos = str.length();

         if(pos != lastPos || !trimEmpty)
            tokens.push_back(typename ContainerT::value_type(str.data()+lastPos,
                  (typename ContainerT::value_type::size_type)(pos - lastPos)));

         break;
      }
      else
      {
         if(pos != lastPos || !trimEmpty)
            tokens.push_back(typename ContainerT::value_type(str.data()+lastPos,
                  (typename ContainerT::value_type::size_type)(pos - lastPos)));
      }

      lastPos = pos + 1;
   }
};
//------------------------------------------------------------------------------------------
double StoD(std::string s, bool &complete)
{
	const char *input = s.c_str();
	char *end;
	double res = strtod(input, &end);
	if (end == input || *end != '\0')
		complete = false;
	else
		complete = true;

	return res;
}
//------------------------------------------------------------------------------------------
double StoD(std::string s)
{
	bool complete;
	double res = StoD(s, complete);
	if (!complete)
		throw Exception("Cannot convert string '" + s + "' to double");

	return res;
}
//------------------------------------------------------------------------------------------
void ReadTwoCols(std::string fname, std::vector<double> &c1, std::vector<double> &c2)	// reads [c1, c2] from 'fname'
{
	std::ifstream sr;
	sr.exceptions(std::ios_base::badbit);
	
	c1.clear();
	c2.clear();
	std::string file_delim = " \t\r\n";

	try
	{
		sr.open(fname);
		if (sr.fail())
			throw Exception("Cannot open file " + fname);

		std::string line;
		while (!sr.eof())
		{
			getline(sr, line);
			std::vector<std::string> SA;
			tokenize(line, SA, file_delim, true);
			if (SA.size() > 0)
			{
				if (SA.size() != 2)
					throw Exception("Wrong file format: 2 values per line expected");

				c1.push_back(StoD(SA[0]));
				c2.push_back(StoD(SA[1]));
			}
		}
		sr.close();
	}
	catch (const std::exception &e)
	{
		if (sr.is_open())
			sr.close();
		throw;
	}
}
//------------------------------------------------------------------------------------------
void CalcSWOF(std::string type, double S, const std::vector<double> &params, double &krw, double &kro)	// params - vector of size 6
{
	if (S < 0 || S > 1)
		throw Exception("'S' should be in [0, 1]");
	if (type != "COREY" && type != "CHIERICI" && type != "LET")
		throw Exception("Type should be COREY, CHIERICI, or LET");

	if (type == "COREY")
	{
		if (params.size() < 2)
			throw Exception("For COREY type two numbers should be specified");
		
		double Nw = params[0];
		double No = params[1];

		kro = pow(1-S, No);
		krw = pow(S, Nw);
	}
	else if (type == "CHIERICI")
	{
		if (params.size() < 4)
			throw Exception("For CHIERICI type four numbers should be specified");

		double Aw = params[0];
		double Lw = params[1];
		double Ao = params[2];
		double Lo = params[3];

		double Rw = S/(1 - S);
		kro = exp(-Ao * pow(Rw, Lo));
		krw = exp(-Aw * pow(Rw, -Lw));
		
		if (S == 0)
			krw = 0;
		if (S == 1)
			kro = 0;
	}
	else if (type == "LET")
	{
		if (params.size() < 6)
			throw Exception("For LET type six numbers should be specified");
		
		double Lw = params[0];
		double Ew = params[1];
		double Tw = params[2];
		double Lo = params[3];
		double Eo = params[4];
		double To = params[5];

		kro = pow(1-S, Lo)/(pow(1-S, Lo) + Eo*pow(S, To));
		krw = pow(S, Lw)/(pow(S, Lw) + Ew*pow(1-S, Tw));
	}
}
//------------------------------------------------------------------------------------------
void CheckSw(const std::vector<double> &Sw)		// throws exception if Sw is not appropriate
{
	if (Sw.size() == 0)
		throw Exception("Saturation array is empty");
	if (Sw[0] != 0)
		throw Exception("First saturation value should be 0");
	if (*--Sw.end() != 1)
		throw Exception("Last saturation value should be 1");
	for (int i = 0; i < int(Sw.size())-1; i++)
		if (Sw[i] >= Sw[i+1])
			throw Exception("Saturation values are not monotonically increasing");
}
//------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	try
	{
		if (argc == 1)
		{
			std::cout << "Calculation of (normalized) relative permeability tables (SWOF or SGOF).\n"
						 "Syntax:\n"
						 "\t./swof.exe Pc.txt COREY Nw No\n"
						 "\t./swof.exe Pc.txt CHIERICI Aw Lw Ao Lo\n"
						 "\t./swof.exe Pc.txt LET Lw Ew Tw Lo Eo To\n"
						 "Pc.txt may be any file name; the file should contain two columns with numbers [Saturation, Pc], with\n"
						 "\t'Saturation' monotonically increasing from 0 to 1, and\n"
						 "\t'Pc' (capillary pressure) - monotonically decreasing.\n"
						 "Values from 'Saturation' will be used in the output table.\n";
			return 0;
		}
		if (argc < 3)
			throw Exception("File name and/or Type are not specified");
		
		std::string fname = argv[1];
		std::string type = argv[2];
		type = ToUpper(type);
		
		std::vector<double> params;
		for (int i = 3; i < argc; i++)
			params.push_back(StoD((std::string)argv[i]));
		
		std::vector<double> Sw, Pc;
		ReadTwoCols(fname, Sw, Pc);
		CheckSw(Sw);
		
		for (size_t i = 0; i < Sw.size(); i++)	// output the final table
		{
			double krw, kro;
			CalcSWOF(type, Sw[i], params, krw, kro);
			printf("%-16.11g\t%-16.11g\t%-16.11g\t%-16.11g\n", Sw[i], krw, kro, Pc[i]);
		}
		printf("/\n");
	}
	catch (std::exception &e)
	{
		std::cout << "ERROR: " << e.what() << "\n";
		return 1;
	}
	
	return 0;
}