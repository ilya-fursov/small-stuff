
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>
#include <cassert>

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
std::vector<double> Autocorr(const std::vector<double> &v);		// autocorrelation for the given series
int Ess(const std::vector<double> &ac, double &res);	// effective sample size by Initial Monotone Sequence Estimator for the given series (res), and lag where it was reached (return)
														// "ac" should be Autocorr of the series in question
double Ess_BW(const std::vector<double> &ac, int M);			// effective sample size by Bartlett Window, with truncatoin at lag "M"
																// "ac" should be Autocorr of the series in question
double mean(const std::vector<double> &v);				// mean of the given series
template <class T>
inline T Min(T x, T y){return (x < y) ? x : y;};
template <class T>
inline T Max(T x, T y){return (x > y) ? x : y;};
//------------------------------------------------------------------------------------------
// creates histogram {xvals, yvals} from 'data'
// X axis interval [xmin, xmax] is split into 'bins' categories
// 'xvals' gets values from the center of each category
// if input xmin >= xmax, then 'xmin' and 'xmax' are found from data
// if only_count == true, then 'yvals' is the count of points in each category,
// otherwise 'yvals' is normalized to make the histogram a PDF (integral = 1)
// function returns the total number of points counted
int Histogram(const std::vector<double> &data, std::vector<double> &xvals, std::vector<double> &yvals, double &xmin, double &xmax, int bins, bool only_count = false)
{
	if (bins <= 0)
		throw Exception("Number of bins <= 0 in Histogram()");
	
	size_t len = data.size();
	if (xmin >= xmax)	// min and max from data
	{
		xmin = std::numeric_limits<double>::max();
		xmax = std::numeric_limits<double>::lowest();
		for (size_t i = 0; i < len; i++)
		{
			if (data[i] < xmin)
				xmin = data[i];
			if (data[i] > xmax)
				xmax = data[i];
		}
		if (xmax == xmin)		// when data.size() == 1
			xmax = xmin + 1;	
	}
	
	double dh = (xmax - xmin)/bins;
	xvals = std::vector<double>(bins);
	yvals = std::vector<double>(bins, 0);
	for (int i = 0; i < bins; i++)
		xvals[i] = xmin + i*dh + dh/2;
	
	// fill yvals
	int count = 0;
	for (size_t i = 0; i < len; i++)
	{
		double x = data[i];
		if (xmin <= x && x <= xmax)
		{
			int ind = (x - xmin)/dh;
			if (ind >= bins)
				ind = bins-1;
			yvals[ind] += 1;
			count++;
		}
	}
	
	// renormalize
	if (!only_count && count > 0)
		for (auto &i : yvals)
			i /= (count * dh);
		
	return count;
}
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
long StoL(std::string s, bool &complete);		// complete = true if whole string is a number
double StoD(std::string s, bool &complete);		// complete = true if whole string is a number
long StoL(std::string s);		// throws exception if whole string is not a number
double StoD(std::string s);		// throws exception if whole string is not a number
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
// MAIN
//------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	size_t nbins = 10;
	size_t col = 1;
	size_t row1 = 1;
	size_t row2 = -1;
	double min = std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::lowest();
	FILE *in = 0, *out = 0;
	
	try
	{
		if (argc <= 2)
		{
			std::cout << "Calculation of histogram (and effective sample size) from input data.\n" \
						 "Syntax: ./Histogram input.txt output.txt [nbins] [col] [row1] [row2] [min] [max]\n"
						 "If [optional parameters] from the right are omitted, default values are used.\n"
						 "input.txt - input file name (up to 131,071 characters per line)\n"
						 "output.txt - output file name\n"
						 "nbins (int) - number of bins (categories), default = 10\n"
						 "col (int) - column in the input file to take data from, default = 1\n"
						 "row1 (int) - take data starting from this row, default = 1\n"
						 "row2 (int) - take data up to this row inclusive (-1 means till the end of file), default = -1\n"
						 "min (double) - minimum bound of the interval (X axis) which is split into categories, default = calculate from data\n"
						 "max (double) - maximum bound of the interval (X axis) which is split into categories, default = calculate from data\n";
			return 1;
		}				 

		
		if (argc > 3)
			nbins = StoL(argv[3]);
		if (argc > 4)
			col = StoL(argv[4]);
		if (argc > 5)
			row1 = StoL(argv[5]);
		if (argc > 6)
			row2 = StoL(argv[6]);
		if (argc > 7)
			min = StoD(argv[7]);
		if (argc > 8)
			max = StoD(argv[8]);	
		
		in = fopen(argv[1], "r");
		if (!in)
			throw Exception("Cannot open file " + std::string(argv[1]));
		const int LINEBUFF = 131072;
		const int SMALLBUFF = 512;
		const int row2_guess = 1000;	// blind guess for capacity of 'res'
		char line[LINEBUFF];
		std::string delim = " \t\r\n";
		std::vector<std::string> toks;
		std::vector<double> res;
		int res_cap = (row2 == size_t(-1))? row2_guess : (row2-row1+1);
		if (res_cap <= 0)
			throw Exception("row2 < row1");
		
		res.reserve(res_cap);
		
		// read the file
		int count = 1;		// line count
		while (!feof(in))		
		{	
			if (!fgets(line, LINEBUFF, in))
				break;				// may happen when empty line with EOF was read
			if (count >= row1 && count <= row2)
			{
				tokenize(line, toks, delim, true);
				if (toks.size() < col)	// check if the necessary column is present
				{
					char msg[SMALLBUFF];
					sprintf(msg, "Line %d only contains %d column(s), target column %d is absent", count, int(toks.size()), col);
					throw Exception(msg);
				}
				res.push_back(StoD(toks[col-1]));		// col is 1-based index
			}
			
			count++;
			if (count > row2)		// reached the end of range
				break;
		}
		
		fclose(in);
		in = 0;
		
		if (res.size() == 0)
			throw Exception("No data were read");
		
		// calculate the histogram
		std::vector<double> xv, yv;
		int c = Histogram(res, xv, yv, min, max, nbins);
		
		// calculate effective sample size
		double ess = 0;
		std::vector<double> ac = Autocorr(res);
		const int lag = Ess(ac, ess);
		
		const int M = Min(Max(5, lag), (int)ac.size()/2);		// truncation_lag = "lag", clipped to [5, N/2] for safety
		double ess_bw = Ess_BW(ac, M);
		
		const int M5 = Min(Max(10, 5*lag), (int)ac.size()/2);	// truncation_lag = "5 x lag", clipped to [10, N/2] for safety
		double ess_bw5 = Ess_BW(ac, M5);
		
		const int M20 = Min(Max(10, 20*lag), (int)ac.size()/2);	// truncation_lag = "20 x lag", clipped to [10, N/2] for safety
		double ess_bw20 = Ess_BW(ac, M20);		
		
		// output
		out = fopen(argv[2], "w");
		for (int i = 0; i < xv.size(); i++)
			fprintf(out, "%-14.10g\t%-14.10g\n", xv[i], yv[i]);
		
		fclose(out);
		out = 0;
		
		std::cout << "Reading lines " << row1 << " to " << count-1 << " in column " << col << "\n";
		std::cout << "Histogram with " << nbins << " bins in [" << min << ", " << max << "], total points counted: " \
				  << c << ", Y-values are showing PDF\n";
		
		printf("ESS = %10.6g (lag %4d)\tESS_BW (M = %4d): %10.6g\tESS_BW (M = %4d): %10.6g\tESS_BW (M = %4d): %10.6g\n", ess, lag, M, ess_bw, M5, ess_bw5, M20, ess_bw20);
		std::cout << "mean = " << mean(res) << "\n";
		
		// output autocorr to file, for reference
		{
			FILE *f = fopen("autocorr_out_file.txt", "w");
			for (size_t i = 0; i < ac.size(); i++)
				fprintf(f, "%g\n", ac[i]);
			fclose(f);		
		}
	}
	catch (std::exception &e)
	{
		fclose(in);
		fclose(out);
		std::cerr << "ERROR: " << e.what() << "\n";
	}
	
	return 0;
}
//------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------
long StoL(std::string s, bool &complete)
{
	const char *input = s.c_str();	// 6.10.2013, using C++98 conversion
	char *end;
	double res = strtol(input, &end, 10);
	if (end == input || *end != '\0')
		complete = false;
	else
		complete = true;

	return res;
}
//------------------------------------------------------------------------------------------
double StoD(std::string s, bool &complete)
{
	const char *input = s.c_str();	// 6.10.2013, using C++98 conversion
	char *end;
	double res = strtod(input, &end);
	if (end == input || *end != '\0')
		complete = false;
	else
		complete = true;

	return res;
}
//------------------------------------------------------------------------------------------
long StoL(std::string s)
{
	bool complete;
	double res = StoL(s, complete);
	if (!complete)
		throw Exception("Cannot convert string '" + s + "' to long int");

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
std::vector<double> Autocorr(const std::vector<double> &v)
{
	size_t len = v.size();
	if (len <= 1)
		throw Exception("Autocorr() should be applied to vectors of size > 1");

	std::vector<double>  res(len, 0.0);
	double *pres = res.data();
	const double *p = v.data();
	double mean = 0;
	for (const auto &i : v)
		mean += i;
	mean /= len;

	for (size_t k = 0; k < len; k++)	// k - lag
	{
		double d = 0;
		for (size_t t = 0; t < len-k; t++)
			d += (p[t] - mean)*(p[t+k] - mean);
		pres[k] = d/(len-1);			// could divide by "len" as well
		if (k > 0)
			pres[k] /= pres[0];			// normalization
	}
	pres[0] /= pres[0];		// normalization

	return res;
}
//------------------------------------------------------------------------------------------
int Ess(const std::vector<double> &ac, double &res)				// effective sample size by Initial Monotone Sequence Estimator for the given series (res), and lag where it was reached (return)
{																// "ac" should be Autocorr of the series in question
	size_t len = ac.size();
	size_t N = len/2;
	const double *pac = ac.data();

    res = 0;
    int lag = 0;
    double G_prev = std::numeric_limits<double>::max();
    for (size_t m = 0; m < N; m++)
    {
    	double G = pac[2*m] + pac[2*m+1];	// Gamma_m
		if (G > 0 && G < G_prev)
		{
			res += G;
			lag = 2*m + 1;					// lag=1  <->  [g0, g1], where gk - autocorrelation with lag k
		}									// lag=3  <->  [g0, g1], [g2, g3]
		else								// etc.
			break;

    	G_prev = G;
    }

    res = (double)len/(2*res - 1);
    return lag;
}
//------------------------------------------------------------------------------------------
double Ess_BW(const std::vector<double> &ac, int M)				// effective sample size by Bartlett Window, with truncatoin at lag "M"
{																// "ac" should be Autocorr of the series in question
	size_t len = ac.size();
	const double *pac = ac.data();
	double res = 0;
	
	assert(M < (int)len);
	for (int i = 1; i <= M; i++)
		res += (1 - double(i)/M)*ac[i];			// Bartlet window is applied here; other windows could also be considered
	
	return double(len)/(1.0 + 2*res);
}
//------------------------------------------------------------------------------------------
double mean(const std::vector<double> &v)
{
	double res = 0;
	for (const auto &x : v)
		res += x;
	
	return res/v.size();
}
//------------------------------------------------------------------------------------------