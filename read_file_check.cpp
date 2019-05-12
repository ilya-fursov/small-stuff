
#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>

//--------------------------------------------------------------------------------------------------
enum KwdType{INTE, REAL, DOUB, CHAR, LOGI};
//--------------------------------------------------------------------------------------------------
class SmryKwd
{
protected:
	KwdType type;		// type of values, e.g. REAL
	size_t len;			// number of values, e.g. 13

	KwdType type_from_str(std::string s);	// convert string CHAR, INTE, REAL, DOUB, LOGI to KwdType
	virtual void read_array(FILE *f, size_t start, size_t count){};		// this function will be non-trivial for SmryKwdData
public:
	std::string name;	// keyword name, e.g. PARAMS

	virtual ~SmryKwd(){};
	void ReadHeader(FILE *f);		// attempts to read the header from the current position; fills 'name', 'type', 'len'
	int ReadNamedHeader(FILE *f);	// searches and reads the nearest header with name = 'name', then fills 'type', 'len'
										// returns 0 on success, 1 on failure (header not found)
	SmryKwd *ReadData(FILE *f);		// reads (from the current position) data of size 'len', type 'type'
										// returns SmryKwdData<type> with 'name', 'type', 'len' as in "this"; 'data' - as read from the file
										// the returned pointer should be *** DELETED *** in the end
	void SkipReadData(FILE *f);		// skips the block of data (from the current position) - to reach the new position in the file
	virtual void cout(){};			// cout all sorts of data -- for debug
};
//--------------------------------------------------------------------------------------------------
template <class T>
class SmryKwdData : public SmryKwd
{
protected:
	virtual void read_array(FILE *f, size_t start, size_t count);		// allocates and reads 'data' [start, start + count) based on "T"	
	
public:
	std::vector<T> data;	
	
	SmryKwdData(const SmryKwd &k) : SmryKwd(k){};
	virtual void cout();
};
//--------------------------------------------------------------------------------------------------
void InitCheckSizes();
inline void ToChar(int32_t x, char *s);			// 0-end should be set for 's' elsewhere
inline int32_t SwapEndian(int32_t x);

template <class T> 
inline T ReadVal(FILE *f);
template <> inline int32_t ReadVal<int32_t>(FILE *f);
template <> inline uint32_t ReadVal<uint32_t>(FILE *f);
template <> inline std::string ReadVal<std::string>(FILE *f);
template <> inline float ReadVal<float>(FILE *f);
template <> inline double ReadVal<double>(FILE *f);
//--------------------------------------------------------------------------------------------------
// *************
// EXAMPLE USAGE
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

int main(int argc, char *argv[])
{
	if (argc <= 1)
	{
		std::cout << "Required arguments: 1 [2]\n1 - file name to read\n2 - KEYWORD to read; if omitted, all keywords are read; if CHAR_INT is specified, simplified output is done";
		return 0;
	}
	
	FILE *f = NULL;
	try
	{
		InitCheckSizes();
		
		f = fopen(argv[1], "rb");
		if (f == NULL)
			throw Exception((std::string)"Cannot open " + argv[1]);
		
		std::string KWD = "";
		if (argc >= 3)
		{
			KWD = argv[2];
			if (KWD.size() < 8)
				KWD += std::string(8 - KWD.size(), ' ');
		}
		
		if (KWD != "CHAR_INT")		
		{
			SmryKwd K;
			while (!feof(f))
			{
				if (KWD == "")
					K.ReadHeader(f);
				else
				{
					K.name = KWD;
					K.ReadNamedHeader(f);
				}
				if (feof(f))
					continue;
				
				SmryKwd *pk = K.ReadData(f);
				pk->cout();
				
				delete pk;
			}
		}
		else	// simplified mode: output every 4-bytes segment
		{
			while (!feof(f))
			{
				int a = ReadVal<int32_t>(f);
				char s[5];
				ToChar(a, s);
				s[4] = 0;
				printf("%8d\t|\t%s\n", a, s);
			}
		}
		
		fclose(f);
	}
	catch (const std::exception &e)
	{
		fclose(f);
		std::cerr << e.what() << "\n";
	}
}
// *************
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
void InitCheckSizes()
{
	assert(sizeof(float) == sizeof(int32_t));
	assert(sizeof(double) == sizeof(int64_t));
}
//--------------------------------------------------------------------------------------------------
inline void ToChar(int32_t x, char *s)			// 0-end should be set for 's' elsewhere
{
	int32_t mask = 0xff000000;

	s[0] = (x >> 24) & 0xFF;
	s[1] = (x >> 16) & 0xFF;
	s[2] = (x >> 8) & 0xFF;
	s[3] = x & 0xFF;
}
//--------------------------------------------------------------------------------------------------
inline int32_t SwapEndian(int32_t x)
{
	int32_t b0, b1, b2, b3;

	b0 = (x & 0x000000ff) << 24;
	b1 = (x & 0x0000ff00) << 8;
	b2 = (x & 0x00ff0000) >> 8;
	b3 = (x & 0xff000000) >> 24;

	return b0 | b1 | b2 | b3;
}
//--------------------------------------------------------------------------------------------------
template <> 
inline int32_t ReadVal<int32_t>(FILE *f)
{
	int32_t res = 0;
	fread(&res, 4, 1, f);
	
	return SwapEndian(res);
}
//--------------------------------------------------------------------------------------------------
template <> 
inline uint32_t ReadVal<uint32_t>(FILE *f)
{
	int32_t res = ReadVal<int32_t>(f);
	return *(uint32_t*)&res;	
}
//--------------------------------------------------------------------------------------------------
template <>
inline std::string ReadVal<std::string>(FILE *f)
{
	char res[9];

	fread(res, 4, 1, f);
	fread(res+4, 4, 1, f);
	res[8] = 0;
	
	return (std::string)res;
}
//--------------------------------------------------------------------------------------------------
template <>
inline float ReadVal<float>(FILE *f)
{
	int32_t res = ReadVal<int32_t>(f);
	return *(float*)&res;
}
//--------------------------------------------------------------------------------------------------
template <>
inline double ReadVal<double>(FILE *f)
{
	int64_t a1 = ReadVal<int32_t>(f);
	int64_t a2 = ReadVal<int32_t>(f);
	int64_t res = (a1 << 32) | a2;
	
	return *(double*)&res;
}
//--------------------------------------------------------------------------------------------------
// class SmryKwd
//--------------------------------------------------------------------------------------------------
KwdType SmryKwd::type_from_str(std::string s)
{
	if (s == "INTE")
		return INTE;
	else if (s == "REAL")
		return REAL;
	else if (s == "DOUB")
		return DOUB;
	else if (s == "CHAR")
		return CHAR;
	else if (s == "LOGI")
		return LOGI;
	else 
		throw Exception("Unknown value type " + s);
}
//--------------------------------------------------------------------------------------------------
void SmryKwd::ReadHeader(FILE *f)
{
	int a = ReadVal<int32_t>(f);
	if (feof(f))
		return;
	assert(a == 16);
	
	name = ReadVal<std::string>(f);
	len = ReadVal<int32_t>(f);
	
	a = ReadVal<int32_t>(f);
	char valtype[5];
	ToChar(a, valtype);
	valtype[4] = 0;
	type = type_from_str(valtype);
}
//--------------------------------------------------------------------------------------------------
int SmryKwd::ReadNamedHeader(FILE *f)
{
	std::string seek_name = name;
	bool finished = false;
	
	while (!finished && !feof(f))
	{
		ReadHeader(f);
		if (name == seek_name)
			finished = true;
		else
			SkipReadData(f);
	}
	
	if (finished)
		return 0;
	else
		return 1;
}
//--------------------------------------------------------------------------------------------------
SmryKwd *SmryKwd::ReadData(FILE *f)
{
	SmryKwd *res = 0;
	if (type == INTE)
		res = new SmryKwdData<int32_t>(*this);
	else if (type == REAL)
		res = new SmryKwdData<float>(*this);
	else if (type == DOUB)
		res = new SmryKwdData<double>(*this);
	else if (type == CHAR)
		res = new SmryKwdData<std::string>(*this);
	else if (type == LOGI)
		res = new SmryKwdData<uint32_t>(*this);
	else
		throw Exception ("Unknown value type " + type);

	int mult = 4;							// sizeof(type)
	if (type == CHAR || type == DOUB)
		mult = 8;
	
	int a = ReadVal<int32_t>(f);
	assert(a == 16);

	size_t total_read = 0;					// number of values currently read
	while (total_read < len)
	{
		int totlen = ReadVal<int32_t>(f);	// number of bytes in current block
		assert(totlen % mult == 0);
		totlen /= mult;						// number of values to read in current block
	
		res->read_array(f, total_read, totlen);
		
		a = ReadVal<int32_t>(f);
		assert(a == totlen*mult);
		total_read += totlen;
	}
	
	return res;
};
//--------------------------------------------------------------------------------------------------
void SmryKwd::SkipReadData(FILE *f)
{
	int mult = 4;							// sizeof(type)
	if (type == CHAR || type == DOUB)
		mult = 8;
	
	int a = ReadVal<int32_t>(f);
	if (feof(f))
		return;
	assert(a == 16);

	size_t total_read = 0;
	while (total_read < len)
	{
		int totlen = ReadVal<int32_t>(f);	// number of bytes in current block
		assert(totlen % mult == 0);
		
		fseek(f, totlen, SEEK_CUR);
		
		a = ReadVal<int32_t>(f);
		assert(a == totlen);
		total_read += totlen/mult;
	}
}
//--------------------------------------------------------------------------------------------------
// class SmryKwdData
//--------------------------------------------------------------------------------------------------
template <class T>
void SmryKwdData<T>::read_array(FILE *f, size_t start, size_t count)
{
	assert(start < len && start + count <= len);
	if (start == 0)
		data = std::vector<T>(len);
	
	for (size_t i = start; i < start + count; i++)
		data[i] = ReadVal<T>(f);
}
//--------------------------------------------------------------------------------------------------
template <class T>
void SmryKwdData<T>::cout()
{
	assert(len == data.size());
	
	std::cout << name << "\t" << len << "\t" << type << "\n";
	for (size_t i = 0; i < len; i++)
	{
		std::cout << data[i];
		if (i < len-1)
			std::cout << "\t";
		else 
			std::cout << "\n";
	}
}
//--------------------------------------------------------------------------------------------------