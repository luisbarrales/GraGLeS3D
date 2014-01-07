#ifndef _io_h_
#define _io_h_

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "applic.h"

#define BUFSIZE 1024
#define GGAPPROACH_MAXCASIZE 64 //20130214, 512 grains is sufficent, current limitation, that could be higher but, 512 grains in 100000 data lines with double precision means 200MB of RAM to cache, only multiples of 8Grains

//#define INCLUDECELLSIZE true  //option to allow for cell size commitment to perform digital line interception with mstoPNG, (set true for mstoGrainSize compatible binaries, in which case you also have to set it to true in mstoPNG's mstoPNG.cpp!!)

#define xWinLoc (30)
#define yWinLoc (50)

#define xWinSize (1000)
#define yWinSize (700)

#define xWinCentre (xWinSize/2)
#define yWinCentre (yWinSize/2)

#define stringsNotEqual(a,b) strcmp(a,b)
#define stringsEqual(a,b) !stringsNotEqual(a,b)

#define READ   0
#define WRITE  1
#define APPEND 2
#define FAILURE 0
#define SUCCESS 1

typedef double Real;

class dataLine;
class dataBlock; 
typedef dataLine * dataLineP;
typedef dataBlock * dataBlockP;
class io;
typedef io * ioP;
class error;
typedef error * errorP;

using namespace std;

typedef struct
{
	union
	{
		float f;
		long i;
		char *s;
	} d;
	char type;

} univData;
typedef univData * univDataP;

class dataLine
{
public:
	dataLine( void );
	dataLine( int size );
	~dataLine( void );

	univDataP dat;
	long dataCount;

	dataLineP next;
	dataLineP prev;
};

class dataBlock
{
public:
	dataBlock( void );
	~dataBlock( void );
	long columnCount;
	long lineCount;
	
	dataLineP first;
	dataLineP last;

	char head[BUFSIZE];
	char name[128];
};

// MPI_IO structs
typedef struct
{
    double time;  //###those values right now only single precision...
    double rx;
    double count;  //should be long ###was float

} profilingData;


typedef struct
{
	double defgrainscnt[GGAPPROACH_MAXCASIZE];  //### this array can in case of 8 grains be too large, better made dynamic later on
	//###all elements which are not necessary are never initialized but also never utilized, bad style I do know, but works for case study
	/*long nucleuscnt;
	long defgrainscnt1;
	long defgrainscnt2;
	long defgrainscnt3;
	long defgrainscnt4;
	long defgrainscnt5;
	long defgrainscnt6;
	long defgrainscnt7;
	long defgrainscnt8;*/
} deformedgrainsGG;



/*MPI_IO structs for GGApproach to avoid oribucket and thousands of file descriptors
typedef struct
{
    float time;
    float rxgrainscnt;
	float deformedgrainscnt[8];  //should be long
} profilingDataGG;*/

typedef struct					//calculation not platform independent!
{
    int OriIndex;               // 4 byte
    float phi1,PHI,phi2;        // 3*4 byte  //Ori Output only single precision ###!
} MPI_IO_OriIndex  ;            //// 16 byte

typedef struct
{
    int MPIRank;              // 4 byte
    int x0, xmax, y0, ymax, z0, zmax;   // 6*4 byte
} MPI_IO_Node;                  //// 28 byte

typedef struct
{
    int OriIndex;               // 4 byte
    char state;                 // 1 byte
} MPI_IO_Cell;                  //// 5 byte

typedef struct
{
	float InfProgress;			// 4 byte  ###coded as 255*(1- c->rxFraction/1) (black is deformed, white is 100% recrystallized)
	char state;					// 1 byte
} MPI_IO_CellInfState;

typedef struct
{
	int OriIndex;				// 4 byte
	float StateAndInfProgress; // 4 byte -1000 for particles, -1 for deformed, 0 to 1 for partially recrystallized

} MPI_IO_CellState;
// ---




class io
{

public:
	io( void  );
	~io( void );

	short open( short rw, const char * filename, FILE ** file );
	short open( const char * filename, ofstream * file );
	short open( const char * filename, ifstream file );
	short write( char * message, FILE * file );
	short write( string message, ofstream * file );

	dataBlockP readDataBlock( const char *name, const char *fileName );
	dataBlockP readDataBlock( const char *name, FILE * file );

	char *newString(const char *s);
	Real getReal( dataLineP line, long column );
	Real getReal2( dataLineP line, long column );
	long getInt( dataLineP line, long column );
	char * getString( dataLineP line, long column );

	Real geTReal( const char *s, dataBlockP db );   
	long geTInt( const char *s, dataBlockP db );
	char * geTString( const char *s, dataBlockP db );

};

class error
{
public:
	
	error( ofstream * log = NULL, int rank = 0 ) :
	  logFile(log), myRank(rank)
	  {
	  }

	~error( void );

	void reportError( char * message, int rank = 0);
	void reportWarning( char * message, int rank = 0);

	ofstream * logFile;
	io ioHdl;
	int myRank;
};

#endif

