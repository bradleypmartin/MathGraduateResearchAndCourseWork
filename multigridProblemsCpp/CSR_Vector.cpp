/*
*	Sample Files created for APPM 6640 - Multigrid Methods
*
*	Written by Toby Jones and Chris Leibs
*	January 4th, 2013
*/

#include "CSR_Vector.h"
#include <stdlib.h>		//For Malloc
#include <time.h>		//For Random
#include <stdio.h>		//For printf,
#include <iostream>		//For file io using fprintf
#include <fstream>		//For ifstream


using namespace std;

//Default Constructor
//Input:
//	
//Output:
//
//Description:
//	Create a vector of length 0 
CSR_Vector::CSR_Vector( )
{
	m_size = 0;
	m_data = 0;
}

//Default Constructor
//Input:
//	Size - Size of Vector
//Output:
//
//Description:
//	Create a vector of length Size ( allocate memory )
CSR_Vector::CSR_Vector( unsigned int Size )
{
	m_size = Size;
	m_data = (double*)malloc( sizeof(double)*m_size);
}

//Default Destructor
//Input:
//	
//Output:
//
//Description:
//	Reset variables
//	free allocated memory
CSR_Vector::~CSR_Vector()
{
	FreeMemory();
	m_size = 0;
}


//FreeMemory
//Input:
//	
//Output:
//
//Description:
//	free allocated memory
void CSR_Vector::FreeMemory()
{
	free( m_data );
	m_data = 0;
}


//Length
//Input:
//	
//Output:
//	Length of Vector
//Description:
//	return the length of the vector
unsigned int CSR_Vector::Length()
{
	return m_size;
}
	
//Set To Zero
//Input:
//	
//Output:
//	
//Description:
//	Set all elements of vector to Zero
void CSR_Vector::SetZero()
{
	for(unsigned int i = 0; i < m_size; i++ )
		m_data[i] = 0.0f;
}
//Set
//Input:
//	CSR_Vector
//Output:
//	none
//Description:
//	copies the data from the input to the vector
void CSR_Vector::Set(CSR_Vector* u )
{
	if( u->Length() != m_size )
	{
		cout << "Mismatched attempt to set from another vector - defaulting to 0" << endl;
		SetZero();
		return;
	}

	//get the data
	for(unsigned int i = 0; i < m_size; i++ )
		m_data[i] = u->GetValue( i );
}

//Set To Constant
//Input:
//	Value 
//Output:
//	
//Description:
//	Set all elements of vector to Value
void CSR_Vector::SetConstant( double Value )
{
	for(unsigned int i = 0; i < m_size; i++ )
		m_data[i] = Value;
}

//Set To Random
//Input:
//	
//Output:
//	
//Description:
//	Set all elements of vector to Random within range 0-1
void CSR_Vector::SetRandom()
{
	//seed random number generator
	srand ( (int)time(NULL) );
	
	//set to random numbers
	for(unsigned int i = 0; i < m_size; i++ )
		m_data[i] = (double)(rand() /(((double)RAND_MAX)));	

}

//GetValue
//Input:
//	Index
//Output:
//	
//Description:
//	Return value of vector at index.
//	If Index is not in range, prints warning and returns 0.0f;
double CSR_Vector::GetValue( unsigned int Index )
{
	if( Index < m_size )
		return m_data[Index];
	printf( "Warning : Tried to index out of allocated memory in CSR_Vector - returning 0.0 as default\n" );
	return 0.0f;
}

//SetValue
//Input:
//	Index and Value
//Output:
//	
//Description:
//	Set the value of the vector to Value at location Index
void CSR_Vector::SetValue( unsigned int Index, double Value )
{
	if( Index < m_size )
		m_data[Index] = Value;
	else
		printf( "Warning : Tried to index out of allocated memory in CSR_Vector - no value has been set\n" );
}

//Add
//Input:
//	CSR_Vector and Factor
//Output:
//	
//Description:
//	Add a vector to data with multiple of Factor
void CSR_Vector::Add( CSR_Vector* Vec, double Factor )
{
	for(unsigned int i = 0; i < m_size; i++ )
		m_data[i] += Factor * Vec->GetValue(i);
}

//Matlab Print
//Input:
//	FileName
//Output:
//	
//Description:
//	Write a text file matlab can load
//		write each value seperated by one space
void CSR_Vector::Matlab_Print( char* FileName )
{
		
	FILE * pFile;

	pFile = fopen( FileName, "w" );
	
	for( unsigned int j = 0; j < m_size; j++ )
	{
		fprintf( pFile,  "%.17g ", m_data[j] );
	}

	fclose( pFile );
}

//Debug_Print
//Input:
//	
//Output:
//	
//Description:
//	Prints the vector to standard out
void CSR_Vector::Debug_Print()
{
	printf("CSR_Vector Print:\n-----\n");
	for(unsigned int i = 0; i < m_size; i++)
		printf("%f\n",m_data[i]);
	printf("-----\n");
}

//Load
//Input:
//	FileName
//Output:
//	
//Description:
//	Read the loadable version of the Vector which is assumed in the format
//		1. Length.
//		2. Each element on a new line.
void CSR_Vector::Load( char* FileName )
{

	//clean up existing
	m_size = 0;
	FreeMemory();

	//open the file
	ifstream inFile(FileName);
	
	//set up temporary buffer for reading the file
	const int MAX_STRING_BUFFER = 256;
	char line[MAX_STRING_BUFFER];

	//load number of elements
	inFile.getline(line, MAX_STRING_BUFFER);
	sscanf(line, "%d", &m_size );
	//load data
	m_data = (double*)malloc( sizeof(double)*m_size);	
	for( unsigned int i = 0; i < m_size; i++) {
		inFile.getline(line, MAX_STRING_BUFFER);
		sscanf(line, "%lf", &(m_data[i]) );
	}
	
	//close file
	inFile.close();
}

//Write
//Input:
//	FileName
//Output:
//	
//Description:
//	Write a loadable version of the Vector
//		1. Write the length.
//		2. Write each element on a new line.
void CSR_Vector::Write( char* FileName )
{
	FILE * pFile;

	pFile = fopen( FileName, "w" );
	
	//Put the vector size on the first line
	fprintf(pFile, "%d\n",m_size);

	for( unsigned int j = 0; j < m_size; j++ )
	{
		fprintf( pFile,  "%.17g\n", m_data[j] );
	}
	fclose( pFile );
}

//GetData
//Input:
//	
//Output:
//	Pointer to the raw data of this vector
//Description:
//	This is useful if you have a routine written designed to work with an
//	array of double from some legavy code
double* CSR_Vector::GetData()
{
	return m_data;
}
