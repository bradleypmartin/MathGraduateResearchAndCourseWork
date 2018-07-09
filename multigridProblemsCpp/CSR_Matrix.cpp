/*
*	Sample Files created for APPM 6640 - Multigrid Methods
*
*	Written by Toby Jones and Chris Leibs
*	January 4th, 2013
*/

#include "CSR_Matrix.h"
#include <stdlib.h>		//For Malloc
#include <time.h>		//For Random
#include <stdio.h>		//For printf,
#include <iostream>		//For file io using fprintf
#include <fstream>		//For ifstream
#include <string.h>
#include <math.h>

using namespace std;

//Default Constructor
//Input:
//	
//Output:
//
//Description:
//	Create an empty CSR matrix
CSR_Matrix::CSR_Matrix()
{
	m_rowData = 0;
	m_colIndices = 0;
	m_rowStartIndices = 0;

	m_numRows = 0;
	m_numCols = 0;

}

//Default Destructor
//Input:
//	
//Output:
//
//Description:
//	Reset variables
//	free allocated memory
CSR_Matrix::~CSR_Matrix()
{
	ClearData();
}

//ClearData
//Input:
//	
//Output:
//
//Description:
//	Reset variables
//	free allocated memory
void CSR_Matrix::ClearData()
{
	if( m_rowData != 0 )
	{
		free( m_rowData );
		m_rowData = 0;
	}
	if( m_colIndices != 0 )
	{
		free( m_colIndices );
		m_colIndices = 0;
	}
	if( m_rowStartIndices != 0 )
	{
		free( m_rowStartIndices );
		m_rowStartIndices = 0;
	}

	m_numRows = 0;
	m_numCols = 0;

}

//Set To Diagonal
//Input:
//	Size : number of rows and columns
//	D : value of diagonal
//Output:
//
//Description:
//	free allocated memory
//  Build a diagonal CSR matrix
void CSR_Matrix::SetToDiagonal( int Size, double D )
{
	//free any memory
	ClearData();

	//square matrix
	m_numRows = Size;
	m_numCols = Size;
	
	//the actual number of non zeros this stores
	m_dataSize = Size;

	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize*sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize*sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	//just set the diagonals to D
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowData[i] = D;
		m_colIndices[i] = i;
		m_rowStartIndices[i] = i;
	}

	//let it know where the last row ends.
	m_rowStartIndices[Size] = m_numRows;
}



//Set To TriDiagonal
//Input:
//	Size : number of rows and columns
//	L : Value below diagonal
//	D : value of diagonal
//	U : Value above diagonal
//Output:
//
//Description:
//	free allocated memory
//  Build a tridiagonal CSR matrix
void CSR_Matrix::SetTridiagonal( int Size,  double L, double D, double U )
{
	//free any memory
	ClearData();

	//square matrix
	m_numRows = Size;
	m_numCols = Size;
	
	//the actual number of non zeros this stores
	m_dataSize = 3*Size -2;

	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize *sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize *sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	int numNonZeros = 0;
	//just set the diagonals to D
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowStartIndices[i] = numNonZeros ;
		if( i == 0 ){
			m_rowData[numNonZeros] = D;
			m_colIndices[numNonZeros] = i;
			numNonZeros++;

			m_rowData[numNonZeros] = U;
			m_colIndices[numNonZeros] = i+1;
			numNonZeros++;
		}
		else if ( i == m_numRows-1)
		{
			m_rowData[numNonZeros] = D;
			m_colIndices[numNonZeros] = i;
			numNonZeros++;

			m_rowData[numNonZeros] = L;
			m_colIndices[numNonZeros] = i-1;
			numNonZeros++;
		}
		else
		{
			m_rowData[numNonZeros] = D;
			m_colIndices[numNonZeros] = i;
			numNonZeros++;

			m_rowData[numNonZeros] = L;
			m_colIndices[numNonZeros] = i-1;
			numNonZeros++;

			m_rowData[numNonZeros] = U;
			m_colIndices[numNonZeros] = i+1;
			numNonZeros++;
			
		}
	}

	//let it know where the last row ends.
	m_rowStartIndices[m_numRows] = numNonZeros;
}

//Set To 2D Crossstencil (5-point)
//
//Input:
//	SideSize : number of interior points in ONE DIRECTION on the 2D grid (square root of the number of rows/cols in operator)
//	Outer: value of matrix on the 4 outside points of the stencil
//	Inner: value of matrix on the diagonal (middle of stencil)
//	
//Output:
//
//Description:
//	free allocated memory
//  Build a 2D 5-point cross stencil matrix (for 2D, 2nd order Laplacian and Romega WJ matrix operators)

void CSR_Matrix::Set2DCrossstencil(int SideSize, double Outer, double Inner)
{
	//free any memory
	ClearData();

	//square matrix
	m_numRows = pow(SideSize,2);
	m_numCols = pow(SideSize,2);
	
	//the actual number of non zeros this stores (first part: tridiagonal "centerblocks." Second: diagonal "sideblocks")
	m_dataSize = SideSize*(3*SideSize -2) + 2*(SideSize)*(SideSize-1);

	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize *sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize *sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	int numNonZeros = 0;
	//just set the diagonals to D
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowStartIndices[i] = numNonZeros;

		if( i%SideSize == 0 )
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize-1);
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize+1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;
			}
		}
		else if ( i%SideSize == (SideSize-1))
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1);
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize-1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;
			}
		}
		else
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1);
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize-1+SideSize*(i/SideSize);
			numNonZeros++;

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize+1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;
			}
		}
	}

	//let it know where the last row ends.
	m_rowStartIndices[m_numRows] = numNonZeros;
}


//Set To 2D VAR Crossstencil (5-point)
//
//Input:
//	SideSize : number of interior points in ONE DIRECTION on the 2D grid (square root of the number of rows/cols in operator)
//	Outer: value of matrix on the 4 outside points of the stencil
//	Inner: value of matrix on the diagonal (middle of stencil)
//	
//Output:
//
//Description:
//	free allocated memory
//  Build a 2D 8-point variational cross stencil matrix (for 2D, 2nd order Laplacian and Romega WJ matrix operators)

void CSR_Matrix::Set2DVARCrossstencil(int SideSize, double Outer, double Inner)
{
	//free any memory
	ClearData();

	//square matrix
	m_numRows = pow(SideSize,2);
	m_numCols = pow(SideSize,2);
	
	//the actual number of non zeros this stores (first part: tridiagonal "centerblocks." Second: diagonal "sideblocks")
	m_dataSize = SideSize*(3*SideSize -2) + 6*(SideSize)*(SideSize-1);

	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize *sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize *sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	int numNonZeros = 0;
	//just set the diagonals to D
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowStartIndices[i] = numNonZeros;

		if( i%SideSize == 0 )
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize-1);
				numNonZeros++;

				///new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize-1)+1;
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize+1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;

				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1)+1;
				numNonZeros++;
			}
		}
		else if ( i%SideSize == (SideSize-1))
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1)-1;
				numNonZeros++;

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1);
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize-1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1)-1;
				numNonZeros++;

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;
			}
		}
		else
		{
			m_rowData[numNonZeros] = Inner;
			m_colIndices[numNonZeros] = i%SideSize+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize >= 1)
			{
				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1)-1;
				numNonZeros++;

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1);
				numNonZeros++;

				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize-1)+1;
				numNonZeros++;
			}

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize-1+SideSize*(i/SideSize);
			numNonZeros++;

			m_rowData[numNonZeros] = Outer;
			m_colIndices[numNonZeros] = i%SideSize+1+SideSize*(i/SideSize);
			numNonZeros++;

			if (i/SideSize < (SideSize-1))
			{
				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1)-1;
				numNonZeros++;

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1);
				numNonZeros++;

				// new below

				m_rowData[numNonZeros] = Outer;
				m_colIndices[numNonZeros] = i%SideSize+(SideSize)*(i/SideSize+1)+1;
				numNonZeros++;
			}
		}
	}

	//let it know where the last row ends.
	m_rowStartIndices[m_numRows] = numNonZeros;
}



//Set2DFWrestrict
//
//Input:
//	SideSize : number of interior points in ONE DIRECTION on the 2D grid (square root of the number of rows/cols in operator)
//	
//Output:
//
//Description:
//	free allocated memory
//  Build a 2D restriction operator (to the next coarser grid) based on full weighting


void CSR_Matrix::Set2DFWrestrict(int SideSize)
{
	//free any memory
	ClearData();

	int finesidesize = SideSize;
	int coarsesidesize = ((SideSize+1)/2)-1;

	m_numRows = pow(coarsesidesize,2);
	m_numCols = pow(finesidesize,2);
	
	m_dataSize = 9*pow(coarsesidesize,2);

	//create memmory to store the non-zero data
	m_rowData = (double*)malloc( m_dataSize *sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize *sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	int numNonZeros = 0;
	
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowStartIndices[i] = numNonZeros;
		
		m_rowData[numNonZeros] = 1.0/16.0;
		m_colIndices[numNonZeros] = 2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/8.0;
		m_colIndices[numNonZeros] = 1+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/16.0;
		m_colIndices[numNonZeros] = 2+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/8.0;
		m_colIndices[numNonZeros] = 2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+finesidesize;
		numNonZeros++;
		
		m_rowData[numNonZeros] = 1.0/4.0;
		m_colIndices[numNonZeros] = 1+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/8.0;
		m_colIndices[numNonZeros] = 2+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/16.0;
		m_colIndices[numNonZeros] = 2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+2*finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/8.0;
		m_colIndices[numNonZeros] = 1+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+2*finesidesize;
		numNonZeros++;

		m_rowData[numNonZeros] = 1.0/16.0;
		m_colIndices[numNonZeros] = 2+2*(i%coarsesidesize)+(i/coarsesidesize)*2*finesidesize+2*finesidesize;
		numNonZeros++;
	}

	//let it know where the last row ends.
	m_rowStartIndices[m_numRows] = numNonZeros;
}



//Set2Dinterpolate
//
//Input:
//	SideSize : number of interior points in ONE DIRECTION on the 2D grid (square root of the number of rows/cols in operator)
//	
//Output:
//
//Description:
//	free allocated memory
//  Build a 2D interpolation operator to the next (finer) grid

void CSR_Matrix::Set2Dinterpolate(int SideSize)
{
	//free any memory
	ClearData();

	int finesidesize = (SideSize+1)*2-1;
	int coarsesidesize = SideSize;

	m_numRows = pow(finesidesize,2);
	m_numCols = pow(coarsesidesize,2);
	
	m_dataSize = 9*pow(coarsesidesize,2);

	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize *sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize *sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	int numNonZeros = 0;
	
	for( int i = 0; i < m_numRows; i++ )
	{
		m_rowStartIndices[i] = numNonZeros;
		
		int blocksubrow = i-(i/finesidesize)*finesidesize;

		if (((i/finesidesize) == 0) || ((i/finesidesize) == finesidesize-1)) // COND A
		{
			if ((i/finesidesize) == 0)
			{
				if ((blocksubrow == 0) || (blocksubrow == finesidesize-1))
				{
					if (blocksubrow == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize+2)/2-1)*coarsesidesize+(blocksubrow+2)/2-1;
						numNonZeros++;
					}
					else
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize+2)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;
					}
				}
				
				else
				{
					if (blocksubrow%2 == 1)
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+2)/2-1)*coarsesidesize+(blocksubrow+1)/2-1;
						numNonZeros++;
					}
					if (blocksubrow%2 == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+2)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+2)/2-1)*coarsesidesize+(blocksubrow)/2;
						numNonZeros++;
					}
				}
			
			}

			else
			{
				if ((blocksubrow == 0) || (blocksubrow == finesidesize-1))
				{
					if (blocksubrow == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow+2)/2-1;
						numNonZeros++;
					}
					else
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;
					}
				}
				
				else
				{
					if (blocksubrow%2 == 1)
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow+1)/2-1;
						numNonZeros++;
					}
					if (blocksubrow%2 == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2;
						numNonZeros++;
					}
				}
			
			}
		}
		
		else
		{
			if ((i/finesidesize)%2 == 1) // COND B
			{
				if ((blocksubrow == 0) || (blocksubrow == finesidesize-1))
				{
					if (blocksubrow == 0)
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] = ((i/finesidesize+1)/2-1)*coarsesidesize+(blocksubrow+2)/2-1;
						numNonZeros++;
					}
					else
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] = ((i/finesidesize+1)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;
					}
				}
				else
				{
					if (blocksubrow%2 == 1)
					{
						m_rowData[numNonZeros] = 1.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+1)/2-1)*coarsesidesize+(blocksubrow+1)/2-1;
						numNonZeros++;
					}
					if (blocksubrow%2 == 0)
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+1)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize+1)/2-1)*coarsesidesize+(blocksubrow)/2;
						numNonZeros++;
					}
				}
			}

			if ((i/finesidesize)%2 == 0) // COND AA
			{
				if ((blocksubrow == 0) || (blocksubrow == finesidesize-1))
				{
					if (blocksubrow == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow+2)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2)*coarsesidesize+(blocksubrow+2)/2-1;
						numNonZeros++;
					}
					else
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] = ((i/finesidesize)/2)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;
					}
				}
				else
				{
					if (blocksubrow%2 == 1)
					{
						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow+1)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/2.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2)*coarsesidesize+(blocksubrow+1)/2-1;
						numNonZeros++;
					}
					if (blocksubrow%2 == 0)
					{
						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2-1)*coarsesidesize+(blocksubrow)/2;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2)*coarsesidesize+(blocksubrow)/2-1;
						numNonZeros++;

						m_rowData[numNonZeros] = 1.0/4.0;
						m_colIndices[numNonZeros] =  ((i/finesidesize)/2)*coarsesidesize+(blocksubrow)/2;
						numNonZeros++;
					}
				}
			}
		}
		
	}

	//let it know where the last row ends.
	m_rowStartIndices[m_numRows] = numNonZeros;
}

//Debug_Print
//Input:
//	
//Output:
//
//Description:
//	Print Dense matrix to std out
void CSR_Matrix::Debug_Print()
{
	
	printf( "\n" );

	//create temp memory
	double* rowData = (double*)malloc( m_numCols*sizeof(double) );
	for( int i = 0; i < m_numRows; i++ )


	{
		//set entire row to zero
		memset( rowData, 0, m_numCols * sizeof(double) );
		
		int start = m_rowStartIndices[i];
		int end = m_rowStartIndices[i+1];
		
		//set all non zeros
		for( int j = start; j < end; j++ )
		{if( m_colIndices[j] >= m_numCols || m_colIndices[j] < 0) 
				cout << "Hmm bad col index " << m_colIndices[j] << "Row " << i << "Col index parameter " << j << "start = " << start << "end = " << end << endl;
			rowData[m_colIndices[j]] = m_rowData[j];
		}
		
		//output all information
		for( int j = 0; j < m_numCols; j++ )
		{
			printf( "%2.5f ", rowData[j] );
		}
		printf( "\n" );
	}

	//clean up temp memory
	free( rowData );
	
}


//Matlab_Print
//Input:
//	FileName
//Output:
//
//Description:
//	Print Dense matrix to file that matlab can load
void CSR_Matrix::Matlab_Print( char* FileName )
{
	
	FILE * pFile;

	pFile = fopen( FileName, "w" );
	
	double* rowData = (double*)malloc( m_numCols*sizeof(double) );
	for( int i = 0; i < m_numRows; i++ )
	{
		memset( rowData, 0, m_numCols * sizeof(double) );
		
		int start = m_rowStartIndices[i];
		int end = m_rowStartIndices[i+1];
		
		for( int j = start; j < end; j++ )
		{
			if( m_colIndices[j] >= m_numCols || m_colIndices[j] < 0) 
				cout << "Hmm bad col index " << m_colIndices[j] << "Row " << i << "Col index parameter " << j << "start = " << start << "end = " << end << endl;
			rowData[m_colIndices[j]] = m_rowData[j];
		}
		
		for( int j = 0; j < m_numCols; j++ )
		{
			fprintf( pFile,  "%.17g ", rowData[j] );
		}

		//new line at end of row
		fprintf( pFile, "\n" );
	}
	free( rowData );
	
	fclose( pFile );
}

//Write
//Input:
//	FileName
//Output:
//	
//Description:
//	Write a loadable version of the Vector
//		1. Write the numRows numCols DataSize.
//		2. Write the row start information
//		3. Write the actual col/data pairings
void CSR_Matrix::Write( char* FileName )
{
	
	FILE * pFile;

	pFile = fopen( FileName, "w" );
	
	//write basic loading information
	
	//Put the vector size on the first line
	fprintf(pFile, "%d %d %d\n",m_numRows, m_numCols, m_dataSize);

	//write the row start indices
	for( int i = 0; i <= m_numRows; i++ )
	{
		fprintf(pFile, "%d\n", m_rowStartIndices[i] );
	}

	//write the data
	for( int i = 0; i < m_dataSize; i++ )
	{
		fprintf(pFile, "%d %.17g\n", m_colIndices[i], m_rowData[i]);
	}

	fclose( pFile );
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
void CSR_Matrix::Load( char* FileName )
{

	//clean up existing
	ClearData();

	//open the file
	ifstream inFile(FileName);
	
	//set up temporary buffer for reading the file
	const int MAX_STRING_BUFFER = 256;
	char line[MAX_STRING_BUFFER];

	//load number of numRows numCols, dataSize
	inFile.getline(line, MAX_STRING_BUFFER);
	sscanf(line, "%d %d %d", &m_numRows, &m_numCols, &m_dataSize);
	
	//create memmory to store the non-zero data / diagonal
	m_rowData = (double*)malloc( m_dataSize*sizeof(double) );
	//create memory to store the column indices
	m_colIndices = (int*)malloc( m_dataSize*sizeof(int) );
	//create memory to store the start of the rows + 1 to know how long last column is
	m_rowStartIndices= (int*)malloc( (m_numRows+1)*sizeof(int) );

	//load the rowStartIndices
	for( int i = 0; i <= m_numRows; i++ )
	{
		inFile.getline(line, MAX_STRING_BUFFER);
		sscanf(line, "%d", &(m_rowStartIndices[i]) );	
	}

	//load the data
	for( int i =0; i < m_dataSize; i++ )
	{
		inFile.getline(line, MAX_STRING_BUFFER);
		sscanf(line, "%d %lf", &(m_colIndices[i]), &(m_rowData[i]) );	

	}
	
	//close file
	inFile.close();
}

void CSR_Matrix::DisplaySize()
{
	cout << endl << m_numRows << endl << m_numCols << endl;
}

//ApplyToVector
//Input:
//	In : The input Vector
//Output:
//	Out : The output Vector
//Description:
//	Out = A*In ( standard vector matrix multiply )
void CSR_Matrix::ApplyToVector( CSR_Vector* In, CSR_Vector* Out )
{

  if( (int)In->Length() != m_numCols )
		printf( "Warning : Size mismatch for Matrix Vector Multiply (In)\n" );
	
  if( (int)Out->Length() != (int)m_numRows )
		printf( "Warning : Size mismatch for Matrix Vector Multiply (Out)\n" );

	for( int i = 0; i < m_numRows; i++ )
	{
		double val = 0;
		for( int j = m_rowStartIndices[i]; j < m_rowStartIndices[i+1]; j++ )
		{
			val += m_rowData[j]*In->GetValue(m_colIndices[j]);
		}
		Out->SetValue( i, val );
	}
}
	


//Get
//Input:
//	Row
//	Col
//Output:
//	The value of A[Row,Col]
//Description:
//	Note if there is not a value stored it must be 0
double CSR_Matrix::Get (int Row, int Col)
{
	if (Row >= m_numRows || Col >= m_numCols)
		printf( "Warning : invalid row or col index for Get in CSR_Matrix\n" );
	
	
	int start = m_rowStartIndices[Row];
	int end = m_rowStartIndices[Row+1];

	for( int i = start; i < end; i++ )
	{
		if( m_colIndices[i] == Col )
			return m_rowData[i];
	}

	return 0.0f;
}



//Scale
//Input:
//	ScaleAmount  : Amount to Scale by
//Output:
//	
//Description:
//	Scales all values of the matrix by ScaleAmount
void CSR_Matrix::Scale( double ScaleAmount )
{
	for( int i =0 ; i < m_dataSize; i++ )
	{
		m_rowData[i] *= ScaleAmount;
	}
}


//ExactSolve
//Input:
//	f  : Right Hand size
//Output:
//	u : Solution to Au = f
//Description:
//	Use Dense Gaussian elimination to solve Au = f
void CSR_Matrix::ExactSolve( CSR_Vector* u, CSR_Vector* f )
{	
	//Create memory to store dense expansion

	//Create pointer for 2D aray
	double**  A = (double**)malloc( sizeof( double* ) * (m_numRows+1) );

	//allocate meory for 2D array
	for( int i =0; i < m_numRows+1; i++ )
	{
		A[i] = (double*)malloc( sizeof( double ) * m_numRows );
		memset( A[i], 0, sizeof( double ) * (m_numRows) );
	}
	
	//For each row set non zeros
	for( int i = 0; i < m_numRows; i++ )
	{
		int start = m_rowStartIndices[i];
		int end = m_rowStartIndices[i+1];
		for( int j = start; j < end; j++ )
		{
			A[m_colIndices[j]][i] = m_rowData[j];
		}
		A[m_numRows][i] = f->GetValue(i);	
	}
	//Debug_PrintMatrix( A, m_numRows, m_numRows +1 );
	GSolve( A, m_numRows, u->GetData() );
	
	for( int i =0; i < m_numRows+1; i++ )
	{
		delete( A[i] );
	}
	
	delete( A );
			
}


//ExactSolve
//Input:
//	a : Stores a dense matrix as well as the right hand side in the last column
//	n : the size of a ( assumed to be n*n )
//Output:
//	x : the solution to Ax = f
//Description:
//	Perform dense guassian elimination
void CSR_Matrix::GSolve(double **a,int n,double *x)
{
   int i,j,k,maxrow;
   double tmp;

   for (i=0;i<n;i++) {

      /* Find the row with the largest first value */
      maxrow = i;
      for (j=i+1;j<n;j++) {
         if (fabs(a[i][j]) > fabs(a[i][maxrow]))
            maxrow = j;
      }

      /* Swap the maxrow and ith row */
      for (k=i;k<n+1;k++) {
         tmp = a[k][i];
         a[k][i] = a[k][maxrow];
         a[k][maxrow] = tmp;
      }

      /* Singular matrix? */
      if (fabs(a[i][i]) < .00000000001)
         printf( "Warning : Singular Matrix" );

      /* Eliminate the ith element of the jth row */
      for (j=i+1;j<n;j++) {
         for (k=n;k>=i;k--) {
            a[k][j] -= a[k][i] * a[i][j] / a[i][i];
         }
      }
   }

   /* Do the back substitution */
   for (j=n-1;j>=0;j--) {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[k][j] * x[k];
      x[j] = (a[n][j] - tmp) / a[j][j];
   }

}




