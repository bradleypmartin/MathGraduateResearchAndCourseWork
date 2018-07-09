/*
*	Sample Files created for APPM 6640 - Multigrid Methods
*
*	Written by Toby Jones and Chris Leibs
*	January 4th, 2013
*/

#ifndef __CSR_MATRIX_H__
#define __CSR_MATRIX_H__

#include "CSR_Vector.h"

class CSR_Matrix
{
public:
	
	//Constructor / Destructor
	CSR_Matrix();
	~CSR_Matrix();
public:

	void SetToDiagonal( int Size, double D );
	void SetTridiagonal( int Size,  double L, double D, double U );

	void Set2DCrossstencil(int SideSize, double Outer, double Inner);
	void Set2DVARCrossstencil(int SideSize, double Outer, double Inner);

	void Set2DFWrestrict(int SideSize);
	void Set2Dinterpolate(int SideSize);

	void Debug_Print();	
	void Matlab_Print( char* FileName  );
	void Write( char* FileName );
	void Load( char* FileName );

	void DisplaySize();
		

	void ApplyToVector( CSR_Vector* In, CSR_Vector* Out );
	
	void Scale( double ScaleAmount );

	void ExactSolve( CSR_Vector* u, CSR_Vector* f );	

	double Get (int Row, int Col);	
	
protected:	
	void ClearData();
	void GSolve(double **a,int n,double *x);
private:

	double* m_rowData;
	int* m_colIndices;
	int* m_rowStartIndices;
	int m_dataSize;
	int m_numRows;
	int m_numCols;			
};


#endif ///__CSR_MATRIX_H__
