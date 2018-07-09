#include "MG_Level.h"
#include "CSR_Vector.h"
#include "CSR_Matrix.h"
#include <math.h>
#include <iostream>
using namespace std;

//Default Constructor
//Input:
//	
//Output:
//
//Description:
//	
MG_Level::MG_Level()
{
	//CSR_Matrix* m_A = new CSR_Matrix();
	
	//CSR_Vector* m_f = new CSR_Vector();
	//CSR_Vector* m_u = new CSR_Vector();

	m_A = new CSR_Matrix();
	m_Romega = new CSR_Matrix();
	m_Restrict = new CSR_Matrix();
	m_Interpolate = new CSR_Matrix();

	m_f = new CSR_Vector();
	m_u = new CSR_Vector();
	m_u_holder = new CSR_Vector();
	m_residual = new CSR_Vector();
	m_correction = new CSR_Vector();

	m_isCoarsestLevel = false;
	m_level = 0;

	m_h = 1.0;
}

MG_Level::MG_Level(int Size)
{
	//CSR_Matrix* m_A = new CSR_Matrix();
	
	//CSR_Vector* m_f = new CSR_Vector(Size);
	//CSR_Vector* m_u = new CSR_Vector(Size);

	m_A = new CSR_Matrix();
	m_Romega = new CSR_Matrix();
	m_Restrict = new CSR_Matrix();
	m_Interpolate = new CSR_Matrix();

	m_f = new CSR_Vector(Size);
	m_u = new CSR_Vector(Size);
	m_u_holder = new CSR_Vector(Size);
	m_residual = new CSR_Vector(Size);
	m_correction = new CSR_Vector(Size);

	m_isCoarsestLevel = false;
	m_level = 0;

	m_h = 1.0;
}

//Default Destructor
//Input:
//	
//Output:
//
//Description:
//
MG_Level::~MG_Level()
{
	delete m_A;
	delete m_Romega;
	delete m_Restrict;
	delete m_Interpolate;

	delete m_f;
	delete m_u;
	delete m_u_holder;
	delete m_residual;
	delete m_correction;

	m_isCoarsestLevel = false;
	m_level = 0;

	m_h = 1.0;
}


//Setup
//Input:
//	Level : The index of this level in the structure
//	NumPoints : The number of points for this level ( NOTE!!! this is the number of interior points )
//	A : The left side of this interval
//	B : The right side of this interval
//	IsCoarsest : True or false
//	(ps with comments - can be done in 9 lines of code )
//Output:
//
//Description:
//
void MG_Level::Setup1D( int Level, int NumPoints, double A, double B, bool IsCoarsest  )
{
	m_level = Level;
	m_h = (B-A)/(NumPoints+1.0);
	m_isCoarsestLevel = IsCoarsest;

	//CSR_Matrix* m_A = new CSR_Matrix();
	m_A -> SetTridiagonal(NumPoints, -1.0/pow(m_h,2.0), 2.0/pow(m_h,2.0), -1.0/pow(m_h,2));

	//CSR_Vector* m_u = new CSR_Vector(NumPoints);
	m_u -> SetZero();

	//CSR_Vector* m_f = new CSR_Vector(NumPoints);
	m_f -> SetZero();

}

//Setup2D
//
//Input:
//	Level : index of the MG level for the MG_Structure to reference
//  NumPoints1Dir: number of points on ONE SIDE of the grid (square root of matrix size)
//  A: lower bound for both x and y (square domain)
//  B: upper bound for both x and y
//  Omega: weighting factor to be used in WJ relaxation
//  Sigma: multiple of "u" present in 1.4 model problem (LHS)
//  IsCoarsest: flag for coarsest level of the MG structure
//	
//Output:
//
//Description:
//	sets up one level of the MG structure (for programming assignment 4)

void MG_Level::Setup2D( int Level, int NumPoints1Dir, double A, double B, double Omega, double Sigma, bool IsCoarsest)
{
	m_level = Level;
	m_h = (B-A)/(NumPoints1Dir+1.0);
	m_isCoarsestLevel = IsCoarsest;

	m_A -> Set2DCrossstencil(NumPoints1Dir,-1.0/pow(m_h,2.0),(4.0+Sigma*pow(m_h,2))/pow(m_h,2));
	m_Romega -> Set2DCrossstencil(NumPoints1Dir,Omega/(4.0+Sigma*pow(m_h,2)),1.0-Omega);

	m_Interpolate -> Set2Dinterpolate(NumPoints1Dir);
	m_Restrict -> Set2DFWrestrict(NumPoints1Dir);

	m_u -> SetZero();
	m_f -> SetZero();
	m_u_holder -> SetZero();
	m_residual -> SetZero();
	m_correction -> SetZero();
}


//Setup2DVAR
//
//Input:
//	Level : index of the MG level for the MG_Structure to reference
//  NumPoints1Dir: number of points on ONE SIDE of the grid (square root of matrix size)
//  A: lower bound for both x and y (square domain)
//  B: upper bound for both x and y
//  Omega: weighting factor to be used in WJ relaxation
//  Sigma: multiple of "u" present in 1.4 model problem (LHS)
//  IsCoarsest: flag for coarsest level of the MG structure
//	
//Output:
//
//Description:
//	sets up one level (with variational operators) of the MG structure (for programming assignment 4)

void MG_Level::Setup2DVAR( int Level, int NumPoints1Dir, double A, double B, double Omega, double Sigma, bool IsCoarsest)
{
	m_level = Level;
	m_h = (B-A)/(NumPoints1Dir+1.0);
	m_isCoarsestLevel = IsCoarsest;

	m_A -> Set2DVARCrossstencil(NumPoints1Dir,(-1.0/3.0)/pow(m_h,2.0),((8.0/3.0)+Sigma*pow(m_h,2.0))/pow(m_h,2.0));
	m_Romega -> Set2DVARCrossstencil(NumPoints1Dir,Omega/((8.0)+3.0*Sigma*pow(m_h,2.0)),1.0-Omega);

	m_Interpolate -> Set2Dinterpolate(NumPoints1Dir);
	m_Restrict -> Set2DFWrestrict(NumPoints1Dir);

	m_u -> SetZero();
	m_f -> SetZero();
	m_u_holder -> SetZero();
	m_residual -> SetZero();
	m_correction -> SetZero();
}

//WJRelax2D
//
//Input:
//	Iterations: number of times to relax
//  Sigma: multiple of "u" present in 1.4 model problem (LHS)
//  IsCoarsest: flag for coarsest level of the MG structure
//	
//Output:
//
//Description:
//	Performs 2D WJ relaxation (for programming assignment 4)

void MG_Level::WJRelax2D(int Iterations, double Omega, double Sigma)
{
	for (int k = 0; k < Iterations; k++)
	{
		m_u_holder -> Set(m_u);
		m_Romega -> ApplyToVector(m_u_holder,m_u);
		m_u -> Add(m_f,Omega*(pow(m_h,2.0)/(4.0+Sigma*pow(m_h,2.0))));
	}

	m_residual->Set(m_f);
	m_A->ApplyToVector(m_u,m_correction);
	m_residual->Add(m_correction,-1.0);
}

void MG_Level::WJRelax2DVAR(int Iterations, double Omega, double Sigma)
{
	for (int k = 0; k < Iterations; k++)
	{
		m_u_holder -> Set(m_u);
		m_Romega -> ApplyToVector(m_u_holder,m_u);
		m_u -> Add(m_f,Omega*(pow(m_h,2.0)/(8.0/3.0+Sigma*pow(m_h,2.0))));
	}

	m_residual->Set(m_f);
	m_A->ApplyToVector(m_u,m_correction);
	m_residual->Add(m_correction,-1.0);
}

//Get A
//Input:
//	
//Output:
//
//Description:
//
CSR_Matrix* MG_Level::GetA()
{
	return m_A;
}
	
//Get Romega
//Input:
//	
//Output:
//
//Description:
//

CSR_Matrix* MG_Level::GetRomega()
{
	return m_Romega;
}

//Get Restriction
//Input:
//	
//Output:
//
//Description:
//

CSR_Matrix* MG_Level::GetRestrict()
{
	return m_Restrict;
}

//Get Interpolation
//Input:
//	
//Output:
//
//Description:
//

CSR_Matrix* MG_Level::GetInterpolation()
{
	return m_Interpolate;
}


//Get F
//Input:
//	
//Output:
//
//Description:
//

CSR_Vector* MG_Level::GetF()
{
	return m_f;
}

//Get U
//Input:
//	
//Output:
//
//Description:
//
CSR_Vector* MG_Level::GetU()
{
	return m_u;
}

//Get Residual
//Input:
//	
//Output:
//
//Description:
//

CSR_Vector* MG_Level::GetResidual()
{
	return m_residual;
}

//Get Correction
//Input:
//	
//Output:
//
//Description:
//

CSR_Vector* MG_Level::GetCorrection()
{
	return m_correction;
}


//Set U
//Input:
//	
//Output:
//
//Description:
//
void MG_Level::SetU( CSR_Vector* U )
{
	m_u->Set(U);
}

//Set U to zero
//Input:
//	
//Output:
//
//Description:
//
void MG_Level::SetUZero()
{
	m_u->SetZero();
}

//Set U to constant
//Input:
//	Constant
//Output:
//
//Description:
//
void MG_Level::SetUConstant( double Const )
{
	m_u->SetConstant( Const );
}

//Set U to Random
//Input:
//	Constant
//Output:
//
//Description:
//
void MG_Level::SetURandom()
{
	m_u->SetRandom();
}
	
//Set F
//Input:
//	
//Output:
//
//Description:
//
void MG_Level::SetF( CSR_Vector* F )
{
	m_f->Set( F );
}

//Set F to zero
//Input:
//	
//Output:
//
//Description:
//
void MG_Level::SetFZero()
{
	m_f->SetZero();
}
	
//Is Coarseset
//Input:
//	
//Output:
//	is coarserest ?
//Description:
//
bool MG_Level::IsCoarsestLevel()
{
	return m_isCoarsestLevel;
}

//Get Size
//Input:
//	
//Output:
//	Number of points for this level
//Description:
//
int MG_Level::GetSize()
{
  return m_u->Length();
}
