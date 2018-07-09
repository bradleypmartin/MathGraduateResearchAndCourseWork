#include "MG_Structure.h"
#include "MG_Level.h"
#include "CSR_Vector.h"
#include "CSR_Matrix.h"

#include <stdlib.h>	//malloc
#include <iostream>	//cout
#include <fstream>  //ofstream

using namespace std;

//Default Constructor
//Input:
//	
//Output:
//
//Description:
//	Create an empty Structure
MG_Structure::MG_Structure()
{
	m_numLevels = 0;
	m_levels = 0;
}

//Default Destructor
//Input:
//	
//Output:
//
//Description:
//	Cleanup memory
MG_Structure::~MG_Structure()
{

	m_numLevels = 0;

	for(int k = 0; k < m_numLevels; k++)
	{
		delete m_levels[k];
	}
	
}

//Setup 1D
//Input:
//	NumPoints - the number of points on the finest level
//	NumLevels - the number of levels you want in your structure
//	A - The start point
//  B - The end point
//Output:
//
//Description:
//	1) Allocate memory for the level pointers
//	2) for each level
//		a) Instantiate a new level
//		b) Call level setup with appropriate parameters
void MG_Structure::Setup1D( int NumPoints, int NumLevels, double A, double B )
{
	//MG_Level** m_levels;
	m_levels = (MG_Level**)malloc(sizeof(MG_Level*)*(NumLevels));
	m_numLevels = NumLevels;

	for (int k = 0; k < NumLevels; k++)
	{
		int levelnumpoints = (NumPoints+1)/pow(2,k)-1;
		bool iscoarsest = false;

		if (k == (NumLevels - 1))
		{
			iscoarsest = true;
		}
		
		m_levels[k] = new MG_Level(levelnumpoints);

		m_levels[k]->Setup1D(k,levelnumpoints,A,B,iscoarsest);
	}
}



//Setup2D
//Input:
// NumPoints1Dir: number of points in 1 direction on FINEST 2D grid (square root of matrix size)
// NumLevels: number of grid levels in MG structure
// A: lower bound for x and y (square domain)
// B: upper bound for x and y
// Omega: weighting factor to be used in WJ relaxation
// Sigma: factor of "u" present from 1.4 model problem (LHS)

//Output:
//
//Description:
//  Allocates memory for and sets up an MG structure (restriction, prolongation, Laplacian operators, etc)
//  prior to performing MG runs for PA 4

void MG_Structure::Setup2D( int NumPoints1Dir, int NumLevels, double A, double B, double Omega, double Sigma)
{
	m_levels = (MG_Level**)malloc(sizeof(MG_Level*)*(NumLevels));
	m_numLevels = NumLevels;

	for (int k = 0; k < NumLevels; k++)
	{
		int levelnumpoints = (NumPoints1Dir+1)/pow(2,k)-1;
		bool iscoarsest = false;

		if (k == (NumLevels - 1))
		{
			iscoarsest = true;
		}
		
		m_levels[k] = new MG_Level(pow(levelnumpoints,2));

		m_levels[k]->Setup2D(k,levelnumpoints,A,B,Omega,Sigma,iscoarsest);
	}
}


//Setup2DVAR
//Input:
// NumPoints1Dir: number of points in 1 direction on FINEST 2D grid (square root of matrix size)
// NumLevels: number of grid levels in MG structure
// A: lower bound for x and y (square domain)
// B: upper bound for x and y
// Omega: weighting factor to be used in WJ relaxation
// Sigma: factor of "u" present from 1.4 model problem (LHS)

//Output:
//
//Description:
//  Allocates memory for and sets up an MG structure (VARIATIONAL restriction, prolongation, Laplacian operators, etc)
//  prior to performing MG runs for PA 4

void MG_Structure::Setup2DVAR( int NumPoints1Dir, int NumLevels, double A, double B, double Omega, double Sigma)
{
	m_levels = (MG_Level**)malloc(sizeof(MG_Level*)*(NumLevels));
	m_numLevels = NumLevels;

	for (int k = 0; k < NumLevels; k++)
	{
		int levelnumpoints = (NumPoints1Dir+1)/pow(2,k)-1;
		bool iscoarsest = false;

		if (k == (NumLevels - 1))
		{
			iscoarsest = true;
		}
		
		m_levels[k] = new MG_Level(pow(levelnumpoints,2));

		m_levels[k]->Setup2DVAR(k,levelnumpoints,A,B,Omega,Sigma,iscoarsest);
	}
}

//VCycle2D
//Input:
// Iterationsdown: number of relaxations to perform on the "downward" sweep of the V-cycle
// Iterationsup: same as above, on the "upward" sweep
// Omega: weighting factor to be used in WJ relaxation
// Sigma: factor of "u" present from 1.4 model problem (LHS)

//Output:
//
//Description:
//  Performs a V-cycle non-recursively for the 2D model problem in PA 4


void MG_Structure::VCycle2D(int Iterationsdown, int Iterationsup, double Omega, double Sigma)
{
	
	for (int k = 0; k < m_numLevels-1; k++)
	{
		if (k > 0)
		{
			m_levels[k]->SetUZero();
			m_levels[k-1]->GetRestrict()->ApplyToVector(m_levels[k-1]->GetResidual(),m_levels[k]->GetF());
		}

		m_levels[k]->WJRelax2D(Iterationsdown, Omega, Sigma);
	}

	m_levels[m_numLevels-1]->SetUZero();
	m_levels[m_numLevels-2]->GetRestrict()->ApplyToVector(m_levels[m_numLevels-2]->GetResidual(),m_levels[m_numLevels-1]->GetF());
	
	DirectSolve(m_levels[m_numLevels-1]->GetA(),m_levels[m_numLevels-1]->GetU(),m_levels[m_numLevels-1]->GetF());

	for (int k = m_numLevels-2; k >= 0; k--)
	{
		m_levels[k+1]->GetInterpolation()->ApplyToVector(m_levels[k+1]->GetU(),m_levels[k]->GetCorrection());
		m_levels[k]->GetU()->Add(m_levels[k]->GetCorrection(),1.0);

		m_levels[k]->WJRelax2D(Iterationsup, Omega, Sigma);
	}
}


void MG_Structure::VCycle2DVAR(int Iterationsdown, int Iterationsup, double Omega, double Sigma)
{
	
	for (int k = 0; k < m_numLevels-1; k++)
	{
		if (k > 0)
		{
			m_levels[k]->SetUZero();
			m_levels[k-1]->GetRestrict()->ApplyToVector(m_levels[k-1]->GetResidual(),m_levels[k]->GetF());
		}

		m_levels[k]->WJRelax2DVAR(Iterationsdown, Omega, Sigma);
	}

	m_levels[m_numLevels-1]->SetUZero();
	m_levels[m_numLevels-2]->GetRestrict()->ApplyToVector(m_levels[m_numLevels-2]->GetResidual(),m_levels[m_numLevels-1]->GetF());
	
	DirectSolve(m_levels[m_numLevels-1]->GetA(),m_levels[m_numLevels-1]->GetU(),m_levels[m_numLevels-1]->GetF());

	for (int k = m_numLevels-2; k >= 0; k--)
	{
		m_levels[k+1]->GetInterpolation()->ApplyToVector(m_levels[k+1]->GetU(),m_levels[k]->GetCorrection());
		m_levels[k]->GetU()->Add(m_levels[k]->GetCorrection(),1.0);

		m_levels[k]->WJRelax2DVAR(Iterationsup, Omega, Sigma);
	}
}
//Direct Solve
//Input:
// A: A matrix in CSR form (this'll be the Laplacian 2nd order 2D operator for PA4)
// u: LHS vector (discrete solution we're looking for)
// f: RHS of the discrete Poisson problem on some grid level

//Output:
//
//Description:
//  Calls exact solve for the coarsest level.  I guess this function isn't really needed
//  (I could call ExactSolve directly, in-line), but oh well.

void MG_Structure::DirectSolve(CSR_Matrix* A, CSR_Vector* u, CSR_Vector* f)
{
	A->ExactSolve(u,f);
}

//RecursiveVCycle2D
//Input:
// StartLevel: Level of the MG structure on which to start the V-Cycle (useful for FMG implementation)
// CurrentLevel: keeps track of the level the V-Cycle is at through the process
// Iterationsdown: number of relaxations to perform on the "downward" sweep of the V-cycle
// Iterationsup: same as above, on the "upward" sweep
// Omega: weighting factor to be used in WJ relaxation
// Sigma: factor of "u" present from 1.4 model problem (LHS)

//Output:
//
//Description:
//  Performs a V-Cycle recursively - to be used in both normal v-cycling AND FMG implementation (PA4)


void MG_Structure::RecursiveVCycle2D(int StartLevel, int CurrentLevel, int Iterationsdown, int Iterationsup, double Omega, double Sigma)
{
	if ((m_levels[CurrentLevel]->IsCoarsestLevel()) == true)
	{
		m_levels[m_numLevels-1]->SetUZero();
		m_levels[m_numLevels-2]->GetRestrict()->ApplyToVector(m_levels[m_numLevels-2]->GetResidual(),m_levels[m_numLevels-1]->GetF());
	
		DirectSolve(m_levels[m_numLevels-1]->GetA(),m_levels[m_numLevels-1]->GetU(),m_levels[m_numLevels-1]->GetF());
	}
	else
	{

		if (CurrentLevel > StartLevel)
		{
			m_levels[CurrentLevel]->SetUZero();
			m_levels[CurrentLevel-1]->GetRestrict()->ApplyToVector(m_levels[CurrentLevel-1]->GetResidual(),m_levels[CurrentLevel]->GetF());
		}

		m_levels[CurrentLevel]->WJRelax2D(Iterationsdown, Omega, Sigma);

		RecursiveVCycle2D(StartLevel, CurrentLevel+1, Iterationsdown, Iterationsup, Omega, Sigma);
		
		m_levels[CurrentLevel+1]->GetInterpolation()->ApplyToVector(m_levels[CurrentLevel+1]->GetU(),m_levels[CurrentLevel]->GetCorrection());
		m_levels[CurrentLevel]->GetU()->Add(m_levels[CurrentLevel]->GetCorrection(),1.0);

		m_levels[CurrentLevel]->WJRelax2D(Iterationsup, Omega, Sigma);;
	}
}


//FMG
//Input:
//
// Iterationsdown: number of relaxations to perform on the "downward" sweep of the V-cycles
// Iterationsup: same as above, on the "upward" sweeps
// Omega: weighting factor to be used in WJ relaxation
// Sigma: factor of "u" present from 1.4 model problem (LHS)

//Output:
//
//Description:
//  Carries out FMG for PA4
//  NOTE: THIS IS ESPECIALLY SET UP FOR THE SPECIFIC TRIAL PROBLEM -Laplace(u) + u = sin(pi*x)*sin(pi*y).
//  Future implementation might accept an array of RHS's for each grid level.

void MG_Structure::FMG(int Iterationsdown, int Iterationsup, double Omega, double Sigma)
{

	int totalpoints = m_levels[m_numLevels-1]->GetSize();
	CSR_Vector* RHS = new CSR_Vector(totalpoints);
	int NumPoints1Dir = pow(totalpoints,0.5);

	const double PI = 4.0*atan(1.0);
	int Index;

	// Setting the RHS for the finest grid level

	for (int m = 0; m < NumPoints1Dir; m++)
	{
		for (int n = 0; n < NumPoints1Dir; n++)
		{
			Index = m*NumPoints1Dir + n;

			RHS->SetValue(Index,sin(PI*(m+1)/(NumPoints1Dir+1))*
				sin(PI*(n+1)/(NumPoints1Dir+1)));
		}

	}
	
	m_levels[m_numLevels-1]->GetF()->Set(RHS);

	delete RHS;

	DirectSolve(m_levels[m_numLevels-1]->GetA(),m_levels[m_numLevels-1]->GetU(),m_levels[m_numLevels-1]->GetF());
	
	for (int Level = m_numLevels-1; Level > 0; Level--)
	{
		
		// Interpolating SOLUTION to next (finer) grid level

		m_levels[Level]->GetInterpolation()->ApplyToVector(m_levels[Level]->GetU(),m_levels[Level-1]->GetU());

		int totalpoints = m_levels[Level-1]->GetSize();
		CSR_Vector* RHS = new CSR_Vector(totalpoints);
		int NumPoints1Dir = pow(totalpoints,0.5);

		// Setting up RHS on next (finer) grid level

		for (int m = 0; m < NumPoints1Dir; m++)
		{
			for (int n = 0; n < NumPoints1Dir; n++)
			{
				Index = m*NumPoints1Dir + n;

				RHS->SetValue(Index,sin(PI*(m+1)/(NumPoints1Dir+1))*
					sin(PI*(n+1)/(NumPoints1Dir+1)));
			}
		}

		m_levels[Level-1]->GetF()->Set(RHS);
		
		delete RHS;

		// Performing one VCycle on appropriate grid level down to the coarsest grid

		RecursiveVCycle2D(Level-1, Level-1, Iterationsdown, Iterationsup, Omega, Sigma);

	}

}



//Get num Levels
//Input:
//	
//Output:
//	Number of levels in this structure
//Description:
//	
int MG_Structure::GetNumLevels()
{
	return m_numLevels;
}

//Get Level
//Input:
//	
//Output:
//	Return a pointer to the request level
//Description:
//	better make sure that it's a valid request!
MG_Level* MG_Structure::GetLevel( int Level )
{
	return m_levels[Level];
}
	
