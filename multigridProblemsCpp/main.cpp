#include "MG_Structure.h"
#include "MG_Level.h"
#include "CSR_Matrix.h"
#include "CSR_Vector.h"

#include <string>	//sprintf
#include <math.h>	//pow
#include <iostream>	//cout
#include <fstream>  //ofstream

using namespace std;


int main(int argc, char** argv )
{

	// Setting up run (used in both V-cycling and FMG implementation

	int maxpow = 7;   // power of 2 (-1) points in one direction on 2D grid

	int NumPoints1Dir = pow(2,maxpow)-1;

	int TotalFinePoints = pow(NumPoints1Dir,2);    // Total number of points on FINEST 2D grid

	int NumLevels = maxpow-2;  // going down to 7 x 7 grid (at the coarsest)

	double A = 0.0;      // Unit square
	double B = 1.0;

	double Omega = 4.0/5.0;  // Ideal for 2D model problem
	double Sigma = 1.0;      // 1.0 picked for a trial problem because 1.0 is awesome

	int VCycles = 5;    // not going to need this many, but what the hell - I'm building tables, yo!

	int IterationsDown = 2;  // number of iterations down and up Vcycles
	int IterationsUp = 2;

	const double PI = 4.0*atan(1.0);

	CSR_Vector* CurrentGuess = new CSR_Vector(TotalFinePoints);      // holds self-explanatory current-guess vectors
	CSR_Vector* CurrentError = new CSR_Vector(TotalFinePoints);
	CSR_Vector* CurrentResidual = new CSR_Vector(TotalFinePoints);
	CSR_Vector* TrueSolution = new CSR_Vector(TotalFinePoints);
	CSR_Vector* RHS = new CSR_Vector(TotalFinePoints);

	CSR_Vector* ErrorNormStore = new CSR_Vector(VCycles);      // These will carry the h-norm of error and residual
	CSR_Vector* ResidualNormStore = new CSR_Vector(VCycles);   // after each v-cycle (or, alternatively, after the FMG implementation)

	ErrorNormStore->SetZero();
	ResidualNormStore->SetZero();

	MG_Structure* MyMultigrid = new MG_Structure();
	
	// TOGGLE VARIATIONAL SETUP HERE >>> AND BELOW, TOO, FOR VCYCLING!

	MyMultigrid->Setup2D(NumPoints1Dir, NumLevels, A, B, Omega, Sigma);
	//MyMultigrid->Setup2DVAR(NumPoints1Dir, NumLevels, A, B, Omega, Sigma);

	int Index;

	double Factor = 1.0/(2*pow(PI,2.0)+Sigma);   // Factor to multiply true solution sin functions by below (when defining true (analytic) solution)


	

	//// IMPLEMENTATION OF VCYCLING ////


	// Setting up finest grid RHS

	for (int m = 0; m < NumPoints1Dir; m++)
	{

		for (int n = 0; n < NumPoints1Dir; n++)
		{
			Index = m*NumPoints1Dir + n;

			TrueSolution->SetValue(Index,
				Factor*sin(PI*(m+1)/(NumPoints1Dir+1))*
				sin(PI*(n+1)/(NumPoints1Dir+1)));

			RHS->SetValue(Index,sin(PI*(m+1)/(NumPoints1Dir+1))*
				sin(PI*(n+1)/(NumPoints1Dir+1)));
		}

	}

	MyMultigrid->GetLevel(0)->GetF()->Set(RHS);

	// Performing V-cycles

	for (int k = 0; k < VCycles; k++)
	{

		// ALSO TOGGLE VARIATIONAL VCYCLING HERE!!!

		//MyMultigrid->VCycle2DVAR(IterationsDown,IterationsUp,Omega,Sigma);
		MyMultigrid->VCycle2D(IterationsDown,IterationsUp,Omega,Sigma);

		CurrentGuess->Set((MyMultigrid->GetLevel(0))->GetU());
		CurrentResidual->Set((MyMultigrid->GetLevel(0))->GetResidual());
		CurrentError->Set(TrueSolution);
		CurrentError->Add(CurrentGuess,-1.0);

		// computing h-norm of residual and error after each v-cycle

		double sumerror2 = 0.0;
		double sumresid2 = 0.0;

		for (int z = 0; z < TotalFinePoints; z++)
		{
			sumerror2 = sumerror2 + (CurrentError->GetValue(z))*(CurrentError->GetValue(z));
			sumresid2 = sumresid2 + (CurrentResidual->GetValue(z))*(CurrentResidual->GetValue(z));
		}

		ErrorNormStore->SetValue(k,pow(sumerror2,0.5)*(1.0/(NumPoints1Dir+1)));
		ResidualNormStore->SetValue(k,pow(sumresid2,0.5)*(1.0/(NumPoints1Dir+1)));
		
	}

	
	//// END VCYCLING IMPLEMENTATION ////

	


	/*

	//// IMPLEMENTATION OF FMG ////

	// Calling FMG (NOTE: RHS's for each grid level are currently built
	// WITHIN the FMG function.  I know this is suboptimal, but I was
	// too lazy at this point to come up with a clever way to pass in
	// each grid level's RHS when it's the new fine grid in FMG.)

	MyMultigrid->FMG(IterationsDown,IterationsUp,Omega,Sigma);

	// setting true (analytic) solution

	for (int m = 0; m < NumPoints1Dir; m++)
	{
		for (int n = 0; n < NumPoints1Dir; n++)
		{
			Index = m*NumPoints1Dir + n;

			TrueSolution->SetValue(Index,
				Factor*sin(PI*(m+1)/(NumPoints1Dir+1))*
				sin(PI*(n+1)/(NumPoints1Dir+1)));
		}
	}

	CurrentGuess->Set((MyMultigrid->GetLevel(0))->GetU());
	CurrentResidual->Set((MyMultigrid->GetLevel(0))->GetResidual());
	CurrentError->Set(TrueSolution);
	CurrentError->Add(CurrentGuess,-1.0);

	// calculating h-norm of error and residual after FMG

	double sumerror2 = 0.0;
	double sumresid2 = 0.0;

	for (int z = 0; z < TotalFinePoints; z++)
	{
		sumerror2 = sumerror2 + (CurrentError->GetValue(z))*(CurrentError->GetValue(z));
		sumresid2 = sumresid2 + (CurrentResidual->GetValue(z))*(CurrentResidual->GetValue(z));
	}

	ErrorNormStore->SetValue(0,pow(sumerror2,0.5)*(1.0/(NumPoints1Dir+1)));
	ResidualNormStore->SetValue(0,pow(sumresid2,0.5)*(1.0/(NumPoints1Dir+1)));

	//// END IMPLEMENTATION OF FMG ////

	*/

	// Printing files for visualization and writeup

	ErrorNormStore->Matlab_Print("errornorms.txt");
	ResidualNormStore->Matlab_Print("residualnorms.txt");
	
	CurrentGuess->Matlab_Print("currentguess.txt");
	CurrentResidual->Matlab_Print("currentresidual.txt");

	// Cleaning house!

	delete CurrentGuess;
	delete CurrentError;
	delete CurrentResidual;
	delete TrueSolution;
	delete RHS;

	delete ErrorNormStore;
	delete ResidualNormStore;

	delete MyMultigrid;

	

	return 0;
}
