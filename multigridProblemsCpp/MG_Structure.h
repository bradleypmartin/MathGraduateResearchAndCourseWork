#ifndef __MG_STRUCTURE_H__
#define __MG_STRUCTURE_H__

#include "MG_Level.h"
#include "CSR_Matrix.h"

class MG_Structure
{
public:
	MG_Structure();
	~MG_Structure();

	void Setup1D( int NumPoints, int NumLevels, double A, double B );
	void Setup2D( int NumPoints1Dir, int NumLevels, double A, double B, double Omega, double Sigma);
	
	void Setup2DVAR( int NumPoints1Dir, int NumLevels, double A, double B, double Omega, double Sigma);
	
	void VCycle2D(int Iterationsdown, int Iterationsup, double Omega, double Sigma);
	void VCycle2DVAR(int Iterationsdown, int Iterationsup, double Omega, double Sigma);

	void RecursiveVCycle2D(int StartLevel, int CurrentLevel, int Iterationsdown, int Iterationsup, double Omega, double Sigma);
	
	void FMG(int Iterationsdown, int Iterationsup, double Omega, double Sigma);
	
	int GetNumLevels();
	MG_Level* GetLevel( int Level );
	
private:

	void DirectSolve(CSR_Matrix* A,CSR_Vector* u, CSR_Vector* f);

private:
	int m_numLevels;
	MG_Level** m_levels;		
};

#endif //__MG_STRUCTURE_H__