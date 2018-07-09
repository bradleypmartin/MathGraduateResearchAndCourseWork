#ifndef __MG_LEVEL_H__
#define __MG_LEVEL_H__

#include "CSR_Matrix.h"


class MG_Level
{
public:
	MG_Level();
	MG_Level(int Size);
	~MG_Level();

	//Setup
	void Setup1D( int Level, int NumPoints, double A, double B, bool IsCoarsest  );
	void Setup2D( int Level, int NumPoints1Dir, double A, double B, double Omega, double Sigma, bool IsCoarsest);
	void Setup2DVAR( int Level, int NumPoints1Dir, double A, double B, double Omega, double Sigma, bool IsCoarsest);
	
	CSR_Matrix* GetA();
	CSR_Matrix* GetRomega();
	CSR_Matrix* GetRestrict();
	CSR_Matrix* GetInterpolation();

	
	CSR_Vector* GetF();
	CSR_Vector* GetU();
	CSR_Vector* GetResidual();
	CSR_Vector* GetCorrection();

	void SetU( CSR_Vector* U );
	void SetUZero();
	void SetUConstant( double Const );
	void SetURandom();
	
	void SetF( CSR_Vector* F );
	void SetFZero();

	void WJRelax2D(int Iterations, double Omega, double Sigma);
	void WJRelax2DVAR(int Iterations, double Omega, double Sigma);
	
	//info about this level
	bool IsCoarsestLevel();
	int GetSize();

	
//helper functions
private:
	
private:
	CSR_Matrix* m_A;
	CSR_Matrix* m_Romega;
	CSR_Matrix* m_Restrict;
	CSR_Matrix* m_Interpolate;
	
	CSR_Vector* m_f;
	CSR_Vector* m_u;
	CSR_Vector* m_u_holder;
	CSR_Vector* m_residual;
	CSR_Vector* m_correction;
	
	bool m_isCoarsestLevel;
	int m_level;

	double m_h;
};

#endif //__MG_LEVEL_1D_H__
