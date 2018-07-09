/*
*	Sample Files created for APPM 6640 - Multigrid Methods
*
*	Written by Toby Jones and Chris Leibs
*	January 4th, 2013
*/

#ifndef __CSR_VECTOR__
#define __CSR_VECTOR__

class CSR_Vector
{
public:
	CSR_Vector( );
	CSR_Vector( unsigned int Size );
	~CSR_Vector();

public:
	unsigned int Length();
	
	void Set(CSR_Vector* u );
	void SetZero();
	void SetConstant( double Value );
	void SetRandom();
	double GetValue( unsigned int Index );
	void SetValue( unsigned int Index, double Value );
	
	void Add( CSR_Vector* Vec, double Factor );
	
	void Matlab_Print( char* FileName );
	void Debug_Print();
	void Load( char* FileName );
	void Write( char* FileName );

	double* GetData();

private:
	void FreeMemory();
private:
	unsigned int m_size;
	double* m_data;
};




#endif ///__CSR_Vector__
