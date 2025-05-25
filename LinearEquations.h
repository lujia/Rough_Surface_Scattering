// LinearEquations.h: interface for the CLinearEquations class.
//
//////////////////////////////////////////////////////////////////////
#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>
#include<vector>
#include "math.h"

using namespace std;


typedef complex <double>  Complex;
typedef vector  <double>  Vector_double;
typedef vector  <Complex> Vector_complex;

class CLinearEquations  
{
	double PI;
	Complex J;
public:
	CLinearEquations();
	virtual ~CLinearEquations();
    ///////////////////////////////////////////////////
	/////  基本的矩阵运算及向量运算（矩阵表达）  //////
	///////////////////////////////////////////////////
	void CLinearEquations::MultiplyMatrixMatrix(Complex **A,Complex **B,Complex **result,int N,int L,int M);//A(N*L),B(L*M),result(N*M)
	void CLinearEquations::MultiplyMatrixVevtor(Complex **A,Complex *vec,Complex *result,int N,int L);//A(N*L),vec(L*1),result(N*1)
	void CLinearEquations::MatrixConjTurn(Complex **A,Complex **result,int N,int L);//A(N*L),result(L*N)
	double CLinearEquations::Norm1(Complex *vec,int N);
	double CLinearEquations::Norm2(Complex *vec,int N);
	///////////////////////////////////////////////////
	/////           基本的运算向量表达            /////
	///////////////////////////////////////////////////
	void CLinearEquations::VectorMultiplyMatrixMatrix(Vector_complex &A,Vector_complex &B,Vector_complex &result,int N,int L,int M);//A(N*L),B(L*M),result(N*M)
	void CLinearEquations::VectorMultiplyMatrixVector(Vector_complex &A,Vector_complex &vec,Vector_complex &result,int N,int L);//A(N*L),vec(N*1),result(N*1)
	void CLinearEquations::VectorMatrixConjTurn(Vector_complex &A,Vector_complex &result,int N,int L);//A(N*L),result(L*N)
    double CLinearEquations::VectorNorm1(Vector_complex &vec,int N);
	double CLinearEquations::VectorNorm2(Vector_complex &vec,int N);
	double CLinearEquations::VectorNorm(Vector_complex &vec,int N);
	void CLinearEquations::VectorMultiplyVector(Vector_complex &vec1,Vector_complex &vec2,Vector_complex &result,int N);//vec1(N),vec2(N),result(N)
	///////////////////////////////////////////////////
	/////   支持实数表达的基本矩阵、向量运算    ///////
	///////////////////////////////////////////////////
	void CLinearEquations::VecMultMatMat(Vector_double &A,Vector_double &B,Vector_double &result,int N,int L,int M);//A(N*L),B(L*M),result(N*M)
	void CLinearEquations::VecMultMatVec(Vector_double &A,Vector_double &vec,Vector_double &result,int N,int L);//A(N*L),vec(N*1),result(N*1)
	void CLinearEquations::VecMatConjTurn(Vector_double &A,Vector_double &result,int N,int L);//A(N*L),result(L*N)
    double CLinearEquations::VecNorm1(Vector_double &vec,int N);
	double CLinearEquations::VecNorm2(Vector_double &vec,int N);
	double CLinearEquations::VecNorm(Vector_double &vec,int N);
	void CLinearEquations::VecMultVec(Vector_double &vec1,Vector_double &vec2,Vector_double &result,int N);//vec1(N),vec2(N),result(N)	
	double CLinearEquations::VecSumVecMultVec(Vector_double &vec1,Vector_double &vec2,int N);
	///////////////////////////////////////////////////
	////               矩阵分解                    ////
	///////////////////////////////////////////////////
    int CLinearEquations::LU(Vector_complex &a,Vector_complex &l,Vector_complex &u,int n);
	int CLinearEquations::QR(double *a,int m,int n,double *q);
	int CLinearEquations::SVD(double *a,int m,int n,double *u,double *v,double eps,int ka);//一般实矩阵的奇异值分解
    void CLinearEquations::sss(double *fg,double *cs);
	void CLinearEquations::ppp(double *a,double *e,double *s,double *v,int m,int n);
	int CLinearEquations::Ginv(double *a,int m,int n,double *aa,double eps,double *u,double *v,int ka);//求广义逆的奇异值分解法
	///////////////////////////////////////////////////
	////          迭代方法解线性方程组             ////
	///////////////////////////////////////////////////
    void CLinearEquations::RCG(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,int Kiter,double eps);
	void CLinearEquations::RLSQR(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,int Kiter,double eps);
	void CLinearEquations::RSVD_LQ(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,double eps);
	void CLinearEquations::RQR_LQ(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N);
	void CLinearEquations::RBI_CG(Vector_double &A,Vector_double &x,Vector_double &b,int N,int Kiter,double eps);
	void CLinearEquations::RBI_CGSTAB(Vector_double &A,Vector_double &x,Vector_double &b,int N,int Kiter,double eps);
	///////////////////////////////////////////////////
	////          非迭代方法解线性方程组           ////
	///////////////////////////////////////////////////
    int CLinearEquations::SVD_LQ(double *a,int m,int n,double *b,double *x,double *aa,double eps,double *u,double *v);//求解线性最小二乘问题的广义逆法
	int CLinearEquations::QR_LQ(double *a,int m,int n,double *b,double *q);
	///////////////////////////////////////////////////
	////            复数同实数的转化               ////
	///////////////////////////////////////////////////
    void CLinearEquations::ExchangeComToReal(Vector_complex &A,Vector_complex &b,int M,int N,Vector_double &At,Vector_double &bt,int Mt,int Nt);
	///////////////////////////////////////////////////
	//////////  贝塞尔函数 汉开尔函数  ////////////////
	///////////////////////////////////////////////////
	Complex CLinearEquations::Bessel1(int nj,Complex x);
	Complex CLinearEquations::Bessel2(int nj,Complex x);
	Complex CLinearEquations::Hankel(int nj,int ku,Complex x);//nj表示阶数，ku表示第一类或第二类汉克尔函数
    //////////////////////////////////////////////
	/////        n可以取任意的整数         ///////
	//////////////////////////////////////////////
	double CLinearEquations::BESSEL1(int n,double x);
	double CLinearEquations::BESSEL2(int n,double x);
	Complex CLinearEquations::HANKEL(int n,int lei,double x);
	//////////////////////////////////////////////////
	//////////     解矩阵(Complex)    ////////////////
	//////////////////////////////////////////////////	    
	void CLinearEquations::LU_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N);
	void CLinearEquations::VSVD_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int M,int N,double eps);
	void CLinearEquations::VQR_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int M,int N);
	int CLinearEquations::GAUSS(Vector_complex &A,Vector_complex &x,Vector_complex &b,int n);
	void CLinearEquations::CGLS(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps);
	void CLinearEquations::CG(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps);
	void CLinearEquations::LSQR(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps);
};

// LinearEquations.cpp: implementation of the CLinearEquations class.
//
//////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CLinearEquations::CLinearEquations()
{
	PI=3.1415926;
	J.real(0);
	J.imag(1.0);

}

CLinearEquations::~CLinearEquations()
{

}
////////////////////////////////////////////
////////  A(N*L),B(L*M),result(N*M)   //////
////////////////////////////////////////////
void CLinearEquations::MultiplyMatrixMatrix(Complex **A,Complex **B,Complex **result,int N,int L,int M)
{
    int n,l,m;
	Complex zero(0,0);
	for(n=0;n<N;n++)
		for(m=0;m<M;m++)
		{
			result[n][m]=zero;
			for(l=0;l<L;l++)
			{
				result[n][m] +=A[n][l]*B[l][m];
			}
			//cout<<result[n][m]<<endl;
		}
}
////////////////////////////////////////////
/////// A(N*L),vec(L*1),result(N*1)   //////
////////////////////////////////////////////
void CLinearEquations::MultiplyMatrixVevtor(Complex **A,Complex *vec,Complex *result,int N,int L)
{
    int n,l;
	Complex zero(0,0);
	for(n=0;n<N;n++)
	{
		result[n]=zero;
		for(l=0;l<L;l++)
		{   
			result[n] +=A[n][l]*vec[l];
		}
		//cout<<result[n]<<endl;
	}
}
////////////////////////////////////////////
/////        A(N*L),result(L*N)        /////
////////////////////////////////////////////
void CLinearEquations::MatrixConjTurn(Complex **A,Complex **result,int N,int L)
{   
	 int n,l;
	 Complex zero(0,0);
	 for(n=0;n<N;n++)
		 for(l=0;l<L;l++)
		 {   
			 result[n][l]=zero;
			 result[n][l]=conj(A[l][n]);
			// cout<<result[n][l]<<endl;
		 }

}
////////////////////////////////////////////
////////////////////////////////////////////
double CLinearEquations::Norm1(Complex *vec,int N)
{
    double result=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += abs(vec[n]);
	}
	return result;  
}
////////////////////////////////////////////
////////////////////////////////////////////
double CLinearEquations::Norm2(Complex *vec,int N)
{
	double result=0;
	double result1=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += norm(vec[n]);
	}
	result1=sqrt(result);
	return result1;
}
////////////////////////////////////////////
/////    A(N*L),B(L*M),result(N*M)     /////
////////////////////////////////////////////
void CLinearEquations::VectorMultiplyMatrixMatrix(Vector_complex &A,Vector_complex &B,Vector_complex &result,int N,int L,int M)
{ 
	int n,l,m;
	Complex zero(0,0);
	for(n=0;n<N;n++)
		for(m=0;m<M;m++)
		{
			result[n*M+m]=zero;
			for(l=0;l<L;l++)
			{
				result[n*M+m] +=A[n*L+l]*B[l*M+m];
			}
			//cout<<result[n][m]<<endl;
		}
}
////////////////////////////////////////////
/////    A(N*L),vec(N*1),result(N*1)   /////
////////////////////////////////////////////
void CLinearEquations::VectorMultiplyMatrixVector(Vector_complex &A,Vector_complex &vec,Vector_complex &result,int N,int L)
{
    int n,l;
	Complex zero(0,0);
	for(n=0;n<N;n++)
	{
		result[n]=zero;
		for(l=0;l<L;l++)
		{   
			result[n] +=A[n*L+l]*vec[l];
		}
		//cout<<result[n]<<endl;
	}
}
////////////////////////////////////////////
///////     A(N*L),result(L*N)      ////////
////////////////////////////////////////////
void CLinearEquations::VectorMatrixConjTurn(Vector_complex &A,Vector_complex &result,int N,int L)
{
   int n,l;
	 Complex zero(0,0);
	 for(n=0;n<N;n++)
		 for(l=0;l<L;l++)
		 {   
			 result[n*L+l]=zero;
			 result[n*L+l]=conj(A[l*N+n]);
			// cout<<result[n][l]<<endl;
		 }
}
////////////////////////////////////////////
////////////////////////////////////////////
double CLinearEquations::VectorNorm1(Vector_complex &vec,int N)
{
	double result=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += abs(vec[n]);
	}
	return result;         

}
////////////////////////////////////////////
////////////////////////////////////////////
double CLinearEquations::VectorNorm2(Vector_complex &vec,int N)
{
    double result=0;
	double result1=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += norm(vec[n]);
	}
	result1=sqrt(result);
	return result1;
}
/////////////////////////////////////////////
/////////////////////////////////////////////
double CLinearEquations::VectorNorm(Vector_complex &vec,int N)
{
    double result=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += norm(vec[n]);
	}
	return result;
}
/////////////////////////////////////////////
/////////////////////////////////////////////
void CLinearEquations::VectorMultiplyVector(Vector_complex &vec1,Vector_complex &vec2,Vector_complex &result,int N)
{
    Complex zero(0,0);
	int n;
	result.assign(N,zero);
	for(n=0;n<N;n++)
	{
		result[n]=vec1[n]*vec2[n];
	}
}
/////////////////////////////////////////////
////    支持实数的基本矩阵、向量运算     //// 
/////////////////////////////////////////////
void CLinearEquations::VecMultMatMat(Vector_double &A,Vector_double &B,Vector_double &result,int N,int L,int M)//A(N*L),B(L*M),result(N*M)
{
	int n,l,m;
	for(n=0;n<N;n++)
		for(m=0;m<M;m++)
		{
			result[n*M+m]=0;
			for(l=0;l<L;l++)
			{
				result[n*M+m] +=A[n*L+l]*B[l*M+m];
			}
			//cout<<result[n][m]<<endl;
		}
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
void CLinearEquations::VecMultMatVec(Vector_double &A,Vector_double &vec,Vector_double &result,int N,int L)//A(N*L),vec(N*1),result(N*1)
{
    int n,l;
	for(n=0;n<N;n++)
	{
		result[n]=0;
		for(l=0;l<L;l++)
		{   
			result[n] +=A[n*L+l]*vec[l];
		}
		//cout<<result[n]<<endl;
	}
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
void CLinearEquations::VecMatConjTurn(Vector_double &A,Vector_double &result,int N,int L)//A(N*L),result(L*N)
{
     int n,l;
	 for(n=0;n<N;n++)
		 for(l=0;l<L;l++)
		 {   
			 result[n*L+l]=0;
			 result[n*L+l]=A[l*N+n];
			// cout<<result[n][l]<<endl;
		 }
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
double CLinearEquations::VecNorm1(Vector_double &vec,int N)
{
	double result=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += fabs(vec[n]);
	}
	return result;  
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
double CLinearEquations::VecNorm2(Vector_double &vec,int N)
{
    double result=0;
	double result1=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += vec[n]*vec[n];
	}
	result1=sqrt(result);
	return result1;
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
double CLinearEquations::VecNorm(Vector_double &vec,int N)
{
    double result=0;
	int n;
	for(n=0;n<N;n++)
	{
		 result += vec[n]*vec[n];
	}
	return result;
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
void CLinearEquations::VecMultVec(Vector_double &vec1,Vector_double &vec2,Vector_double &result,int N)//vec1(N),vec2(N),result(N)
{
	int n;
	result.assign(N,0);
	for(n=0;n<N;n++)
	{
		result[n]=vec1[n]*vec2[n];
	}
}
/////////////////////////////////////////////
////
/////////////////////////////////////////////
double CLinearEquations::VecSumVecMultVec(Vector_double &vec1,Vector_double &vec2,int N)
{
	double result=0;
	for(int n=0;n<N;n++)
	{
		result +=vec1[n]*vec2[n];
	}
	return result;
}
//-----------------------------------------------------------------------------
//                  复数贝塞尔及汉克尔函数
//-----------------------------------------------------------------------------
Complex CLinearEquations::Bessel1(int nj,Complex x)
{
	static double a[6]={2.2499997,1.2656208,0.3163866,0.0444479,0.0039444,0.0002100};
	static double c[6]={0.00000077,0.00552740,0.00009512,0.00137237,0.00072805,0.00014476};
	static double d[6]={0.04166397,0.00003954,0.00262573,0.00054125,0.00029333,0.00013558};
	static double e[6]={0.56249985,0.21093573,0.03954289,0.00443319,0.00031761,0.00001109};
	static double g[6]={0.00000156,0.01659667,0.00017105,0.00249511,0.00113653,0.00020033};
	static double h[6]={0.12499612,0.00005650,0.00637879,0.00074348,0.00079824,0.00029166};
	Complex re;
	Complex uv,x1,f0,c0;
	double Module_x=sqrt(x.real()*(x.real())+x.imag()*(x.imag()));
	if(nj == 0)
	{
		
		if(Module_x <= 3.0)
		{
			re=1.0-a[0]*pow(x/3.0,2)+a[1]*pow(x/3.0,4)-a[2]*pow(x/3.0,6)+a[3]*pow(x/3.0,8)
				-a[4]*pow(x/3.0,10)+a[5]*pow(x/3.0,12);
		}
		else
		{
			uv=3.0/x;
			x1=1.0/sqrt(x);
			f0=0.79788456-c[0]*uv-c[1]*pow(uv,2)-c[2]*pow(uv,3)+c[3]*pow(uv,4)-c[4]*pow(uv,5)+c[5]*pow(uv,6);
			c0=x-0.78539816-d[0]*uv-d[1]*pow(uv,2)+d[2]*pow(uv,3)-d[3]*pow(uv,4)-d[4]*pow(uv,5)+d[5]*pow(uv,6);
			re=x1*f0*cos(c0);
		}
	}
	else if(nj == 1)//1阶的第一类贝塞尔函数
	{
		if(Module_x <= 3.0)
		{
			re=0.5-e[0]*pow(x/3.0,2)+e[1]*pow(x/3.0,4)-e[2]*pow(x/3.0,6)+e[3]*pow(x/3.0,8)
				-e[4]*pow(x/3.0,10)+e[5]*pow(x/3.0,12);
			re=x*re;
		}
		else
		{
			uv=3.0/x;
			x1=1.0/sqrt(x);
			f0=0.79788456+g[0]*uv+g[1]*pow(uv,2)+g[2]*pow(uv,3)-g[3]*pow(uv,4)+g[4]*pow(uv,5)-g[5]*pow(uv,6);
			c0=x-2.35619449+h[0]*uv+h[1]*pow(uv,2)-h[2]*pow(uv,3)+h[3]*pow(uv,4)-h[4]*pow(uv,5)-h[5]*pow(uv,6);
			re=x1*f0*cos(c0);
		}
	}
	
	else
	{
		cout<<"nj must be 0 or 1"<<endl;
	}
	return re;
}
//-------------------------------------------------------------------------------
Complex CLinearEquations::Bessel2(int nj,Complex x)
{
	static double b[6]={0.60559366,0.74350384,0.25300117,0.04261214,0.00427916,0.00024846};
	static double c[6]={0.00000077,0.00552740,0.00009512,0.00137237,0.00072805,0.00014476};
	static double d[6]={0.04166397,0.00003954,0.00262573,0.00054125,0.00029333,0.00013558};
	static double f[6]={0.2212091,2.1682709,1.3164827,0.3123951,0.0400976,0.0027873};
	static double g[6]={0.00000156,0.01659667,0.00017105,0.00249511,0.00113653,0.00020033};
	static double h[6]={0.12499612,0.00005650,0.00637879,0.00074348,0.00079824,0.00029166};
	Complex re;
	Complex uv,x1,f0,c0;
	double Module_x=sqrt(x.real()*(x.real())+x.imag()*(x.imag()));
	if(nj == 0)
	{
		if(Module_x <= 3.0)
		{
			re=(2.0/PI)*log(x/2.0)*Bessel1(0,x)+0.36746691+b[0]*pow(x/3.0,2)-b[1]*pow(x/3.0,4)
				+b[2]*pow(x/3.0,6)-b[3]*pow(x/3.0,8)+b[4]*pow(x/3.0,10)-b[5]*pow(x/3.0,12);//
		}
		else
		{
			uv=3.0/x;
			x1=1.0/sqrt(x);
			f0=0.79788456-c[0]*uv-c[1]*pow(uv,2)-c[2]*pow(uv,3)+c[3]*pow(uv,4)-c[4]*pow(uv,5)+c[5]*pow(uv,6);
			c0=x-0.78539816-d[0]*uv-d[1]*pow(uv,2)+d[2]*pow(uv,3)-d[3]*pow(uv,4)-d[4]*pow(uv,5)+d[5]*pow(uv,6);
			re=x1*f0*sin(c0);		    
		}
	}
	else if(nj == 1)//1阶的第二类贝塞尔函数
	{
		if(Module_x<= 3.0)
		{          
			re=(2.0/PI)*x*log(x/2.0)*Bessel1(1,x)-0.6366198+f[0]*pow(x/3.0,2)+f[1]*pow(x/3.0,4)
				-f[2]*pow(x/3.0,6)+f[3]*pow(x/3.0,8)-f[4]*pow(x/3.0,10)+f[5]*pow(x/3.0,12);
			re=re/x;
		}
		else
		{
			uv=3.0/x;
			x1=1.0/sqrt(x);
			f0=0.79788456+g[0]*uv+g[1]*pow(uv,2)+g[2]*pow(uv,3)-g[3]*pow(uv,4)+g[4]*pow(uv,5)-g[5]*pow(uv,6);
			c0=x-2.35619449+h[0]*uv+h[1]*pow(uv,2)-h[2]*pow(uv,3)+h[3]*pow(uv,4)-h[4]*pow(uv,5)-h[5]*pow(uv,6);
			re=x1*f0*sin(c0);
		}
	}
	
	else
	{
		cout<<"nj must be 0 or 1"<<endl;
	}
	return re;
}
//-----------------------------------------------------------------------------------
Complex CLinearEquations::Hankel(int nj,int ku,Complex x)//nj表示阶数，ku表示第一类或第二类汉克尔函数
{
	Complex re;
	Complex a,b;
	a=Bessel1(nj,x);
	b=Bessel2(nj,x);
	switch(ku) 
	{
	case 1: {
		
		re=a+J*b;
		break;
			}
	case 2: {
				
		re=a-J*b;
		break;
			}
	}	
	return re;
}
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//////////////////////////////////////////////
//////////////////////////////////////////////
void CLinearEquations::CGLS(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps)
{
  Complex zero(0,0);
  Vector_complex vector_d,vector_p,vector_s,vector_q,matrix_A_turn;
  vector_d.assign(N,zero);   vector_p.assign(N,zero);
  vector_s.assign(N,zero);   vector_q.assign(N,zero);
  x.assign(N,zero);          matrix_A_turn.assign(N*N,zero);
  int k;
//////////////INITIAL
  for(k=0;k<N;k++)
  {
	  vector_d[k]=b[k];
  }
  VectorMatrixConjTurn(A,matrix_A_turn,N,N);
  VectorMultiplyMatrixVector(matrix_A_turn,vector_d,vector_p,N,N);
  for(k=0;k<N;k++)
  {
	  vector_s[k]=vector_p[k];
  }
//////////////MAIN LOOP
  for(int iter=0;iter<Kiter;iter++)
  {
	  VectorMultiplyMatrixVector(A,vector_p,vector_q,N,N);
	  double temp_s_pre=VectorNorm(vector_s,N);
	  double temp_q=VectorNorm(vector_q,N);
	  double alpha=temp_s_pre/temp_q;
	  /////////////
	  for(k=0;k<N;k++)
	  {
         x[k]  +=alpha*vector_p[k];
		 vector_d[k] -=alpha*vector_q[k];
	  }
	  /////////////
	  VectorMultiplyMatrixVector(matrix_A_turn,vector_d,vector_s,N,N);
	  double temp_s_lat=VectorNorm(vector_s,N);
	  double beta=temp_s_lat/temp_s_pre;
	  for(k=0;k<N;k++)
	  {
		  vector_p[k]=vector_s[k]+beta*vector_p[k];
	  }
	  double temp=VectorNorm2(vector_d,N)/VectorNorm2(b,N);
      cout<<"rcs["<<iter<<"]--------:"<<temp<<endl;
	  if(temp<eps)
	  {
		  break;
	  }	  
  }

}
//////////////////////////////////////////////
//////////////////////////////////////////////
void CLinearEquations::CG(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps)
{
	Complex zero(0,0);
	Vector_complex vector_r,vector_p,vector_q,vector_Ar,matrix_A_turn;
	vector_r.assign(N,zero);   
	vector_p.assign(N,zero);   vector_q.assign(N,zero);
	vector_Ar.assign(N,zero);  matrix_A_turn.assign(N*N,zero);
///////////////INITIAL
	x.assign(N,zero);
	int m;
	for(m=0;m<N;m++)
	{
		vector_r[m]=b[m];
	}
	VectorMatrixConjTurn(A,matrix_A_turn,N,N);
	VectorMultiplyMatrixVector(matrix_A_turn,vector_r,vector_p,N,N);
//////////////MAIN LOOP
	for(int iter=0;iter<Kiter;iter++)
	{
		double temp_r_pre=VectorNorm(vector_r,N);
		double temp_p=VectorNorm(vector_p,N);
		double alpha=temp_r_pre/temp_p;
		//////////////alpha
		VectorMultiplyMatrixVector(A,vector_p,vector_q,N,N);
		for(m=0;m<N;m++)
		{
			x[m]  +=alpha*vector_p[m];
			vector_r[m] -=alpha*vector_q[m];
		}
		//////////////x[N],r[N]
		double temp_r_lat=VectorNorm(vector_r,N);
		double beta=temp_r_lat/temp_r_pre;
		VectorMultiplyMatrixVector(matrix_A_turn,vector_r,vector_Ar,N,N);
		for(m=0;m<N;m++)
		{
			vector_p[m]=vector_Ar[m]+beta*vector_p[m];
		}
		double temp=VectorNorm2(vector_r,N)/VectorNorm2(b,N);
		cout<<"rcs["<<iter<<"]-------:"<<temp<<endl;
		if(temp<eps)
		{
			break;
		}
	}
}
//////////////////////////////////////////////
//////////////////////////////////////////////
void CLinearEquations::LSQR(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N,int Kiter,double eps)
{
   Complex zero(0,0);
   Vector_complex vector_u,vector_v,vector_w,matrix_A_turn,vector_temp,vector_delta_x;
   Vector_complex vector_u_temp,vector_v_temp,vector_r;
   vector_u.assign(N,zero);     vector_v.assign(N,zero);
   vector_w.assign(N,zero);     matrix_A_turn.assign(N*N,zero);
   vector_temp.assign(N,zero);  vector_delta_x.assign(N,zero);
   vector_u_temp.assign(N,zero);vector_v_temp.assign(N,zero);
   vector_r.assign(N,zero);
/////////////////INITIAL
   x.assign(N,zero);
   double beta=VectorNorm2(b,N);
   int i;
   for(i=0;i<N;i++)
   {
       vector_u[i]=b[i]/beta;
   }
   VectorMatrixConjTurn(A,matrix_A_turn,N,N);
   VectorMultiplyMatrixVector(matrix_A_turn,vector_u,vector_v_temp,N,N);
   double alpha=VectorNorm2(vector_v_temp,N);
   for(i=0;i<N;i++)
   {
	   vector_v[i]=vector_v_temp[i]/alpha;
	   vector_w[i]=vector_v[i];
   }
   double phi_1=beta;
   double p_1=alpha;
//////////////////MAIN LOOP
   for(int iter=0;iter<Kiter;iter++)
   {
	   VectorMultiplyMatrixVector(A,vector_v,vector_temp,N,N);
	   for(i=0;i<N;i++)
	   {
		   vector_u_temp[i]=vector_temp[i]-alpha*vector_u[i];
	   }
	   beta=VectorNorm2(vector_u_temp,N);
	   for(i=0;i<N;i++)
	   {
		   vector_u[i]=vector_u_temp[i]/beta;
	   }
	   ////////////////u
	   VectorMultiplyMatrixVector(matrix_A_turn,vector_u,vector_temp,N,N);
	   for(i=0;i<N;i++)
	   {
		   vector_v_temp[i]=vector_temp[i]-beta*vector_v[i];
	   }
	   alpha=VectorNorm2(vector_v_temp,N);
	   for(i=0;i<N;i++)
	   {
		   vector_v[i]=vector_v_temp[i]/alpha;
	   }
	   ///////////////v
	   double p=sqrt(p_1*p_1+beta*beta);
	   double c=p_1/p;
	   double s=beta/p;
	   double theta=s*alpha;
	          p_1=-c*alpha;
	   double phi=c*phi_1;
	          phi_1=s*phi_1;
	   for(i=0;i<N;i++)
	   {
		   x[i]  +=(phi/p)*vector_w[i];
		   vector_w[i]=vector_v[i]-(theta/p)*vector_w[i];
		   vector_delta_x[i]=(phi/p)*vector_w[i];
	   }
	   VectorMultiplyMatrixVector(A,x,vector_temp,N,N);
       for(i=0;i<N;i++)
	   {
		   vector_r[i]=b[i]-vector_temp[i];
	   }
	   double temp=VectorNorm2(vector_r,N)/VectorNorm2(b,N);
	   cout<<"rcs["<<iter<<"]-------:"<<temp<<endl;
	   if(temp<eps)
	   {
		   break;
	   }	   
   }
}
//////////////////////////////////////////////
//////////////////////////////////////////////
double CLinearEquations::BESSEL1(int n,double x)
{
    int i,m;
    double t,y,z,p,q,s,b0,b1;
    static double a[6]={ 57568490574.0,-13362590354.0,
             651619640.7,-11214424.18,77392.33017,
            -184.9052456};
    static double b[6]={ 57568490411.0,1029532985.0,
             9494680.718,59272.64853,267.8532712,1.0};
    static double c[6]={ 72362614232.0,-7895059235.0,
             242396853.1,-2972611.439,15704.4826,
             -30.16036606};
    static double d[6]={ 144725228443.0,2300535178.0,
             18583304.74,99447.43394,376.9991397,1.0};
    static double e[5]={ 1.0,-0.1098628627e-02,
             0.2734510407e-04,-0.2073370639e-05,
             0.2093887211e-06};
    static double f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
    static double g[5]={ 1.0,0.183105e-02,
             -0.3516396496e-04,0.2457520174e-05,
             -0.240337019e-06};
    static double h[5]={ 0.4687499995e-01,
             -0.2002690873e-03,0.8449199096e-05,
             -0.88228987e-06,0.105787412e-06};
    t=fabs(x);
    if (n<0) n=-n;
    if (n!=1)
      { if (t<8.0)
          { y=t*t; p=a[5]; q=b[5];
	    for (i=4; i>=0; i--)
              { p=p*y+a[i]; q=q*y+b[i];}
            p=p/q;
          }
        else
          { z=8.0/t; y=z*z;
            p=e[4]; q=f[4];
            for (i=3; i>=0; i--)
              { p=p*y+e[i]; q=q*y+f[i];}
            s=t-0.785398164;
            p=p*cos(s)-z*q*sin(s);
            p=p*sqrt(0.636619772/t);
          }
      }
    if (n==0) return(p);
    b0=p;
    if (t<8.0)
      { y=t*t; p=c[5]; q=d[5];
        for (i=4; i>=0; i--)
          { p=p*y+c[i]; q=q*y+d[i];}
        p=x*p/q;
      }
    else
      { z=8.0/t; y=z*z;
        p=g[4]; q=h[4];
        for (i=3; i>=0; i--)
          { p=p*y+g[i]; q=q*y+h[i];}
        s=t-2.356194491;
        p=p*cos(s)-z*q*sin(s);
        p=p*x*sqrt(0.636619772/t)/t;
      }
    if (n==1) return(p);
    b1=p;
    if (x==0.0) return(0.0);
    s=2.0/t;
    if (t>1.0*n)
      { if (x<0.0) b1=-b1;
        for (i=1; i<=n-1; i++)
          { p=s*i*b1-b0; b0=b1; b1=p;}
      }
    else
      { m=(n+(int)sqrt(40.0*n))/2;
        m=2*m;
        p=0.0; q=0.0; b0=1.0; b1=0.0;
        for (i=m-1; i>=0; i--)
          { t=s*(i+1)*b0-b1;
            b1=b0; b0=t;
            if (fabs(b0)>1.0e+10)
              { b0=b0*1.0e-10; b1=b1*1.0e-10;
                p=p*1.0e-10; q=q*1.0e-10;
              }
            if ((i+2)%2==0) q=q+b0;
            if ((i+1)==n) p=b1;
          }
        q=2.0*q-b0; p=p/q;
      }
    if ((x<0.0)&&(n%2==1)) p=-p;
    return(p);
}
//////////////////////////////////////////////
//////////////////////////////////////////////
double CLinearEquations::BESSEL2(int n,double x)
{
    int i;
    double y,z,p,q,s,b0,b1;
    static double a[6]={ -2.957821389e+9,7.062834065e+9,
             -5.123598036e+8,1.087988129e+7,-8.632792757e+4,
             2.284622733e+2};
    static double b[6]={ 4.0076544269e+10,7.452499648e+8,
           7.189466438e+6,4.74472647e+4,2.261030244e+2,1.0};
    static double c[6]={ -4.900604943e+12,1.27527439e+12,
            -5.153438139e+10,7.349264551e+8,-4.237922726e+6,
             8.511937935e+3};
    static double d[7]={ 2.49958057e+13,4.244419664e+11,
            3.733650367e+9,2.245904002e+7,1.02042605e+5,
            3.549632885e+2,1.0};
    static double e[5]={ 1.0,-0.1098628627e-02,
             0.2734510407e-04,-0.2073370639e-05,
             0.2093887211e-06};
    static double f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
    static double g[5]={ 1.0,0.183105e-02,
             -0.3516396496e-04,0.2457520174e-05,
             -0.240337019e-06};
    static double h[5]={ 0.4687499995e-01,
             -0.2002690873e-03,0.8449199096e-05,
             -0.88228987e-06,0.105787412e-06};
    if (n<0) n=-n;
    if (x<0.0) x=-x;
    if (x==0.0) return(-1.0e+70);
    if (n!=1)
      { if (x<8.0)
          { y=x*x; p=a[5]; q=b[5];
	    for (i=4; i>=0; i--)
              { p=p*y+a[i]; q=q*y+b[i];}
            p=p/q+0.636619772*BESSEL1(0,x)*log(x);
          }
        else
          { z=8.0/x; y=z*z;
            p=e[4]; q=f[4];
            for (i=3; i>=0; i--)
              { p=p*y+e[i]; q=q*y+f[i];}
            s=x-0.785398164;
            p=p*sin(s)+z*q*cos(s);
            p=p*sqrt(0.636619772/x);
          }
      }
    if (n==0) return(p);
    b0=p;
    if (x<8.0)
      { y=x*x; p=c[5]; q=d[6];
        for (i=4; i>=0; i--)
          { p=p*y+c[i]; q=q*y+d[i+1];}
        q=q*y+d[0];
        p=x*p/q+0.636619772*(BESSEL1(1,x)*log(x)-1.0/x);;
      }
    else
      { z=8.0/x; y=z*z;
        p=g[4]; q=h[4];
        for (i=3; i>=0; i--)
          { p=p*y+g[i]; q=q*y+h[i];}
        s=x-2.356194491;
        p=p*sin(s)+z*q*cos(s);
        p=p*sqrt(0.636619772/x);
      }
    if (n==1) return(p);
    b1=p;
    s=2.0/x;
    for (i=1; i<=n-1; i++)
      { p=s*i*b1-b0; b0=b1; b1=p;}
    return(p);
}
//////////////////////////////////////////////
//////////////////////////////////////////////
Complex CLinearEquations::HANKEL(int n,int lei,double x)
{
	Complex re;
	Complex a,b;
	a=BESSEL1(n,x);
	b=BESSEL2(n,x);
	switch(lei) 
	{
	case 1: {
		
		re=a+J*b;
		break;
			}
	case 2: {
				
		re=a-J*b;
		break;
			}
	}	
	return re;
}
//////////////////////////////////////////////
//////////////////////////////////////////////
int CLinearEquations::GAUSS(Vector_complex &A,Vector_complex &x,Vector_complex &b,int n)
 // double ar[],ai[],br[],bi[];
{
	double *ar,*ai,*br,*bi;
	ar=new double[n*n];  ai=new double[n*n];
	br=new double[n];    bi=new double[n];
	for(int m=0;m<n*n;m++)
	{
		ar[m]=A[m].real();
		ai[m]=A[m].imag();
	}
	for( m=0;m<n;m++)
	{
		br[m]=b[m].real();
		bi[m]=b[m].imag();
	}
	int *js,l,k,i,j,is,u,v;
    double p,q,s,d;
    //js=malloc(n*sizeof(int));
	js=new int[n];
    for (k=0;k<=n-2;k++)
      { d=0.0;
        for (i=k;i<=n-1;i++)
        for (j=k;j<=n-1;j++)
          { u=i*n+j;
            p=ar[u]*ar[u]+ai[u]*ai[u];
            if (p>d) {d=p;js[k]=j;is=i;}
          }
        if (d+1.0==1.0)
          { 
//			delete[js];
			cout<<"err**fail1"<<endl;
			//free(js); printf("err**fail\n");
            return(0);
          }
        if (is!=k)
          { for (j=k;j<=n-1;j++)
              { u=k*n+j; v=is*n+j;
                p=ar[u]; ar[u]=ar[v]; ar[v]=p;
                p=ai[u]; ai[u]=ai[v]; ai[v]=p;
              }
            p=br[k]; br[k]=br[is]; br[is]=p;
            p=bi[k]; bi[k]=bi[is]; bi[is]=p;
          }
        if (js[k]!=k)
          for (i=0;i<=n-1;i++)
            { u=i*n+k; v=i*n+js[k];
              p=ar[u]; ar[u]=ar[v]; ar[v]=p;
              p=ai[u]; ai[u]=ai[v]; ai[v]=p;
            }
        v=k*n+k;
        for (j=k+1;j<=n-1;j++)
          { u=k*n+j;
            p=ar[u]*ar[v]; q=-ai[u]*ai[v];
            s=(ar[v]-ai[v])*(ar[u]+ai[u]);
            ar[u]=(p-q)/d; ai[u]=(s-p-q)/d;
          }
        p=br[k]*ar[v]; q=-bi[k]*ai[v];
        s=(ar[v]-ai[v])*(br[k]+bi[k]);
        br[k]=(p-q)/d; bi[k]=(s-p-q)/d;
        for (i=k+1;i<=n-1;i++)
          { u=i*n+k;
            for (j=k+1;j<=n-1;j++)
              { v=k*n+j; l=i*n+j;
                p=ar[u]*ar[v]; q=ai[u]*ai[v];
                s=(ar[u]+ai[u])*(ar[v]+ai[v]);
                ar[l]=ar[l]-p+q;
                ai[l]=ai[l]-s+p+q;
              }
            p=ar[u]*br[k]; q=ai[u]*bi[k];
            s=(ar[u]+ai[u])*(br[k]+bi[k]);
            br[i]=br[i]-p+q; bi[i]=bi[i]-s+p+q;
          }
      }
    u=(n-1)*n+n-1;
    d=ar[u]*ar[u]+ai[u]*ai[u];
    if (d+1.0==1.0)
      { 
//		delete[js];
	    cout<<"err**fail2"<<endl;
		//free(js); printf("err**fail\n");
        return(0);
      }
    p=ar[u]*br[n-1]; q=-ai[u]*bi[n-1];
    s=(ar[u]-ai[u])*(br[n-1]+bi[n-1]);
    br[n-1]=(p-q)/d; bi[n-1]=(s-p-q)/d;
    for (i=n-2;i>=0;i--)
    for (j=i+1;j<=n-1;j++)
      { u=i*n+j;
        p=ar[u]*br[j]; q=ai[u]*bi[j];
        s=(ar[u]+ai[u])*(br[j]+bi[j]);
        br[i]=br[i]-p+q;
        bi[i]=bi[i]-s+p+q;
      }
    js[n-1]=n-1;
    for (k=n-1;k>=0;k--)
      if (js[k]!=k)
        { p=br[k]; br[k]=br[js[k]]; br[js[k]]=p;
          p=bi[k]; bi[k]=bi[js[k]]; bi[js[k]]=p;
        }
//	  delete[js];
    //free(js);
   for(m=0;m<n;m++)
   {
	   x[m]=br[m]+J*bi[m];
   }
    return(1);
}
/////////////////////////////////////////
/////////////////////////////////////////
int CLinearEquations::LU(Vector_complex &a,Vector_complex &l,Vector_complex &u,int n)
//  int n;
//  double a[],l[],u[];
{ int i,j,k,w,v,ll;
    for (k=0; k<=n-2; k++)
      { ll=k*n+k;
	if (
		abs(a[ll])+1.0==1.0)
          { cout<<"fail"<<endl;
	        return(0);}
        for (i=k+1; i<=n-1; i++)
	  { w=i*n+k; a[w]=a[w]/a[ll];}
        for (i=k+1; i<=n-1; i++)
          { w=i*n+k;
            for (j=k+1; j<=n-1; j++)
              { v=i*n+j;
                a[v]=a[v]-a[w]*a[k*n+j];
              }
          }
      }
    for (i=0; i<=n-1; i++)
      { for (j=0; j<i; j++)
          { w=i*n+j; l[w]=a[w]; u[w]=0.0;}
        w=i*n+i;
        l[w]=1.0; u[w]=a[w];
        for (j=i+1; j<=n-1; j++)
          { w=i*n+j; l[w]=0.0; u[w]=a[w];}
      }
    return(1);
}
//////////////////////////////////////
//////////////////////////////////////
void CLinearEquations::LU_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int N)
{
  ////////////////
  //INITIAL
  ////////////////
  Complex  zero(0,0); 
  Vector_complex L,U,y;
  L.assign(N*N,zero);   U.assign(N*N,zero);   y.assign(N,zero);
  ////////////////
  //LU分解
  ////////////////
  LU(A,L,U,N);
  ////////////////
  //  解Ly=b
  ////////////////
  /* for(int i=0;i<N;i++)
   {
	  for(int j=0;j<N;j++)
	  {
		  cout<<U[i*N+j]<<"  ";
	  }
	  cout<<"   "<<endl;
   }*/
  y[0]=b[0];
  for(int i=1;i<N;i++)
  {   
	  Complex temp_ly=zero;
	  for(int k=0;k<i;k++)
	  {
		   temp_ly +=L[i*N+k]*y[k];
	  }
	  y[i]=b[i]-temp_ly;
  }
  ////////////////
  //  解Ux=y
  ////////////////
  x[N-1]=y[N-1]/U[(N-1)*N+N-1];
  for(int k=N-2;k>=0;k--)
  {
	  Complex temp_c=y[k];
	  for(int j=k+1;j<N;j++)
	  {
           temp_c=temp_c-U[k*N+j]*x[j];
	  }
	  x[k]=temp_c/U[k*N+k];
  }
  //y.assign(N,zero);
  //VectorMultiplyMatrixVector(A,x,y,N,N);
  /*for(i=0;i<N;i++)
  {
	  cout<<"y["<<i<<"]----:"<<y[i]<<endl;
  }*/
}
//////////////////////////////////////
//////////////////////////////////////
void CLinearEquations::RCG(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,int Kiter,double eps)
{
   //////////////
   // INITIAL
   //////////////
   Vector_double Mat_AT,Vec_d,Vec_Ax,Vec_p,Vec_s,Vec_q,Vec_rre;
   x.assign(N,0);         Vec_rre.assign(M,0);
   Mat_AT.assign(M*N,0);  Vec_p.assign(N,0);
   Vec_d.assign(M,0);     Vec_s.assign(N,0);
   Vec_Ax.assign(M,0);    Vec_q.assign(M,0); 
   for(int m=0;m<M;m++)
   {
	   Vec_d[m]=b[m];
   }
   VecMatConjTurn(A,Mat_AT,M,N);
   VecMultMatVec(Mat_AT,Vec_d,Vec_p,N,M);
   VecMultMatVec(Mat_AT,Vec_d,Vec_s,N,M);
   //////////////
   // MAIN LOOP
   //////////////
   for(int k=0;k<Kiter;k++)
   {
       /////////////////
	   //
	   /////////////////
       VecMultMatVec(A,Vec_p,Vec_q,M,N);
	   double temp_s_pre=VecNorm(Vec_s,N);
	   double temp_q=VecNorm(Vec_q,M);
	   double alpha=temp_s_pre/temp_q;
	   /////////////////
	   //  更新变量
	   /////////////////
	   for(int n=0;n<N;n++)
	   {
		   x[n]=x[n]+alpha*Vec_p[n];
	   }
	   /////////////////
	   //  中断条件
	   /////////////////
	   double temp_x=VecNorm2(b,M);
	   VecMultMatVec(A,x,Vec_Ax,M,N);
	   for(int m=0;m<M;m++)
	   {
		   Vec_rre[m]=b[m]-Vec_Ax[m];
	   }
	   double temp_rre=VecNorm2(Vec_rre,M);
	   double RRE=temp_rre/temp_x;
	   cout<<" RRE["<<k<<"]--------:"<<RRE<<endl;
	   if(RRE<eps)
	   {
		   break;
	   }
	   for(m=0;m<M;m++)
	   {
		   Vec_d[m]=Vec_d[m]-alpha*Vec_q[m];
	   }
	   /////////////////
	   // 更新迭代变量
	   /////////////////
	   VecMultMatVec(Mat_AT,Vec_d,Vec_s,N,M);
	   /////////////////
	   double temp_s_lat=VecNorm(Vec_s,N);
	   double beta=temp_s_lat/temp_s_pre;
	   /////////////////
	   for(n=0;n<N;n++)
	   {
		   Vec_p[n]=Vec_s[n]+beta*Vec_p[n];
	   }
   }
}
//////////////////////////////////////
////
//////////////////////////////////////
void CLinearEquations::ExchangeComToReal(Vector_complex &A,Vector_complex &b,int M,int N,Vector_double &At,Vector_double &bt,int Mt,int Nt)
{
   ////////////////
   //  INITIAL
   ////////////////
   Mt=2*M;
   Nt=2*N;
   Vector_double Ar,Ai,br,bi;
   Ar.assign(M*N,0); Ai.assign(M*N,0);
   br.assign(M,0);   bi.assign(M,0);
   int m,n;
   for(m=0;m<M;m++)
   {
	   for(n=0;n<N;n++)
	   {
		   Ar[m*N+n]=A[m*N+n].real();
		   Ai[m*N+n]=A[m*N+n].imag();
	   }
	   br[m]=b[m].real();
	   bi[m]=b[m].imag();
   }
   /////////////////
   //  MAIN LOOP
   /////////////////
   //////////////////////////
   //      Ar Ai  At
   //////////////////////////
   for(m=0;m<Mt;m++)
   {
	   if(m<M)
	   {
		   for(n=0;n<N;n++)
		   {
			   At[m*Nt+n]=Ar[m*N+n];
		   }
		   for(n=N;n<Nt;n++)
		   {
			   At[m*Nt+n]=-Ai[m*N+n-N];
		   }
	   }
	   else
	   {
		   for(n=0;n<N;n++)
		   {
			   At[m*Nt+n]=Ai[(m-M)*N+n];
		   }
		   for(n=N;n<Nt;n++)
		   {
			   At[m*Nt+n]=Ar[(m-M)*N+n-N];
		   }
	   }
   }
   ////////////////////////
   //     bt br bi
   ////////////////////////
   for(m=0;m<Mt;m++)
   {
	   if(m<M)
	   {
		   bt[m]=br[m];
	   }
	   else
	   {
		   bt[m]=bi[m-M];
	   }
   }
}
//////////////////////////////////////
////
//////////////////////////////////////
void CLinearEquations::RLSQR(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,int Kiter,double eps)
{
	///////////////
	// INITIAL
	///////////////
    Vector_double Mat_AT,Vec_u,Vec_v,Vec_w,Vec_temp_v,Vec_temp_u,Vec_Av,Vec_ATu,Vec_Ax,Vec_rre;
	x.assign(N,0);          Vec_rre.assign(M,0);       Vec_Ax.assign(M,0);
	Mat_AT.assign(M*N,0);   Vec_u.assign(M,0);         Vec_v.assign(N,0);
	Vec_w.assign(N,0);      Vec_temp_u.assign(M,0);    Vec_temp_v.assign(N,0);
	Vec_Av.assign(M,0);     Vec_ATu.assign(N,0);
	int m,n,k;
	double beta=VecNorm2(b,M);
	for(m=0;m<M;m++)
	{
		Vec_u[m]=b[m]/beta;
	}
	VecMatConjTurn(A,Mat_AT,M,N);
	VecMultMatVec(Mat_AT,Vec_u,Vec_temp_v,N,M);
	double alpha=VecNorm2(Vec_temp_v,N);
	for(n=0;n<N;n++)
	{
		Vec_w[n]=Vec_temp_v[n]/alpha;
		Vec_v[n]=Vec_temp_v[n]/alpha;
	}
	double temp_Phi=beta;
	double temp_p=alpha;
	////////////////
	// MAIN LOOP
	////////////////
	for(k=0;k<Kiter;k++)
	{
		////////////
		//  u,v
		////////////
		VecMultMatVec(A,Vec_v,Vec_Av,M,N);
		for(m=0;m<M;m++)
		{
			Vec_temp_u[m]=Vec_Av[m]-alpha*Vec_u[m];
		}
		beta=VecNorm2(Vec_temp_u,M);
		for(m=0;m<M;m++)
		{
			Vec_u[m]=Vec_temp_u[m]/beta;
		}
		////////////
		VecMultMatVec(Mat_AT,Vec_u,Vec_ATu,N,M);
		for(n=0;n<N;n++)
		{
			Vec_temp_v[n]=Vec_ATu[n]-beta*Vec_v[n];
		}
		alpha=VecNorm2(Vec_temp_v,N);
		for(n=0;n<N;n++)
		{
			Vec_v[n]=Vec_temp_v[n]/alpha;
		}
		//////////////
		//
		//////////////
		double p=sqrt(temp_p*temp_p+beta*beta);
		double c=temp_p/p;
		double s=beta/p;
		double theta=s*alpha;
		       temp_p=-c*alpha;
		double Phi=c*temp_Phi;
		       temp_Phi=s*temp_Phi;
		for(n=0;n<N;n++)
		{
			x[n]=x[n]-(Phi/p)*Vec_w[n];
			Vec_w[n]=Vec_v[n]-(theta/p)*Vec_w[n];
		}
		//////////////
		//
		//////////////
		double temp_rre=fabs(Phi/p)*VecNorm2(Vec_w,N);
		cout<<"RRE["<<k<<"]---------:"<<temp_rre<<endl;
		if(temp_rre<eps)
		{
			break;
		}
	}
}
//////////////////////////////////////
////
//////////////////////////////////////
void CLinearEquations::RBI_CG(Vector_double &A,Vector_double &x,Vector_double &b,int N,int Kiter,double eps)
{
    ////////////
	// INITIAL
	////////////
	Vector_double Vector_r,Vector_w,Vector_p,Vector_z;
	Vector_double Vector_temp_r,Vector_temp_w,Vector_temp_p,Vector_temp_z;
	Vector_double K,AT,Vector_Ax,Vector_rre;             //预处理矩阵
    Vector_r.assign(N,0);  Vector_temp_r.assign(N,0);
	Vector_w.assign(N,0);  Vector_temp_w.assign(N,0);
	Vector_p.assign(N,0);  Vector_temp_p.assign(N,0);
	Vector_z.assign(N,0);  Vector_temp_z.assign(N,0);
	x.assign(N,0);         K.assign(N,0);      AT.assign(N*N,0);   
	Vector_Ax.assign(N,0); Vector_rre.assign(N,0);
	int n,k;
	for(n=0;n<N;n++)
	{
		Vector_r[n]=b[n];
		Vector_temp_r[n]=b[n];
		K[n]=A[n*N+n];
	}
	double temp_p_pre=0;
	//////////////
	//    K
	//////////////
	////////////////////
	//    MAIN LOOP
	////////////////////
	for(k=1;k<Kiter;k++)
	{
		//////////////
		//
		//////////////
		for(n=0;n<N;n++)
		{
			Vector_w[n]=Vector_r[n]/K[n];
			Vector_temp_w[n]=Vector_temp_r[n]/K[n];
		}
		double temp_p_lat=VecSumVecMultVec(Vector_w,Vector_temp_w,N);
		if(temp_p_lat==0)
		{
			cout<<"problem occuring!"<<endl;
		}
		//////////////
		//
		//////////////
		if(k==1)
		{
			for(n=1;n<N;n++)
			{
				Vector_p[n]=Vector_w[n];
				Vector_temp_p[n]=Vector_temp_w[n];
			}
		}
		else
		{
		    double beta=temp_p_lat/temp_p_pre;
			for(n=0;n<N;n++)
			{
				Vector_p[n]=Vector_w[n]+beta*Vector_p[n];
				Vector_temp_p[n]=Vector_temp_w[n]+beta*Vector_temp_p[n];
			}
		}
		//////////////
		//
		//////////////
		temp_p_pre=temp_p_lat;
		VecMultMatVec(A,Vector_p,Vector_z,N,N);
		VecMultMatVec(AT,Vector_temp_p,Vector_temp_z,N,N);
        double alpha=temp_p_lat/VecSumVecMultVec(Vector_temp_p,Vector_z,N);
		for(n=0;n<N;n++)
		{
            x[n]=x[n]+alpha*Vector_p[n];
			Vector_r[n]=Vector_r[n]-alpha*Vector_z[n];
			Vector_temp_r[n]  -=alpha*Vector_temp_z[n];
		}
		//////////////
		//  break
		//////////////
		double temp_b=VecNorm2(b,N);
		VecMultMatVec(A,x,Vector_Ax,N,N);
		for(n=0;n<N;n++)
		{
			Vector_rre[n]=b[n]-Vector_Ax[n];
		}
		double temp_RRE=VecNorm2(Vector_rre,N);
		double RRE=temp_RRE/temp_b;
		cout<<"RRE["<<k<<"]------:"<<RRE<<endl;
		if(RRE<eps)
		{
			break;
		}
	}    
}
//////////////////////////////////////
////
//////////////////////////////////////
void CLinearEquations::RBI_CGSTAB(Vector_double &A,Vector_double &x,Vector_double &b,int N,int Kiter,double eps)
{
   ////////////////
   //
   ////////////////
   Vector_double Vector_r,Vector_p,Vector_v,Vector_s,Vector_t;
   Vector_double Vector_temp_r,Vector_Ax,Vector_rre;
   Vector_r.assign(N,0);     Vector_temp_r.assign(N,0);
   Vector_p.assign(N,0);     Vector_v.assign(N,0);
   Vector_s.assign(N,0);     Vector_t.assign(N,0);
   Vector_rre.assign(N,0);   Vector_Ax.assign(N,0);
   double temp_p_pre=1.0;
   double temp_p_lat=1.0;
   double alpha=1.0;
   double w=1.0;
   int n,k;
   for(n=0;n<N;n++)
   {
	   Vector_r[n]=b[n];
	   Vector_temp_r[n]=b[n];
   }
   /////////////////
   //
   /////////////////
   for(k=1;k<Kiter;k++)
   {
	   temp_p_lat=VecSumVecMultVec(Vector_r,Vector_temp_r,N);
	   double temp_p=temp_p_lat/temp_p_pre;
	   double beta=temp_p*alpha/w;
	   for(n=0;n<N;n++)
	   {
		   Vector_p[n]=Vector_r[n]+beta*(Vector_p[n]-w*Vector_v[n]);
	   }
	   ///////////////////
	   VecMultMatVec(A,Vector_p,Vector_v,N,N);
	   alpha=temp_p_lat/VecSumVecMultVec(Vector_temp_r,Vector_v,N);
	   for(n=0;n<N;n++)
	   {
		   Vector_s[n]=Vector_r[n]-alpha*Vector_v[n];
	   }
	   VecMultMatVec(A,Vector_s,Vector_t,N,N);
	   w=VecSumVecMultVec(Vector_t,Vector_s,N)/VecNorm(Vector_t,N);
	   for(n=0;n<N;n++)
	   {
		   x[n]=x[n]+alpha*Vector_p[n]+w*Vector_s[n];
	   }
	   /////////////////
	   //
	   /////////////////
	    double temp_b=VecNorm2(b,N);
        VecMultMatVec(A,x,Vector_Ax,N,N);
		for(n=0;n<N;n++)
		{
			Vector_rre[n]=b[n]-Vector_Ax[n];
		}
		double temp_RRE=VecNorm2(Vector_rre,N);
		double RRE=temp_RRE/temp_b;
		cout<<"RRE["<<k<<"]------:"<<RRE<<endl;
		if(RRE<eps)
		{
			break;
		}
		/////////////////
		//
		/////////////////
		temp_p_pre=temp_p_lat;
		for(n=0;n<N;n++)
		{
			Vector_r[n]=Vector_s[n]-w*Vector_t[n];
		}
   }
}
//////////////////////////////////////
////A[M][N] U[M][M] V[N][N] KA=max(M,N)+1
//////////////////////////////////////  
int CLinearEquations::SVD(double *a,int m,int n,double *u,double *v,double eps,int ka)
{ 
	int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh;
    double *s,*e,*w,*fg,*cs;
    s=new double[ka];   fg=new double[2];
	e=new double[ka];   cs=new double[2];
	w=new double[ka];
    it=60; k=n;
    if (m-1<n) k=m-1;
    l=m;
    if (n-2<m) l=n-2;
    if (l<0) l=0;
    ll=k;
    if (l>k) ll=l;
    if (ll>=1)
      { 
		for (kk=1; kk<=ll; kk++)
          { 
			if (kk<=k)
              { 
				d=0.0;
                for (i=kk; i<=m; i++)
                  { 
					ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];
				}
                s[kk-1]=sqrt(d);
                if (s[kk-1]!=0.0)
                  { 
					ix=(kk-1)*n+kk-1;
                    if (a[ix]!=0.0)
                      { 
						s[kk-1]=fabs(s[kk-1]);
                        if (a[ix]<0.0) s[kk-1]=-s[kk-1];
                      }
                    for (i=kk; i<=m; i++)
                      {
						iy=(i-1)*n+kk-1;
                        a[iy]=a[iy]/s[kk-1];
                      }
                    a[ix]=1.0+a[ix];
                  }
                s[kk-1]=-s[kk-1];
              }
            if (n>=kk+1)
              { 
				for (j=kk+1; j<=n; j++)
                  { 
					if ((kk<=k)&&(s[kk-1]!=0.0))
                      { 
						d=0.0;
                        for (i=kk; i<=m; i++)
                          { 
							ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+a[ix]*a[iy];
                          }
                        d=-d/a[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++)
                          { 
							ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            a[ix]=a[ix]+d*a[iy];
                          }
                      }
                    e[j-1]=a[(kk-1)*n+j-1];
                  }
              }
            if (kk<=k)
              { 
				for (i=kk; i<=m; i++)
                  {
					ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                  }
              }
            if (kk<=l)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                  { 
					if (e[kk]!=0.0)
                      {
						e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                  }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { 
					for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        { 
						  ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                  }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      {
		for (j=k+1; j<=nn; j++)
          { 
			for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      {
		for (ll=1; ll<=k; ll++)
          {
			kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { 
				if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    {
					  d=0.0;
                      for (i=kk; i<=m; i++)
                        {
						  ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { 
						  ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { 
					  ix=(i-1)*m+kk-1; 
					  u[ix]=-u[ix];
				  }
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { 
				for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { 
		kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { 
			for (j=kk+1; j<=n; j++)
              { 
				d=0.0;
                for (i=kk+1; i<=n; i++)
                  {
					ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  {
					ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
      { 
		if (mm==0)
          { 
			ppp(a,e,s,v,m,n);
            delete(s); delete(e); delete(w); return(1);
          }
        if (it==0)
          { 
			ppp(a,e,s,v,m,n);
            delete(s); delete(e); delete(w); return(-1);
          }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          {
		    d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
          }
        if (kk==mm-1)
          { 
			kk=kk+1;
            if (s[kk-1]<0.0)
              {
				s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { 
					ix=(i-1)*n+kk-1; 
					v[ix]=-v[ix];
				}
              }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
              {
				d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { 
					  ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { 
					  ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
              }
            it=60;
            mm=mm-1;
          }
        else
          {
			ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              {
				d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { 
				kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { 
					shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  {
					sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { 
						  ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { 
							ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { 
						i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          {
							fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            {
							  ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { 
						fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { 
							  ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
    return(1);
  }
void CLinearEquations::ppp(double *a,double *e,double *s,double *v,int m,int n)
{ 
	int i,j,p,q;
    double d;
    if (m>=n) i=n;
    else i=m;
    for (j=1; j<=i-1; j++)
      {
		a[(j-1)*n+j-1]=s[j-1];
        a[(j-1)*n+j]=e[j-1];
      }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      {
		p=(i-1)*n+j-1; q=(j-1)*n+i-1;
        d=v[p]; v[p]=v[q]; v[q]=d;
      }
    return;
}
//////////////////////////////////////////////////////////////////////////
void CLinearEquations::sss(double *fg,double *cs)
{   
	double r,d;
    if ((fabs(fg[0])+fabs(fg[1]))==0.0)
      { 
		cs[0]=1.0; cs[1]=0.0; d=0.0;}
    else 
      {
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
          { d=fabs(d);
            if (fg[0]<0.0) d=-d;
          }
        if (fabs(fg[1])>=fabs(fg[0]))
          { 
			d=fabs(d);
            if (fg[1]<0.0) d=-d;
          }
        cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
    else
      if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
}
//////////////////////////////////////////////////////////////////////////
//a[M][N],AA[N][M](广异逆矩阵),U[M][M],V[N][N],ka=max(M,N)+1
//////////////////////////////////////////////////////////////////////////
int CLinearEquations::Ginv(double *a,int m,int n,double *aa,double eps,double *u,double *v,int ka)
{ 
	int i,j,k,l,t,p,q,f;
    i=SVD(a,m,n,u,v,eps,ka);
    if (i<0) return(-1);
    j=n;
    if (m<n) j=m;
    j=j-1;
    k=0;
    while ((k<=j)&&(a[k*n+k]!=0.0)) k=k+1;
    k=k-1;
    for (i=0; i<=n-1; i++)
    for (j=0; j<=m-1; j++)
      { t=i*m+j; aa[t]=0.0;
        for (l=0; l<=k; l++)
          { f=l*n+i; p=j*m+l; q=l*n+l;
            aa[t]=aa[t]+v[f]*u[p]/a[q];
          }
      }
    return(1);
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
int CLinearEquations::SVD_LQ(double *a,int m,int n,double *b,double *x,double *aa,double eps,double *u,double *v)
{ 
	int ka;
	if(m>n)
	{
		ka=m+1;
	}
	else
	{
		ka=n+1;
	}
	int i,j;
    i=Ginv(a,m,n,aa,eps,u,v,ka);
    if (i<0) return(-1);
    for (i=0; i<=n-1; i++)
      { 
		x[i]=0.0;
        for (j=0; j<=m-1; j++)
          x[i]=x[i]+aa[i*m+j]*b[j];
      }
    return(1);
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
void CLinearEquations::VSVD_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int M,int N,double eps)
{
   Vector_double At,bt,xt;
   int Mt=2*M;
   int Nt=2*N;
   At.assign((2*N)*(2*M),0);
   bt.assign((2*M),0);
   xt.assign(2*N,0);
   ExchangeComToReal(A,b,M,N,At,bt,Mt,Nt);
   RSVD_LQ(At,xt,bt,Mt,Nt,eps);
   for(int i=0;i<N;i++)
   {
	   x[i]=xt[i]+J*xt[i+N];
   }
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
void CLinearEquations::RSVD_LQ(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N,double eps)
{
   double *a,*aa,*u,*v,*bs,*xs;
   a=new double[M*N];
   aa=new double[M*N];
   u=new double[M*M];
   v=new double[N*N];
   bs=new double[M];
   xs=new double[N];
   int i,j;
   for(i=0;i<M;i++)
   {
	   for(j=0;j<N;j++)
	   {
		   a[i*N+j]=A[i*N+j];
	   }
	   bs[i]=b[i];
   }
   SVD_LQ(a,M,N,bs,xs,aa,eps,u,v);
   for(i=0;i<N;i++)
   {
	   x[i]=xs[i];
   }
   delete a;
   delete aa;
   delete u;
   delete v;
   delete bs;
   delete xs;
}
//////////////////////////////////////////////////////////////////////////
//A[M][N],Q[M][M],R[M][N],R=A
/////////////////////////////////////////////////////////////////////////
int CLinearEquations::QR(double *a,int m,int n,double *q)
{
	int i,j,k,l,nn,p,jj;
    double u,alpha,w,t;
    if (m<n)
      { 
		cout<<"fail\n";
		return(0);}
    for (i=0; i<=m-1; i++)
    for (j=0; j<=m-1; j++)
      { l=i*m+j; q[l]=0.0;
        if (i==j) q[l]=1.0;
      }
    nn=n;
    if (m==n) nn=m-1;
    for (k=0; k<=nn-1; k++)
      { u=0.0; l=k*n+k;
        for (i=k; i<=m-1; i++)
          { w=fabs(a[i*n+k]);
            if (w>u) u=w;
          }
        alpha=0.0;
        for (i=k; i<=m-1; i++)
          { t=a[i*n+k]/u; alpha=alpha+t*t;}
        if (a[l]>0.0) u=-u;
        alpha=u*sqrt(alpha);
        if (fabs(alpha)+1.0==1.0)
          { cout<<"fail\n"; return(0);}
        u=sqrt(2.0*alpha*(alpha-a[l]));
        if ((u+1.0)!=1.0)
          { a[l]=(a[l]-alpha)/u;
            for (i=k+1; i<=m-1; i++)
              { p=i*n+k; a[p]=a[p]/u;}
            for (j=0; j<=m-1; j++)
              { t=0.0;
                for (jj=k; jj<=m-1; jj++)
                  t=t+a[jj*n+k]*q[jj*m+j];
                for (i=k; i<=m-1; i++)
                  { p=i*m+j; q[p]=q[p]-2.0*t*a[i*n+k];}
              }
            for (j=k+1; j<=n-1; j++)
              { t=0.0;
                for (jj=k; jj<=m-1; jj++)
                  t=t+a[jj*n+k]*a[jj*n+j];
                for (i=k; i<=m-1; i++)
                  { p=i*n+j; a[p]=a[p]-2.0*t*a[i*n+k];}
              }
            a[l]=alpha;
            for (i=k+1; i<=m-1; i++)
              a[i*n+k]=0.0;
          }
      }
    for (i=0; i<=m-2; i++)
    for (j=i+1; j<=m-1;j++)
      { p=i*m+j; l=j*m+i;
        t=q[p]; q[p]=q[l]; q[l]=t;
      }
    return(1);
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
int CLinearEquations::QR_LQ(double *a,int m,int n,double *b,double *q)
{ 
	int i,j;
    double d,*c;
	c=new double[n];
    i=QR(a,m,n,q);
    if (i==0) { delete(c); return(0);}
    for (i=0; i<=n-1; i++)
      { d=0.0;
        for (j=0; j<=m-1; j++)
          d=d+q[j*m+i]*b[j];
        c[i]=d;
      }
    b[n-1]=c[n-1]/a[n*n-1];
    for (i=n-2; i>=0; i--)
      { d=0.0;
        for (j=i+1; j<=n-1; j++)
          d=d+a[i*n+j]*b[j];
        b[i]=(c[i]-d)/a[i*n+i];
      }
    delete(c); return(1);
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
void CLinearEquations::VQR_LQ(Vector_complex &A,Vector_complex &x,Vector_complex &b,int M,int N)
{
   Vector_double At,bt,xt;
   int Mt=2*M;
   int Nt=2*N;
   At.assign((2*N)*(2*M),0);
   bt.assign((2*M),0);
   xt.assign(2*N,0);
   ExchangeComToReal(A,b,M,N,At,bt,Mt,Nt);
   RQR_LQ(At,xt,bt,Mt,Nt);
   for(int i=0;i<N;i++)
   {
	   x[i]=xt[i]+J*xt[i+N];
   }
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
void CLinearEquations::RQR_LQ(Vector_double &A,Vector_double &x,Vector_double &b,int M,int N)
{
   double *a,*bqr,*q;
   a=new double[M*N];
   bqr=new double[M];
   q=new double[M*M];
   int i,j;
   for(i=0;i<M;i++)
   {
	   for(j=0;j<N;j++)
	   {
		   a[i*N+j]=A[i*N+j];
	   }
	   bqr[i]=b[i];
   }
   QR_LQ(a,M,N,bqr,q);
   for(i=0;i<N;i++)
   {
	   x[i]=bqr[i];
   }
   delete a;
   delete bqr;
   delete q;
}
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////







