
#pragma warning(disable:4786)
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <CString>
#include <math.h>
#include <process.h>
#include <utility>
#include <time.h>

typedef std::complex<double> Complex;
using namespace std;

typedef vector <Complex> Vec_Com;
typedef vector <double>  Vec_dou;
 
const double PI = 3.1415926535897932384626433832795;
const Complex i_(0.0,1.0);
//真空中的介电常数和磁导率
const double ebslun0 = 1./(36.*PI)*1.e-9;
const double miu0 = (4.*PI)*1.e-7;

#define  maxIndex  10000L

//#define  maxPrimeFactor        1007//最大素因子经过修改比较大
#define  maxPrimeFactor        37//最大素因子经过修改比较大
#define  maxPrimeFactorDiv2    (maxPrimeFactor+1)/2
#define  maxFactorCount        20

static double  c3_1 = -1.5000000000000E+00;  /*  c3_1 = cos(2*PI/3)-1;          */
static double  c3_2 =  8.6602540378444E-01;  /*  c3_2 = sin(2*PI/3);            */
                                          
static double  u5   =  1.2566370614359E+00;  /*  u5   = 2*PI/5;                 */
static double  c5_1 = -1.2500000000000E+00;  /*  c5_1 = (cos(u5)+cos(2*u5))/2-1;*/
static double  c5_2 =  5.5901699437495E-01;  /*  c5_2 = (cos(u5)-cos(2*u5))/2;  */
static double  c5_3 = -9.5105651629515E-01;  /*  c5_3 = -sin(u5);               */
static double  c5_4 = -1.5388417685876E+00;  /*  c5_4 = -(sin(u5)+sin(2*u5));   */
static double  c5_5 =  3.6327126400268E-01;  /*  c5_5 = (sin(u5)-sin(2*u5));    */
static double  c8   =  7.0710678118655E-01;  /*  c8 = 1/sqrt(2);    */

static int      groupOffset,dataOffset,blockOffset,adr;
static int      groupNo,dataNo,blockNo,twNo;
static double   omega, tw_re,tw_im;
static double   twiddleRe[maxPrimeFactor], twiddleIm[maxPrimeFactor],
                trigRe[maxPrimeFactor], trigIm[maxPrimeFactor],
                zRe[maxPrimeFactor], zIm[maxPrimeFactor];
static double   vRe[maxPrimeFactorDiv2], vIm[maxPrimeFactorDiv2];
static double   wRe[maxPrimeFactorDiv2], wIm[maxPrimeFactorDiv2];

ofstream outresult("Test_Result.txt",ios::trunc);

//存储二维粗糙面矩阵元素的类
class Zmn_Cell
{
public:
	//int Zs_bw;//半带宽
	int Zmn_i; //行数
	int Zmn_j; //列数
	Complex Zmn_z;//矩阵元素
/*	Zmn_Cell()
	{
		Zmn_i = 1;
		Zmn_j = 1;
		Zmn_z = 1.;
	}*/
};


typedef pair<int,int> Pair;//pair
typedef vector<double> Vec_Dou;//元素为双整型实数的向量
typedef vector<Complex> Vec_Com;//元素为复数的向量
typedef vector<Zmn_Cell> Vec_Zmn_Cell;//元素为Zmn_Cell的向量
typedef map<Pair,Vec_Com* > MapStruct;//map
typedef map<Pair,Vec_Com* >::value_type ValType;//map的计数器


class Mathe_Functions
{
public:
	//////////////Complex////////////////////////////////
	//向量的2-范数的平方，如果严格取2-范数的话，对结果开方就可以了
	double L2_Norm(Vec_Com &x);
	double L2_Norm(Complex *x,int n);
	double L2_Norm(Vec_Com x[],int m);
	double L2_Norm_Inf(Vec_Com &x);//向量的无穷范数
	double L2_Norm_Inf(Complex *x,int n);


	//将向量b的值赋给向量a
    void Assign_Vecor(Vec_Com &a,Vec_Com &b);
	void Assign_Vecor(Complex *a,Complex *b,int m);
	void Assign_Vecor(Vec_Com a[6],Vec_Com b[6]);

	//复数贝塞尔函数
    Complex IBessel1(int,Complex);//复数的第一类贝塞尔函数
    Complex IBessel2(int,Complex);//复数的第二类贝塞尔函数
    //第一、二类汉克尔函数  
    Complex Hankel(int,int,Complex);
    Complex Hankel(int,int,double);
	Complex *Unite_Matrix(Complex *a,Complex *b,int n);//将两个n维的矩阵合并成一个2*n维的矩阵
	Complex *Plus_Matrix(Complex *a,Complex *b,int n);//将两个n维的矩阵相加得到一个n维的矩阵
	
	//向量VT的共轭与向量V的点积
	Complex Multi_vt_v(Vec_Com &vt,Vec_Com &v);
	Complex Multi_vt_v(Complex *vt,Complex *v,int m);
	Complex Multi_vt_v(Vec_Com vt[6],Vec_Com v[6]);

	//阶乘:factorial
	double factor(int);
	//复数的开方
	Complex Isqrt(Complex x);
	//以复数格式输出复数
	void display(Complex re);
	
	//几种插值函数
	Complex *Linear_Interpolate(Complex *M,double *XC,double *XD,int ncg,int ndg);//线形插值
	Complex Lagrange(double *x,Complex *y,int n,double t);//复数的Lagrange插值
	Complex *Lagrange_Interpolate(Complex *y,double *x,double *t,int n,int nt);//n个复数的Lagrange插值
	double Lagrange(double *x,double *y,int n,double t);//实数的Lagrange插值
    double *Lagrange_Interpolate(double *y,double *x,double *t,int n,int nt);//n个实数的Lagrange插值
	Complex lg3(Vec_Dou &x,Complex *y,int n,double t);//一元三点插值
    Complex *Lg3_Interpolate(Complex *y,Vec_Dou &x,Vec_Dou &t,int n,int nt);
	Complex pqs(double *x,Complex *y,int n,double t);//连分式插值
	Complex *pqs_Interpolate(Complex *y,double *x,double *t,int n,int nt);
    Complex hmt(double *x,Complex *y,double *dy,int n,double t);//埃尔米特插值
	Complex *hmt_Interpolate(Complex *y,double *x,double *dy,double *t,int n,int nt);
    Complex atk(double *x,Complex *y,int n,double t,double eps);//埃特金逐步插值
	Complex *atk_Interpolate(Complex *y,double *x,double *t,int n,int nt,double eps);
	
	//n个数平均为ncg个数的向量
	Complex *Average(Complex *M,int n1,int ncg);
	
	//快速付立叶变换及其应用
	//任意维数的1-D FFT程序，当Bool为-1时为逆FFT，为其它值时为FFT
	void fft(int n, Complex x[], Complex y[], int Bool);
    
	void transTableSetup(int sofar[], int actual[], int remain[], int *nFact, int *nPoints);
    void permute(int nPoint, int nFact, int fact[], int remain[], double xRe[], double xIm[], double yRe[], double yIm[]);
    void twiddleTransf(int sofarRadix, int radix, int remainRadix, double yRe[], double yIm[]);
    void factorize(int n, int *nFact, int fact[]);
    void initTrig(int radix);
    void fft_4(double aRe[], double aIm[]);
    void fft_5(double aRe[], double aIm[]);
    void fft_8();
    void fft_10();
    void fft_odd(int radix);
	
	//1-D FFT 的应用
	//应用1-D FFT 快速计算对称TOEPLITZ矩阵向量积，TOEPLITZ矩阵只需存储一行或一列
	Complex *fft_toematrix_vector(Complex *g,Complex *x,int n2);
	
	//任意维数的2D-FFT：x，y都为n*n维
	void fft2(int n, Complex x[], Complex y[], int Bool);
	//应用2-D FFT 快速计算二维对称TOEPLITZ矩阵向量积，TOEPLITZ矩阵只需存储扩展后的一行或一列
	void Multi_BTTB_Vec(Vec_Com &Zfs,Vec_Com &vec,int nxy,Vec_Com &result);
	//x，y都为nx*ny维，将x变换为nx*ny维的2D-FFT
	void fft2(int nx, int ny, Complex x[], Complex y[], int Bool);
};
Complex Mathe_Functions::Isqrt(Complex x)
{
	Complex result;
	double re,im;
	re=sqrt((sqrt(x.real()*x.real()+x.imag()*x.imag())+x.real())/2.0);
	if(x.imag() >= 0.0)
		im=sqrt((sqrt(x.real()*x.real()+x.imag()*x.imag())-x.real())/2.0);
	else
		im=-sqrt((sqrt(x.real()*x.real()+x.imag()*x.imag())-x.real())/2.0);
	result = re + i_*im;
	return result;
}
void Mathe_Functions::display(Complex re)
{
	outresult<<re.real();
	if(re.imag() >= 0.0)
	{
	    outresult<<"+"<<re.imag()<<"i";
	}
	else
	{
	    outresult<<re.imag()<<"i";
	}
}
double Mathe_Functions::factor(int mn)
{
	double re;
	int i;
	long double Value_fact;
	Value_fact = 1;
	if (mn == 0)
	{
		re = 1.0;
	}
	else
	{
		for ( i=1; i<=mn; i++)
		{
			Value_fact = Value_fact * i;
		}
		re = double(Value_fact);
	}
	
	return re;
}
double Mathe_Functions::L2_Norm_Inf(Vec_Com &x)
{//向量的无穷范数：向量中绝对值最大的元素
	double result,tem;
	int i;
	int n=x.size();
	result = abs(x[0]);	
	for( i=1; i<n; i++ )
	{
		tem = abs(x[i]);
		if ( result < tem )
		{
			result = tem;
		}
	}
	return result;
}
double Mathe_Functions::L2_Norm_Inf(Complex *x,int n)
{//向量的无穷范数：向量中绝对值最大的元素
	double result,tem;
	int i;
	result = abs(x[0]);
	for( i=1; i<n; i++ )
	{
		tem = abs(x[i]);
		if ( result < tem )
		{
			result = tem;
		}
	}
	return result;
}
double Mathe_Functions::L2_Norm(Vec_Com &x)
{//向量的2范数是平方和的开方
	double result;
	int i;
	int n=x.size();
	result=0.0;
	for(i=0;i<n;i++)
	{
		result=result+norm(x[i]);//
	}
	return result;
}
double Mathe_Functions::L2_Norm(Vec_Com x[],int m)
{//向量的2范数是平方和的开方
	double result;
	int k,i;
	int n = x[0].size();
	result=0.0;
	for ( k=0; k<m; k++ )
	{
		for(i=0;i<n;i++)
		{
			result += norm(x[k][i]);//
		}
	}	
	return result;
}

double Mathe_Functions::L2_Norm(Complex *x,int n)
{//向量的2范数是平方和的开方
	double result;
	int i;
	result=0.0;
	for(i=0;i<n;i++)
	{
		result=result+norm(x[i]);//
	}
	return result;
}
Complex *Mathe_Functions::Unite_Matrix(Complex *M1,Complex *M2,int n)
{
	Complex *result=new Complex[2*n];
	for(int i=0;i<n;i++)
	{
		result[i]=M1[i];
		result[i+n]=M2[i];
	}
	return result;
}
Complex *Mathe_Functions::Plus_Matrix(Complex *M1,Complex *M2,int n)
{
	Complex *result=new Complex[n];
	for(int i=0;i<n;i++)
	{
		result[i]=M1[i]+M2[i];
	}
	return result;
}
void Mathe_Functions::Assign_Vecor(Vec_Com &a,Vec_Com &b)
{
	int i;
	int m=a.size();
	for(i=0; i<m; i++)
		b[i]=a[i];
}
void Mathe_Functions::Assign_Vecor(Vec_Com a[6],Vec_Com b[6])
{
	int k,i;
	int m=a[0].size();
	for ( k=0; k<6; k++ )
	{
		for(i=0; i<m; i++)
			b[k][i]=a[k][i];
	}	
}

void Mathe_Functions::Assign_Vecor(Complex *a,Complex *b,int m)
{
	int i;
	for(i=0; i<m; i++)
		b[i]=a[i];
}
Complex Mathe_Functions::Multi_vt_v(Vec_Com &vt,Vec_Com &v)
{
	Complex re;
	re=0.0;
	int m=vt.size();
	for (int i=0; i<m; i++)
	{
		re+=conj(vt[i])*v[i];
	}
	return re;
	//cout<<re<<endl;
}
Complex Mathe_Functions::Multi_vt_v(Vec_Com vt[6],Vec_Com v[6])
{
	Complex re;
	int m=vt[0].size();
	for ( int k=0; k<6; k++ )
	{
		for (int i=0; i<m; i++)
		{
			re += conj(vt[k][i])*v[k][i];
		}
	}	
	return re;
}

Complex Mathe_Functions::Multi_vt_v(Complex *vt,Complex *v,int m)
{
	Complex re;
	re=0.0;
	for (int i=0; i<m; i++)
	{
		re+=conj(vt[i])*v[i];
	}
	return re;
}
//-----------------------------------------------------------------------------
//                  复数贝塞尔及汉克尔函数
//-----------------------------------------------------------------------------
Complex Mathe_Functions::IBessel1(int nj,Complex x)
{
	static double a[6]={2.2499997,1.2656208,0.3163866,0.0444479,0.0039444,0.0002100};
	static double c[6]={0.00000077,0.00552740,0.00009512,0.00137237,0.00072805,0.00014476};
	static double d[6]={0.04166397,0.00003954,0.00262573,0.00054125,0.00029333,0.00013558};
	static double e[6]={0.56249985,0.21093573,0.03954289,0.00443319,0.00031761,0.00001109};
	static double g[6]={0.00000156,0.01659667,0.00017105,0.00249511,0.00113653,0.00020033};
	static double h[6]={0.12499612,0.00005650,0.00637879,0.00074348,0.00079824,0.00029166};
	Complex re;
	Complex uv,x1,f0,c0;
	if(nj == 0)
	{
		if(abs(x) <= 3.0)
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
		if(abs(x) <= 3.0)
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
Complex Mathe_Functions::IBessel2(int nj,Complex x)
{
	static double b[6]={0.60559366,0.74350384,0.25300117,0.04261214,0.00427916,0.00024846};
	static double c[6]={0.00000077,0.00552740,0.00009512,0.00137237,0.00072805,0.00014476};
	static double d[6]={0.04166397,0.00003954,0.00262573,0.00054125,0.00029333,0.00013558};
	static double f[6]={0.2212091,2.1682709,1.3164827,0.3123951,0.0400976,0.0027873};
	static double g[6]={0.00000156,0.01659667,0.00017105,0.00249511,0.00113653,0.00020033};
	static double h[6]={0.12499612,0.00005650,0.00637879,0.00074348,0.00079824,0.00029166};
	Complex re;
	Complex uv,x1,f0,c0;
	if(nj == 0)
	{
		if(abs(x) <= 3.0)
		{
			
			re=(2.0/PI)*log(x/2.0)*IBessel1(0,x)+0.36746691+b[0]*pow(x/3.0,2)-b[1]*pow(x/3.0,4)
				+b[2]*pow(x/3.0,6)-b[3]*pow(x/3.0,8)+b[4]*pow(x/3.0,10)-b[5]*pow(x/3.0,12);
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
		if(abs(x) <= 3.0)
		{
			re=(2.0/PI)*x*log(x/2.0)*IBessel1(1,x)-0.6366198+f[0]*pow(x/3.0,2)+f[1]*pow(x/3.0,4)
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
Complex Mathe_Functions::Hankel(int nj,int ku,Complex x)
//nj表示阶数，ku表示第一类或第二类汉克尔函数
{
	Complex re;
	Complex a,b;
	a=IBessel1(nj,x);
	b=IBessel2(nj,x);
	switch(ku) 
	{
	case 1: {
		re=a+i_*b;
		//re.Real=a.Real-b.Imag;
	    //re.Imag=a.Imag+b.Real;
		break;
			}
	case 2: {
		re=a-i_*b;
		//re.Real=a.Real+b.Imag;
	    //re.Imag=a.Imag-b.Real;
		break;
			}
	}

	return re;
}
Complex Mathe_Functions::Hankel(int nj,int ku,double x)
//nj表示阶数，ku表示第一类或第二类汉克尔函数
{
	Complex re;
	switch(ku) 
	{
	case 1: {
		if(nj == 0)
		{
			re=_j0(x)+i_*_y0(x);
		}
		else if(nj == 1)
		{
			re=_j1(x)+i_*_y1(x);
		}
		else
		{
			re=_jn(nj,x)+i_*_yn(nj,x);
		}
		break;
			}
	case 2: {
		if(nj == 0)
		{
			re=_j0(x)-i_*_y0(x);
		}
		else if(nj == 1)
		{
			re=_j1(x)-i_*_y1(x);
		}
		else
		{
			re=_jn(nj,x)-i_*_yn(nj,x);
		}
		break;
			}
	}
	return re;
}
//------------------------------------------------------------------------------------
//====================================================================================
Complex *Mathe_Functions::Linear_Interpolate(Complex *M,double *XC,double *XD,int ncg,int ndg)
{
	//该程序就适用于n1为奇数，即3，5，7...
	Complex *re=new Complex[ndg];
	int tn;
	int n1=ndg/ncg;
	//vector<double> a;
	//map<double >
	for(int i=0;i<ndg;i++)
	{
		if(i <= (n1-1)/2)
			re[i]=M[0]+(M[1]-M[0])/(XC[1]-XC[0])*(XD[i]-XC[0]);
		else if(i > ndg-1-(n1-1)/2)
			re[i]=M[ncg-2]+(M[ncg-1]-M[ncg-2])/(XC[ncg-1]-XC[ncg-2])*(XD[i]-XC[ncg-2]);
		else
		{
			tn=(i+1)/n1;
			re[i]=M[tn-1]+(M[tn]-M[tn-1])/(XC[tn]-XC[tn-1])*(XD[i]-XC[tn-1]);
		}

	}
	return re;
}
Complex Mathe_Functions::Lagrange(double *x,Complex *y,int n,double t)
{ 
	int i,j,k,m;
    Complex z,s;
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0];return(z);}
    if (n==2)
      { z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
        return(z);
      }
    i=0;
    while ((x[i]<t)&&(i<n)) i=i+1;
    k=i-4;
    if (k<0) k=0;
    m=i+3;
    if (m>n-1) m=n-1;
    for (i=k;i<=m;i++)
      { s=1.0;
        for (j=k;j<=m;j++)
          if (j!=i) s=s*(t-x[j])/(x[i]-x[j]);
        z=z+s*y[i];
      }
    return(z);
}
Complex *Mathe_Functions::Lagrange_Interpolate(Complex *y,double *x,double *t,int n,int nt)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	Complex *re=new Complex[nt];
	for(int i=0;i<nt;i++)
		re[i]=Lagrange(x,y,n,t[i]);
	return re;
}
double Mathe_Functions::Lagrange(double *x,double *y,int n,double t)
{ 
	int i,j,k,m;
    double z,s;
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0];return(z);}
    if (n==2)
      { z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
        return(z);
      }
    i=0;
    while ((x[i]<t)&&(i<n)) i=i+1;
    k=i-4;
    if (k<0) k=0;
    m=i+3;
    if (m>n-1) m=n-1;
    for (i=k;i<=m;i++)
      { s=1.0;
        for (j=k;j<=m;j++)
          if (j!=i) s=s*(t-x[j])/(x[i]-x[j]);
        z=z+s*y[i];
      }
    return(z);
}
double *Mathe_Functions::Lagrange_Interpolate(double *y,double *x,double *t,int n,int nt)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	double *re=new double[nt];
	for(int i=0;i<nt;i++)
		re[i]=Lagrange(x,y,n,t[i]);
	return re;
}
Complex Mathe_Functions::lg3(Vec_Dou &x,Complex *y,int n,double t)
{ 
	int i,j,k,m;
    Complex z,s;
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0]; return(z);}
    if (n==2)
      { z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
        return(z);
      }
    if (t<=x[1]) { k=0; m=2;}
    else if (t>=x[n-2]) { k=n-3; m=n-1;}
    else
      { k=1; m=n;
        while (m-k!=1)
          { i=(k+m)/2;
            if (t<x[i-1]) m=i;
            else k=i;
          }
        k=k-1; m=m-1;
        if (fabs(t-x[k])<fabs(t-x[m])) k=k-1;
        else m=m+1;
      }
    z=0.0;
    for (i=k;i<=m;i++)
      { s=1.0;
        for (j=k;j<=m;j++)
          if (j!=i) s=s*(t-x[j])/(x[i]-x[j]);
        z=z+s*y[i];
      }
    return(z);
}
Complex *Mathe_Functions::Lg3_Interpolate(Complex *y,Vec_Dou &x,Vec_Dou &t,int n,int nt)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	Complex *re=new Complex[nt];
	for(int i=0;i<nt;i++)
	{
		re[i]=lg3(x,y,n,t[i]);
		//cout<<t[i]<<endl;
	}
	return re;
}
Complex Mathe_Functions::pqs(double *x,Complex *y,int n,double t)//连分式插值
{ 
	int i,j,k,m,l;
	Complex z,h,b[8];
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0]; return(z);}
    if (n<=8) { k=0; m=n;}
    else if (t<x[4]) { k=0; m=8;}
    else if (t>x[n-5]) { k=n-8; m=8;}
    else
      { k=1; j=n;
        while (j-k!=1)
          { i=(k+j)/2;
            if (t<x[i-1]) j=i;
            else k=i;
          }
        k=k-4; m=8;
      }
    b[0]=y[k];
    for (i=2;i<=m;i++)
      { h=y[i+k-1]; l=0; j=1;
        while ((l==0)&&(j<=i-1))
            { if (sqrt(norm(h-b[j-1]))+1.0==1.0) l=1;//??
              else h=(x[i+k-1]-x[j+k-1])/(h-b[j-1]);
              j=j+1;
            }
        b[i-1]=h;
        if (l!=0) b[i-1]=1.0e+35;
      }
    z=b[m-1];
    for (i=m-1;i>=1;i--) z=b[i-1]+(t-x[i+k-1])/z;
    return(z);
}
Complex *Mathe_Functions::pqs_Interpolate(Complex *y,double *x,double *t,int n,int nt)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	Complex *re=new Complex[nt];
	for(int i=0;i<nt;i++)
		re[i]=pqs(x,y,n,t[i]);
	return re;
}
Complex Mathe_Functions::hmt(double *x,Complex *y,double *dy,int n,double t)//埃尔米特插值
{
	int i,j;
	Complex z,p,q,s;
    z=0.0;
    for (i=1;i<=n;i++)
      { s=1.0;
        for (j=1;j<=n;j++)
          if (j!=i) s=s*(t-x[j-1])/(x[i-1]-x[j-1]);
        s=s*s;
        p=0.0;
        for (j=1;j<=n;j++)
          if (j!=i) p=p+1.0/(x[i-1]-x[j-1]);
        q=y[i-1]+(t-x[i-1])*(dy[i-1]-2.0*y[i-1]*p);
        z=z+q*s;
      }
    return(z);
}
Complex *Mathe_Functions::hmt_Interpolate(Complex *y,double *x,double *dy,double *t,int n,int nt)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	Complex *re=new Complex[nt];
	for(int i=0;i<nt;i++)
		re[i]=hmt(x,y,dy,n,t[i]);
	return re;
}
Complex Mathe_Functions::atk(double *x,Complex *y,int n,double t,double eps)//埃特金逐步插值
{
	int i,j,k,m,l;
	Complex z,xx[10],yy[10];
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0]; return(z);}
    m=10;
    if (m>n) m=n;
    if (t<=x[0]) k=1;
    else if (t>=x[n-1]) k=n;
    else
      { k=1; j=n;
        while ((k-j!=1)&&(k-j!=-1))
          { l=(k+j)/2;
            if (t<x[l-1]) j=l;
            else k=l;
          }
        if (fabs(t-x[l-1])>fabs(t-x[j-1])) k=j;
      }
    j=1; l=0;
    for (i=1;i<=m;i++)
      { k=k+j*l;
        if ((k<1)||(k>n))
          { l=l+1; j=-j; k=k+j*l;}
        xx[i-1]=x[k-1]; yy[i-1]=y[k-1];
        l=l+1; j=-j;
      }
    i=0;
    do
      { i=i+1; z=yy[i];
        for (j=0;j<=i-1;j++)
          z=yy[j]+(t-xx[j])*(yy[j]-z)/(xx[j]-xx[i]);
        yy[i]=z;
      }
    while ((i!=m-1)&&(sqrt(norm(yy[i]-yy[i-1]))>eps));
    return(z);
}
Complex *Mathe_Functions::atk_Interpolate(Complex *y,double *x,double *t,int n,int nt,double eps)
{//y[]为已知向量，x[]为已知，t[]为未知，n为已知长度，nt为未知长度,nt可为任意长度
	Complex *re=new Complex[nt];
	for(int i=0;i<nt;i++)
		re[i]=atk(x,y,n,t[i],eps);
	return re;
}
Complex *Mathe_Functions::Average(Complex *M,int n1,int ncg)
{
	int i,tt;
	Complex *re=new Complex[ncg];
	
	i=0;
    for(tt=0;tt<ncg;tt++)
	{	
		while(i/n1 == tt)
		{
			re[tt]+=M[i];
			i++;
		}
		re[tt]=re[tt]/double(n1);
	}
	return re;
}
Complex *Mathe_Functions::fft_toematrix_vector(Complex *g,Complex *x,int n2)
{
	//Toeplitz 矩阵为对称矩阵，则只取一行元素即可
	int i;
	Complex *x1=new Complex[2*n2];
	Complex *x2=new Complex[2*n2];
	Complex *yk=new Complex[2*n2];
	Complex *re=new Complex[n2];
	for(i=0;i<n2;i++)//右端向量的零填充
	{
		x1[i]=x[i];
	}
	for(i=0;i<2*n2;i++)//Toeplitz矩阵的变换
	{
		if(i < n2)
		{
			x2[i]=g[i];
		}
		else if(i == n2) {}
		else
		{
			x2[i]=g[n2-(i-n2)];
		}
	}
	fft(2*n2,x1,x1,1);
	fft(2*n2,x2,x2,1);

	for(i=0;i<2*n2;i++)
	{
		yk[i]=x1[i]*x2[i];
	}
	//逆离散傅里叶变换，从1到n2，而只需要1到n2/2
	fft(2*n2,yk,yk,-1);
    
	for(i=0;i<n2;i++)
	{
		re[i]=yk[i];
	}
	delete [] x1;
	delete [] x2;
	delete [] yk;
	return re;
/*	//L. Tsang 书中的算法
    int i;
	Complex *x1=new Complex[2*n2];
	Complex *x2=new Complex[2*n2];
	Complex *yk=new Complex[2*n2];
	Complex *re=new Complex[n2];
	for(i=0;i<2*n2;i++)//右端向量的零填充
	{
		if(i < n2)
		{
			x1[i]=x[i];
		}
		else
		{
			x1[i]=0.0;
		}
	}
	for(i=0;i<2*n2;i++)//Toeplitz矩阵的变换
	{
		if(i < n2)
		{
			x2[i]=g[i];
		}
		else if(i == n2)
		{
			x2[i]=0.0;
		}
		else
		{
			x2[i]=g[2*n2-i];
		}
	}
	fft(2*n2,x1,x1,1);
	fft(2*n2,x2,x2,1);
	//fft(x1,2*n2,1);
	//fft(x2,2*n2,1);
	for(i=0;i<2*n2;i++)
	{
		yk[i]=x1[i]*x2[i];
		//yk[i].Real=x1[i].Real*x2[i].Real-x1[i].Imag*x2[i].Imag;
		//yk[i].Imag=x1[i].Real*x2[i].Imag+x2[i].Real*x1[i].Imag;
	}
	fft(2*n2,yk,yk,-1);
    //fft(yk,2*n2,-1);//逆离散傅里叶变换，从1到n2，而只需要1到n2/2
	for(i=0;i<n2;i++)
	{
		re[i]=yk[i];
	}
	delete [] x1;
	delete [] x2;
	delete [] yk;
	return re;*/
}
void Mathe_Functions::Multi_BTTB_Vec(Vec_Com &Zfs,Vec_Com &vec,int nxy,Vec_Com &result)
{
	int i,j,u;
	//Construct C_Zfs_First_Column and vec
	Complex *M_Ch=new Complex[4*nxy*nxy];
	Complex *M_Cv=new Complex[4*nxy*nxy];
	result.assign(nxy*nxy);
	for (i=0; i<2*nxy; i++)
	{
		if (i < nxy)
		{
			for (j=0; j<2*nxy; j++)
			{
				u=i*2*nxy+j;
				if (j < nxy)
				{
					M_Ch[u] = Zfs[i*nxy+j];
				}
				else if (j == nxy){}
				else
				{
					M_Ch[u] = M_Ch[i*2*nxy+nxy-(j-nxy)];
				}
			}
		}
		else if (i == nxy){}
		else
		{
			for (j=0; j<2*nxy; j++)
			{
				u=i*2*nxy+j;
				M_Ch[u] = M_Ch[(nxy-(i-nxy))*2*nxy+j];
			}			
		}
	}
	
	for (i=0; i<nxy; i++)
	{
		for (j=0; j<nxy; j++)
		{
			u=i*2*nxy+j;
			M_Cv[u] = vec[i*nxy+j];
		}
	}
	
/*	for (i=0; i<2*Nxy; i++)
	{
		for (j=0; j<2*Nxy; j++)
		{
			u=i*2*Nxy+j;
			display(M_Ch[u]);
		}
		outzm<<";"<<endl;
	}*/
	//2D FFT
	//检验2D-FFT
		
	fft2(2*nxy,M_Ch,M_Ch,1);
	fft2(2*nxy,M_Cv,M_Cv,1);
	
	for (i=0; i<2*nxy; i++)
	{
		for (j=0; j<2*nxy; j++)
		{
			u=i*2*nxy+j;
			M_Ch[u] = M_Ch[u] * M_Cv[u];
		}
	}

	fft2(2*nxy,M_Ch,M_Ch,-1);

	for (i=0; i<nxy; i++)
	{
		for (j=0; j<nxy; j++)
		{
			u=i*2*nxy+j;
			result[i*nxy+j] = M_Ch[u];
			//cout<<result[i*nxy+j]<<endl;
		}
	}

	delete [] M_Ch;
	delete [] M_Cv;
}
void Mathe_Functions::fft2(int n, Complex x[], Complex y[], int Bool)
{//二维复序列2D-FFT的行列算法	
	int i,j;

	Complex *M_fx=new Complex[n];
	Complex *M_Fy=new Complex[n];
	//行的一维FFT
	for (i=0; i<n;i++)
	{
		for (j=0; j<n; j++)
		{
			M_fx[j] = x[i*n+j];
		}
		fft(n,M_fx,M_Fy,Bool);
		for (j=0; j<n; j++)
		{
			y[i*n+j] = M_Fy[j];
		}
	}
	//列的一维FFT
	for (j=0; j<n; j++)
	{
		for (i=0; i<n; i++)
		{
			M_fx[i] = y[i*n+j];
		}
		fft(n,M_fx,M_Fy,Bool);
		for (i=0; i<n; i++)
		{
			y[i*n+j] = M_Fy[i];
		}
	}

	delete [] M_fx;
	delete [] M_Fy;
}
/*
void Mathe_Functions::fft2(int nx, int ny, Complex x[], Complex y[], int Bool)
{//二维复序列2D-FFT的行列算法，将x变换为nx*ny维	
	int i,j;

	int n = nx * ny;
	Complex *M_fx=new Complex[nx];
	Complex *M_Fy=new Complex[ny];
	//行的一维FFT
	for (i=0; i<nx;i++)
	{
		for (j=0; j<ny; j++)
		{
			M_fx[j] = x[i*ny+j];
		}
		fft(ny,M_fx,M_Fy,Bool);
		for (j=0; j<ny; j++)
		{
			y[i*ny+j] = M_Fy[j];
		}
	}
	//列的一维FFT
	for (j=0; j<ny; j++)
	{
		for (i=0; i<nx; i++)
		{
			M_fx[i] = y[i*ny+j];
		}
		fft(nx,M_fx,M_Fy,Bool);
		for (i=0; i<ny; i++)
		{
			y[i*ny+j] = M_Fy[i];
		}
	}

	delete [] M_fx;
	delete [] M_Fy;
}
*/
void Mathe_Functions::fft2(int nx, int ny, Complex x[], Complex y[], int Bool)
{//二维复序列2D-FFT的行列算法，将x变换为nx*ny维	
	int i,j;

	int n = nx * ny;
	Complex *M_fx=new Complex[nx];
	Complex *M_Fy=new Complex[ny];
	//行的一维FFT
	for (i=0; i<nx;i++)
	{
		for (j=0; j<ny; j++)
		{
			M_Fy[j] = x[i*ny+j];
		}
		fft(ny,M_Fy,M_Fy,Bool);
		for (j=0; j<ny; j++)
		{
			y[i*ny+j] = M_Fy[j];
		}
	}
	//列的一维FFT
	for (j=0; j<ny; j++)
	{
		for (i=0; i<nx; i++)
		{
			M_fx[i] = y[i*ny+j];
		}
		fft(nx,M_fx,M_fx,Bool);
		for (i=0; i<nx; i++)
		{
			y[i*ny+j] = M_fx[i];
		}
	}

	delete [] M_fx;
	delete [] M_Fy;
}
void Mathe_Functions::fft(int n, Complex x[], Complex y[], int Bool)
{
	int i;
	int   sofarRadix[maxFactorCount], 
          actualRadix[maxFactorCount], 
          remainRadix[maxFactorCount];
    int   nFactor;
    int   count;
	double *xRe=new double[n];
	double *xIm=new double[n];
	double *yRe=new double[n];
	double *yIm=new double[n];
	if(Bool == -1)
	{
		for(i=0;i<n;i++)
		{
			xRe[i]=real(x[i]);
			xIm[i]=-imag(x[i]);
		}
	}
	else
	{
		for(i=0;i<n;i++)
		{
			xRe[i]=real(x[i]);
			xIm[i]=imag(x[i]);
		}
	}
    

    //PI = 4*atan(1);    

    transTableSetup(sofarRadix, actualRadix, remainRadix, &nFactor, &n);
    permute(n, nFactor, actualRadix, remainRadix, xRe, xIm, yRe, yIm);

    for (count=1; count<=nFactor; count++)
      twiddleTransf(sofarRadix[count], actualRadix[count], remainRadix[count], 
                    yRe, yIm);
	if(Bool == -1)
	{
		for(i=0;i<n;i++)
		{
			y[i]=(yRe[i]-i_*yIm[i])/(1.0*n);
			//y[i].Real=yRe[i]/(1.0*n);
			//y[i].Imag=-yIm[i]/(1.0*n);
		}
	}
	else
	{
		for(i=0;i<n;i++)
		{
			y[i]=yRe[i]+i_*yIm[i];
			//y[i].Real=yRe[i];
			//y[i].Imag=yIm[i];
		}
	}
	delete [] xRe;
	delete [] xIm;
	delete [] yRe;
	delete [] yIm;
}   /* fft */

void Mathe_Functions::transTableSetup(int sofar[], int actual[], int remain[], int *nFact, int *nPoints)
{
    int i;

    factorize(*nPoints, nFact, actual);
    if (actual[1] > maxPrimeFactor)
    {
        cout<<"\nPrime factor of FFT length too large :"<<actual[1];
        cout<<"\nPlease modify the value of maxPrimeFactor in mixfft.c";
        exit(1);
    }
    remain[0]=*nPoints;
    sofar[1]=1;
    remain[1]=*nPoints / actual[1];
    for (i=2; i<=*nFact; i++)
    {
        sofar[i]=sofar[i-1]*actual[i-1];
        remain[i]=remain[i-1] / actual[i];
    }
}   /* transTableSetup */

void Mathe_Functions::permute(int nPoint, int nFact, int fact[], int remain[], double xRe[], double xIm[], double yRe[], double yIm[])
{
    int i,j,k;
    int count[maxFactorCount]; 

    for (i=1; i<=nFact; i++) count[i]=0;
    k=0;
    for (i=0; i<=nPoint-2; i++)
    {
        yRe[i] = xRe[k];
        yIm[i] = xIm[k];
        j=1;
        k=k+remain[j];
        count[1] = count[1]+1;
        while (count[j] >= fact[j])
        {
            count[j]=0;
            k=k-remain[j-1]+remain[j+1];
            j=j+1;
            count[j]=count[j]+1;
        }
    }
    yRe[nPoint-1]=xRe[nPoint-1];
    yIm[nPoint-1]=xIm[nPoint-1];
}   /* permute */


void Mathe_Functions::twiddleTransf(int sofarRadix, int radix, int remainRadix, double yRe[], double yIm[])
{   /* twiddleTransf */ 
    double  cosw, sinw, gem;
    double  t1_re,t1_im, t2_re,t2_im, t3_re,t3_im;
    double  t4_re,t4_im, t5_re,t5_im;
    double  m2_re,m2_im, m3_re,m3_im, m4_re,m4_im;
    double  m1_re,m1_im, m5_re,m5_im;
    double  s1_re,s1_im, s2_re,s2_im, s3_re,s3_im;
    double  s4_re,s4_im, s5_re,s5_im;


    initTrig(radix);
    omega = 2*PI/(double)(sofarRadix*radix);
    cosw =  cos(omega);
    sinw = -sin(omega);
    tw_re = 1.0;
    tw_im = 0;
    dataOffset=0;
    groupOffset=dataOffset;
    adr=groupOffset;
    for (dataNo=0; dataNo<sofarRadix; dataNo++)
    {
        if (sofarRadix>1)
        {
            twiddleRe[0] = 1.0; 
            twiddleIm[0] = 0.0;
            twiddleRe[1] = tw_re;
            twiddleIm[1] = tw_im;
            for (twNo=2; twNo<radix; twNo++)
            {
                twiddleRe[twNo]=tw_re*twiddleRe[twNo-1]
                               - tw_im*twiddleIm[twNo-1];
                twiddleIm[twNo]=tw_im*twiddleRe[twNo-1]
                               + tw_re*twiddleIm[twNo-1];
            }
            gem   = cosw*tw_re - sinw*tw_im;
            tw_im = sinw*tw_re + cosw*tw_im;
            tw_re = gem;                      
        }
        for (groupNo=0; groupNo<remainRadix; groupNo++)
        {
            if ((sofarRadix>1) && (dataNo > 0))
            {
                zRe[0]=yRe[adr];
                zIm[0]=yIm[adr];
                blockNo=1;
                do {
                    adr = adr + sofarRadix;
                    zRe[blockNo]=  twiddleRe[blockNo] * yRe[adr]
                                 - twiddleIm[blockNo] * yIm[adr];
                    zIm[blockNo]=  twiddleRe[blockNo] * yIm[adr]
                                 + twiddleIm[blockNo] * yRe[adr]; 
                    
                    blockNo++;
                } while (blockNo < radix);
            }
            else
                for (blockNo=0; blockNo<radix; blockNo++)
                {
                   zRe[blockNo]=yRe[adr];
                   zIm[blockNo]=yIm[adr];
                   adr=adr+sofarRadix;
                }
            switch(radix) {
              case  2  : gem=zRe[0] + zRe[1];
                         zRe[1]=zRe[0] -  zRe[1]; zRe[0]=gem;
                         gem=zIm[0] + zIm[1];
                         zIm[1]=zIm[0] - zIm[1]; zIm[0]=gem;
                         break;
              case  3  : t1_re=zRe[1] + zRe[2]; t1_im=zIm[1] + zIm[2];
                         zRe[0]=zRe[0] + t1_re; zIm[0]=zIm[0] + t1_im;
                         m1_re=c3_1*t1_re; m1_im=c3_1*t1_im;
                         m2_re=c3_2*(zIm[1] - zIm[2]); 
                         m2_im=c3_2*(zRe[2] -  zRe[1]);
                         s1_re=zRe[0] + m1_re; s1_im=zIm[0] + m1_im;
                         zRe[1]=s1_re + m2_re; zIm[1]=s1_im + m2_im;
                         zRe[2]=s1_re - m2_re; zIm[2]=s1_im - m2_im;
                         break;
              case  4  : t1_re=zRe[0] + zRe[2]; t1_im=zIm[0] + zIm[2];
                         t2_re=zRe[1] + zRe[3]; t2_im=zIm[1] + zIm[3];

                         m2_re=zRe[0] - zRe[2]; m2_im=zIm[0] - zIm[2];
                         m3_re=zIm[1] - zIm[3]; m3_im=zRe[3] - zRe[1];

                         zRe[0]=t1_re + t2_re; zIm[0]=t1_im + t2_im;
                         zRe[2]=t1_re - t2_re; zIm[2]=t1_im - t2_im;
                         zRe[1]=m2_re + m3_re; zIm[1]=m2_im + m3_im;
                         zRe[3]=m2_re - m3_re; zIm[3]=m2_im - m3_im;
                         break;
              case  5  : t1_re=zRe[1] + zRe[4]; t1_im=zIm[1] + zIm[4];
                         t2_re=zRe[2] + zRe[3]; t2_im=zIm[2] + zIm[3];
                         t3_re=zRe[1] - zRe[4]; t3_im=zIm[1] - zIm[4];
                         t4_re=zRe[3] - zRe[2]; t4_im=zIm[3] - zIm[2];
                         t5_re=t1_re + t2_re; t5_im=t1_im + t2_im;
                         zRe[0]=zRe[0] + t5_re; zIm[0]=zIm[0] + t5_im;
                         m1_re=c5_1*t5_re; m1_im=c5_1*t5_im;
                         m2_re=c5_2*(t1_re - t2_re); 
                         m2_im=c5_2*(t1_im - t2_im);

                         m3_re=-c5_3*(t3_im + t4_im); 
                         m3_im=c5_3*(t3_re + t4_re);
                         m4_re=-c5_4*t4_im; m4_im=c5_4*t4_re;
                         m5_re=-c5_5*t3_im; m5_im=c5_5*t3_re;

                         s3_re=m3_re - m4_re; s3_im=m3_im - m4_im;
                         s5_re=m3_re + m5_re; s5_im=m3_im + m5_im;
                         s1_re=zRe[0] + m1_re; s1_im=zIm[0] + m1_im;
                         s2_re=s1_re + m2_re; s2_im=s1_im + m2_im;
                         s4_re=s1_re - m2_re; s4_im=s1_im - m2_im;

                         zRe[1]=s2_re + s3_re; zIm[1]=s2_im + s3_im;
                         zRe[2]=s4_re + s5_re; zIm[2]=s4_im + s5_im;
                         zRe[3]=s4_re - s5_re; zIm[3]=s4_im - s5_im;
                         zRe[4]=s2_re - s3_re; zIm[4]=s2_im - s3_im;
                         break;
              case  8  : fft_8(); break;
              case 10  : fft_10(); break;
              default  : fft_odd(radix); break;
            }
            adr=groupOffset;
            for (blockNo=0; blockNo<radix; blockNo++)
            {
                yRe[adr]=zRe[blockNo]; yIm[adr]=zIm[blockNo];
                adr=adr+sofarRadix;
            }
            groupOffset=groupOffset+sofarRadix*radix;
            adr=groupOffset;
        }
        dataOffset=dataOffset+1;
        groupOffset=dataOffset;
        adr=groupOffset;
    }
}   /* twiddleTransf */

void Mathe_Functions::factorize(int n, int *nFact, int fact[])
{
    int i,j,k;
    int nRadix;
    int radices[7];
    int factors[maxFactorCount];

    nRadix    =  6;  
    radices[1]=  2;
    radices[2]=  3;
    radices[3]=  4;
    radices[4]=  5;
    radices[5]=  8;
    radices[6]= 10;

    if (n==1)
    {
        j=1;
        factors[1]=1;
    }
    else j=0;
    i=nRadix;
    while ((n>1) && (i>0))
    {
      if ((n % radices[i]) == 0)
      {
        n=n / radices[i];
        j=j+1;
        factors[j]=radices[i];
      }
      else  i=i-1;
    }
    if (factors[j] == 2)   /*substitute factors 2*8 with 4*4 */
    {   
      i = j-1;
      while ((i>0) && (factors[i] != 8)) i--;
      if (i>0)
      {
        factors[j] = 4;
        factors[i] = 4;
      }
    }
    if (n>1)
    {
        for (k=2; k<sqrt(n)+1; k++)
            while ((n % k) == 0)
            {
                n=n / k;
                j=j+1;
                factors[j]=k;
            }
        if (n>1)
        {
            j=j+1;
            factors[j]=n;
        }
    }               
    for (i=1; i<=j; i++)         
    {
      fact[i] = factors[j-i+1];
    }
    *nFact=j;
}   /* factorize */

void Mathe_Functions::initTrig(int radix)
{
    int i;
    double w,xre,xim;

    w=2*PI/radix;
    trigRe[0]=1; trigIm[0]=0;
    xre=cos(w); 
    xim=-sin(w);
    trigRe[1]=xre; trigIm[1]=xim;
    for (i=2; i<radix; i++)
    {
        trigRe[i]=xre*trigRe[i-1] - xim*trigIm[i-1];
        trigIm[i]=xim*trigRe[i-1] + xre*trigIm[i-1];
    }
}   /* initTrig */

void Mathe_Functions::fft_4(double aRe[], double aIm[])
{
    double  t1_re,t1_im, t2_re,t2_im;
    double  m2_re,m2_im, m3_re,m3_im;

    t1_re=aRe[0] + aRe[2]; t1_im=aIm[0] + aIm[2];
    t2_re=aRe[1] + aRe[3]; t2_im=aIm[1] + aIm[3];

    m2_re=aRe[0] - aRe[2]; m2_im=aIm[0] - aIm[2];
    m3_re=aIm[1] - aIm[3]; m3_im=aRe[3] - aRe[1];

    aRe[0]=t1_re + t2_re; aIm[0]=t1_im + t2_im;
    aRe[2]=t1_re - t2_re; aIm[2]=t1_im - t2_im;
    aRe[1]=m2_re + m3_re; aIm[1]=m2_im + m3_im;
    aRe[3]=m2_re - m3_re; aIm[3]=m2_im - m3_im;
}   /* fft_4 */


void Mathe_Functions::fft_5(double aRe[], double aIm[])
{    
    double  t1_re,t1_im, t2_re,t2_im, t3_re,t3_im;
    double  t4_re,t4_im, t5_re,t5_im;
    double  m2_re,m2_im, m3_re,m3_im, m4_re,m4_im;
    double  m1_re,m1_im, m5_re,m5_im;
    double  s1_re,s1_im, s2_re,s2_im, s3_re,s3_im;
    double  s4_re,s4_im, s5_re,s5_im;

    t1_re=aRe[1] + aRe[4]; t1_im=aIm[1] + aIm[4];
    t2_re=aRe[2] + aRe[3]; t2_im=aIm[2] + aIm[3];
    t3_re=aRe[1] - aRe[4]; t3_im=aIm[1] - aIm[4];
    t4_re=aRe[3] - aRe[2]; t4_im=aIm[3] - aIm[2];
    t5_re=t1_re + t2_re; t5_im=t1_im + t2_im;
    aRe[0]=aRe[0] + t5_re; aIm[0]=aIm[0] + t5_im;
    m1_re=c5_1*t5_re; m1_im=c5_1*t5_im;
    m2_re=c5_2*(t1_re - t2_re); m2_im=c5_2*(t1_im - t2_im);

    m3_re=-c5_3*(t3_im + t4_im); m3_im=c5_3*(t3_re + t4_re);
    m4_re=-c5_4*t4_im; m4_im=c5_4*t4_re;
    m5_re=-c5_5*t3_im; m5_im=c5_5*t3_re;

    s3_re=m3_re - m4_re; s3_im=m3_im - m4_im;
    s5_re=m3_re + m5_re; s5_im=m3_im + m5_im;
    s1_re=aRe[0] + m1_re; s1_im=aIm[0] + m1_im;
    s2_re=s1_re + m2_re; s2_im=s1_im + m2_im;
    s4_re=s1_re - m2_re; s4_im=s1_im - m2_im;

    aRe[1]=s2_re + s3_re; aIm[1]=s2_im + s3_im;
    aRe[2]=s4_re + s5_re; aIm[2]=s4_im + s5_im;
    aRe[3]=s4_re - s5_re; aIm[3]=s4_im - s5_im;
    aRe[4]=s2_re - s3_re; aIm[4]=s2_im - s3_im;
}   /* fft_5 */

void Mathe_Functions::fft_8()
{
    double  aRe[4], aIm[4], bRe[4], bIm[4], gem;

    aRe[0] = zRe[0];    bRe[0] = zRe[1];
    aRe[1] = zRe[2];    bRe[1] = zRe[3];
    aRe[2] = zRe[4];    bRe[2] = zRe[5];
    aRe[3] = zRe[6];    bRe[3] = zRe[7];

    aIm[0] = zIm[0];    bIm[0] = zIm[1];
    aIm[1] = zIm[2];    bIm[1] = zIm[3];
    aIm[2] = zIm[4];    bIm[2] = zIm[5];
    aIm[3] = zIm[6];    bIm[3] = zIm[7];

    fft_4(aRe, aIm); fft_4(bRe, bIm);

    gem    = c8*(bRe[1] + bIm[1]);
    bIm[1] = c8*(bIm[1] - bRe[1]);
    bRe[1] = gem;
    gem    = bIm[2];
    bIm[2] =-bRe[2];
    bRe[2] = gem;
    gem    = c8*(bIm[3] - bRe[3]);
    bIm[3] =-c8*(bRe[3] + bIm[3]);
    bRe[3] = gem;
    
    zRe[0] = aRe[0] + bRe[0]; zRe[4] = aRe[0] - bRe[0];
    zRe[1] = aRe[1] + bRe[1]; zRe[5] = aRe[1] - bRe[1];
    zRe[2] = aRe[2] + bRe[2]; zRe[6] = aRe[2] - bRe[2];
    zRe[3] = aRe[3] + bRe[3]; zRe[7] = aRe[3] - bRe[3];

    zIm[0] = aIm[0] + bIm[0]; zIm[4] = aIm[0] - bIm[0];
    zIm[1] = aIm[1] + bIm[1]; zIm[5] = aIm[1] - bIm[1];
    zIm[2] = aIm[2] + bIm[2]; zIm[6] = aIm[2] - bIm[2];
    zIm[3] = aIm[3] + bIm[3]; zIm[7] = aIm[3] - bIm[3];
}   /* fft_8 */

void Mathe_Functions::fft_10()
{
    double  aRe[5], aIm[5], bRe[5], bIm[5];

    aRe[0] = zRe[0];    bRe[0] = zRe[5];
    aRe[1] = zRe[2];    bRe[1] = zRe[7];
    aRe[2] = zRe[4];    bRe[2] = zRe[9];
    aRe[3] = zRe[6];    bRe[3] = zRe[1];
    aRe[4] = zRe[8];    bRe[4] = zRe[3];

    aIm[0] = zIm[0];    bIm[0] = zIm[5];
    aIm[1] = zIm[2];    bIm[1] = zIm[7];
    aIm[2] = zIm[4];    bIm[2] = zIm[9];
    aIm[3] = zIm[6];    bIm[3] = zIm[1];
    aIm[4] = zIm[8];    bIm[4] = zIm[3];

    fft_5(aRe, aIm); fft_5(bRe, bIm);

    zRe[0] = aRe[0] + bRe[0]; zRe[5] = aRe[0] - bRe[0];
    zRe[6] = aRe[1] + bRe[1]; zRe[1] = aRe[1] - bRe[1];
    zRe[2] = aRe[2] + bRe[2]; zRe[7] = aRe[2] - bRe[2];
    zRe[8] = aRe[3] + bRe[3]; zRe[3] = aRe[3] - bRe[3];
    zRe[4] = aRe[4] + bRe[4]; zRe[9] = aRe[4] - bRe[4];

    zIm[0] = aIm[0] + bIm[0]; zIm[5] = aIm[0] - bIm[0];
    zIm[6] = aIm[1] + bIm[1]; zIm[1] = aIm[1] - bIm[1];
    zIm[2] = aIm[2] + bIm[2]; zIm[7] = aIm[2] - bIm[2];
    zIm[8] = aIm[3] + bIm[3]; zIm[3] = aIm[3] - bIm[3];
    zIm[4] = aIm[4] + bIm[4]; zIm[9] = aIm[4] - bIm[4];
}   /* fft_10 */

void Mathe_Functions::fft_odd(int radix)
{
    double  rere, reim, imre, imim;
    int     i,j,k,n,max;

    n = radix;
    max = (n + 1)/2;
    for (j=1; j < max; j++)
    {
      vRe[j] = zRe[j] + zRe[n-j];
      vIm[j] = zIm[j] - zIm[n-j];
      wRe[j] = zRe[j] - zRe[n-j];
      wIm[j] = zIm[j] + zIm[n-j];
    }

    for (j=1; j < max; j++)
    {
        zRe[j]=zRe[0]; 
        zIm[j]=zIm[0];
        zRe[n-j]=zRe[0]; 
        zIm[n-j]=zIm[0];
        k=j;
        for (i=1; i < max; i++)
        {
            rere = trigRe[k] * vRe[i];
            imim = trigIm[k] * vIm[i];
            reim = trigRe[k] * wIm[i];
            imre = trigIm[k] * wRe[i];
            
            zRe[n-j] += rere + imim;
            zIm[n-j] += reim - imre;
            zRe[j]   += rere - imim;
            zIm[j]   += reim + imre;

            k = k + j;
            if (k >= n)  k = k - n;
        }
    }
    for (j=1; j < max; j++)
    {
        zRe[0]=zRe[0] + vRe[j]; 
        zIm[0]=zIm[0] + wIm[j];
    }
}   /* fft_odd */