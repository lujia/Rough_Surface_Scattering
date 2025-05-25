#include "LinearEquations.h"
//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
const int  M=360;
//const double  PI = 3.1415926535897932384626433832795;
const Complex J(0,1.0);
const Complex zero(0.,0.);

//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
	
class CylinderQuestion 
{
public:
    double Yita;
    double Gama;
    double Lamda;
    double R;
    double H;
    double ds;
    double K;
    double Sita0;
	double dx_t;

   
    
	Vector_complex ZT;
    Vector_complex J0;
    Vector_complex Ei_T;
    Vector_double xtt;
    Vector_double ztt;
	Vector_double ds_Muti;
	Vector_double df_z;
	Vector_double ddf_z;
	double Sita;
public:
	CylinderQuestion();
	~CylinderQuestion(); 
	void Cy_PoufenMatix();
	void Cy_PoufenMatix_above_target();
	void Cy_PoufenMatix_buried_target();
	void Cy_PutEi();
	void Cy_ImplenceMatix();
	void Cy_ResultofJ();
	void Muti_target_poufen(double xpp,Vector_double &x_number,Vector_double &z_number,Vector_double &ds_num,double aa,double bb,int n);
	void Muti_target_ImplenceMatix(Vector_double &x_number,Vector_double &z_number,Vector_complex &ZT_num,Vector_double &ds_num,int n);

};


CylinderQuestion::CylinderQuestion()
{
	ZT.assign(M*M);
	Ei_T.assign(M);
	J0.assign(M);
	xtt.assign(M);
	ztt.assign(M);
	df_z.assign(M);
	ddf_z.assign(M);
	ds_Muti.assign(M);
    Sita=0;

}


CylinderQuestion::~CylinderQuestion()
{



}


void CylinderQuestion::Cy_PoufenMatix()
{
	cout<<"进行目标的剖分......------>>"<<endl;
	
	Lamda=3e8/3e9; 
    R=5.*Lamda;
		
	H=10.*Lamda;
	ds=2.*PI*R/M;
	Sita0=30.*PI/180.;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;


	for(int i=0;i<=M-1;i++)
	{
		Sita=i*2.*PI/M ;
		xtt[i]=R*cos(Sita) ;
		ztt[i]=R*sin(Sita)+H;		
	}

	df_z[0] = (ztt[1]-ztt[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
        df_z[i] = (ztt[i+1]-ztt[i-1])/(2*dx_t);
	}
	df_z[M-1] = (ztt[0]-ztt[M-2])/(2*dx_t);

	ddf_z[0] = (df_z[1]-df_z[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
      ddf_z[i] = (df_z[i+1]-df_z[i-1])/(2*dx_t);
	}
	ddf_z[M-1] = (df_z[0]-df_z[M-2])/(2*dx_t);
}


void CylinderQuestion::Cy_PoufenMatix_buried_target()
{
	cout<<"进行目标的剖分......------>>"<<endl;
	
	Lamda=3e8/3e9; 
    R=1.*Lamda;;
		
	H=-5.*Lamda;
	ds=2.*PI*R/M;
	Sita0=30.*PI/180.;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;

	double xp = 0.*Lamda;


	for(int i=0;i<=M-1;i++)
	{
		Sita=i*2.*PI/M ;
		xtt[i]=R*cos(Sita)+xp ;
		ztt[i]=R*sin(Sita)+H;		
	}

	df_z[0] = (ztt[1]-ztt[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
        df_z[i] = (ztt[i+1]-ztt[i-1])/(2*dx_t);
	}
	df_z[M-1] = (ztt[0]-ztt[M-2])/(2*dx_t);

	ddf_z[0] = (df_z[1]-df_z[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
      ddf_z[i] = (df_z[i+1]-df_z[i-1])/(2*dx_t);
	}
	ddf_z[M-1] = (df_z[0]-df_z[M-2])/(2*dx_t);
}


void CylinderQuestion::Cy_PoufenMatix_above_target()
{
	cout<<"进行目标的剖分......------>>"<<endl;
	
	Lamda=3e8/1.5e9; 
    R=2.*Lamda;;
		
	H=5.*Lamda;
	ds=2.*PI*R/M;
	Sita0=30.*PI/180.;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;

	double xp = 0.*Lamda;


	for(int i=0;i<=M-1;i++)
	{
		Sita=i*2.*PI/M ;
		xtt[i]=R*cos(Sita)+xp ;
		ztt[i]=R*sin(Sita)+H;		
	}

	df_z[0] = (ztt[1]-ztt[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
        df_z[i] = (ztt[i+1]-ztt[i-1])/(2*dx_t);
	}
	df_z[M-1] = (ztt[0]-ztt[M-2])/(2*dx_t);

	ddf_z[0] = (df_z[1]-df_z[M-1])/(2*dx_t);
	for (i=1;i<M-1;i++)
	{
      ddf_z[i] = (df_z[i+1]-df_z[i-1])/(2*dx_t);
	}
	ddf_z[M-1] = (df_z[0]-df_z[M-2])/(2*dx_t);
}


void CylinderQuestion::Cy_PutEi()
{

	
    double K=2.*PI/Lamda;
    double Sita0=PI/6.;
	for(int i=0;i<=M-1;i++)
	{
		Ei_T[i]= exp(J*K*(xtt[i]*sin(Sita0)-ztt[i]*cos(Sita0)));
		
	}
}


void CylinderQuestion::Cy_ImplenceMatix()
{

	cout<<"进行目标阻抗矩阵的填充......------>>"<<endl;
	
    Yita=120.*PI;
    Gama=1.781072418;   
    ds=2.*PI*R/M;
    K=2.*PI/Lamda;
    Sita0=PI/6.;

	
	CLinearEquations linear;
	for(int i=0;i<=M-1;i++)
	{
		for(int j=0;j<=M-1;j++)
		{
			if(i!=j)
			{
				double kk=K*sqrt((xtt[i]-xtt[j])*(xtt[i]-xtt[j])+(ztt[i]-ztt[j])*(ztt[i]-ztt[j]));
				ZT[i*M+j]=(K*Yita/4.)*ds*linear.Hankel(0,1,kk);
			}
			else
			{
                ZT[i*M+j]=(K*Yita/4.)*ds*(1.+2.*J/PI*(log(Gama*K*ds/4.)-1.));	
			}
		}
	}
}


void CylinderQuestion::Cy_ResultofJ()
{

	

    Cy_PoufenMatix();
    Cy_PutEi();
    Cy_ImplenceMatix();

	CLinearEquations linear;

	cout<<"用CG法求目标表面电流......------>>"<<endl;

	linear.CG(ZT,J0,Ei_T,M,M,0.01);
	ofstream Out1;
	Out1.open("J0.txt");

	for(int m=0;m<M;m++)
	{
//		Out1<<J1[m].real()<<"   "<<J1[m].imag()<<"   "<<abs(J1[m])<<endl;
//		cout<<J1[m].real()<<"   "<<J1[m].imag()<<"   "<<abs(J1[m])<<endl;
		Out1<<m*360./M<<" "<<real(J0[m])<<"  "<<imag(J0[m])<<" "<<abs(J0[m])<<endl;
	}
}


void CylinderQuestion::Muti_target_poufen(double xpp,Vector_double &x_number,Vector_double &z_number,Vector_double &ds_num,double aa,double bb,int n)
	{
		cout<<"进行目标的剖分......------>>"<<endl;
	

	    Lamda=3e8/3e9; 
        R=5.*Lamda;
	    H=10.*Lamda;
	    x_number.assign(n);
	    z_number.assign(n);


	    for(int i=0;i<M;i++)
		{
		   Sita=i*2.*PI/M ;
		   x_number[i]=aa*cos(Sita)+xpp ;
		   z_number[i]=bb*sin(Sita)+H;	
		   ds_num[i] = sqrt(pow(aa*sin(Sita),2.0)+pow(bb*cos(Sita),2.0))*2.*PI/M;	
		}



	}

void CylinderQuestion::Muti_target_ImplenceMatix(Vector_double &x_number,Vector_double &z_number,Vector_complex &ZT_num,Vector_double &ds_num ,int n)
{
	cout<<"进行目标阻抗矩阵的填充......------>>"<<endl;
	
    Yita=120.*PI;
    Gama=1.781072418;   
    K=2.*PI/Lamda;
    Sita0=PI/6.;

	ZT_num.assign(n*n);

	
	CLinearEquations linear;
	for(int i=0;i<=M-1;i++)
	{
		for(int j=0;j<=M-1;j++)
		{
			if(i!=j)
			{
				double kk=K*sqrt((x_number[i]-x_number[j])*(x_number[i]-x_number[j])+(z_number[i]-z_number[j])*(z_number[i]-z_number[j]));
				ZT_num[i*M+j]=(K*Yita/4.)*ds_num[j]*linear.Hankel(0,1,kk);
			}
			else
			{
                ZT_num[i*M+j]=(K*Yita/4.)*ds_num[j]*(1.+2.*J/PI*(log(Gama*K*ds_num[j]/4.)-1.));	
			}
		}
	}
}






