#include "resolve_matrix.h"

ofstream outrs("rough_surface.txt",ios::trunc);
ofstream outx("xx.txt",ios::trunc);
ofstream outy("yy.txt",ios::trunc);
ofstream outz("zz.txt",ios::trunc);
ofstream outdfx("dfx.txt",ios::trunc);
ofstream outdfy("dfy.txt",ios::trunc);
ofstream outddfx("ddfx.txt",ios::trunc);
ofstream outddfy("ddfy.txt",ios::trunc);

class ThreeD//建立复矢量三维坐标 
{ public:
     Complex TD_x;
	 Complex TD_y;
	 Complex TD_z;

  public:
	  ThreeD()
	  {
		  TD_x = 0.0;
		  TD_y = 0.0;
		  TD_z = 0.0;
	  };
	  ThreeD(Complex x1,Complex y1,Complex z1)
	  {
		  TD_x = x1;
		  TD_y = y1;
		  TD_z = z1;
	  };
	  ~ThreeD(){};
	  void Assign(Complex,Complex,Complex);
	  ThreeD operator +(ThreeD t);/////////操作符重载
	  ThreeD operator -(ThreeD t);
	  ThreeD operator =(ThreeD t);
	  //矢量与任意单位向量做内积
	  Complex operator *( ThreeD t);// 内积操作
	  ThreeD operator *(Complex t);
	  ThreeD operator /(Complex t);
	  double Module(ThreeD t);
	  int operator==(ThreeD);
	  Complex Point_Multiply(ThreeD t1,ThreeD t2); //内积函数
	  ThreeD Xc(ThreeD t1,ThreeD t2); //外积函数
	  //重载流插入和析取运算符必须为非成员函数，即为友元函数
	  friend ostream& operator << (ostream&,ThreeD&);
	  friend istream& operator >> (istream&,ThreeD&);
};

class RoughSurface: public resolve_matrix
{
public:
	RoughSurface();
	~RoughSurface();
	//一维粗糙面
	void rand_n();
	void gauss_surface();
	double wk(double);
	void rand_fn();
	void fractal_surface();
	void ocean_surface();
	void ocean_surface_PM(double U195);
	double owk(double);
	double owk_PM(double kx,double U195);
	//二维粗糙面
	//高斯
	void Two_D_rand_n(int nxy,double mmean,double vvari);
    void Two_D_Gauss_Surface(int nxy,double mm,double vv,double hh,double llx,double lly,double rl);
	void TD_Gauss_Surface(int nxy,double mm,double vv,double hh,double llx,double lly,double rl);
    double Two_D_wk(double kx,double ky);
	//分形
	void Two_D_Frctal_Rand_n(int nnf,double mmean,double vvari);
	void Two_D_Fractal_surface(int nxy,int nnf,double hh,double mm,double vv,double rrl,double kkl,double kku,double fra_s);
	void show_rs();
public:
	//n:number of sampling points on surface,seed:seed for random number generator 
	int N;//
	double seed;//随机数
	//double *x,*z,*df_x,*ddf_x;//
	Vec_Dou r_n,x,y,z,df_x,ddf_x,df_y,ddf_y,ddf_xy;
	//h:rms height,lc: correlation length,rl:rough surface length 
	double h,lc,RL,mean,vari;//u为均值，g为方差
	double lx,ly;
	int Nxy;
	double kl,ku,Fractal_S;
	int nf;
	double us;
	double wave_len;
	double dx_Suf;
};
//////////////////////////////////////////////////////////////////////////
//////////////矢量的四则运算/////////////////
void ThreeD::Assign (Complex x1,Complex y1,Complex z1)
{
	TD_x=x1;
	TD_y=y1;
	TD_z=z1;
}
ThreeD ThreeD::operator + (ThreeD t)
{
	ThreeD temp;
	temp.TD_x=TD_x+t.TD_x;
	temp.TD_y=TD_y+t.TD_y;
	temp.TD_z=TD_z+t.TD_z;
	return temp;
}
ThreeD ThreeD::operator - (ThreeD t)
{
	ThreeD temp;
	temp.TD_x=TD_x-t.TD_x;
	temp.TD_y=TD_y-t.TD_y;
	temp.TD_z=TD_z-t.TD_z;
	return temp;
}
ThreeD ThreeD::operator = (ThreeD t)
{
	TD_x=t.TD_x ;
	TD_y=t.TD_y ;
	TD_z=t.TD_z ;
	return *this;
}
Complex ThreeD::operator * (ThreeD t)
{ 
	Complex temp;
	if (t.TD_x != 0.0 && t.TD_y == 0.0 && t.TD_z == 0.0)
	{
		temp = TD_x*t.TD_x;
	}
	if (t.TD_x == 0.0 && t.TD_y != 0.0 && t.TD_z == 0.0)
	{
		temp = TD_y*t.TD_y;
	}
	if (t.TD_x == 0.0 && t.TD_y == 0.0 && t.TD_z != 0.0)
	{
		temp = TD_z*t.TD_z;
	}

	return temp;
}
ThreeD ThreeD::operator *(Complex t)
{
	ThreeD temp;
	temp.TD_x=TD_x*t;
	temp.TD_y=TD_y*t;
	temp.TD_z=TD_z*t;
	return temp;
}
ThreeD ThreeD::operator /(Complex t)
{
	ThreeD temp;
	temp.TD_x=TD_x/t;
	temp.TD_y=TD_y/t;
	temp.TD_z=TD_z/t;
	return temp;
}
double ThreeD::Module (ThreeD t)
{
	return abs(sqrt(t.TD_x *t.TD_x +t.TD_y *t.TD_y +t.TD_z *t.TD_z));
}
Complex ThreeD::Point_Multiply(ThreeD t1,ThreeD t2)
{
	Complex temp;
	temp = t1.TD_x*t2.TD_x + t1.TD_y*t2.TD_y + t1.TD_z*t2.TD_z;
	return temp;
}
ThreeD ThreeD::Xc(ThreeD t1,ThreeD t2)//外积函数
{
	ThreeD temp;
	temp.TD_x=t1.TD_y*t2.TD_z-t1.TD_z*t2.TD_y;
	temp.TD_y=t1.TD_z*t2.TD_x-t1.TD_x*t2.TD_z;
	temp.TD_z=t1.TD_x*t2.TD_y-t1.TD_y*t2.TD_x;
	return temp;
}
int ThreeD::operator==(ThreeD t)
{
	if(TD_x==t.TD_x && TD_y==t.TD_y && TD_z==t.TD_z)
		return 1;
	else
		return 0;
}
ostream& operator << (ostream& osobject,ThreeD& x)
{
	osobject<<x.TD_x<<"  "<<x.TD_y<<"  "<<x.TD_z<<endl;
	
	return osobject;
}

//////////////////////////矢量四则运算结束////////////////////////////////
//////////////////////////////////////////////////////////////////////////

RoughSurface::RoughSurface()
{
	N=0;
	seed=0;
	h=0.0;
	lc=0.0;
	RL=0.0;
	mean=0.0;
	vari=0.0;
	kl=0.0;
	ku=0.0;
	Fractal_S=0.0;
	us=0.0;
}
RoughSurface::~RoughSurface()
{
	r_n.clear();
	x.clear();
	y.clear();
	z.clear();
	df_x.clear();
	ddf_x.clear();
	df_y.clear();
	ddf_y.clear();
}
void RoughSurface::rand_n()
{
    int i,k,m;
	r_n.assign(N);
	double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for(k=0;k<N;k++)
	{
		t=0.0;
		//cout<<"第i个random number"<<endl;
        for(i=1;i<=12;i++)
        {
			seed=seed*w+v;
			m=(int)(seed/s);
            seed=seed-m*s;
			t=t+seed/s;

		}
        r_n[k]=mean+vari*(t-6.0);

		//cout<<r_n[k]<<endl;
    }
	seed++;
}
void RoughSurface::gauss_surface()
{
	rand_n();


	int i;
	//double r_n[10]={0.6787,0.1629,-0.5375,-0.7282,0.5821,-0.4467,0.2605,1.4576,-1.2747,-1.3908};
	Complex *b=new Complex[N];
	double *kx=new double[N];
	x.assign(N);
	for(i=-N/2+1;i<=N/2;i++)
	{
		x[i+N/2-1]=i*RL/N;
      
	}
	for(i=-N/2+1;i<=N/2;i++)
	{
		kx[i+N/2-1]=2*PI*i/RL;

	}
	b[N/2-1]=wk(kx[N/2-1])*r_n[N-2];//b(0)
	b[N-1]=wk(kx[N-1])*r_n[N-1];//b(N/2)
	for(i=0;i<N/2-1;i++)
	{
		b[i]=wk(kx[i])*(r_n[2*i]+i_*r_n[2*i+1])/sqrt(2);//b(-N/2+1)到b(-1)
		//b[i].Imag=wk(kx[i])*r_n[2*i+1]/sqrt(2);

	}
	for(i=N/2;i<N-1;i++)
	{
		b[i]=conj(b[i-2*(i+1-N/2)]);//b(1)到b(N/2-1)
		//b[i].Imag=-b[i-2*(i+1-N/2)].Imag;
	}
	Complex *xt=new Complex[N];
	for(i=0;i<N;i++)
	{
		if(i<=N/2)
		{
			xt[i]=b[i+N/2-1];
		}
		else
		{
			xt[i]=b[i-N/2-1];
		}
		
	}
	Complex *ft=new Complex[N];
	z.assign(N);
	df_x.assign(N);
	ddf_x.assign(N);
//	fft(N,xt,ft,-1);

	fft2(N,1,xt,ft,-1);
	
	for(i=0;i<N;i++)
	{
		if(i<=N/2-2)
		{
			z[i]=real(ft[i+N/2+1])*N/RL;
		}
		else
		{
			z[i]=real(ft[i-N/2+1])*N/RL;
		}
			  
		//z[i] = 0.0;//平板（直线）的模型
	}
	double dx=RL/N;
	df_x[0]=(z[1]-z[N-1])/(2*dx);
	df_x[N-1]=(z[0]-z[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		df_x[i]=(z[i+1]-z[i-1])/(2*dx);
	}
	////////////
	ddf_x[0]=(df_x[1]-df_x[N-1])/(2*dx);
	ddf_x[N-1]=(df_x[0]-df_x[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		ddf_x[i]=(df_x[i+1]-df_x[i-1])/(2*dx);
	}
	/////////////
	delete [] b;
	delete [] kx;
	delete [] xt;
	delete [] ft;
	show_rs();

}
double RoughSurface::wk(double kx)
{
	double z;
	//wk为Gaussian spectrum;
	z=sqrt(sqrt(PI)*RL*pow(h,2)*lc*exp(-pow(kx*lc*0.5,2)));

	return z;
}
void RoughSurface::rand_fn()
{
    int i,k,m;
	r_n.assign(nf,0);
	double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for(k=0;k<nf;k++)
	{
		t=0.0;
		//cout<<"第i个random number"<<endl;
        for(i=1;i<=12;i++)
        {
			seed=seed*w+v;
			m=(int)(seed/s);
            seed=seed-m*s;
			t=t+seed/s;
		}
        r_n[k]=mean+vari*(t-6.0);
		r_n[k]=2*PI*r_n[k];
		//cout<<r_n[k]<<endl;
    }
}
void RoughSurface::fractal_surface()
{
	//N=num;
	//RL=rrl;
	//us=uus;
	//kl=kkl;
	//ku=kku;
	///随机面的产生
	rand_fn();
   	int i,k;
	x.assign(N,0);
	for(i=0;i<N;i++)
	{
		x[i]=(i+1/2)*RL/N-RL/2;
	}
	z.assign(N);
	df_x.assign(N);
	ddf_x.assign(N);
	double b=pow(ku/kl,1.0/(nf-1));
	double cn=sqrt(2*(1-pow(b,2*(Fractal_S-2)))/(1-pow(b,2*(Fractal_S-2)*nf)));
	for(i=0;i<N;i++)
	{
		z[i]=0.0;
		for(k=0;k<nf;k++)
		{
			z[i]=z[i]+pow(b,(Fractal_S-2)*k)*sin(kl*pow(b,k)*x[i]+r_n[k]);
		}
        z[i]=h*cn*z[i];
	}
	double dx=RL/N;
	////////////////////////////////////////////////////////////
	//计算一阶导数的第一种公式：用有限差分
	df_x[0]=(z[1]-z[N-1])/(2*dx);
	df_x[N-1]=(z[0]-z[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		df_x[i]=(z[i+1]-z[i-1])/(2*dx);
	}
	////////////////////////////////////////////////////////////
    //计算二阶导数的第一种公式：用有限差分
	ddf_x[0]=(df_x[1]-df_x[N-1])/(2*dx);
	ddf_x[N-1]=(df_x[0]-df_x[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		ddf_x[i]=(df_x[i+1]-df_x[i-1])/(2*dx);
	}/**/
	//cout<<rn<<endl;
	/////////////////////////////////////////////////////////////
	//计算二阶导数的第二种公式：书中公式（4.1.110）
	//用这种方法计算二阶导数较好
	///////////////////////////////////////////////////////////////////
	//display_sur();
}

void RoughSurface::ocean_surface()
{
	rand_n();
	int i;
	double *kx=new double[N];
	Complex *b=new Complex[N];
	
	x.assign(N);
	for(i=-N/2+1;i<N/2+1;i++)
	{
		x[i+N/2-1]=i*RL/N;
	}
	for(i=-N/2+1;i<=N/2;i++)
	{
		kx[i+N/2-1]=2*PI*i/RL;
	}
	b[N/2-1]=owk(kx[N/2-1])*r_n[N-2];
	//b[N/2-1].Imag=0.0;
	b[N-1]=owk(kx[N-1])*r_n[N-1];
    //b[N-1].Imag=0.0;
	for(i=0;i<N/2-1;i++)
	{
		b[i]=owk(kx[i])*(r_n[2*i]+i_*r_n[2*i+1])/sqrt(2);//b(-N/2+1)到b(-1)
		//b[i].Imag=owk(kx[i])*r_n[2*i+1]/sqrt(2);
	}
	for(i=N/2;i<N-1;i++)
	{
		//b[i]=b[i-2*(i+1-N/2)];
		b[i]=conj(b[i-2*(i+1-N/2)]);//b(1)到b(N/2-1)

		//b[i].Imag=-b[i-2*(i+1-N/2)].Imag;//复共轭
	}
	Complex *xt=new Complex[N];
	for(i=0;i<N;i++)
	{
		if(i<=N/2)
		{
			xt[i]=b[i+N/2-1];

		}
		else
		{
			xt[i]=b[i-N/2-1];
		}
	}
	Complex *ft=new Complex[N];
	z.assign(N);
	df_x.assign(N);
	ddf_x.assign(N);
	fft(N,xt,ft,-1);
	
	for(i=0;i<N;i++)
	{
		if(i<=N/2-2)
		{
			z[i]=real(ft[i+N/2+1])*N/RL;
		}
		else
		{
			z[i]=real(ft[i-N/2+1])*N/RL;
		}
		//cout<<x[i]<<"     "<<f_r[i]<<endl;
	}
	double dx=RL/N;
	df_x[0]=(z[1]-z[N-1])/(2*dx);
	df_x[N-1]=(z[0]-z[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		df_x[i]=(z[i+1]-z[i-1])/(2*dx);
	}
	//cur_draft();
    //计算二阶导数的第一种公式：用有限差分
	ddf_x[0]=(df_x[1]-df_x[N-1])/(2*dx);
	ddf_x[N-1]=(df_x[0]-df_x[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		ddf_x[i]=(df_x[i+1]-df_x[i-1])/(2*dx);
	}
	//计算二阶导数的第二种公式：书中公式（4.1.110）
	delete [] b;
	delete [] kx;
	delete [] xt;
	delete [] ft;
	show_rs();

}

void RoughSurface::ocean_surface_PM(double U195)
{
	rand_n();
	int i;
	double *kx=new double[N];
	Complex *b=new Complex[N];
	
	x.assign(N);
	for(i=-N/2+1;i<N/2+1;i++)
	{
		x[i+N/2-1]=i*RL/N;
	}
	for(i=-N/2+1;i<=N/2;i++)
	{
		kx[i+N/2-1]=2*PI*i/RL;
	}
	b[N/2-1]=owk_PM(kx[N/2-1],U195)*r_n[N-2];
	//b[N/2-1].Imag=0.0;
	b[N-1]=owk_PM(kx[N-1],U195)*r_n[N-1];
    //b[N-1].Imag=0.0;
	for(i=0;i<N/2-1;i++)
	{
		b[i]=owk_PM(kx[i],U195)*(r_n[2*i]+i_*r_n[2*i+1])/sqrt(2);//b(-N/2+1)到b(-1)
		//b[i].Imag=owk(kx[i])*r_n[2*i+1]/sqrt(2);
	}
	for(i=N/2;i<N-1;i++)
	{
		//b[i]=b[i-2*(i+1-N/2)];
		b[i]=conj(b[i-2*(i+1-N/2)]);//b(1)到b(N/2-1)
		//b[i].Imag=-b[i-2*(i+1-N/2)].Imag;//复共轭
	}
	Complex *xt=new Complex[N];
	for(i=0;i<N;i++)
	{
		if(i<=N/2)
		{
			xt[i]=b[i+N/2-1];
		}
		else
		{
			xt[i]=b[i-N/2-1];
		}
	}
	Complex *ft=new Complex[N];
	z.assign(N);
	df_x.assign(N);
	ddf_x.assign(N);
	fft(N,xt,ft,-1);
	
	for(i=0;i<N;i++)
	{
		if(i<=N/2-2)
		{
			z[i]=real(ft[i+N/2+1])*N/RL;
		}
		else
		{
			z[i]=real(ft[i-N/2+1])*N/RL;
		}
		//cout<<x[i]<<"     "<<f_r[i]<<endl;
	}
	double dx=RL/N;
	df_x[0]=(z[1]-z[N-1])/(2*dx);
	df_x[N-1]=(z[0]-z[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		df_x[i]=(z[i+1]-z[i-1])/(2*dx);
	}
	//cur_draft();
    //计算二阶导数的第一种公式：用有限差分
	ddf_x[0]=(df_x[1]-df_x[N-1])/(2*dx);
	ddf_x[N-1]=(df_x[0]-df_x[N-2])/(2*dx);
	for(i=1;i<N-1;i++)
	{
		ddf_x[i]=(df_x[i+1]-df_x[i-1])/(2*dx);
	}
	//计算二阶导数的第二种公式：书中公式（4.1.110）
	delete [] b;
	delete [] kx;
	delete [] xt;
	delete [] ft;
	show_rs();

}

double RoughSurface::owk_PM(double kx,double U195)
{
     double wk;
	 double a0 = 0.0081;
	 double bita = 0.74;
	 double g = 9.81;

    if (kx==0)
    {
		wk=0;
    }
     else
	 {
	   wk=sqrt(2.0*PI*RL*a0/(4.0*pow(abs(kx),3.0)*exp(-bita*g*g/(kx*kx*pow(U195,4.0)))));

	 }

	 return wk;



}
double RoughSurface::owk(double kx)
{
	double wk;
	//wk为ocean spectrum;
	kx=fabs(kx);
	if(kx>kl && kx<ku)
	{
//	cout<<ku<<"  "<<kl<<endl;

		double kr=fabs(kx);
    	double a0=0.004;
    	double a=0.225;
    	double b=1.25;
    	double kj=2.0;
        double gs;
    	gs=9.81+(7.25e-5)*pow(kx,2);
    	double z0;
    	z0=(6.84e-5)/us+(4.28e-3)*pow(us,2)-4.43e-4;
     	double U195;
     	U195=us/0.4*log(19.5/z0);
//		cout<<"U195="<<U195;
//		U195=5.;
    	double kc;
    	kc=9.81/pow(U195,2);
    	if(kr>kj)
		{
			wk=sqrt(2.0*PI*RL*(a0/pow(kr,3))*pow(b*kr*pow(us,2)/gs,a*log10(kr/kj)));
		}
    	else
		{
			wk=sqrt(2.0*PI*RL*(a0/pow(kr,3))*exp(-0.74*pow(kc/kr,2)));
		}
	}

	else
	{
		wk=0.0;
	}
									

	return wk;
}
void RoughSurface::Two_D_rand_n(int nxy,double mmean,double vvari)
{
	int i,k,m;
	r_n.assign(nxy*nxy);
	double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for(k=0;k<nxy*nxy;k++)
	{
		t=0.0;
		//cout<<"第i个random number"<<endl;
        for(i=1;i<=12;i++)
        {
			seed=seed*w+v;
			m=(int)(seed/s);
            seed=seed-m*s;
			t=t+seed/s;
		}
        r_n[k]=mmean+vvari*(t-6.0);
		r_n[k]=2.0*PI*r_n[k];
		//cout<<r_n[k]<<endl;
    }
}
void RoughSurface::Two_D_Gauss_Surface(int nxy,double mm,double vv,double hh,double llx,double lly,double rl)
{	
	Two_D_rand_n(nxy,mean,vari);
	int i,j,u,v;
	Vec_Com bri;
	bri.assign(nxy*nxy);
	Vec_Dou kx,ky;
	kx.assign(nxy);
	ky.assign(nxy);
	x.assign(nxy*nxy);
	y.assign(nxy*nxy);
	z.assign(nxy*nxy);
	df_x.assign(nxy*nxy);
	ddf_x.assign(nxy*nxy);
	df_y.assign(nxy*nxy);
	ddf_y.assign(nxy*nxy);
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			x[i*nxy+j]=-RL/2.0+i*RL/nxy+RL/nxy/2.0;
		    y[i*nxy+j]=-RL/2.0+j*RL/nxy+RL/nxy/2.0;
		}
	}
	for(i=-nxy/2+1;i<=nxy/2;i++)
		kx[i+nxy/2-1]=2.0*PI*i/RL;
	for(i=-nxy/2+1;i<=nxy/2;i++)
		ky[i+nxy/2-1]=2.0*PI*i/RL;
	/////////////////////////////////////////////////
	//系数b_r,b_i
	int count;
	count=0;
	for(i=-nxy/2+1; i<0; i++)
	{
		for(j=-nxy/2+1; j<0; j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])/sqrt(2.0)*(r_n[count]+i_*r_n[count+1]);
			count+=2;
		}
	}
	for(i=1;i<nxy/2;i++)
	{
		for(j=1;j<nxy/2;j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			v=(-i+nxy/2-1)*nxy+(-j+nxy/2-1);
			bri[u]=conj(bri[v]);
		}
	}
	for(i=-nxy/2+1;i<0;i++)
	{
		for(j=1;j<nxy/2;j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])/sqrt(2.0)*(r_n[count]+i_*r_n[count+1]);
			count+=2;
		}
	}
	for(i=1;i<nxy/2;i++)
	{
		for(j=-nxy/2+1;j<0;j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			v=(-i+nxy/2-1)*nxy+(-j+nxy/2-1);
			bri[u]=conj(bri[v]);
		}
	}
	////////
	i=0;
	for(j=-nxy/2+1;j<nxy/2;j++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	i=nxy/2;
	for(j=-nxy/2+1;j<nxy/2;j++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	/////////
	j=nxy/2;
	for(i=-nxy/2+1;i<=nxy/2;i++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	j=0;
	for(i=-nxy/2+1;i<=nxy/2;i++)
	{
		if(i!=0 && i!=nxy/2)
		{
			
	    	u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
			count++;
		}
	}
	//cout<<count<<endl;
	////////////对系数的处理（包括周期扩展和变换到逆离散傅里叶变换的系数）//////////////
	Vec_Com xt1,xt2;
	xt1.assign(nxy*nxy);
	xt2.assign(nxy*nxy);
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(i<=nxy/2)
			{
				xt1[u]=bri[(i+nxy/2-1)*nxy+j];
			}
			else
			{
				xt1[u]=bri[(i-nxy/2-1)*nxy+j];
			}
			//cout<<xt1[u]<<"    "<<u<<endl;
		}
	}
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(j<=nxy/2)
			{
				xt2[u]=xt1[i*nxy+j+nxy/2-1];
			}
			else
			{
				xt2[u]=xt1[i*nxy+j-nxy/2-1];
			}
			//cout<<xt2[u]<<"    "<<u<<endl;
		}
	}
	/////////////////////////////////////////////////
	//fft2(nxy,xt2,xt2,-1);
    double temp;
	Vec_Com zt;
	zt.assign(nxy*nxy);
	int k,w;
	//time_t ltime3,ltime4,ltime5;
	//time( &ltime3 );
	i=0;
	while( i < nxy*nxy )
	{
	    for(k=0; k<nxy; k++)
		{
			for(w=0; w<nxy; w++)
			{
				v=k*nxy+w;
				temp=2.0*PI*(k*x[i]/RL+w*y[i]/RL);
				zt[i] = zt[i] + xt2[v]*exp(i_*temp)/(RL*RL);
				//zt[i]=zt[i]+real(xt2[v])*cos(temp)-imag(xt2[v])*sin(temp);
			}
		}
		//z[i]=real(zt[i]);
		//cout<<zt[i]<<"    "<<i<<endl;
		i++;
	}
	//time( &ltime4 );
	//ltime5=ltime4-ltime3;
	//cout<<ltime5<<" seconds"<<endl;
	//////////////对DFT结果的处理（周期扩展）
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(i<=nxy/2-2)
			{
				xt1[u]=zt[(i+nxy/2+1)*nxy+j];
			}
	    	else
			{
				xt1[u]=zt[(i-nxy/2+1)*nxy+j];
			}
		}
	}
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(j<=nxy/2-2)
			{
				z[u]=real(xt1[i*nxy+j+nxy/2+1]);
			}
	    	else
			{
				z[u]=real(xt1[i*nxy+j-nxy/2+1]);
			}
		}
	}
	///////////用有限差分法计算偏导数/////////////
	//dfx
	double dxy;
	dxy=RL/nxy;
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(i == 0)
				df_x[u]=(z[1*nxy+j]-z[(nxy-1)*nxy+j])/(2*dxy);
			else if(i == (nxy-1))
				df_x[u]=(z[0*nxy+j]-z[(nxy-2)*nxy+j])/(2*dxy);
			else
				df_x[u]=(z[(i+1)*nxy+j]-z[(i-1)*nxy+j])/(2*dxy);
		}
	}
	//dfy
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(j == 0)
				df_y[u]=(z[i*nxy+1]-z[i*nxy+(nxy-1)])/(2*dxy);
			else if(j == (nxy-1))
				df_y[u]=(z[i*nxy+0]-z[i*nxy+(nxy-2)])/(2*dxy);
			else
				df_y[u]=(z[i*nxy+j+1]-z[i*nxy+j-1])/(2*dxy);
		}
	}
	//ddfx
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(i == 0)
				ddf_x[u]=(df_x[1*nxy+j]-df_x[(nxy-1)*nxy+j])/(2*dxy);
			else if(i == (nxy-1))
				ddf_x[u]=(df_x[0*nxy+j]-df_x[(nxy-2)*nxy+j])/(2*dxy);
			else
				ddf_x[u]=(df_x[(i+1)*nxy+j]-df_x[(i-1)*nxy+j])/(2*dxy);
		}
	}
	//ddfy
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			if(j == 0)
				ddf_y[u]=(df_x[i*nxy+1]-df_x[i*nxy+(nxy-1)])/(2*dxy);
			else if(j == (nxy-1))
				ddf_y[u]=(df_x[i*nxy+0]-df_x[i*nxy+(nxy-2)])/(2*dxy);
			else
				ddf_y[u]=(df_x[i*nxy+j+1]-df_x[i*nxy+j-1])/(2*dxy);
		}
	}
	for(i=0;i<nxy*nxy;i++)
	{
		outx<<x[i]<<endl;
		outy<<y[i]<<endl;
		outz<<z[i]<<endl;
		outdfx<<df_x[i]<<endl;
		outdfy<<df_y[i]<<endl;
		outddfx<<ddf_x[i]<<endl;
		outddfy<<ddf_y[i]<<endl;
	}

	bri.clear();
	kx.clear();
	ky.clear();
	xt1.clear();
	xt2.clear();
	zt.clear();
	
}
void RoughSurface::TD_Gauss_Surface(int nxy,double mm,double vv,double hh,double llx,double lly,double rl)
{
	Two_D_rand_n(nxy,mean,vari);
	
	int i,j,u,v;
/*	for (i=0; i<nxy*nxy; i++)
	{
		outrs<<r_n[i]<<"  "<<i<<endl;
	}*/
	
	Vec_Com bri;
	bri.assign(nxy*nxy);
	Vec_Dou kx,ky;
	kx.assign(nxy);
	ky.assign(nxy);
	x.assign(nxy*nxy);
	y.assign(nxy*nxy);
	z.assign(nxy*nxy);
	df_x.assign(nxy*nxy);
	ddf_x.assign(nxy*nxy);
	df_y.assign(nxy*nxy);
	ddf_y.assign(nxy*nxy);
	ddf_xy.assign(nxy*nxy);
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			x[i*nxy+j]=-RL/2.0+i*RL/nxy+RL/nxy/2.0;
		    y[i*nxy+j]=-RL/2.0+j*RL/nxy+RL/nxy/2.0;
		}
	}
	for(i=-nxy/2+1;i<=nxy/2;i++)
		kx[i+nxy/2-1]=2*PI*i/RL;
	for(i=-nxy/2+1;i<=nxy/2;i++)
		ky[i+nxy/2-1]=2*PI*i/RL;
	//////////////////////////////////////////////////////////////////////////
	//输出kx,ky
/*	for (i=-nxy/2+1; i<=nxy/2; i++)
	{
		for (j=-nxy/2+1; j<=nxy/2; j++)
		{
			outresult<<ky[i+nxy/2-1]<<"  ";
		}
		outresult<<endl;
	}*/
	/////////////////////////////////////////////////
	//系数bri
	int count;
	count=0;
	for(i=-nxy/2+1; i<0; i++)
	{
		for(j=-nxy/2+1; j<0; j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])/sqrt(2.0)*(r_n[count]+i_*r_n[count+1]);
			count+=2;
		}
	}
	for(i=1; i<nxy/2; i++)
	{
		for(j=1; j<nxy/2; j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			v=(-i+nxy/2-1)*nxy+(-j+nxy/2-1);
			bri[u]=conj(bri[v]);
		}
	}
	for(i=-nxy/2+1; i<0; i++)
	{
		for(j=1; j<nxy/2; j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])/sqrt(2.0)*(r_n[count]+i_*r_n[count+1]);
			count+=2;
		}
	}
	for(i=1; i<nxy/2; i++)
	{
		for(j=-nxy/2+1; j<0; j++)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			v=(-i+nxy/2-1)*nxy+(-j+nxy/2-1);
			bri[u]=conj(bri[v]);
		}
	}
	////////
	i=0;
	for(j=-nxy/2+1; j<nxy/2; j++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	i=nxy/2;
	for(j=-nxy/2+1; j<nxy/2; j++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	/////////
	j=nxy/2;
	for(i=-nxy/2+1; i<=nxy/2; i++)
	{
		u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
		bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
		count++;
	}
	j=0;
	for(i=-nxy/2+1; i<nxy/2; i++)
	{
		if(i!=0)
		{
			u=(i+nxy/2-1)*nxy+(j+nxy/2-1);
			bri[u]=Two_D_wk(kx[i+nxy/2-1],ky[j+nxy/2-1])*r_n[count]/sqrt(2.0);
			count++;
		}
	}
	//cout<<count<<endl;
	////////////对系数的处理（包括周期扩展和变换到逆离散傅里叶变换的系数）//////////////
	Vec_Com xt1;
	xt1.assign(nxy*nxy);
	Complex *xt2=new Complex[nxy*nxy];
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i <= nxy/2)
			{
				xt1[u]=bri[(i+nxy/2-1)*nxy+j];
			}
			else
			{
				xt1[u]=bri[(i-nxy/2-1)*nxy+j];
			}
			//cout<<xt1[u]<<"    "<<u<<endl;
		}
	}
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j <= nxy/2)
			{
				xt2[u]=xt1[i*nxy+j+nxy/2-1];
			}
			else
			{
				xt2[u]=xt1[i*nxy+j-nxy/2-1];
			}
			//cout<<xt2[u]<<"    "<<u<<endl;
		}
	}
	
	fft2(nxy,xt2,xt2,-1);
	
	//////////////对DFT结果的处理（周期扩展）
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i <= nxy/2-2)
			{
				xt1[u]=xt2[(i+nxy/2+1)*nxy+j];
			}
	    	else
			{
				xt1[u]=xt2[(i-nxy/2+1)*nxy+j];
			}
		}
	}
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j <= nxy/2-2)
			{
				z[u]=real(xt1[i*nxy+j+nxy/2+1])*(nxy*nxy)/(RL*RL);
			}
	    	else
			{
				z[u]=real(xt1[i*nxy+j-nxy/2+1])*(nxy*nxy)/(RL*RL);
			}
		}
	}
	///////////用有限差分法计算偏导数/////////////
	//dfx
	double dxy;
	dxy=RL/nxy;
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i == 0)
				df_x[u]=(z[1*nxy+j]-z[(nxy-1)*nxy+j])/(2.0*dxy);
			else if(i == (nxy-1))
				df_x[u]=(z[0*nxy+j]-z[(nxy-2)*nxy+j])/(2.0*dxy);
			else
				df_x[u]=(z[(i+1)*nxy+j]-z[(i-1)*nxy+j])/(2.0*dxy);
		}
	}
	//dfy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				df_y[u]=(z[i*nxy+1]-z[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				df_y[u]=(z[i*nxy+0]-z[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				df_y[u]=(z[i*nxy+j+1]-z[i*nxy+j-1])/(2.0*dxy);
		}
	}
	//ddfx
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i == 0)
				ddf_x[u]=(df_x[1*nxy+j]-df_x[(nxy-1)*nxy+j])/(2.0*dxy);
			else if(i == (nxy-1))
				ddf_x[u]=(df_x[0*nxy+j]-df_x[(nxy-2)*nxy+j])/(2.0*dxy);
			else
				ddf_x[u]=(df_x[(i+1)*nxy+j]-df_x[(i-1)*nxy+j])/(2.0*dxy);
		}
	}
	//ddfy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				ddf_y[u]=(df_y[i*nxy+1]-df_y[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				ddf_y[u]=(df_y[i*nxy+0]-df_y[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				ddf_y[u]=(df_y[i*nxy+j+1]-df_y[i*nxy+j-1])/(2.0*dxy);
		}
	}
	//ddfxy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				ddf_xy[u]=(df_x[i*nxy+1]-df_x[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				ddf_xy[u]=(df_x[i*nxy+0]-df_x[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				ddf_xy[u]=(df_x[i*nxy+j+1]-df_x[i*nxy+j-1])/(2.0*dxy);
		}
	}
/*	for(i=0; i<nxy*nxy; i++)
	{
		outx<<x[i]<<endl;
		outy<<y[i]<<endl;
		outz<<z[i]<<endl;
    	outrs<<x[i]<<" "<<y[i]<<" "<<z[i]<<endl;
//      outdfx<<df_x[i]<<endl;
//		outdfy<<df_y[i]<<endl;
//		outddfx<<ddf_x[i]<<endl;
//		outddfy<<ddf_y[i]<<endl;
	}*/

	bri.clear();
	kx.clear();
	ky.clear();
	xt1.clear();
	delete [] xt2;
}
double RoughSurface::Two_D_wk(double kx,double ky)
{
	double wk;
	//wk为Gaussian spectrum;
	wk=2.0*PI*sqrt(RL*RL*lx*ly*h*h/4.0/PI*exp(-kx*kx*lx*lx/4.0-ky*ky*ly*ly/4.0));
	return wk;
}
//////////////////////////////////////////////////////////////
/////////////////////fractal rough surface////////////////////
void RoughSurface::Two_D_Frctal_Rand_n(int nnf,double mmean,double vvari)
{
    int i,k,m;
	r_n.assign(3*nnf);
	double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for(k=0;k<3*nnf;k++)
	{
		t=0.0;
		//cout<<"第i个random number"<<endl;
        for(i=1;i<=12;i++)
        {
			seed=seed*w+v;
			m=(int)(seed/s);
            seed=seed-m*s;
			t=t+seed/s;
		}
        r_n[k]=mmean+vvari*(t-6.0);
		r_n[k]=2.0*PI*r_n[k];
		//cout<<r_n[k]<<endl;
    }
}
void RoughSurface::Two_D_Fractal_surface(int nxy,int nnf,double hh,double mm,double vv,double rrl,double kkl,double kku,double fra_s)
{
	
	///随机数的产生
	Two_D_Frctal_Rand_n(nnf,mm,vv);
   	int i,j,u,k;
	x.assign(nxy*nxy);
	y.assign(nxy*nxy);
	for(i=0;i<nxy;i++)
	{
		for(j=0;j<nxy;j++)
		{
			u=i*nxy+j;
			x[u]=-rrl/2+i*rrl/nxy+RL/nxy/2;
		    y[u]=-rrl/2+j*rrl/nxy+rrl/nxy/2;
		}
	}
	z.assign(nxy*nxy);
	df_x.assign(nxy*nxy);
	ddf_x.assign(nxy*nxy);
	df_y.assign(nxy*nxy);
	ddf_y.assign(nxy*nxy);
	ddf_xy.assign(nxy*nxy);
	double b=pow(kku/kkl,1.0/(nnf-1.0));
	double cn=sqrt(2.0*(1.0-pow(b,2.0*(fra_s-2.0)))/(1.0-pow(b,2.0*(fra_s-2.0)*nnf)));
	for(i=0;i<nxy*nxy;i++)
	{
		z[i]=0.0;
		for(k=0;k<nnf;k++)
		{
			z[i] += r_n[3*k]*pow(b,(fra_s-2)*k)*sin(kl*pow(b,k)*(x[i]*cos(r_n[3*k+1])+y[i]*sin(r_n[3*k+1]))+r_n[3*k+2]);
		}
        z[i] = hh*cn*z[i];
	}

	///////////用有限差分法计算偏导数/////////////
	//dfx
	double dxy;
	dxy=rrl/nxy;
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i == 0)
				df_x[u]=(z[1*nxy+j]-z[(nxy-1)*nxy+j])/(2.0*dxy);
			else if(i == (nxy-1))
				df_x[u]=(z[0*nxy+j]-z[(nxy-2)*nxy+j])/(2.0*dxy);
			else
				df_x[u]=(z[(i+1)*nxy+j]-z[(i-1)*nxy+j])/(2.0*dxy);
		}
	}
	//dfy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				df_y[u]=(z[i*nxy+1]-z[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				df_y[u]=(z[i*nxy+0]-z[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				df_y[u]=(z[i*nxy+j+1]-z[i*nxy+j-1])/(2.0*dxy);
		}
	}
	//ddfx
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i == 0)
				ddf_x[u]=(df_x[1*nxy+j]-df_x[(nxy-1)*nxy+j])/(2.0*dxy);
			else if(i == (nxy-1))
				ddf_x[u]=(df_x[0*nxy+j]-df_x[(nxy-2)*nxy+j])/(2.0*dxy);
			else
				ddf_x[u]=(df_x[(i+1)*nxy+j]-df_x[(i-1)*nxy+j])/(2.0*dxy);
		}
	}
	//ddfy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				ddf_y[u]=(df_y[i*nxy+1]-df_y[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				ddf_y[u]=(df_y[i*nxy+0]-df_y[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				ddf_y[u]=(df_y[i*nxy+j+1]-df_y[i*nxy+j-1])/(2.0*dxy);
		}
	}
	//ddfxy
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(j == 0)
				ddf_xy[u]=(df_x[i*nxy+1]-df_x[i*nxy+(nxy-1)])/(2.0*dxy);
			else if(j == (nxy-1))
				ddf_xy[u]=(df_x[i*nxy+0]-df_x[i*nxy+(nxy-2)])/(2.0*dxy);
			else
				ddf_xy[u]=(df_x[i*nxy+j+1]-df_x[i*nxy+j-1])/(2.0*dxy);
		}
	}

/*	for(i=0; i<nxy*nxy; i++)
	{
		outx<<x[i]<<endl;
		outy<<y[i]<<endl;
		outz<<z[i]<<endl;
//		outrs<<x[i]/wave_len<<" "<<y[i]/wave_len<<" "<<z[i]/wave_len<<endl;
//      outdfx<<df_x[i]<<endl;
//		outdfy<<df_y[i]<<endl;
//		outddfx<<ddf_x[i]<<endl;
//		outddfy<<ddf_y[i]<<endl;
	}
	for(i=0;i<nxy*nxy;i++)
	{
		outrs<<x[i]<<"  "<<y[i]<<"   "<<z[i]<<endl;
	}*/

}
void RoughSurface::show_rs()
{
	int i;
	for(i=0;i<N;i++)
	{
		outrs<<x[i]/wave_len<<"     "<<z[i]<<endl;
		//cout<<z[i]<<endl;
	}
}