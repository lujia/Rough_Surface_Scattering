
#define H_post_computer

#include "matrix.h"


class post_cp: public matrix
{
public:
    post_cp();
	~post_cp();
	//一维问题
    void One_D_Conductor_Ocean_Js0(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);
    void One_D_Conductor_Ocean_Es0(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);
    void One_D_Conductor_Gauss_Es0(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);
    void One_D_Conductor_Gauss_Js0(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);
    void One_D_Conductor_Gauss_Jsd(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);
    void One_D_Conductor_Gauss(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps);

public:
	
	int n_sc;
	int N_Realization;
	double *sc_an;//scattering angles in degree 
	double *sc_p;//bistatic scattering coef	
	vector <double> sc_hl;

};

post_cp::post_cp()
{
	sc_an=NULL;
	sc_p=NULL;
	n_sc=0;
}
post_cp::~post_cp()
{
	delete [] sc_an;
	delete [] sc_p;
	sc_hl.clear();
}


void post_cp::One_D_Conductor_Ocean_Js0(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
	mean=0.;              //均值
	vari=1.;              //方差
	N=num;              //total number of sample points 
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;             //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;  //wave number
//	g_p=RL/4.;         //tapering parameter for incident wave 
	g_p = ggp;

	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;

		ocean_surface();
		One_D_Conductor_BMIA_CAG_Matrix_Js0();

}
void post_cp::One_D_Conductor_Ocean_Es0(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
	mean=0.;              //均值
	vari=1.;              //方差
	N=num;              //total number of sample points 
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;             //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;  //wave number
//	g_p=RL/4.;         //tapering parameter for incident wave
	g_p = ggp;
	
	
	ocean_surface();
	One_D_Conductor_BMIA_CAG_Matrix_Es0();


	
}

void post_cp::One_D_Conductor_Gauss_Js0(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
	mean=0.;              //均值
	vari=1.;              //方差
	N=num;                  //total number of sample points 
	h = hh;
	lc = llc;
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;                  //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;      //incident angle in degree
	N_Realization=nn_r;      //total number of surface realizations 
	n_sc=nn_sc;              //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;       //wave number
//	g_p=RL/4.;               //tapering parameter for incident wave
	g_p = ggp;
	
	gauss_surface();
	One_D_Conductor_BMIA_CAG_Matrix_Js0();


}


void post_cp::One_D_Conductor_Gauss_Es0(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
	mean=0.;              //均值
	vari=1.;              //方差
	N=num;              //total number of sample points 
	h = hh;
	lc = llc;
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;             //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;  //wave number
//	g_p=RL/4.;         //tapering parameter for incident wave
	g_p = ggp;
	
	

	gauss_surface();
	One_D_Conductor_BMIA_CAG_Matrix_Es0();

}


void post_cp::One_D_Conductor_Gauss_Jsd(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
	mean=0.;              //均值
	vari=1.;              //方差
	N=num;              //total number of sample points 
	h = hh;
	lc = llc;
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;             //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;  //wave number
//	g_p=RL/4.;         //tapering parameter for incident wave
	g_p = ggp;

	Vec_Com vss;
	  vss.assign(N);

	double d=10.*wave_len;
	
	double xsz,ww,kg,fac;
    for(int i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+d*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		vss[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		
	}
	
	

	gauss_surface();
	One_D_Conductor_BMIA_CAG_Matrix_Jsd(vss,Jsdd);

}


void post_cp::One_D_Conductor_Gauss(int num,double hh,double llc,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double ggp,double eeps)
{
    mean=0.;              //均值
	vari=1.;              //方差
	N=num;              //total number of sample points 
	h = hh;
	lc = llc;
	seed=123456.;            //seed for random number generator 
	wave_len=wwave_len;
	RL=rrl;             //rough surface length
	us=uss;
	kl=kkl;
	ku=kku;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	k_w=2*PI/wave_len;  //wave number
//	g_p=RL/4.;         //tapering parameter for incident wave
	g_p = ggp;
	
	

	gauss_surface();
}