#ifndef H_post_computer
#define H_post_computer
//#include "FIPWA_OD_Cond.h"
//#include "FIPWA_OD_Diel.h"
//#include "FIPWA_TD_Cond.h"
//#include "FIPWA_TWO_D_Cond.h"
#include "matrix.h"

ofstream out("scat_para.txt",ios::trunc);
ofstream outdb("scat_paradb.txt",ios::trunc);
ofstream outrcs("rcs.txt",ios::trunc);
ofstream outrcsdb("rcsdb.txt",ios::trunc);

ofstream out_inci_bt("inci_bt.txt",ios::trunc);
ofstream out_current_mag("current_mag.txt",ios::trunc);
ofstream outre_hh("scat_para_hh.txt",ios::trunc);
ofstream outredb_hh("scat_paradb_hh.txt",ios::trunc);
ofstream outre_vh("scat_para_vh.txt",ios::trunc);
ofstream outredb_vh("scat_paradb_vh.txt",ios::trunc);
//ofstream out_Reflectivity("Reflectivity.txt",ios::trunc);

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



	void One_D_Conductor_Scatter_P(int num,double hh,double llc,double rrl,double wave_le,double ii_an,int nn_r,int nn_sc,int bbw,double eeps,double crit);
	void One_D_Conductor_Ocean_Scatter(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double eeps);
	void One_D_Medium_Scatter_P(double gg_p,int num,double hh,double llc,double rrl,double wave_le,
						double ii_an,int nn_r,int nn_sc,int rrd_bw,int nn1,Complex rr_em,Complex permit,double eeps,double crit);
	void Inci_Br_Temp(double gg_p,int num,double hh,double llc,double rrl,double wave_le,
						double i_an_down,double i_an_up,int nn_r,int nn_sc,int rrd_bw,int nn1,Complex rr_em,Complex permit,double eeps,double crit);
   	double Bright_Temperature();//u,g,n,seed,h,lc,rl,wave_len,i_an,n_r,n_sc,bw,r_em,permittivity
	double Reflectivity_Calcu(double *Sc,int m);
	double Bright_Temperature(double *Sc,int m);
	void show_result();
	//FIPWA
	void One_D_Conductor_Gauss_Scatter_FIPWA(int num,double hh,double llc,double rrl,double wwave_len,double ii_an,double eeps);
	void One_D_Dielectric_Gauss_Scatter_FIPWA(int num,double hh,double llc,double rrl,double wwave_len,double ii_an,Complex rr_em,Complex permit,double eeps);
    void One_D_Dielectric_Ocean_Scatter_FIPWA(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,Complex rr_em,Complex permit,double eeps);
	

	//二维问题
	//标量波入射二维导体粗糙面
    void Two_D_Conductor_Gauss_Scatter_P(int nxy,double rrl,double rrd,int Para_Num,double hh,double llx,double lly,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit);
   	void Two_D_Conductor_Fractal_Scatter_P(int nxy,int nnf,double rrl,double rrd,int Para_Num,double hh,double kkl,double kku,double fra_s,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit);

	//矢量波入射二维导体粗糙面
	void Two_D_Conductor_Vector_Gauss_Scatter_P(int nxy,double rrl,double rrd,int Para_Num,double hh,double llx,double lly,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit,int po_h,int po_v);

	//分形粗糙面散射系数的计算
	void Two_D_Conductor_Vector_Fractal_Scatter_P(int nxy,int nnf,double rrl,double rrd,int Para_Num,double hh,double kkl,double kku,double fra_s,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit,int po_h,int po_v);

	//从文本文档读取数据
	void Read_From_Text(char *path,int N_total);
public:
	double *sc_an;//scattering angles in degree 
	double *sc_p;//bistatic scattering coef
	vector <double> RCS;
	int n_sc;
	int N_Realization;
	vector <double> sc_hl;
	double sc_az;
	vector <double> sc_p_hh,sc_p_vh;
};
#endif
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

/*void post_cp::One_D_Conductor_Gauss_Scatter_FIPWA(int num,double hh,double llc,double rrl,double wwave_len,double ii_an,double eeps)
{
	mean = 0.0;           //均值
	vari = 1.0;           //方差
	seed = 123456.0;      //seed for random number generator 
	N = num;
	h = hh;
	lc = llc;
	RL = rrl;
	wave_len = wwave_len;
	//入射角：在FIPWA中是如何表示的???
	i_an = ii_an*PI/180.0;
	eps = eeps;
	k_w = 2.0 * PI / wave_len;
	//对于一维的情况，g用叶红霞结果不正确
	g_p=RL/4.0;         //tapering parameter for incident wave
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果
	//g_p=6.0/pow(cos(i_an),1.5);
	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;
	cout<<"一维标量锥形波的入射功率= "<<OD_Taper_Pow()<<endl;

	Get_list();
	CalcZijMatrix();
	CalcExci();
	CG();
	CalcRCS();
}
void post_cp::One_D_Dielectric_Gauss_Scatter_FIPWA(int num,double hh,double llc,double rrl,double wwave_len,double ii_an,Complex rr_em,Complex permit,double eeps)
{
	mean = 0.0;           //均值
	vari = 1.0;           //方差
	seed = 123456.0;      //seed for random number generator 
	N = num;
	h = hh;
	lc = llc;
	RL = rrl;
	wave_len = wwave_len;
	fre=3.0e8/wave_len;
	//入射角：在FIPWA中是如何表示的???
	i_an = ii_an*PI/180.0;
	r_em=rr_em;
	permittivity=permit;
	eps = eeps;
	k_w = 2.0 * PI / wave_len;
	k_w1=2.0*PI*fre*sqrt(((4.0*PI)*1.0e-7)*(1.0/(36.0*PI*1.0e9))*permittivity);//区域1的介电常数
	//对于一维的情况，g用叶红霞结果不正确
	g_p=RL/4.0;         //tapering parameter for incident wave
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果
	//g_p=6.0/pow(cos(i_an),1.5);
	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;
	cout<<"一维标量锥形波的入射功率= "<<OD_Taper_Pow()<<endl;

	Get_list();
	CalcZijMatrix();
	CalcExci();
	CG();
	CalcRCS();
}
void post_cp::One_D_Dielectric_Ocean_Scatter_FIPWA(int num,double uss,double kkl,double kku,double rrl,double wwave_len,
												   double ii_an,Complex rr_em,Complex permit,double eeps)
{
	mean = 0.0;           //均值
	vari = 1.0;           //方差
	seed = 123456.0;      //seed for random number generator 
	N = num;
	us = uss;
	kl = kkl;
	ku = kku;
	RL = rrl;
	wave_len = wwave_len;
	fre=3.0e8/wave_len;
	
	i_an = ii_an*PI/180.0;//入射角
	r_em=rr_em;
	permittivity=permit;
	eps = eeps;
	k_w = 2.0 * PI / wave_len;
	k_w1=2.0*PI*fre*sqrt(((4.0*PI)*1.0e-7)*(1.0/(36.0*PI*1.0e9))*permittivity);//区域1的介电常数
	//对于一维的情况，g用叶红霞结果不正确
	g_p=RL/4.0;         //tapering parameter for incident wave
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果
	//g_p=6.0/pow(cos(i_an),1.5);
	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;
	cout<<"一维标量锥形波的入射功率= "<<OD_Taper_Pow()<<endl;

	Get_list();
	CalcZijMatrix();
	CalcExci();
	CG();
	CalcRCS();
}*/
void post_cp::Read_From_Text(char *path,int N_total)
{
	int i;
	Complex current;
	Vec_Com Result_Total;
	ifstream filein;
	filein.open(path);
	for( i=0; i<N_total; i++ )
	{
		filein>>current;
		Result_Total.push_back(current);
	}
	filein.close();
	for (i=0; i<N_total; i++)
	{
		outresult<<Result_Total[i]<<"   "<<i<<endl;
	}
}
void post_cp::Two_D_Conductor_Gauss_Scatter_P(int nxy,double rrl,double rrd,int Para_Num,double hh,double llx,double lly,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit)
{
	mean=0.0;              //均值
	vari=1.0;              //方差
	seed=123456.0;         //seed for random number generator 
	Nxy=nxy;
	RL=rrl;
	RD = rrd;
	h=hh;
	lx=llx;
	ly=lly;
	wave_len=wwave_len;

	Inci_HL=in_hl*PI/180.0;//incident angle in degree
	Inci_AZ=in_az*PI/180.0;
	sc_az = ssc_az*PI/180.0;
	N_Realization = nn_r;
	n_sc = nn_sc;
	eps = eeps;
	criterion = crit;
	k_w = 2.0*PI/wave_len;  //wave number
	g_p = RL/3.0;         //tapering parameter for incident wave

	int i,j,k;
	double denom;
	//锥形波入射功率
	denom=1.0+pow(cos(Inci_HL),2)+2.0*pow(tan(Inci_HL),2);
	denom=denom/(2.0*pow(k_w*g_p*cos(Inci_HL),2));
	denom=8.0*pow(PI,3)*pow(g_p,2)*cos(Inci_HL)*(1.0-denom);
	cout<<denom<<endl;
	double dan=180.0/n_sc;
	sc_hl.assign(n_sc);
	sc_p=new double[n_sc];
	for(i=0; i<n_sc; i++)
	{
		sc_hl[i]=-90.0+(i+0.5)*dan;
	}
	double san;
	double tsc_az;
	double dxy=RL/Nxy;
	Complex integ;
	Complex ksc;
	double temp;
	int Num = Nxy * Nxy;
	for(i=0; i<N_Realization; i++)
	{
		TD_Gauss_Surface(Nxy,mean,vari,h,lx,ly,RL);
		
		cout<<"The "<<i+1<<"th random rough surface;"<<endl;
		Two_D_Conductor_Scalar_SMFIA_CAG_Matrix(Para_Num);//SMFSIA
		//Two_D_Conductor_Scalar_MoM_Matrix();//点配匹矩量法
		
		for(j=0; j<n_sc; j++)
		{
			san=sc_hl[j]*PI/180.0;
			if(san >= 0)
			{
				tsc_az = sc_az;
				integ = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					ksc = -x[k]*sin(san)*cos(tsc_az)-y[k]*sin(san)*sin(tsc_az)-z[k]*cos(san);
					integ += dxy*dxy*Uxy[k]*exp(i_*k_w*ksc);
				}
			}
			else
			{
				san = -san;//对于三维散射问题，散射角为与z轴的夹角，则在上半空间，散射角范围为0~90度
				tsc_az = PI + sc_az;
				integ = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					ksc = -x[k]*sin(san)*cos(tsc_az)-y[k]*sin(san)*sin(tsc_az)-z[k]*cos(san);
					integ += dxy*dxy*Uxy[k]*exp(i_*k_w*ksc);
				}
			}
			temp = pow(abs(integ),2);
			temp = temp/denom;//用入射波功率规范化
			sc_p[j] = (temp+i*sc_p[j])/(i+1.0);
		}
	}

	for(i=0; i<n_sc; i++)
	{
		outdb<<sc_hl[i]<<" "<<10*log10(sc_p[i])<<endl;
		out<<sc_hl[i]<<" "<<sc_p[i]<<endl;
	}
}
void post_cp::Two_D_Conductor_Vector_Gauss_Scatter_P(int nxy,double rrl,double rrd,int Para_Num,double hh,double llx,double lly,double wwave_len,double in_hl,
											   double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit,int po_h,int po_v)
{
	mean = 0.0;              //均值
	vari = 1.0;              //方差
	seed = 123456.0;         //seed for random number generator 
	Nxy = nxy;
	RL = rrl;
	RD = rrd;
	h = hh;
	lx = llx;
	ly = lly;
	wave_len = wwave_len;

	Inci_HL = in_hl*PI/180.0;//incident angle in degree
	Inci_AZ = in_az*PI/180.0;
	sc_az = ssc_az*PI/180.0;
	N_Realization = nn_r;
	n_sc = nn_sc;
	eps = eeps;
	criterion = crit;
	Po_H = po_h;
	Po_V = po_v;
	k_w = 2.0*PI/wave_len;  //wave number
	g_p = RL/3.0;         //tapering parameter for incident wave

	int i,j,k;
	double denom;
	//矢量锥形波入射功率：到底哪一种入射功率正确？？？
	denom = Incidence_Vector_Wave_Pinc_h(k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
	cout<<"第一种矢量锥形波入射功率= "<<denom<<endl;
	//平面波功率
//	denom = cos(Inci_HL)*RL*RL/(2.0*120.0*PI);
//	cout<<denom<<endl;
	double dan = 180.0/n_sc;
	sc_hl.assign(n_sc);
	sc_p_hh.assign(n_sc);
	sc_p_vh.assign(n_sc);
	for(i=0; i<n_sc; i++)
	{
		sc_hl[i]=-90.0+(i+0.5)*dan;
	}
	double san;
	double tsc_az;
	double dxy = RL/Nxy;
	int Num = Nxy*Nxy;
	Complex integ_hh,integ_vh;
	double ksc;
	double temp_hh,temp_vh;
	double Ph_inc;
	for(i=0; i<N_Realization; i++)
	{
		TD_Gauss_Surface(Nxy,mean,vari,h,lx,ly,RL);
	    cout<<"The "<<i+1<<"th random rough surface;"<<endl;
		Ph_inc = Two_D_Conductor_Vector_Ph_inc(nxy,k_w,Inci_HL,Inci_AZ,g_p);
		cout<<"第二种矢量锥形波入射功率= "<<Ph_inc<<endl;
		
		Two_D_Conductor_Vector_SMFIA_CAG_Matrix(Para_Num);
		//Two_D_Conductor_Vector_SMFIA_CAG_Matrix_Integ();
		//Two_D_Conductor_Vector_MoM_Matrix();//点配匹矩量法
		//Two_D_Conductor_Vector_MoM_Matrix_Integ();//矩阵元素采用积分的点配匹矩量法
		
		for(j=0; j<n_sc; j++)
		{
			san = sc_hl[j]*PI/180.0;
			if(san >= 0)
			{
				tsc_az = sc_az;//散射方位角为零的平面上
				integ_hh = 0.0;
				integ_vh = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					//水平极化的散射系数
					ksc = x[k]*sin(san)*cos(tsc_az) + y[k]*sin(san)*sin(tsc_az) + z[k]*cos(san);
	 		       	integ_hh = integ_hh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*sin(tsc_az) - Uxy[k+Num]*cos(tsc_az));	
					//垂直极化的散射系数
					integ_vh = integ_vh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*(df_x[k]*sin(san) - cos(san)*cos(tsc_az))
						 + Uxy[k+Num]*(df_y[k]*sin(san) - cos(san)*sin(tsc_az)));
				}
			}
			else
			{
				san = -san;
				tsc_az = PI + sc_az;
				integ_hh = 0.0;
				integ_vh = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					//水平极化的散射系数
					ksc = x[k]*sin(san)*cos(tsc_az) + y[k]*sin(san)*sin(tsc_az) + z[k]*cos(san);
	 		       	integ_hh = integ_hh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*sin(tsc_az) - Uxy[k+Num]*cos(tsc_az));	
					//垂直极化的散射系数
					integ_vh = integ_vh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*(df_x[k]*sin(san) - cos(san)*cos(tsc_az))
						 + Uxy[k+Num]*(df_y[k]*sin(san) - cos(san)*sin(tsc_az)));
				}
			}
			//水平极化的散射系数
			integ_hh = integ_hh*(120.0*PI)*i_*k_w/(4.0*PI);
			temp_hh = pow(abs(integ_hh),2.0) / Ph_inc;
		//	temp_hh = pow(abs(integ_hh),2.0) / denom;
		//	temp_hh = pow(abs(integ_hh),2.0);//没有归一化
			sc_p_hh[j] = (temp_hh+i*sc_p_hh[j])/(i+1);
            //垂直极化的散射系数
			integ_vh = integ_vh*(120.0*PI)*i_*k_w/(4.0*PI);
			temp_vh = pow(abs(integ_vh),2.0) / Ph_inc;
		//	temp_vh = pow(abs(integ_vh),2.0) / denom;
		//	temp_vh = pow(abs(integ_vh),2.0);//没有归一化
			sc_p_vh[j] = (temp_vh+i*sc_p_vh[j])/(i+1);
		}
	}

	for(i=0; i<n_sc; i++)
	{
		//outredb_hh<<sc_hl[i]<<" "<<10*log10(sc_p_hh[i])<<endl;
		outre_hh<<sc_hl[i]<<" "<<sc_p_hh[i]<<endl;
		//outredb_vh<<sc_hl[i]<<" "<<10*log10(sc_p_vh[i])<<endl;
		outre_vh<<sc_hl[i]<<" "<<sc_p_vh[i]<<endl;
	}
}
void post_cp::Two_D_Conductor_Vector_Fractal_Scatter_P(int nxy,int nnf,double rrl,double rrd,int Para_Num,double hh,double kkl,double kku,double fra_s,double wwave_len,double in_hl,
													   double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit,int po_h,int po_v)
{
	mean = 0.0;              //均值
	vari = 1.0;              //方差
	seed = 112233.0;         //seed for random number generator 
	Nxy = nxy;
	nf = nnf;
	RL = rrl;
	RD = rrd;
	h = hh;
	kl = kkl;
	ku = kku;
	Fractal_S = fra_s;
	wave_len = wwave_len;

	Inci_HL = in_hl*PI/180.0;//incident angle in degree
	Inci_AZ = in_az*PI/180.0;
	sc_az = ssc_az*PI/180.0;
	N_Realization = nn_r;
	n_sc = nn_sc;
	eps = eeps;
	criterion = crit;
	Po_H = po_h;
	Po_V = po_v;
	k_w = 2.0*PI/wave_len;  //wave number
	g_p = RL/3.0;         //tapering parameter for incident wave

	int i,j,k;
	double denom;
	//矢量锥形波入射功率
	denom = Incidence_Vector_Wave_Pinc_h(k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
	cout<<"第一种矢量锥形波入射功率= "<<denom<<endl;
	//平面波功率
//	denom = cos(Inci_HL)*RL*RL/(2.0*120.0*PI);
//	cout<<denom<<endl;
	double dan = 180.0/n_sc;
	sc_hl.assign(n_sc);
	sc_p_hh.assign(n_sc);
	sc_p_vh.assign(n_sc);
	for(i=0; i<n_sc; i++)
	{
		sc_hl[i]=-90.0+(i+0.5)*dan;
	}
	double san;
	double tsc_az;
	double dxy = RL/Nxy;
	int Num = Nxy*Nxy;
	Complex integ_hh,integ_vh;
	double ksc;
	double temp_hh,temp_vh;
	double Ph_inc;
	for(i=0; i<N_Realization; i++)
	{
		Two_D_Fractal_surface(Nxy,nf,h,mean,vari,RL,kl,ku,Fractal_S);
	    cout<<"The "<<i+1<<"th random rough surface;"<<endl;
		Ph_inc = Two_D_Conductor_Vector_Ph_inc(nxy,k_w,Inci_HL,Inci_AZ,g_p);
		cout<<"第二种矢量锥形波入射功率= "<<Ph_inc<<endl;
		Two_D_Conductor_Vector_SMFIA_CAG_Matrix(Para_Num);
		//Two_D_Conductor_Vector_SMFIA_CAG_Matrix_Integ();
		//Two_D_Conductor_Vector_MoM_Matrix();//点配匹矩量法
		//Two_D_Conductor_Vector_MoM_Matrix_Integ();//矩阵元素采用积分的点配匹矩量法
		
		for(j=0; j<n_sc; j++)
		{
			san = sc_hl[j]*PI/180.0;
			if(san >= 0)
			{
				tsc_az = sc_az;//散射方位角为零的平面上
				integ_hh = 0.0;
				integ_vh = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					//水平极化的散射系数
					ksc = x[k]*sin(san)*cos(tsc_az) + y[k]*sin(san)*sin(tsc_az) + z[k]*cos(san);
	 		       	integ_hh = integ_hh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*sin(tsc_az) - Uxy[k+Num]*cos(tsc_az));	
					//垂直极化的散射系数
					integ_vh = integ_vh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*(df_x[k]*sin(san) - cos(san)*cos(tsc_az))
						 + Uxy[k+Num]*(df_y[k]*sin(san) - cos(san)*sin(tsc_az)));
				}
			}
			else
			{
				san = -san;
				tsc_az = PI + sc_az;
				integ_hh = 0.0;
				integ_vh = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					//水平极化的散射系数
					ksc = x[k]*sin(san)*cos(tsc_az) + y[k]*sin(san)*sin(tsc_az) + z[k]*cos(san);
	 		       	integ_hh = integ_hh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*sin(tsc_az) - Uxy[k+Num]*cos(tsc_az));	
					//垂直极化的散射系数
					integ_vh = integ_vh + dxy*dxy*exp(-i_*k_w*ksc)*(Uxy[k]*(df_x[k]*sin(san) - cos(san)*cos(tsc_az))
						 + Uxy[k+Num]*(df_y[k]*sin(san) - cos(san)*sin(tsc_az)));
				}
			}
			//水平极化的散射系数
			integ_hh = integ_hh*(120.0*PI)*i_*k_w/(4.0*PI);
			temp_hh = pow(abs(integ_hh),2.0) / Ph_inc;
		//	temp_hh = pow(abs(integ_hh),2.0) / denom;
			//temp_hh = pow(abs(integ_hh),2.0);//没有归一化
			sc_p_hh[j] = (temp_hh+i*sc_p_hh[j])/(i+1);
            //垂直极化的散射系数
			integ_vh = integ_vh*(120.0*PI)*i_*k_w/(4.0*PI);
			temp_vh = pow(abs(integ_vh),2.0) / Ph_inc;
		//	temp_vh = pow(abs(integ_vh),2.0) / denom;
			//temp_vh = pow(abs(integ_vh),2.0);//没有归一化
			sc_p_vh[j] = (temp_vh+i*sc_p_vh[j])/(i+1);
		}
	}

	for(i=0; i<n_sc; i++)
	{
		//outredb_hh<<sc_hl[i]<<" "<<10*log10(sc_p_hh[i])<<endl;
		outre_hh<<sc_hl[i]<<" "<<sc_p_hh[i]<<endl;
		//outredb_vh<<sc_hl[i]<<" "<<10*log10(sc_p_vh[i])<<endl;
		outre_vh<<sc_hl[i]<<" "<<sc_p_vh[i]<<endl;
	}
}
void post_cp::Two_D_Conductor_Fractal_Scatter_P(int nxy,int nnf,double rrl,double rrd,int Para_Num,double hh,double kkl,double kku,double fra_s,double wwave_len,double in_hl,double in_az,double ssc_az,int nn_r,int nn_sc,double eeps,double crit)
{
	mean=0.0;              //均值
	vari=1.0;              //方差
	seed=112233.0;         //seed for random number generator 
	Nxy=nxy;
	nf = nnf;
	RL=rrl;
	RD = rrd;
	h=hh;
	kl = kkl;
	ku = kku;
	Fractal_S = fra_s;
	wave_len=wwave_len;
	Inci_HL=in_hl*PI/180.0;//incident angle in degree
	Inci_AZ=in_az*PI/180.0;
	sc_az = ssc_az*PI/180.0;
	N_Realization = nn_r;
	n_sc = nn_sc;
	eps = eeps;
	criterion = crit;
	k_w = 2.0*PI/wave_len;  //wave number
	g_p = RL/3.0;         //tapering parameter for incident wave

	int i,j,k;
	double denom;
	//锥形波入射功率
	denom=1.0+pow(cos(Inci_HL),2)+2.0*pow(tan(Inci_HL),2);
	denom=denom/(2.0*pow(k_w*g_p*cos(Inci_HL),2));
	denom=8.0*pow(PI,3)*pow(g_p,2)*cos(Inci_HL)*(1.0-denom);
	cout<<denom<<endl;
	double dan=180.0/n_sc;
	sc_hl.assign(n_sc);
	sc_p=new double[n_sc];
	for(i=0; i<n_sc; i++)
	{
		sc_hl[i]=-90.0+(i+0.5)*dan;
	}
	double san;
	double tsc_az;
	double dxy=RL/Nxy;
	Complex integ;
	Complex ksc;
	double temp;
	int Num = Nxy * Nxy;
	for(i=0; i<N_Realization; i++)
	{
		Two_D_Fractal_surface(Nxy,nf,h,mean,vari,RL,kl,ku,Fractal_S);
		
		cout<<"The "<<i+1<<"th random rough surface;"<<endl;
		Two_D_Conductor_Scalar_SMFIA_CAG_Matrix(Para_Num);//SMFSIA
		//Two_D_Conductor_Scalar_MoM_Matrix();//点配匹矩量法
				
		for(j=0; j<n_sc; j++)
		{
			san=sc_hl[j]*PI/180.0;
			if(san >= 0)
			{
				tsc_az = sc_az;
				integ = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					ksc = -x[k]*sin(san)*cos(tsc_az)-y[k]*sin(san)*sin(tsc_az)-z[k]*cos(san);
					integ += dxy*dxy*Uxy[k]*exp(i_*k_w*ksc);
				}
			}
			else
			{
				san = -san;//对于三维散射问题，散射角为与z轴的夹角，则在上半空间，散射角范围为0~90度
				tsc_az = PI + sc_az;
				integ = 0.0;
	  	    	for(k=0; k<Num; k++)
				{
					ksc = -x[k]*sin(san)*cos(tsc_az)-y[k]*sin(san)*sin(tsc_az)-z[k]*cos(san);
					integ += dxy*dxy*Uxy[k]*exp(i_*k_w*ksc);
				}
			}
			temp = pow(abs(integ),2);
			temp = temp/denom;//用入射波功率规范化
			//temp = temp;//没有规范化
			sc_p[j] = (temp+i*sc_p[j])/(i+1.0);
		}
	}

	for(i=0; i<n_sc; i++)
	{
		outdb<<sc_hl[i]<<" "<<10*log10(sc_p[i])<<endl;
		out<<sc_hl[i]<<" "<<sc_p[i]<<endl;
	}
	out_Reflectivity<<Reflectivity_Calcu(sc_p,n_sc);
}


void post_cp::One_D_Conductor_Ocean_Scatter(int num,double uss,double kkl,double kku,double rrl,double wwave_len,double ii_an,int nn_r,int nn_sc,int bbw,double eeps)
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
	g_p=RL/4.;         //tapering parameter for incident wave 

	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;
	
	int i,j,k;
	double denom;
	denom=1.0+2.0*pow(tan(i_an),2);
	denom=denom/(2.0*pow(k_w*g_p*cos(i_an),2));
	denom=8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-denom);
	double dan=180.0/n_sc/2.0;
	sc_an=new double[n_sc];
	sc_p=new double[n_sc];
	RCS.assign(n_sc);
	for(i=0;i<n_sc;i++)
	{
		sc_an[i]=-90.0+(2*(i+1)-1)*dan;
	}
	double san;
	double dx=RL/N;
	Complex integ;
	double ksc;
	double temp;
	for(i=0;i<N_Realization;i++)
	{
		cout<<"The "<<i+1<<" random rough surface"<<endl;

		ocean_surface();
		One_D_Conductor_BMIA_CAG_Matrix();
//		One_D_Conductor_MoM_Matrix();
		for(j=0;j<n_sc;j++)
		{
			san=sc_an[j]*PI/180.0;
			integ=0.0;
			for(k=0;k<N;k++)
			{
				ksc=-k_w*(x[k]*sin(san)+z[k]*cos(san));
				integ=integ+dx*rex[k]*exp(i_*ksc);
			}
			//散射系数
			temp = norm(integ)/denom;
			sc_p[j] = (temp+i*sc_p[j])/(i+1);
			//RCS
			temp = norm(integ) / (4.0 * k_w);
			RCS[j] = (temp + i * RCS[j]) / (i + 1);
		}
	}
}
void post_cp::One_D_Conductor_Scatter_P(int num,double hh,double llc,double rrl,double wave_le,double ii_an,int nn_r,int nn_sc,int bbw,double eeps,double crit)
{
	mean=0.0;              //均值
	vari=1.0;              //方差
	N=num;              //total number of sample points 
	seed=123456.0;            //seed for random number generator 
	h=hh;               //rms height，均方根高度
	lc=llc;             //correlation length 
	RL=rrl;             //rough surface length 
	wave_len=wave_le;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=bbw;
	eps = eeps;
	criterion = crit;
	k_w=2.0*PI/wave_len;  //wave number
	//对于一维的情况，g用叶红霞结果不正确
	g_p=RL/4.0;         //tapering parameter for incident wave
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果
	//g_p=6.0/pow(cos(i_an),1.5);
	cout<<"锥形参数g= "<<g_p<<" ;"<<endl;
	
	int i,j,k;
	double denom;
	denom=1.0+2.0*pow(tan(i_an),2);
	denom=denom/(2.0*pow(k_w*g_p*cos(i_an),2));
	denom=8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-denom);
	double dan=180.0/n_sc/2.0;
	sc_an=new double[n_sc];
	sc_p=new double[n_sc];
	RCS.assign(n_sc);
	for(i=0;i<n_sc;i++)
	{
		sc_an[i]=-90.0+(2*(i+1)-1)*dan;
	}
	double san;
	double dx=RL/N;
	Complex integ;
	double ksc;
	double temp;
	for(i=0;i<N_Realization;i++)
	{
		cout<<"The "<<i+1<<" random rough surface"<<endl;

//		ocean_surface();
		gauss_surface();
		One_D_Conductor_BMIA_CAG_Matrix();
//		One_D_Conductor_MoM_Matrix();
		for(j=0;j<n_sc;j++)
		{
			san=sc_an[j]*PI/180.0;
			integ=0.0;
			for(k=0;k<N;k++)
			{
				ksc=-k_w*(x[k]*sin(san)+z[k]*cos(san));
				integ=integ+dx*rex[k]*exp(i_*ksc);
			}
			//散射系数
			temp = norm(integ)/denom;
			sc_p[j] = (temp+i*sc_p[j])/(i+1);
			//RCS
			temp = norm(integ) / (4.0 * k_w);
			RCS[j] = (temp + i * RCS[j]) / (i + 1);
		}
	}

/*	for(i=0; i<n_sc; i++)
	{
		outdb<<sc_an[i]<<" "<<10*log10(sc_p[i])<<endl;
		out<<sc_an[i]<<" "<<sc_p[i]<<endl;
	}*/
}
void post_cp::One_D_Medium_Scatter_P(double gg_p,int num,double hh,double llc,double rrl,double wave_le,
						double ii_an,int nn_r,int nn_sc,int rrd_bw,int nn1,Complex rr_em,Complex permit,double eeps,double crit)
{
	mean=0.0;              //均值
	vari=1.0;              //方差
	N=num;              //total number of sample points 
	seed=123456.0;            //seed for random number generator 
	h=hh;               //rms height，均方根高度
	lc=llc;             //correlation length 
	RL=rrl;             //rough surface length 
	wave_len=wave_le;
	i_an=ii_an*PI/180.0;//incident angle in degree
	N_Realization=nn_r;           //total number of surface realizations 
	n_sc=nn_sc;         //number of scattered angles from -90 deg to 90 deg 
	rd_bw=rrd_bw;
	n1=nn1;
	r_em=rr_em;
	permittivity=permit;
	eps=eeps;
	criterion=crit;
	k_w=2.0*PI/wave_len;  //wave number
	//g_p=RL/4.0;         //tapering parameter for incident wave 
	g_p=gg_p;
	
	fre=3e8/wave_len;
	int i,j,k;
	double denom;
	denom=1.0+2.0*pow(tan(i_an),2.0);
	denom=denom/(2.0*pow(k_w*g_p*cos(i_an),2.0));
	denom=8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-denom);
	double dan=180.0/n_sc/2.0;
	sc_an=new double[n_sc];
	sc_p=new double[n_sc];
	for(i=0;i<n_sc;i++)
	{
		sc_an[i]=-90.0+(2*(i+1)-1)*dan;
	}
	double san;
	double dx=RL/N;
	Complex integ;
	double ksc,fsc;
	double temp;
	
	for(i=0; i < N_Realization; i++)
	{
		cout<<"The "<<i+1<<" random rough surface"<<endl;
		
		gauss_surface();
		//cout<<"gauss_surface"<<endl;
		One_D_Medium_PBTG_BMIA_CAG_Matrix();
		//One_D_Medium_MoM_Matrix();
		for(j=0;j<n_sc;j++)
		{
			san=sc_an[j]*PI/180.0;
			integ=0.0;
			for(k=0;k<N;k++)
			{
				ksc=-k_w*(x[k]*sin(san)+z[k]*cos(san));
				fsc=k_w*(df_x[k]*sin(san)-cos(san));
				integ=integ+dx*(-rex[k]+rex[k+N]*i_*fsc)*exp(i_*ksc);
			}
			
			temp=norm(integ)/denom;
			sc_p[j]=(temp+i*sc_p[j])/(i+1);
		}
	}

}
void post_cp::show_result()
{
	int i;
	double Reflectivity,Emissivity,Bright_T;
	Reflectivity=0.0;
	for(i=0;i<n_sc;i++)
	{
		outdb<<sc_an[i]<<"       "<<10*log10(sc_p[i])<<endl;
		out<<sc_an[i]<<"       "<<sc_p[i]<<endl;
		//outrcsdb<<sc_an[i]<<"       "<<10*log10(RCS[i])<<endl;
		//outrcs<<sc_an[i]<<"       "<<RCS[i]<<endl;
		Reflectivity=Reflectivity+PI/n_sc*sc_p[i];
	}
	Emissivity=1.0-Reflectivity;
	Bright_T=Emissivity*300.0;//300为物理温度
	out_Reflectivity<<"Reflectivity= "<<Reflectivity<<endl;
	out_Reflectivity<<"Emissivity= "<<Emissivity<<endl;
	out_Reflectivity<<"Bright_T= "<<Bright_T<<endl;
	if (fabs(Reflectivity - 1.0) > 0.01)
	{
		out_Reflectivity<<"误差大于0.01，不满足能量守恒，即精度要求不够"<<endl;
	}
	for(i=0;i<N;i++)
	{
		out_current_mag<<x[i]/wave_len<<"  "<<sqrt(norm(rex[i]))<<endl;
	}
}
double post_cp::Bright_Temperature()
{
	double Bright_T;
	double Reflectivity,Absorptivity;
	Reflectivity=0.0;
	for(int i=0;i<n_sc;i++)
	{
		Reflectivity=Reflectivity+PI/n_sc*sc_p[i];
	}
	Absorptivity=1-Reflectivity;
	Bright_T=Absorptivity*300.0;//300为物理温度
	out_Reflectivity<<"Reflectivity= "<<Reflectivity<<endl;
	out_Reflectivity<<"Absorptivity= "<<Absorptivity<<endl;
	out_Reflectivity<<"Bright_T= "<<Bright_T<<endl;
	return Bright_T;
}
void post_cp::Inci_Br_Temp(double gg_p,int num,double hh,double llc,double rrl,double wave_le,
						double i_an_down,double i_an_up,int nn_r,int nn_sc,int rrd_bw,int nn1,Complex rr_em,Complex permit,double eeps,double crit)
{
	double incident,B_T;
	for(incident=i_an_down;incident<=i_an_up;incident+=10.0)
	{
		One_D_Medium_Scatter_P(gg_p,num,hh,llc,rrl,wave_le,incident,nn_r,nn_sc,rrd_bw,nn1,rr_em,permit,eeps,crit);
		B_T=Bright_Temperature(sc_p,nn_sc);
		cout<<incident<<"      "<<B_T<<endl;
		out_inci_bt<<incident<<"      "<<B_T<<endl;
	}
}
double post_cp::Reflectivity_Calcu(double *Sc,int m)
{
	double re;
	re=0.0;
	for(int i=0;i<m;i++)
	{
		re=re+PI/m*Sc[i];
	}
	return re;
}
double post_cp::Bright_Temperature(double *Sc,int m)
{
	double Bright_T;
	double Refle,Absorptivity;
	Refle=Reflectivity_Calcu(Sc,m);
	Absorptivity=1-Refle;
	Bright_T=Absorptivity*300.0;//300K为物理温度
	
	return Bright_T;
}