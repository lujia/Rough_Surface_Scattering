#include "Rough_Surface.h"

//一维良导体粗糙面的阻抗矩阵不是Toeplitz矩阵，但是对称矩阵，因此只存一半元素就可以了(没实现)
ofstream outstep("criterion.txt",ios::trunc);
ofstream outvm("vm.txt",ios::trunc);
ofstream outzm("zm.txt",ios::trunc);
ofstream outJs0("Js0.txt",ios::trunc);
ofstream outJsd("Jsd.txt",ios::trunc);
ofstream outEs0("Es0.txt",ios::trunc);
ofstream out_Current("current.txt",ios::trunc);
ofstream out_current_mag1("current_mag1.txt",ios::trunc);
ofstream out_current_mag2("current_mag2.txt",ios::trunc);
ofstream out_Reflectivity("Reflectivity.txt",ios::trunc);

class matrix: public RoughSurface
{


public:
	matrix();
	//一维
    void One_D_Conductor_BMIA_CAG_Matrix_Js0();
	void One_D_Conductor_BMIA_CAG_Matrix_Es0();
    void One_D_Conductor_BMIA_CAG_Matrix_Jsd(Vec_Com &ETF_I,Vec_Com &Jsd_I);

	void One_D_Conductor_MoM_Matrix();
	void One_D_Conductor_BMIA_CAG_Matrix();
	void One_D_Medium_MoM_Matrix();
	void One_D_Medium_PBTG_BMIA_CAG_Matrix();
	//一维弱相关矩阵的计算
	Complex para_a_te(int Bool,Complex k,double x);
	Complex *computer_ym_te(Complex k_w,double dx,Complex *x,Vec_Dou &fmn,int nn,int bw,int n1);
	Complex para_a_tm(int Bool,Complex k,double x);
	Complex *computer_ym_tm(Complex k_w,double dx,Complex *rx,Vec_Dou &xmn,Vec_Dou &dfx,Vec_Dou &fmn,int nn,int bw,int n1);
	Complex *computer_ym(Complex k_w,double dx,Complex *rx,Vec_Dou &Cx,Vec_Dou &x,Vec_Dou &dfx,Vec_Dou &fmn,int ncg,int ndg,int rd_bw,Complex r_em,int n1);

	//二维标量波入射导体粗糙面
	void Two_D_Conductor_Scalar_SMFIA_CAG_Matrix(int Para_Num);
	void Two_D_Conductor_Scalar_MoM_Matrix();
	
	//二维标量波入射中的积分
	Complex Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(double dxy,double xm,double ym,double xn,double yn,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex kk);//精确计算非奇异矩阵元素
    Complex Two_D_Conductor_Scalar_Integ_Matrix_Element_Zfs(double dxy,double xm,double ym,double xn,double yn,Complex kk);
	Complex Two_D_Conductor_Scalar_Func_G_0_Zs(Complex x[2],double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex k);
    Complex Two_D_Conductor_Scalar_Func_G_0_Zfs(Complex x[2],double xm,double ym,Complex k);
	//格林函数奇异时G0的左端积分
	Complex Two_D_Conductor_Scalar_Integ_Sing_Left(double dxy,double xm,double ym,double dfx,double dfy,double ddfxx,double ddfxy,double ddfyy,Complex kk);
   	Complex Two_D_Conductor_Scalar_Integ_Sing_Func_Left(Complex x[2],double xm,double ym,double dfx,double dfy,double ddfxx,double ddfxy,double ddfyy,Complex k);
	//格林函数奇异时G0的奇异点积分
	double Two_D_Conductor_Scalar_Integ_Sing(double dxy,double dfxm,double dfym);
    double Two_D_Conductor_Scalar_Integ_Sing_Function(double x[1],double dxy,double dfxm,double dfym);
	void Two_D_Conductor_Scalar_Integ_Up_Down(double dxy,int j,double y[2]);
    void Two_D_Conductor_Scalar_Integ_Up_Down(double dxy,int j,double xn,double yn,double y[2]);
	//后来修改的奇异积分
	double Two_D_Conductor_Scalar_Integ_Sing_Function(double x[1],double ym,double dxy,double dfxm,double dfym);
    double Two_D_Conductor_Scalar_Integ_Sing(double ym,double dxy,double dfx,double dfy);
	//二维矢量波入射导体粗糙面
	void Two_D_Conductor_Vector_SMFIA_CAG_Matrix(int Para_Num);
	void Two_D_Conductor_Vector_SMFIA_CAG_Matrix_Integ(int Para_Num);
	void Two_D_Conductor_Vector_MoM_Matrix();
	void Two_D_Conductor_Vector_MoM_Matrix_Integ();
    Complex G_R(double R,Complex k);
	Complex G_R(Complex R,Complex k);
	//二维矢量波入射的奇异性处理
	Complex Two_D_Con_Vec_Wave_Integ_Sing_Left_No(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk);
	Complex Two_D_Con_Vec_Wave_Integ_Sing_Left_No_Function(char a,Complex x[2],double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex kk);
	Complex Two_D_Con_Vec_Wave_Integ_Sing_Right(char xyz,double ym,double dfxm,double dfym,double dxy,Complex kk);
	Complex Two_D_Con_Vec_Wave_Integ_Sing_Right_Function(char xyz,Complex x[1],double ym,double dfxm,double dfym,double dxy,Complex kk);
	void Two_D_Conductor_Vector_Integ_Up_Down(double dxy,int j,double yn,double y[2]);
	Complex Two_D_Con_Vec_Wave_Integ_Sing_Unite_G(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk);
	Complex Two_D_Con_Vec_Wave_Integ_Sing(string LRUD,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk);
	//二维矢量波入射的积分计算矩阵元素
	//Zs
	Complex Two_D_Conductor_Vector_Integ_Function_Zs(string LRUD,Complex x[2],double xm,double ym,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num);
	Complex Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(string LRUD,double dxy,double xm,double xn,double ym,double yn,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num);
	void Two_D_Conductor_Vector_Integ_Up_Down(double dxy,int j,double xn,double yn,double y[2]);
	//Zfs
	Complex Two_D_Conductor_Vector_Integ_Function_Zfs(string LRUD,Complex x[2],double xm,double ym,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num);
	Complex Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(string LRUD,double dxy,double xm,double xn,double ym,double yn,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num);

	//二维矢量锥形波入射
	Complex Incidence_Vector_Wave(char xyz,double xm,double ym,double zm,double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g);
	void Up_Down_Incidence_Vector_Wave(int j,double y[2],Complex xup[2],double rl,int nxy,double wave_num,double in_hl);
	Complex Incidence_Vector_Wave_Integral(char xyz,Complex x[4],double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g);
	void Up_Down_Incidence_Vector_Wave_Test(int j,double y[2],double rl,int nxy,double Value_UD);
	Complex Incidence_Vector_Wave_Test(char xyz,double xm,double ym,double zm,double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g,double Value_UD);
	//徐鹏方法
	ThreeD Hi_xyz(double xi,double yi,double zi,int nxy,double wave_num,double In_HL,double In_AZ,double Value_g);
	ThreeD Ei_xyz(double xi,double yi,double zi,int nxy,double wave_num,double In_HL,double In_AZ,double Value_g);
	double Two_D_Conductor_Vector_Ph_inc(int nxy,double wave_num,double In_HL,double In_AZ,double Value_g);
	//入射波功率
	double Incidence_Vector_Integ_Pinc_h(double xup[4],double wave_num,double In_HL,double In_AZ,double Value_g);
	void Up_Down_Incidence_Vector_Pinc_h(int j,double y[2],double xup[2],double rl,int nxy,double wave_num);
	double Incidence_Vector_Wave_Pinc_h(double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g);
	//锥形波近似
	Complex Incidence_Vector_Wave_Approximation(char xyz,double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g);
	Complex Incidence_Vector_Wave_Approximation_Hi_All_Polarization(char xyz,double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g,int Pola_H,int Pola_V);
    ThreeD ThreeD_hkp(double wave_num,double kx,double ky,double in_hl,double in_az);
	ThreeD ThreeD_vkp(double wave_num,double kx,double ky,double in_hl,double in_az);
	//平面波入射
	//TE极化平面波：分别为磁场、电场
	Complex Incidence_Vector_Plane_Wave_Hi(char xyz,double xm,double ym,double zm,double wave_num,double In_HL,double In_AZ);
	Complex Incidence_Vector_Plane_Wave_Ei(char xyz,double xm,double ym,double zm,double wave_num,double In_HL,double In_AZ);
	
	//二维介质粗糙面的电磁散射仿真
	void Two_D_Dielectric_SMCG_Matrix(int Para_a_Num,int Para_b_Num);
	void Two_D_Dielectric_MoM_Matrix();
    Complex Two_D_Die_Integ_Sing_Unite_G(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk);
    Complex Two_D_Die_Integ_Sing_Unite_g(double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk);
	Complex Diel_G_R(double R,Complex k);
	Complex Diel_g_R(double R,Complex k);

public:
	//g_p为tapering parameter for incident wave 
	double k_w,g_p;
	Complex k_w1,permittivity;
	//w_num1,w_num2分别为媒介1和媒介2中的波数
	double w_num1;
	Complex w_num2;
	//w_Im1,w_Im2分别为媒介1和媒介2中的波阻抗
	double w_Im1;
	Complex w_Im2;
	double i_an;//入射角和波长
	//r_em为一个决定参数：为1时为TE wave;and为 e1/e for TM wave.
	Complex r_em;
	//带状矩阵的半带宽
	int rd_bw;
	double fre;
	//粗格数、细格数及倍数,且n1必须是奇数
	int Ndg,Ncg,n1;
	//强相关距离
	double RD;
	double Inci_HL;
	double Inci_AZ;
	Vec_Com Uxy;
	Vec_Com Diel_Uxy[6];
	//极化方向选择：Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
	int Po_H;
	int Po_V;

	Vec_Com Js0; //粗糙面感应电流；
	Vec_Com Jsdd;
	Vec_Com Ei_S;
//	Vec_Com ETF;
	
};

matrix::matrix()
{
	k_w=0.0;
	g_p=0.0;
	i_an=0.0;
	wave_len=0.0;
	k_w1=0.0;
	permittivity=0.0;
	r_em=0.0;
	rd_bw=0;
	fre=0.0;
	Ndg=0;
	Ncg=0;
	n1=0;
}
//////////////////////////////////////////////////////////////////////////
//一维导体粗糙面
void matrix::One_D_Conductor_BMIA_CAG_Matrix_Js0()
{
	int i,j,u;
	int mbw=2*rd_bw+1;
	Complex *zm,*vm;
	zm=new Complex[N*mbw];
	vm=new Complex[N];
	k_w=2*PI/wave_len;
	double dx=RL/N;
	
	double euler=1.78107;
	double dl;
	double abx;
	double Yita=120*PI;
	double d=10.*wave_len;  //目标距离粗糙面的高度；
//    double d=10.*wave_len; 
    Vec_Com Es0;
		Es0.assign(N);
	Vec_dou XX;
	     XX.assign(N);
	
	     Js0.assign(N);
	Vec_Com zss;
	      zss.assign(N*N);
	Vec_Com vss;
	      vss.assign(N);
		  
		  Ei_S.assign(N);
/*	//采用压缩矩阵存储带状矩阵:强矩阵
	for(i=0;i<N;i++)
	{
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= rd_bw)
			{
				if(abs(i-j) <= rd_bw)
				{
					if(i == j)
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rd_bw)) <= rd_bw && (j+i-rd_bw) < N)
				{
					if(i == (j+i-rd_bw))
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j+i-rd_bw],2)+pow(z[i]-z[j+i-rd_bw],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			//cout<<zm[u]<<"   "<<u<<endl;
		}
	}*/
	//
	for(i=0;i<N;i++)
	{
		
		for(j=0;j<N;j++)
		{
			dl=dx*sqrt(1.0+pow(df_x[j],2));
			u=i*N+j;
			if(i == j)
			{
						
			  zss[u]=(k_w*Yita/4.)*dl*(1.+2.*i_/PI*(log(euler*k_w*dl/4.)-1.));
			}
		    else
			{
		       abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		       zss[u]=(k_w*Yita/4.)*dl*Hankel(0,1,abx);
			}
		}
	}
	
	//v矩阵的填充
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+d*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		vss[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		Ei_S[i] = vss[i];
		
	}

    resolve_FBM(zss,vss,Js0,N);
    
	for(i=0;i<N;i++)
	{
	 XX[i]=x[i]/wave_len;
	 outJs0<<XX[i]<<" "<<real(Js0[i])<<"  "<<imag(Js0[i])<<"  "<<abs(Js0[i])<<endl;
	}
	

/*	
	//迭代法解矩阵方程
	//resolve_m(zm,N,vm);
	//resolve_bm(zm,vm,n,rd_bw);
	resolve_cgm(zm,vm,N,rd_bw);
	//for(i=0;i<n;i++)
	//	outzm<<rex[i]<<endl;
	/////////////////////////////////////////////
	Complex *ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
    //Complex *ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);

	Complex *C,*r,*C0;
	C=new Complex[N];
	r=new Complex[N];
	C0=new Complex[N];
	
	for(i=0;i<N;i++)
	{
		C[i]=vm[i]-ym[i];
		r[i]=-ym[i];
		
	}
	//cout<<Norm2(ym,n)<<endl;
	double Nom=L2_Norm(vm,N);
	double te;
	//te=sqrt(Norm2(r,N)/Nom);
	cout<<sqrt(L2_Norm(r,N)/Nom)<<endl;
	while(sqrt(L2_Norm(r,N)/Nom) > 0.01)
	{
		te=sqrt(L2_Norm(r,N)/Nom);
		cout<<"BMIA迭代标准＝"<<te<<endl;
		for(i=0;i<N;i++)
		{
			C0[i]=C[i];
			//C0[i].display();
		}
		//C0=C;//??
		//resolve_m(zm,N,C);
		//resolve_bm(zm,C,N,rd_bw);
		resolve_cgm(zm,C,N,rd_bw);
		//ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);
		ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
		for(i=0;i<N;i++)
		{
			C[i]=vm[i]-ym[i];
			r[i]=C[i]-C0[i];
			//C[i].display();
		}
		te=sqrt(L2_Norm(r,N)/Nom);
	}
	//////////////////////////////////////////////////
	delete [] zm;
	delete [] vm;
	delete [] ym;


	for(i=0;i<N;i++)
	{
//     Js0[i]=rex[i]*i_/(k_w*Yita)/sqrt(1+pow(df_x[i],2));
	 Js0[i]=rex[i];
	 XX[i]=x[i]/wave_len;
	 outJs0<<XX[i]<<" "<<real(Js0[i])<<"  "<<imag(Js0[i])<<"  "<<abs(Js0[i])<<endl;
	}
*/
}

void matrix::One_D_Conductor_BMIA_CAG_Matrix_Jsd(Vec_Com &ETF_I, Vec_Com &Jsd_I)
{

	
	int i,j,u;
	int mbw=2*rd_bw+1;
	Complex *zm,*vm;
	zm=new Complex[N*mbw];
	vm=new Complex[N];
	k_w=2*PI/wave_len;
	double dx=RL/N;
	
	double euler=1.78107;
	double dl;
	double abx;
	double Yita=120*PI;
	double d=10.*wave_len;  //目标距离粗糙面的高度；

    Vec_Com Es0;
		Es0.assign(N);
	Vec_dou XX;
	     XX.assign(N);
	Vec_Com zss;
	      zss.assign(N*N);
	Vec_Com vss;
	      vss.assign(N);
    Vec_Com Jsd;
	     Jsd.assign(N);


		
for(i=0;i<N;i++)
	{
		
		for(j=0;j<N;j++)
		{
			dl=dx*sqrt(1.0+pow(df_x[j],2));
			u=i*N+j;
			if(i == j)
			{
						
			  zss[u]=(k_w*Yita/4.)*dl*(1.+2.*i_/PI*(log(euler*k_w*dl/4.)-1.));
			}
		    else
			{
		       abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		       zss[u]=(k_w*Yita/4.)*dl*Hankel(0,1,abx);
			}
		}
	}
	
	//v矩阵的填充
    for(i=0;i<N;i++)
	{
		vss[i] = ETF_I[i]; 
	}

    resolve_FBM(zss,vss,Jsd,N);

	for (i=0;i<N;i++)
	{
	  Jsd_I[i]=Jsd[i];
	}

/*    
	for(i=0;i<N;i++)
	{
	 XX[i]=x[i]/wave_len;
	 outJsd<<XX[i]<<" "<<real(Jsd[i])<<"  "<<imag(Jsd[i])<<"  "<<abs(Jsd[i])<<endl;
	}
/*	//采用压缩矩阵存储带状矩阵:强矩阵
	for(i=0;i<N;i++)
	{
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= rd_bw)
			{
				if(abs(i-j) <= rd_bw)
				{
					if(i == j)
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rd_bw)) <= rd_bw && (j+i-rd_bw) < N)
				{
					if(i == (j+i-rd_bw))
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j+i-rd_bw],2)+pow(z[i]-z[j+i-rd_bw],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			//cout<<zm[u]<<"   "<<u<<endl;
		}
	}
*/	//
/*	for(i=0;i<N;i++)
	{
		dl=dx*sqrt(1.0+pow(df_x[i],2));
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= rd_bw)
			{
				
				if(abs(i-j) <= rd_bw)
				{
					if(i == j)
					{
						
			         	zm[u]=(k_w*Yita/4.)*dl*(1.+2.*i_/PI*(log(euler*k_w*dl/4.)-1.));
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        		zm[u]=(k_w*Yita/4.)*dl*Hankel(0,1,abx);
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rd_bw)) <= rd_bw && (j+i-rd_bw) < N)
				{
					if(i == (j+i-rd_bw))
					{
						
			         	zm[u]=(k_w*Yita/4.)*dl*(1.+2.*i_/PI*(log(euler*k_w*dl/4.)-1.));
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j+i-rd_bw],2)+pow(z[i]-z[j+i-rd_bw],2));
		        		zm[u]=(k_w*Yita/4.)*dl*Hankel(0,1,abx);
					}
				}
			}
			//cout<<zm[u]<<"   "<<u<<endl;
		}
	}
	//FFT计算弱矩阵
	
	//v矩阵的填充
//	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		vm[i] = ETF_I[i]; 
	}
	
	//迭代法解矩阵方程
	//resolve_m(zm,N,vm);
	//resolve_bm(zm,vm,n,rd_bw);
	resolve_cgm(zm,vm,N,rd_bw);
	//for(i=0;i<n;i++)
	//	outzm<<rex[i]<<endl;
	/////////////////////////////////////////////
	Complex *ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
    //Complex *ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);

	Complex *C,*r,*C0;
	C=new Complex[N];
	r=new Complex[N];
	C0=new Complex[N];
	
	for(i=0;i<N;i++)
	{
		C[i]=vm[i]-ym[i];
		r[i]=-ym[i];
		
	}
	//cout<<Norm2(ym,n)<<endl;
	double Nom=L2_Norm(vm,N);
	double te;
	//te=sqrt(Norm2(r,N)/Nom);
	cout<<sqrt(L2_Norm(r,N)/Nom)<<endl;
	while(sqrt(L2_Norm(r,N)/Nom) > 0.01)
	{
		te=sqrt(L2_Norm(r,N)/Nom);
		cout<<"BMIA迭代标准＝"<<te<<endl;
		for(i=0;i<N;i++)
		{
			C0[i]=C[i];
			//C0[i].display();
		}
		//C0=C;//??
		//resolve_m(zm,N,C);
		//resolve_bm(zm,C,N,rd_bw);
		resolve_cgm(zm,C,N,rd_bw);
		//ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);
		ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
		for(i=0;i<N;i++)
		{
			C[i]=vm[i]-ym[i];
			r[i]=C[i]-C0[i];
			//C[i].display();
		}
		te=sqrt(L2_Norm(r,N)/Nom);
	}
	//////////////////////////////////////////////////
	delete [] zm;
	delete [] vm;
	delete [] ym;
	for(i=0;i<N;i++)
	{
//    Jsd[i]=rex[i]*i_/(k_w*Yita)/sqrt(1+pow(df_x[i],2));
      Jsd[i]=rex[i];
//	 XX[i]=x[i]/wave_len;
//	 outJsd<<XX[i]<<"   "<<abs(Jsd[i])<<endl;
	}
*/
}



void matrix::One_D_Conductor_BMIA_CAG_Matrix_Es0()
{
int i,j,u;
	int mbw=2*rd_bw+1;
	Complex *zm,*vm;
	zm=new Complex[N*mbw];
	vm=new Complex[N];
	k_w=2*PI/wave_len;
	double dx=RL/N;
	
	double euler=1.78107;
	double dl;
	double abx;
	double Yita=120*PI;
	double d=10.*wave_len;  //目标距离粗糙面的高度；

    Vec_Com Es0;
		Es0.assign(N);
	Vec_dou XX;
	     XX.assign(N);
	
	     Js0.assign(N);
	Vec_Com zss;
	      zss.assign(N*N);
	Vec_Com vss;
	      vss.assign(N);
/*	//采用压缩矩阵存储带状矩阵:强矩阵
	for(i=0;i<N;i++)
	{
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= rd_bw)
			{
				if(abs(i-j) <= rd_bw)
				{
					if(i == j)
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rd_bw)) <= rd_bw && (j+i-rd_bw) < N)
				{
					if(i == (j+i-rd_bw))
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j+i-rd_bw],2)+pow(z[i]-z[j+i-rd_bw],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			//cout<<zm[u]<<"   "<<u<<endl;
		}
	}*/
	//
	for(i=0;i<N;i++)
	{
		dl=dx*sqrt(1.0+pow(df_x[i],2));
		for(j=0;j<N;j++)
		{
			u=i*N+j;
			if(i == j)
			{
						
			  zss[u]=(k_w*Yita/4.)*dl*(1.+2.*i_/PI*(log(euler*k_w*dl/4.)-1.));
			}
		    else
			{
		       abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		       zss[u]=(k_w*Yita/4.)*dl*Hankel(0,1,abx);
			}
		}
	}
	
	//v矩阵的填充
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+d*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		vss[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		
	}

    resolve_FBM(zss,vss,Js0,N);
    
	for(i=0;i<N;i++)
	{
	 XX[i]=x[i]/wave_len;
	 outJs0<<XX[i]<<" "<<real(Js0[i])<<"  "<<imag(Js0[i])<<"  "<<abs(Js0[i])<<endl;
	}

////// 离水平面10倍波长的高度的散射场///////////////////////////
	double abx1;

	for( i=0;i<N;i++)
	{
		for( j=0;j<N;j++)
		{
			dl=dx*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(x[i]-x[j],2)+pow(d-z[j],2));
			Es0[i]=Es0[i]-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Js0[j];
		}
	}
	for(i=0;i<N;i++)
	{
	  outEs0<<XX[i]<<" "<<real(Es0[i])<<"  "<<imag(Es0[i])<<"  "<<abs(Es0[i])<<endl;
	}	
///////////////////////////////////////////////////////////////////

	
}
//////////////////////////////////////////////////////////////////////////
//一维介质粗糙面
//////////////////////////////////////////////////////////////////////////
//二维导体粗糙面
//////////////////////////////////////////////////////////////////////////
//二维介质粗糙面

Complex matrix::Two_D_Die_Integ_Sing_Unite_G(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk)
{//由二维介质粗糙面的积分方程中可以看出：G1,2(R)的奇异性只与二维导体粗糙面中的G(R)相差一个“-”号
	Complex re;
	re = - Two_D_Con_Vec_Wave_Integ_Sing_Unite_G(xyz,xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk);
	return re;
}
Complex matrix::Two_D_Die_Integ_Sing_Unite_g(double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk)
{
	Complex re;
	re = Two_D_Conductor_Scalar_Integ_Sing_Left(dxy,xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,kk)
		+Two_D_Conductor_Scalar_Integ_Sing(ym,dxy,dfxm,dfym);
	return re;
}
void matrix::Two_D_Dielectric_MoM_Matrix()
{
	int i,j,u;
	mean = 0.;
	vari = 1.;
	seed = 112233.;
	fre = 3.e9;
	wave_len = 3.e8/fre;
	Nxy = 8;
	nf=100;
	h=8.1202e-4;         //rms height 
//	h = 0.5*wave_len;
	kl=100;            //lower wavenumber cutoff 
	ku=4000;           //upper wavenumber cutoff 
	Fractal_S=1.5;      
	RL = ( Nxy / 4.0 )*wave_len;
	RD = 0.5*wave_len;//强相关距离
	Inci_HL = 10.*PI/180.0;//incident angle in degree
	Inci_AZ = 0.*PI/180.0;
	g_p = RL/3.0;         //tapering parameter for incident wave
		
	permittivity = 25. + i_*1.;
	w_num1 = 2.*PI/wave_len;
	w_num2 = 2.*PI*fre*sqrt(1.*miu0*ebslun0*permittivity);//区域2的介电常数
	w_Im1 = sqrt(miu0/ebslun0);//区域1中的波阻抗
	w_Im2 = sqrt(1.*miu0/(ebslun0*permittivity));//区域2中的波阻抗

	eps = 1.0e-2;//迭代法的收敛准则
	criterion = 1.0e-2;//SMFSIA的收敛准则

	Two_D_Fractal_surface(Nxy,nf,h,mean,vari,RL,kl,ku,Fractal_S);

	//////////////////////////////////////////////////////////////////////////
	//
	int Num = Nxy*Nxy;//二维表面总的节点个数
	double dxy = RL/Nxy;
	cout<<"单位间隔为 dxy= "<<dxy/wave_len<<"*Wavelength ;"<<endl;
	
	cout<<"强弱相关距离 RD= "<<RD/wave_len<<"*Wavelength ;"<<endl;
	if ( RD > RL)
	{
		cout<<"RD must be smaller than RL;"<<endl;
		exit(0);//程序终止
	}

	//入射波:x,y,z三个分量
	ThreeD Hi,Ei;
	ThreeD the;
	Vec_Com Vmn;
	Vmn.assign(6*Num);
	
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		Hi = Hi_xyz(x[i],y[i],z[i],Nxy,w_num1,Inci_HL,Inci_AZ,g_p);
		Ei = Ei_xyz(x[i],y[i],z[i],Nxy,w_num1,Inci_HL,Inci_AZ,g_p);
		the.Assign(0.,-1.,-df_y[i]);
		Vmn[i] = the.Point_Multiply(the,Hi);
		the.Assign(1.,0.,df_x[i]);
		Vmn[i+Num] = the.Point_Multiply(the,Hi);
		the.Assign(-df_x[i],-df_y[i],1.);
		Vmn[i+2*Num] = the.Point_Multiply(the,Ei);
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Vmn[0][i])+norm(Vmn[1][i])+norm(Vmn[2][i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Ei.TD_x)+norm(Ei.TD_y)+norm(Ei.TD_z))<<endl;
	}
//	cout<<L2_Norm(Vmn)<<endl;

	cout<<"入射波计算完毕; "<<endl;

	//强矩阵Zs
	int RD_num = int(RD/dxy);//整数倍
	cout<<"强弱相关距离相对单位间隔的整数倍 RD_num= "<<RD_num<<endl;
	//////////////////////////////////////////////////////////////////////////
/*  Vec_Com zfs_00;
	double rou;
	zfs_00.assign(Num*Num);
	for ( i=0; i<Num; i++ )
	{
		for ( j=0; j<Num; j++ )
		{
			u = i*Num + j;
			rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
			rou=rou/dxy;
			if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
			{
			}
			else
			{
				rou = rou * dxy;
				zfs_00[u] = dxy*dxy*Diel_G_R(rou,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
			}
			display(zfs_00[u]);
			outresult<<" ";
		}
		outresult<<";"<<endl;
	}*/
/*	Vec_Zmn_Cell Zfs,zfs_tr;
	Vec_Com Zfs_full;
	Zfs_full.assign(6*Num*6*Num);
	Zmn_Cell zsmn_t;
	double R_S,rou;
	Complex temp;
	for ( i=0; i<6*Num; i++ )
	{
		for ( j=0; j<6*Num; j++ )
		{
			u = i*6*Num + j;
			if ( i < Num && j < Num )
			{
				rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i < Num && j >= Num && j<2*Num)
			{
				rou=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(-df_y[i]*(x[i]-x[j-Num])+df_y[j-Num]*(x[i]-x[j-Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i < Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(df_y[i]*(z[i]-z[j-2*Num])+(y[i]-y[j-2*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i < Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(rou,w_num1)*(df_y[i]*df_x[j-3*Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i < Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(rou,w_num1)*(df_y[i]*df_y[j-4*Num]+1.);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i < Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i]-x[j-5*Num],2.0)+pow(y[i]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-5*Num )
					{					
					}
					else
					{	
					}
				}
			}
			
			else if ( i >= Num && i<2*Num && j<Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(df_x[j]*(y[i-Num]-y[j])-df_x[i-Num]*(y[i-Num]-y[j]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= Num && i<2*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(-(z[i-Num]-z[j-Num])+df_y[j-Num]*(y[i-Num]-y[j-Num])+df_x[i-Num]*(x[i-Num]-x[j-Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= Num && i<2*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(-df_x[i-Num]*(z[i-Num]-z[j-2*Num])-(x[i-Num]-x[j-2*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= Num && i<2*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(rou,w_num1)*(-1.-df_x[j-3*Num]-df_x[i-Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= Num && i<2*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(rou,w_num1)*(df_y[j-4*Num]*df_x[i-Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= Num && i<2*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-5*Num],2.0)+pow(y[i-Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}

			else if ( i >= 2*Num && i<3*Num && j<Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num1)*(df_x[i-2*Num]-df_x[j]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num1)*(df_y[i-2*Num]-df_y[j-Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-2*Num],2.0)+pow(y[i-2*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(-df_x[i-2*Num]*df_x[j-3*Num]*(y[i-2*Num]-y[j-3*Num])+df_y[i-2*Num]*(df_x[j-3*Num]*(x[i-2*Num]-x[j-3*Num])-(z[i-2*Num]-z[j-3*Num]))-(y[i-2*Num]-y[j-3*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(df_x[i-2*Num]*((z[i-2*Num]-z[j-4*Num])-df_y[j-4*Num]*(y[i-2*Num]-y[j-4*Num]))+df_y[i-2*Num]*df_y[j-4*Num]*(x[i-2*Num]-x[j-4*Num])+(x[i-2*Num]-x[j-4*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num1)*(df_x[i-2*Num]*(x[i-2*Num]-x[j-5*Num])+df_y[i-2*Num]*(y[i-2*Num]-y[j-5*Num])-(z[i-2*Num]-z[j-5*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}

			else if ( i >= 3*Num && i<4*Num && j<Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num2)*(df_y[i-3*Num]*df_x[j]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num2)*(df_y[i-3*Num]*df_y[j-Num]+1.);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-2*Num],2.0)+pow(y[i-3*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(df_y[i-3*Num]*(y[i-3*Num]-y[j-3*Num])+df_x[j-3*Num]*(x[i-3*Num]-x[j-3*Num])-(z[i-3*Num]-z[j-3*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(-df_y[i-3*Num]*(x[i-3*Num]-x[j-4*Num])+df_y[j-4*Num]*(x[i-3*Num]-x[j-4*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(1./permittivity)*Diel_G_R(rou,w_num2)*(df_y[i-3*Num]*(z[i-3*Num]-z[j-5*Num])+(y[i-3*Num]-y[j-5*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}

			else if ( i >= 4*Num && i<5*Num && j<Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num2)*(-1.-df_x[i-4*Num]*df_x[j]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(rou,w_num2)*(-df_x[i-4*Num]*df_y[j-Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-2*Num],2.0)+pow(y[i-4*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-2*Num )
					{					
					}
					else
					{	
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(df_x[j-3*Num]*(y[i-4*Num]-y[j-3*Num])-df_x[i-4*Num]*(y[i-4*Num]-y[j-3*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(-(z[i-4*Num]-z[j-4*Num])+df_y[j-4*Num]*(y[i-4*Num]-y[j-4*Num])+df_x[i-4*Num]*(x[i-4*Num]-x[j-4*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(1./permittivity)*Diel_G_R(rou,w_num2)*(-df_x[i-4*Num]*(z[i-4*Num]-z[j-5*Num])+(x[i-4*Num]-x[j-5*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}

			else if ( i >= 5*Num && i<6*Num && j<Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(-df_x[i-5*Num]*df_x[j]*(y[i-5*Num]-y[j])+df_y[i-5*Num]*(df_x[j]*(x[i-5*Num]-x[j])-(z[i-5*Num]-z[j]))-(y[i-5*Num]-y[j]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(df_x[i-5*Num]*(z[i-5*Num]-z[j-Num])+df_y[j-Num]*df_y[i-5*Num]*(x[i-5*Num]-x[j-Num])-df_y[j-Num]*df_x[i-5*Num]*(y[i-5*Num]-y[j-Num])+(x[i-5*Num]-x[j-Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*Diel_G_R(rou,w_num2)*(df_x[i-5*Num]*(x[i-5*Num]-x[j-2*Num])+df_y[i-5*Num]*(y[i-5*Num]-y[j-2*Num])-(z[i-5*Num]-z[j-2*Num]));
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(-i_*(w_num1*w_Im1*(1./permittivity)))*Diel_g_R(rou,w_num2)*(df_x[i-5*Num]-df_x[j-3*Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
				}
				else
				{
					rou = rou * dxy;
					temp = dxy*dxy*(-i_*(w_num1*w_Im1*(1./permittivity)))*Diel_g_R(rou,w_num2)*(df_y[i-5*Num]-df_y[j-4*Num]);
					zsmn_t.Zmn_i = i;
					zsmn_t.Zmn_j = j;
					zsmn_t.Zmn_z = temp;
					Zfs.push_back(zsmn_t);
					Zfs_full[u] = temp;
				}
			}
			else// if ( i >= 5*Num && i<6*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-5*Num],2.0)+pow(y[i-5*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}
			//display(Zfs_full[u]);
			//outresult<<" ";
		}
		//outresult<<";"<<endl;
	}

	Vec_Zmn_Cell Zfs_tran;
	for ( i=0; i<Zfs.size(); i++ )
	{
		zsmn_t.Zmn_i = Zfs[i].Zmn_j;
		zsmn_t.Zmn_j = Zfs[i].Zmn_i;
		zsmn_t.Zmn_z = conj(Zfs[i].Zmn_z);
		Zfs_tran.push_back(zsmn_t);
	}*/
/*	Vec_Com Zfs_full_tran;
	Zfs_full_tran.assign(6*Num*6*Num);
	for ( i=0; i<6*Num; i++ )
	{
		for ( j=0; j<6*Num; j++ )
		{
			Zfs_full_tran[i*6*Num+j] = conj(Zfs_full[j*6*Num+i]);
			display(Zfs_full_tran[i*6*Num+j]);
			outresult<<" ";
		}
		outresult<<";"<<endl;
	}
	Vec_Com tt;
	Two_D_Diel_Matrix_Vector_Multiply(' ',Zfs_tran,Vmn,tt);
//	Multiply_Full_Matrix_Vector(' ',Zfs_full_tran,Vmn,tt);
	for ( i=0; i<6*Num; i++ )
	{
		outresult<<tt[i]<<"   "<<i<<endl;
	}

*/
	//////////////////////////////////////////////////////////////////////////
	//Zs
/*	Vec_Zmn_Cell Zsmn;
	Zmn_Cell zsmn_t;
	double R_S,rou;
	Complex temp;
	for ( i=0; i<6*Num; i++ )
	{
		for ( j=0; j<6*Num; j++ )
		{
			if ( i < Num && j < Num )
			{
				rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j )
					{
						temp = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i < Num && j >= Num && j<2*Num)
			{
				rou=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-Num )
					{					
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j-Num])+df_y[j-Num]*(x[i]-x[j-Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i < Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-2*Num )
					{
						temp = df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0)+pow(z[i]-z[j-2*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j-2*Num])+(y[i]-y[j-2*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i < Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-3*Num )
					{
						temp = -i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0)+pow(z[i]-z[j-3*Num],2.0));
						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j-3*Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i < Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-4*Num )
					{
						temp = -i_*(w_num1/w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0)+pow(z[i]-z[j-4*Num],2.0));
						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j-4*Num]+1.);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i < Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i]-x[j-5*Num],2.0)+pow(y[i]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-5*Num )
					{					
					}
					else
					{	
					}
				}
			}
			
			else if ( i >= Num && i<2*Num && j<Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j )
					{					
					}
					else
					{	
						R_S = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i-Num]-y[j])-df_x[i-Num]*(y[i-Num]-y[j]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-Num )
					{
						temp = 0.5 + df_y[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) + df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i-Num]-z[j-Num])+df_y[j-Num]*(y[i-Num]-y[j-Num])+df_x[i-Num]*(x[i-Num]-x[j-Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-2*Num )
					{
						temp = -df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0)+pow(z[i-Num]-z[j-2*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-Num]*(z[i-Num]-z[j-2*Num])-(x[i-Num]-x[j-2*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-3*Num )
					{
						temp = -i_*(w_num1/w_Im1)*(-1.-df_x[i-Num]-df_x[i-Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0)+pow(z[i-Num]-z[j-3*Num],2.0));
						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j-3*Num]-df_x[i-Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-4*Num )
					{
						temp = i_*(w_num1/w_Im1)*df_y[i-Num]*df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0)+pow(z[i-Num]-z[j-4*Num],2.0));
						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j-4*Num]*df_x[i-Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-5*Num],2.0)+pow(y[i-Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}

			else if ( i >= 2*Num && i<3*Num && j<Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0)+pow(z[i-2*Num]-z[j],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i-2*Num]-df_x[j]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0)+pow(z[i-2*Num]-z[j-Num],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i-2*Num]-df_y[j-Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-2*Num],2.0)+pow(y[i-2*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-3*Num )
					{
						temp = (-df_x[i-2*Num]*df_x[i-2*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0)+pow(z[i-2*Num]-z[j-3*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-2*Num]*df_x[j-3*Num]*(y[i-2*Num]-y[j-3*Num])+df_y[i-2*Num]*(df_x[j-3*Num]*(x[i-2*Num]-x[j-3*Num])-(z[i-2*Num]-z[j-3*Num]))-(y[i-2*Num]-y[j-3*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-4*Num )
					{
						temp = df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_x[i-2*Num]*df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + (df_y[i-2*Num]*df_y[i-2*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{	
						R_S = sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0)+pow(z[i-2*Num]-z[j-4*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*((z[i-2*Num]-z[j-4*Num])-df_y[j-4*Num]*(y[i-2*Num]-y[j-4*Num]))+df_y[i-2*Num]*df_y[j-4*Num]*(x[i-2*Num]-x[j-4*Num])+(x[i-2*Num]-x[j-4*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-5*Num )
					{
						temp = 0.5 + df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0)+pow(z[i-2*Num]-z[j-5*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*(x[i-2*Num]-x[j-5*Num])+df_y[i-2*Num]*(y[i-2*Num]-y[j-5*Num])-(z[i-2*Num]-z[j-5*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}

			else if ( i >= 3*Num && i<4*Num && j<Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j )
					{
						temp = i_*(w_num1*w_Im1)*df_y[i-3*Num]*df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0)+pow(z[i-3*Num]-z[j],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_x[j]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-Num )
					{
						temp = i_*(w_num1*w_Im1)*(df_y[i-3*Num]*df_y[i-3*Num]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0)+pow(z[i-3*Num]-z[j-Num],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_y[j-Num]+1.);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-2*Num],2.0)+pow(y[i-3*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-3*Num )
					{
						temp = -0.5 + df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0)+pow(z[i-3*Num]-z[j-3*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(y[i-3*Num]-y[j-3*Num])+df_x[j-3*Num]*(x[i-3*Num]-x[j-3*Num])-(z[i-3*Num]-z[j-3*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-4*Num )
					{			
					}
					else
					{			
						R_S = sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0)+pow(z[i-3*Num]-z[j-4*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i-3*Num]*(x[i-3*Num]-x[j-4*Num])+df_y[j-4*Num]*(x[i-3*Num]-x[j-4*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-5*Num )
					{
						temp = (1./permittivity)*df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{			
						R_S = sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0)+pow(z[i-3*Num]-z[j-5*Num],2.0));
						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(z[i-3*Num]-z[j-5*Num])+(y[i-3*Num]-y[j-5*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}

			else if ( i >= 4*Num && i<5*Num && j<Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j )
					{
						temp = i_*(w_num1*w_Im1)*(-1.-df_x[i-4*Num]*df_x[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0)+pow(z[i-4*Num]-z[j],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i-4*Num]*df_x[j]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-Num )
					{
						temp = i_*(w_num1*w_Im1)*(-df_x[i-4*Num]*df_y[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0)+pow(z[i-4*Num]-z[j-Num],2.0));
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i-4*Num]*df_y[j-Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-2*Num],2.0)+pow(y[i-4*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-2*Num )
					{					
					}
					else
					{	
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-3*Num )
					{					
					}
					else
					{	
						R_S = sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0)+pow(z[i-4*Num]-z[j-3*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j-3*Num]*(y[i-4*Num]-y[j-3*Num])-df_x[i-4*Num]*(y[i-4*Num]-y[j-3*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-4*Num )
					{
						temp = -0.5 - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_y[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{	
						R_S = sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0)+pow(z[i-4*Num]-z[j-4*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i-4*Num]-z[j-4*Num])+df_y[j-4*Num]*(y[i-4*Num]-y[j-4*Num])+df_x[i-4*Num]*(x[i-4*Num]-x[j-4*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-5*Num )
					{
						temp = -(1./permittivity)*df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0)+pow(z[i-4*Num]-z[j-5*Num],2.0));
						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i-4*Num]*(z[i-4*Num]-z[j-5*Num])+(x[i-4*Num]-x[j-5*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}

			else if ( i >= 5*Num && i<6*Num && j<Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j )
					{
						temp = (-df_x[i-5*Num]*df_x[i-5*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_x[i-5*Num])*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0)+pow(z[i-5*Num]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i-5*Num]*df_x[j]*(y[i-5*Num]-y[j])+df_y[i-5*Num]*(df_x[j]*(x[i-5*Num]-x[j])-(z[i-5*Num]-z[j]))-(y[i-5*Num]-y[j]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-Num )
					{
						temp = df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_y[i-5*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{	
						R_S = sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0)+pow(z[i-5*Num]-z[j-Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(z[i-5*Num]-z[j-Num])+df_y[j-Num]*df_y[i-5*Num]*(x[i-5*Num]-x[j-Num])-df_y[j-Num]*df_x[i-5*Num]*(y[i-5*Num]-y[j-Num])+(x[i-5*Num]-x[j-Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-2*Num )
					{
						temp = -0.5 + df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
					else
					{	
						R_S = sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0)+pow(z[i-5*Num]-z[j-2*Num],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(x[i-5*Num]-x[j-2*Num])+df_y[i-5*Num]*(y[i-5*Num]-y[j-2*Num])-(z[i-5*Num]-z[j-2*Num]));
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-3*Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0)+pow(z[i-5*Num]-z[j-3*Num],2.0));
						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i-5*Num]-df_x[j-3*Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-4*Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0)+pow(z[i-5*Num]-z[j-4*Num],2.0));
						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i-5*Num]-df_y[j-4*Num]);
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zsmn.push_back(zsmn_t);
					}
				}
			}
			else// if ( i >= 5*Num && i<6*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-5*Num],2.0)+pow(y[i-5*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}
		}
	}
	
	Vec_Zmn_Cell Zsmn_tran;
	for ( i=0; i<Zsmn.size(); i++ )
	{
		zsmn_t.Zmn_i = Zsmn[i].Zmn_j;
		zsmn_t.Zmn_j = Zsmn[i].Zmn_i;
		zsmn_t.Zmn_z = conj(Zsmn[i].Zmn_z);
		Zsmn_tran.push_back(zsmn_t);
	}

	Vec_Com t0;
//	Two_D_Diel_Matrix_Vector_Multiply(' ',Zsmn_tran,Vmn,t0);
	Two_D_Diel_Matrix_Vector_Multiply('*',Zsmn,Vmn,t0);
	for ( i=0; i<6*Num; i++ )
	{
		outresult<<t0[i]<<"   "<<i<<endl;
	}
*/

    //////////////////////////////////////////////////////////////////////////
/*	Vec_Com Zsmn;
	Zsmn.assign(6*Num*6*Num);
	double rou,R_S;
	for ( i=0; i<6*Num; i++ )
	{
		for ( j=0; j<6*Num; j++ )
		{
			u = i*6*Num + j;
			if ( i < Num && j < Num )
			{
				rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j )
					{
						Zsmn[u] = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
					}
				}
			}
			else if ( i < Num && j >= Num && j<2*Num)
			{
				rou=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-Num )
					{					
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j-Num])+df_y[j-Num]*(x[i]-x[j-Num]));
					}
				}
			}
			else if ( i < Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-2*Num )
					{
						Zsmn[u] = df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0)+pow(z[i]-z[j-2*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j-2*Num])+(y[i]-y[j-2*Num]));
					}
				}
			}
			else if ( i < Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-3*Num )
					{
						Zsmn[u] = -i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
					}
					else
					{
						R_S = sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0)+pow(z[i]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j-3*Num]);
					}
				}
			}
			else if ( i < Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-4*Num )
					{
						Zsmn[u] = -i_*(w_num1/w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0)+pow(z[i]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j-4*Num]+1.);
					}
				}
			}
			else if ( i < Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i]-x[j-5*Num],2.0)+pow(y[i]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i == j-5*Num )
					{					
					}
					else
					{	
					}
				}
			}
			
			else if ( i >= Num && i<2*Num && j<Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j )
					{					
					}
					else
					{	
						R_S = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i-Num]-y[j])-df_x[i-Num]*(y[i-Num]-y[j]));
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-Num )
					{
						Zsmn[u] = 0.5 + df_y[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) + df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i-Num]-z[j-Num])+df_y[j-Num]*(y[i-Num]-y[j-Num])+df_x[i-Num]*(x[i-Num]-x[j-Num]));
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-2*Num )
					{
						Zsmn[u] = -df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0)+pow(z[i-Num]-z[j-2*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-Num]*(z[i-Num]-z[j-2*Num])-(x[i-Num]-x[j-2*Num]));
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-3*Num )
					{
						Zsmn[u] = -i_*(w_num1/w_Im1)*(-1.-df_x[i-Num]-df_x[i-Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0)+pow(z[i-Num]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j-3*Num]-df_x[i-Num]);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-4*Num )
					{
						Zsmn[u] = i_*(w_num1/w_Im1)*df_y[i-Num]*df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0)+pow(z[i-Num]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j-4*Num]*df_x[i-Num]);
					}
				}
			}
			else if ( i >= Num && i<2*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-Num]-x[j-5*Num],2.0)+pow(y[i-Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}

			else if ( i >= 2*Num && i<3*Num && j<Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0)+pow(z[i-2*Num]-z[j],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i-2*Num]-df_x[j]);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0)+pow(z[i-2*Num]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i-2*Num]-df_y[j-Num]);
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-2*Num],2.0)+pow(y[i-2*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-3*Num )
					{
						Zsmn[u] = (-df_x[i-2*Num]*df_x[i-2*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0)+pow(z[i-2*Num]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-2*Num]*df_x[j-3*Num]*(y[i-2*Num]-y[j-3*Num])+df_y[i-2*Num]*(df_x[j-3*Num]*(x[i-2*Num]-x[j-3*Num])-(z[i-2*Num]-z[j-3*Num]))-(y[i-2*Num]-y[j-3*Num]));
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-4*Num )
					{
						Zsmn[u] = df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_x[i-2*Num]*df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + (df_y[i-2*Num]*df_y[i-2*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
					}
					else
					{	
						R_S = sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0)+pow(z[i-2*Num]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*((z[i-2*Num]-z[j-4*Num])-df_y[j-4*Num]*(y[i-2*Num]-y[j-4*Num]))+df_y[i-2*Num]*df_y[j-4*Num]*(x[i-2*Num]-x[j-4*Num])+(x[i-2*Num]-x[j-4*Num]));
					}
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-2*Num == j-5*Num )
					{
						Zsmn[u] = 0.5 + df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
					}
					else
					{		
						R_S = sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0)+pow(z[i-2*Num]-z[j-5*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*(x[i-2*Num]-x[j-5*Num])+df_y[i-2*Num]*(y[i-2*Num]-y[j-5*Num])-(z[i-2*Num]-z[j-5*Num]));
					}
				}
			}

			else if ( i >= 3*Num && i<4*Num && j<Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j )
					{
						Zsmn[u] = i_*(w_num1*w_Im1)*df_y[i-3*Num]*df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0)+pow(z[i-3*Num]-z[j],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_x[j]);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-Num )
					{
						Zsmn[u] = i_*(w_num1*w_Im1)*(df_y[i-3*Num]*df_y[i-3*Num]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0)+pow(z[i-3*Num]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_y[j-Num]+1.);
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-2*Num],2.0)+pow(y[i-3*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-2*Num )
					{					
					}
					else
					{		
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-3*Num )
					{
						Zsmn[u] = -0.5 + df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0)+pow(z[i-3*Num]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(y[i-3*Num]-y[j-3*Num])+df_x[j-3*Num]*(x[i-3*Num]-x[j-3*Num])-(z[i-3*Num]-z[j-3*Num]));
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-4*Num )
					{			
					}
					else
					{			
						R_S = sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0)+pow(z[i-3*Num]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i-3*Num]*(x[i-3*Num]-x[j-4*Num])+df_y[j-4*Num]*(x[i-3*Num]-x[j-4*Num]));
					}
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-3*Num == j-5*Num )
					{
						Zsmn[u] = (1./permittivity)*df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
					}
					else
					{			
						R_S = sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0)+pow(z[i-3*Num]-z[j-5*Num],2.0));
						Zsmn[u] = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(z[i-3*Num]-z[j-5*Num])+(y[i-3*Num]-y[j-5*Num]));
					}
				}
			}

			else if ( i >= 4*Num && i<5*Num && j<Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j )
					{
						Zsmn[u] = i_*(w_num1*w_Im1)*(-1.-df_x[i-4*Num]*df_x[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0)+pow(z[i-4*Num]-z[j],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i-4*Num]*df_x[j]);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-Num )
					{
						Zsmn[u] = i_*(w_num1*w_Im1)*(-df_x[i-4*Num]*df_y[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0)+pow(z[i-4*Num]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i-4*Num]*df_y[j-Num]);
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-2*Num],2.0)+pow(y[i-4*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-2*Num )
					{					
					}
					else
					{	
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-3*Num )
					{					
					}
					else
					{	
						R_S = sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0)+pow(z[i-4*Num]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j-3*Num]*(y[i-4*Num]-y[j-3*Num])-df_x[i-4*Num]*(y[i-4*Num]-y[j-3*Num]));
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-4*Num )
					{
						Zsmn[u] = -0.5 - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_y[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
					}
					else
					{	
						R_S = sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0)+pow(z[i-4*Num]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i-4*Num]-z[j-4*Num])+df_y[j-4*Num]*(y[i-4*Num]-y[j-4*Num])+df_x[i-4*Num]*(x[i-4*Num]-x[j-4*Num]));
					}
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-4*Num == j-5*Num )
					{
						Zsmn[u] = -(1./permittivity)*df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0)+pow(z[i-4*Num]-z[j-5*Num],2.0));
						Zsmn[u] = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i-4*Num]*(z[i-4*Num]-z[j-5*Num])+(x[i-4*Num]-x[j-5*Num]));
					}
				}
			}

			else if ( i >= 5*Num && i<6*Num && j<Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j )
					{
						Zsmn[u] = (-df_x[i-5*Num]*df_x[i-5*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_x[i-5*Num])*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0)+pow(z[i-5*Num]-z[j],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i-5*Num]*df_x[j]*(y[i-5*Num]-y[j])+df_y[i-5*Num]*(df_x[j]*(x[i-5*Num]-x[j])-(z[i-5*Num]-z[j]))-(y[i-5*Num]-y[j]));
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=Num && j<2*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-Num )
					{
						Zsmn[u] = df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_y[i-5*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
					}
					else
					{	
						R_S = sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0)+pow(z[i-5*Num]-z[j-Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(z[i-5*Num]-z[j-Num])+df_y[j-Num]*df_y[i-5*Num]*(x[i-5*Num]-x[j-Num])-df_y[j-Num]*df_x[i-5*Num]*(y[i-5*Num]-y[j-Num])+(x[i-5*Num]-x[j-Num]));
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=2*Num && j<3*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-2*Num )
					{
						Zsmn[u] = -0.5 + df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
					}
					else
					{	
						R_S = sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0)+pow(z[i-5*Num]-z[j-2*Num],2.0));
						Zsmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(x[i-5*Num]-x[j-2*Num])+df_y[i-5*Num]*(y[i-5*Num]-y[j-2*Num])-(z[i-5*Num]-z[j-2*Num]));
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=3*Num && j<4*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-3*Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0)+pow(z[i-5*Num]-z[j-3*Num],2.0));
						Zsmn[u] = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i-5*Num]-df_x[j-3*Num]);
					}
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=4*Num && j<5*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-4*Num )
					{					
					}
					else
					{		
						R_S = sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0)+pow(z[i-5*Num]-z[j-4*Num],2.0));
						Zsmn[u] = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i-5*Num]-df_y[j-4*Num]);
					}
				}
			}
			else// if ( i >= 5*Num && i<6*Num && j>=5*Num && j<6*Num )
			{
				rou=sqrt(pow(x[i-5*Num]-x[j-5*Num],2.0)+pow(y[i-5*Num]-y[j-5*Num],2.0));
				rou=rou/dxy;
				if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
				{
					if ( i-5*Num == j-5*Num )
					{					
					}
					else
					{		
					}
				}
			}
			display(Zsmn[u]);
			outresult<<" ";
		}
		outresult<<";"<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
    double R_S;
	Vec_Com Zmn;
	Zmn.assign(6*Num*6*Num);
	for ( i=0; i<6*Num; i++ )
	{
		for ( j=0; j<6*Num; j++ )
		{
			u = i*6*Num + j;
			if ( i < Num && j < Num )
			{
				if ( i == j )
				{
					Zmn[u] = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
				}
				else
				{
					R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
				}
			}
			else if ( i < Num && j >= Num && j<2*Num)
			{
				if ( i == j-Num )
				{
				}
				else
				{
					R_S = sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j-Num])+df_y[j-Num]*(x[i]-x[j-Num]));
				}
			}
			else if ( i < Num && j>=2*Num && j<3*Num )
			{
				if ( i == j-2*Num )
				{
					Zmn[u] = df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
				}
				else
				{
					R_S = sqrt(pow(x[i]-x[j-2*Num],2.0)+pow(y[i]-y[j-2*Num],2.0)+pow(z[i]-z[j-2*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j-2*Num])+(y[i]-y[j-2*Num]));
				}
			}
			else if ( i < Num && j>=3*Num && j<4*Num )
			{
				if ( i == j-3*Num )
				{
					Zmn[u] = -i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
				}
				else
				{
					R_S = sqrt(pow(x[i]-x[j-3*Num],2.0)+pow(y[i]-y[j-3*Num],2.0)+pow(z[i]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j-3*Num]);
				}
			}
			else if ( i < Num && j>=4*Num && j<5*Num )
			{
				if ( i == j-4*Num )
				{
					Zmn[u] = -i_*(w_num1/w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i]-x[j-4*Num],2.0)+pow(y[i]-y[j-4*Num],2.0)+pow(z[i]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j-4*Num]+1.);
				}
			}
			else if ( i < Num && j>=5*Num && j<6*Num )
			{				
			}
			
			else if ( i >= Num && i<2*Num && j<Num )
			{
				if ( i-Num == j )
				{					
				}
				else
				{	
					R_S = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i-Num]-y[j])-df_x[i-Num]*(y[i-Num]-y[j]));
				}
			}
			else if ( i >= Num && i<2*Num && j>=Num && j<2*Num )
			{
				if ( i-Num == j-Num )
				{
					Zmn[u] = 0.5 + df_y[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) + df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i-Num]-z[j-Num])+df_y[j-Num]*(y[i-Num]-y[j-Num])+df_x[i-Num]*(x[i-Num]-x[j-Num]));
				}
			}
			else if ( i >= Num && i<2*Num && j>=2*Num && j<3*Num )
			{
				if ( i-Num == j-2*Num )
				{
					Zmn[u] = -df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('X',x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-Num]-x[j-2*Num],2.0)+pow(y[i-Num]-y[j-2*Num],2.0)+pow(z[i-Num]-z[j-2*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-Num]*(z[i-Num]-z[j-2*Num])-(x[i-Num]-x[j-2*Num]));
				}
			}
			else if ( i >= Num && i<2*Num && j>=3*Num && j<4*Num )
			{
				if ( i-Num == j-3*Num )
				{
					Zmn[u] = -i_*(w_num1/w_Im1)*(-1.-df_x[i-Num]-df_x[i-Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-Num]-x[j-3*Num],2.0)+pow(y[i-Num]-y[j-3*Num],2.0)+pow(z[i-Num]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j-3*Num]-df_x[i-Num]);
				}
			}
			else if ( i >= Num && i<2*Num && j>=4*Num && j<5*Num )
			{
				if ( i-Num == j-4*Num )
				{
					Zmn[u] = i_*(w_num1/w_Im1)*df_y[i-Num]*df_x[i-Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-Num],y[i-Num],df_x[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-Num]-x[j-4*Num],2.0)+pow(y[i-Num]-y[j-4*Num],2.0)+pow(z[i-Num]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j-4*Num]*df_x[i-Num]);
				}
			}
			else if ( i >= Num && i<2*Num && j>=5*Num && j<6*Num )
			{				
			}

			else if ( i >= 2*Num && i<3*Num && j<Num )
			{
				if ( i-2*Num == j )
				{					
				}
				else
				{		
					R_S = sqrt(pow(x[i-2*Num]-x[j],2.0)+pow(y[i-2*Num]-y[j],2.0)+pow(z[i-2*Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i-2*Num]-df_x[j]);
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=Num && j<2*Num )
			{
				if ( i-2*Num == j-Num )
				{					
				}
				else
				{	
					R_S = sqrt(pow(x[i-2*Num]-x[j-Num],2.0)+pow(y[i-2*Num]-y[j-Num],2.0)+pow(z[i-2*Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i-2*Num]-df_y[j-Num]);
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=2*Num && j<3*Num )
			{				
			}
			else if ( i >= 2*Num && i<3*Num && j>=3*Num && j<4*Num )
			{
				if ( i-2*Num == j-3*Num )
				{
					Zmn[u] = (-df_x[i-2*Num]*df_x[i-2*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-2*Num]-x[j-3*Num],2.0)+pow(y[i-2*Num]-y[j-3*Num],2.0)+pow(z[i-2*Num]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i-2*Num]*df_x[j-3*Num]*(y[i-2*Num]-y[j-3*Num])+df_y[i-2*Num]*(df_x[j-3*Num]*(x[i-2*Num]-x[j-3*Num])-(z[i-2*Num]-z[j-3*Num]))-(y[i-2*Num]-y[j-3*Num]));
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=4*Num && j<5*Num )
			{
				if ( i-2*Num == j-4*Num )
				{
					Zmn[u] = df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - df_x[i-2*Num]*df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + (df_y[i-2*Num]*df_y[i-2*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
				}
				else
				{	
					R_S = sqrt(pow(x[i-2*Num]-x[j-4*Num],2.0)+pow(y[i-2*Num]-y[j-4*Num],2.0)+pow(z[i-2*Num]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*((z[i-2*Num]-z[j-4*Num])-df_y[j-4*Num]*(y[i-2*Num]-y[j-4*Num]))+df_y[i-2*Num]*df_y[j-4*Num]*(x[i-2*Num]-x[j-4*Num])+(x[i-2*Num]-x[j-4*Num]));
				}
			}
			else if ( i >= 2*Num && i<3*Num && j>=5*Num && j<6*Num )
			{
				if ( i-2*Num == j-5*Num )
				{
					Zmn[u] = 0.5 + df_x[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) + df_y[i-2*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-2*Num],y[i-2*Num],df_x[i-2*Num],df_y[i-2*Num],ddf_x[i-2*Num],ddf_xy[i-2*Num],ddf_y[i-2*Num],dxy,w_num1);
				}
				else
				{		
					R_S = sqrt(pow(x[i-2*Num]-x[j-5*Num],2.0)+pow(y[i-2*Num]-y[j-5*Num],2.0)+pow(z[i-2*Num]-z[j-5*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i-2*Num]*(x[i-2*Num]-x[j-5*Num])+df_y[i-2*Num]*(y[i-2*Num]-y[j-5*Num])-(z[i-2*Num]-z[j-5*Num]));
				}
			}

			else if ( i >= 3*Num && i<4*Num && j<Num )
			{
				if ( i-3*Num == j )
				{
					Zmn[u] = i_*(w_num1*w_Im1)*df_y[i-3*Num]*df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-3*Num]-x[j],2.0)+pow(y[i-3*Num]-y[j],2.0)+pow(z[i-3*Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_x[j]);
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=Num && j<2*Num )
			{
				if ( i-3*Num == j-Num )
				{
					Zmn[u] = i_*(w_num1*w_Im1)*(df_y[i-3*Num]*df_y[i-3*Num]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-3*Num]-x[j-Num],2.0)+pow(y[i-3*Num]-y[j-Num],2.0)+pow(z[i-3*Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i-3*Num]*df_y[j-Num]+1.);
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=2*Num && j<3*Num )
			{				
			}
			else if ( i >= 3*Num && i<4*Num && j>=3*Num && j<4*Num )
			{
				if ( i-3*Num == j-3*Num )
				{
					Zmn[u] = -0.5 + df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + df_x[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-3*Num]-x[j-3*Num],2.0)+pow(y[i-3*Num]-y[j-3*Num],2.0)+pow(z[i-3*Num]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(y[i-3*Num]-y[j-3*Num])+df_x[j-3*Num]*(x[i-3*Num]-x[j-3*Num])-(z[i-3*Num]-z[j-3*Num]));
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=4*Num && j<5*Num )
			{
				if ( i-3*Num == j-4*Num )
				{			
				}
				else
				{			
					R_S = sqrt(pow(x[i-3*Num]-x[j-4*Num],2.0)+pow(y[i-3*Num]-y[j-4*Num],2.0)+pow(z[i-3*Num]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i-3*Num]*(x[i-3*Num]-x[j-4*Num])+df_y[j-4*Num]*(x[i-3*Num]-x[j-4*Num]));
				}
			}
			else if ( i >= 3*Num && i<4*Num && j>=5*Num && j<6*Num )
			{
				if ( i-3*Num == j-5*Num )
				{
					Zmn[u] = (1./permittivity)*df_y[i-3*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-3*Num],y[i-3*Num],df_x[i-3*Num],df_y[i-3*Num],ddf_x[i-3*Num],ddf_xy[i-3*Num],ddf_y[i-3*Num],dxy,w_num2);
				}
				else
				{			
					R_S = sqrt(pow(x[i-3*Num]-x[j-5*Num],2.0)+pow(y[i-3*Num]-y[j-5*Num],2.0)+pow(z[i-3*Num]-z[j-5*Num],2.0));
					Zmn[u] = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i-3*Num]*(z[i-3*Num]-z[j-5*Num])+(y[i-3*Num]-y[j-5*Num]));
				}
			}

			else if ( i >= 4*Num && i<5*Num && j<Num )
			{
				if ( i-4*Num == j )
				{
					Zmn[u] = i_*(w_num1*w_Im1)*(-1.-df_x[i-4*Num]*df_x[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-4*Num]-x[j],2.0)+pow(y[i-4*Num]-y[j],2.0)+pow(z[i-4*Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i-4*Num]*df_x[j]);
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=Num && j<2*Num )
			{
				if ( i-4*Num == j-Num )
				{
					Zmn[u] = i_*(w_num1*w_Im1)*(-df_x[i-4*Num]*df_y[i-4*Num])*Two_D_Die_Integ_Sing_Unite_g(x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-4*Num]-x[j-Num],2.0)+pow(y[i-4*Num]-y[j-Num],2.0)+pow(z[i-4*Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i-4*Num]*df_y[j-Num]);
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=2*Num && j<3*Num )
			{				
			}
			else if ( i >= 4*Num && i<5*Num && j>=3*Num && j<4*Num )
			{
				if ( i-4*Num == j-3*Num )
				{					
				}
				else
				{	
					R_S = sqrt(pow(x[i-4*Num]-x[j-3*Num],2.0)+pow(y[i-4*Num]-y[j-3*Num],2.0)+pow(z[i-4*Num]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j-3*Num]*(y[i-4*Num]-y[j-3*Num])-df_x[i-4*Num]*(y[i-4*Num]-y[j-3*Num]));
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=4*Num && j<5*Num )
			{
				if ( i-4*Num == j-4*Num )
				{
					Zmn[u] = -0.5 - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_y[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
				}
				else
				{	
					R_S = sqrt(pow(x[i-4*Num]-x[j-4*Num],2.0)+pow(y[i-4*Num]-y[j-4*Num],2.0)+pow(z[i-4*Num]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i-4*Num]-z[j-4*Num])+df_y[j-4*Num]*(y[i-4*Num]-y[j-4*Num])+df_x[i-4*Num]*(x[i-4*Num]-x[j-4*Num]));
				}
			}
			else if ( i >= 4*Num && i<5*Num && j>=5*Num && j<6*Num )
			{
				if ( i-4*Num == j-5*Num )
				{
					Zmn[u] = -(1./permittivity)*df_x[i-4*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-4*Num],y[i-4*Num],df_x[i-4*Num],df_y[i-4*Num],ddf_x[i-4*Num],ddf_xy[i-4*Num],ddf_y[i-4*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-4*Num]-x[j-5*Num],2.0)+pow(y[i-4*Num]-y[j-5*Num],2.0)+pow(z[i-4*Num]-z[j-5*Num],2.0));
					Zmn[u] = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i-4*Num]*(z[i-4*Num]-z[j-5*Num])+(x[i-4*Num]-x[j-5*Num]));
				}
			}

			else if ( i >= 5*Num && i<6*Num && j<Num )
			{
				if ( i-5*Num == j )
				{
					Zmn[u] = (-df_x[i-5*Num]*df_x[i-5*Num]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_x[i-5*Num])*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
				}
				else
				{		
					R_S = sqrt(pow(x[i-5*Num]-x[j],2.0)+pow(y[i-5*Num]-y[j],2.0)+pow(z[i-5*Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i-5*Num]*df_x[j]*(y[i-5*Num]-y[j])+df_y[i-5*Num]*(df_x[j]*(x[i-5*Num]-x[j])-(z[i-5*Num]-z[j]))-(y[i-5*Num]-y[j]));
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=Num && j<2*Num )
			{
				if ( i-5*Num == j-Num )
				{
					Zmn[u] = df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + (df_y[i-5*Num]*df_y[i-5*Num]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - df_y[i-5*Num]*df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
				}
				else
				{	
					R_S = sqrt(pow(x[i-5*Num]-x[j-Num],2.0)+pow(y[i-5*Num]-y[j-Num],2.0)+pow(z[i-5*Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(z[i-5*Num]-z[j-Num])+df_y[j-Num]*df_y[i-5*Num]*(x[i-5*Num]-x[j-Num])-df_y[j-Num]*df_x[i-5*Num]*(y[i-5*Num]-y[j-Num])+(x[i-5*Num]-x[j-Num]));
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=2*Num && j<3*Num )
			{
				if ( i-5*Num == j-2*Num )
				{
					Zmn[u] = -0.5 + df_x[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('X',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) + df_y[i-5*Num]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i-5*Num],y[i-5*Num],df_x[i-5*Num],df_y[i-5*Num],ddf_x[i-5*Num],ddf_xy[i-5*Num],ddf_y[i-5*Num],dxy,w_num2);
				}
				else
				{	
					R_S = sqrt(pow(x[i-5*Num]-x[j-2*Num],2.0)+pow(y[i-5*Num]-y[j-2*Num],2.0)+pow(z[i-5*Num]-z[j-2*Num],2.0));
					Zmn[u] = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i-5*Num]*(x[i-5*Num]-x[j-2*Num])+df_y[i-5*Num]*(y[i-5*Num]-y[j-2*Num])-(z[i-5*Num]-z[j-2*Num]));
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=3*Num && j<4*Num )
			{
				if ( i-5*Num == j-3*Num )
				{					
				}
				else
				{		
					R_S = sqrt(pow(x[i-5*Num]-x[j-3*Num],2.0)+pow(y[i-5*Num]-y[j-3*Num],2.0)+pow(z[i-5*Num]-z[j-3*Num],2.0));
					Zmn[u] = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i-5*Num]-df_x[j-3*Num]);
				}
			}
			else if ( i >= 5*Num && i<6*Num && j>=4*Num && j<5*Num )
			{
				if ( i-5*Num == j-4*Num )
				{					
				}
				else
				{		
					R_S = sqrt(pow(x[i-5*Num]-x[j-4*Num],2.0)+pow(y[i-5*Num]-y[j-4*Num],2.0)+pow(z[i-5*Num]-z[j-4*Num],2.0));
					Zmn[u] = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i-5*Num]-df_y[j-4*Num]);
				}
			}
			else// if ( i >= 5*Num && i<6*Num && j>=5*Num && j<6*Num )
			{				
			}
		}
	}

	cout<<"矩阵填充完毕！"<<endl;
	Vec_Com re;
	resolve_m(Zmn,Vmn,re);
	for ( i=0; i<6*Num; i++ )
	{
		out_Current<<re[i]<<"  "<<i<<endl;
	}

}
void matrix::Two_D_Dielectric_SMCG_Matrix(int Para_a_Num,int Para_b_Num)
{
	int i,j,k;	
	mean = 0.;
	vari = 1.;
	seed = 112233.;
	fre = 3.e9;
	wave_len = 3.e8/fre;
	Nxy = 8;
	nf=100;
	h=8.1202e-4;         //rms height 
//	h = 0.5*wave_len;
	kl=100;            //lower wavenumber cutoff 
	ku=4000;           //upper wavenumber cutoff 
	Fractal_S=1.5;      
	RL = ( Nxy / 4.0 )*wave_len;
	RD = 0.5*wave_len;//强相关距离
	Inci_HL = 10.*PI/180.0;//incident angle in degree
	Inci_AZ = 0.*PI/180.0;
	g_p = RL/3.0;         //tapering parameter for incident wave
		
	permittivity = 25. + i_*1.;
	w_num1 = 2.*PI/wave_len;
	w_num2 = 2.*PI*fre*sqrt(1.*miu0*ebslun0*permittivity);//区域2的介电常数
	w_Im1 = sqrt(miu0/ebslun0);//区域1中的波阻抗
	w_Im2 = sqrt(1.*miu0/(ebslun0*permittivity));//区域2中的波阻抗

	eps = 1.0e-2;//迭代法的收敛准则
	criterion = 1.0e-2;//SMFSIA的收敛准则

	Two_D_Fractal_surface(Nxy,nf,h,mean,vari,RL,kl,ku,Fractal_S);

	//////////////////////////////////////////////////////////////////////////
	//
	int Num = Nxy*Nxy;//二维表面总的节点个数
	double dxy = RL/Nxy;
	cout<<"单位间隔为 dxy= "<<dxy/wave_len<<"*Wavelength ;"<<endl;
	
	cout<<"强弱相关距离 RD= "<<RD/wave_len<<"*Wavelength ;"<<endl;
	if ( RD > RL)
	{
		cout<<"RD must be smaller than RL;"<<endl;
		exit(0);//程序终止
	}

/*	Vec_Com tv;
	tv.assign(Num);
	for ( i=0; i<Num; i++ )
	{
		tv[i] = i+0.5;
	}
//	tv.clear();
//	tv.erase(tv.begin(),tv.end());
	cout<<tv.size()<<endl;
	Vec_Com::iterator ii;
	for ( ii=tv.begin(); ii<tv.end(); ii++ )
	{
		cout<<*ii<<"   "<<ii<<endl;
	}*/

	//入射波:x,y,z三个分量
	ThreeD Hi,Ei;
	ThreeD the;
	Vec_Com Vmn[6];
	for ( i=0; i<6; i++ )
	{
		Vmn[i].assign(Num);
	}
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		Hi = Hi_xyz(x[i],y[i],z[i],Nxy,w_num1,Inci_HL,Inci_AZ,g_p);
		Ei = Ei_xyz(x[i],y[i],z[i],Nxy,w_num1,Inci_HL,Inci_AZ,g_p);
		the.Assign(0.,-1.,-df_y[i]);
		Vmn[0][i] = the.Point_Multiply(the,Hi);
		the.Assign(1.,0.,df_x[i]);
		Vmn[1][i] = the.Point_Multiply(the,Hi);
		the.Assign(-df_x[i],-df_y[i],1.);
		Vmn[2][i] = the.Point_Multiply(the,Ei);
//		cout<<Vmn[2][i]<<"   "<<i<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Vmn[0][i])+norm(Vmn[1][i])+norm(Vmn[2][i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Ei.TD_x)+norm(Ei.TD_y)+norm(Ei.TD_z))<<endl;
	}
//	cout<<L2_Norm(Vmn,6)<<endl;

	cout<<"入射波计算完毕; "<<endl;

	//强矩阵Zs
	int RD_num = int(RD/dxy);//整数倍
	cout<<"强弱相关距离相对单位间隔的整数倍 RD_num= "<<RD_num<<endl;
	//////////////////////////////////////////////////////////////////////////
	//
	Vec_Zmn_Cell Zs[6][6];
	Vec_Com Zs_Diag[6];
	int Band_Ser,Band_Num;//带状矩阵的序号和子带矩阵的数目
	Zmn_Cell zsmn_t;
	Complex temp;
	int Zs_bw;
	double R_S;
	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";
	
	//////////////////////////////////////////////////////////////////////////
	//第零带
	Band_Ser = 0;
	for( Band_Num = 0; Band_Num < Nxy; Band_Num++ )
	{
		Zs_bw = RD_num;
		for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
		{
			if (i-Band_Num*Nxy <= Zs_bw)
			{
				for ( j=Band_Num*Nxy; j <= Band_Num*Nxy+Zs_bw+(i-Band_Num*Nxy); j++ )
				{
					if ( i == j )
					{
						//第一个积分方程
						temp = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);
						
						Zs_Diag[0].push_back(temp);

//						test_zmn_cell(Zs[0][0][0],zsmn_t);
//						cout<<Zs[0][0][0].Zmn_i<<" "<<Zs[0][0][0].Zmn_j<<" "<<Zs[0][0][0].Zmn_z<<endl;
//						cout<<zsmn_t.Zmn_i<<" "<<zsmn_t.Zmn_j<<" "<<zsmn_t.Zmn_z<<endl;

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[0][1].push_back(zsmn_t);

						temp = df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][2].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][3].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][4].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[0][5].push_back(zsmn_t);

						//第二个积分方程
						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[1][0].push_back(zsmn_t);

						temp = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][1].push_back(zsmn_t);

						Zs_Diag[1].push_back(temp);

						temp = -df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][2].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*(-1.-df_x[i]-df_x[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][3].push_back(zsmn_t);

						temp = i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][4].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[1][5].push_back(zsmn_t);

						//第三个积分方程
						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][0].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][1].push_back(zsmn_t);

						temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][2].push_back(zsmn_t);

						Zs_Diag[2].push_back(temp);

						temp = (-df_x[i]*df_x[i]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][3].push_back(zsmn_t);

						temp = df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - df_x[i]*df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + (df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][4].push_back(zsmn_t);

						temp = 0.5 + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = i_*(w_num1*w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][0].push_back(zsmn_t);

						temp = i_*(w_num1*w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][1].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[3][2].push_back(zsmn_t);
						
						temp = -0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][3].push_back(zsmn_t);

						Zs_Diag[3].push_back(temp);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[3][4].push_back(zsmn_t);

						temp = (1./permittivity)*df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = i_*(w_num1*w_Im1)*(-1.-df_x[i]*df_x[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][0].push_back(zsmn_t);

						temp = i_*(w_num1*w_Im1)*(-df_x[i]*df_y[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][1].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[4][2].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[4][3].push_back(zsmn_t);

						temp = -0.5 - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][4].push_back(zsmn_t);

						Zs_Diag[4].push_back(temp);

						temp = -(1./permittivity)*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = (-df_x[i]*df_x[i]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (df_y[i]*df_x[i])*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][0].push_back(zsmn_t);

						temp = df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][1].push_back(zsmn_t);

						temp = -0.5 + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][2].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][3].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][4].push_back(zsmn_t);

						temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][5].push_back(zsmn_t);

						Zs_Diag[5].push_back(temp);
					}
					else
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
			}
			else
			{
				for ( j=Band_Num*Nxy+(i-Band_Num*Nxy)-Zs_bw; j < Band_Num*Nxy+((i-Band_Num*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Nxy; j++ )
				{
					if ( i == j )
					{
						//第一个积分方程
						temp = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);
						
						Zs_Diag[0].push_back(temp);
						
						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[0][1].push_back(zsmn_t);

						temp = df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][2].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][3].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[0][4].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[0][5].push_back(zsmn_t);

						//第二个积分方程
						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[1][0].push_back(zsmn_t);

						temp = 0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][1].push_back(zsmn_t);

						Zs_Diag[1].push_back(temp);

						temp = -df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][2].push_back(zsmn_t);

						temp = -i_*(w_num1/w_Im1)*(-1.-df_x[i]-df_x[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][3].push_back(zsmn_t);

						temp = i_*(w_num1/w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[1][4].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[1][5].push_back(zsmn_t);

						//第三个积分方程
						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][0].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][1].push_back(zsmn_t);

						temp = 0.0;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[2][2].push_back(zsmn_t);

						Zs_Diag[2].push_back(temp);

						temp = (-df_x[i]*df_x[i]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][3].push_back(zsmn_t);

						temp = df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - df_x[i]*df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + (df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][4].push_back(zsmn_t);

						temp = 0.5 + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num1);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = i_*(w_num1*w_Im1)*df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][0].push_back(zsmn_t);

						temp = i_*(w_num1*w_Im1)*(df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][1].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[3][2].push_back(zsmn_t);
						
						temp = -0.5 + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][3].push_back(zsmn_t);

						Zs_Diag[3].push_back(temp);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[3][4].push_back(zsmn_t);

						temp = (1./permittivity)*df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = i_*(w_num1*w_Im1)*(-1.-df_x[i]*df_x[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][0].push_back(zsmn_t);

						temp = i_*(w_num1*w_Im1)*(-df_x[i]*df_y[i])*Two_D_Die_Integ_Sing_Unite_g(x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][1].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[4][2].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[4][3].push_back(zsmn_t);

						temp = -0.5 - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][4].push_back(zsmn_t);

						Zs_Diag[4].push_back(temp);

						temp = -(1./permittivity)*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (1./permittivity)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = (-df_x[i]*df_x[i]-1.)*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (df_y[i]*df_x[i])*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][0].push_back(zsmn_t);

						temp = df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + (df_y[i]*df_y[i]+1.)*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - df_y[i]*df_x[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][1].push_back(zsmn_t);

						temp = -0.5 + df_x[i]*Two_D_Die_Integ_Sing_Unite_G('X',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) + df_y[i]*Two_D_Die_Integ_Sing_Unite_G('Y',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2) - Two_D_Die_Integ_Sing_Unite_G('Z',x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,w_num2);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs[5][2].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][3].push_back(zsmn_t);

						//为零不需要存储
						//temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][4].push_back(zsmn_t);

						temp = 0.;
						//zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						//Zs[5][5].push_back(zsmn_t);

						Zs_Diag[5].push_back(temp);
					}
					else
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
			}
		}
	}
	//从左到右的带
	for( Band_Ser = 1; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
			{
				if (i-Band_Num*Nxy <= Zs_bw)
				{
					for ( j=Band_Num*Nxy+Band_Ser*Nxy; j <= Band_Num*Nxy+Band_Ser*Nxy+Zs_bw+(i-Band_Num*Nxy); j++ )
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
				else
				{
					for ( j=Band_Num*Nxy+Band_Ser*Nxy+(i-Band_Num*Nxy)-Zs_bw; j < Band_Num*Nxy+Band_Ser*Nxy+((i-Band_Num*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Band_Ser*Nxy+Nxy; j++ )
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
			}				
		}
	}

	//当Band_Ser=RD_num时；矩阵为一条线
	for( i=0; i < Num-RD_num*Nxy; i++ )
	{
		//第一个积分方程
		j = i+RD_num*Nxy;
		R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_i = i;
		zsmn_t.Zmn_j = j;
		zsmn_t.Zmn_z = temp;
		Zs[0][0].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][1].push_back(zsmn_t);
						
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][2].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][3].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][4].push_back(zsmn_t);

		//第二个积分方程
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][0].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][2].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][3].push_back(zsmn_t);

		temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][4].push_back(zsmn_t);
						
		//第三个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][3].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][4].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][5].push_back(zsmn_t);

		//第四个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][3].push_back(zsmn_t);

		//为零
		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][4].push_back(zsmn_t);

		temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][5].push_back(zsmn_t);

		//第五个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][3].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][4].push_back(zsmn_t);

		temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][5].push_back(zsmn_t);

		//第六个积分方程
		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][0].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][2].push_back(zsmn_t);
		
		temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][3].push_back(zsmn_t);
		
		temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][4].push_back(zsmn_t);		
	}
	//从上到下的带
	for( Band_Ser = 1; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			for ( i=Band_Num*Nxy+Band_Ser*Nxy; i < Nxy+Band_Num*Nxy+Band_Ser*Nxy; i++ )
			{
				if (i-Band_Num*Nxy-Band_Ser*Nxy <= Zs_bw)
				{
					for ( j=Band_Num*Nxy; j <= Band_Num*Nxy+Zs_bw+(i-Band_Num*Nxy-Band_Ser*Nxy); j++ )
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
				else
				{
					for ( j=Band_Num*Nxy+(i-Band_Num*Nxy-Band_Ser*Nxy)-Zs_bw; j < Band_Num*Nxy+((i-Band_Num*Nxy-Band_Ser*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Nxy; j++ )
					{
						//第一个积分方程
						R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs[0][0].push_back(zsmn_t);	
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][1].push_back(zsmn_t);
						
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[0][4].push_back(zsmn_t);

						//第二个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][3].push_back(zsmn_t);

						temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[1][4].push_back(zsmn_t);
						
						//第三个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][4].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[2][5].push_back(zsmn_t);

						//第四个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][3].push_back(zsmn_t);

						//为零
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[3][5].push_back(zsmn_t);

						//第五个积分方程
						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][0].push_back(zsmn_t);

						temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][3].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][4].push_back(zsmn_t);

						temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[4][5].push_back(zsmn_t);

						//第六个积分方程
						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][0].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][1].push_back(zsmn_t);

						temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][2].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][3].push_back(zsmn_t);

						temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
						//cout<<temp<<endl;
						zsmn_t.Zmn_z = temp;
						Zs[5][4].push_back(zsmn_t);
					}
				}
			}				
		}
	}

	//当Band_Ser=RD_num时；矩阵为一条线
	for( j=0; j < Num-RD_num*Nxy; j++ )
	{
		//第一个积分方程
		i = j+RD_num*Nxy;
		R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_i = i;
		zsmn_t.Zmn_j = j;
		zsmn_t.Zmn_z = temp;
		Zs[0][0].push_back(zsmn_t);	
						
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][1].push_back(zsmn_t);
						
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][2].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][3].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[i]*df_y[j]+1.);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[0][4].push_back(zsmn_t);

		//第二个积分方程
		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][0].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*(z[i]-z[j])-(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][2].push_back(zsmn_t);

		temp = dxy*dxy*(-i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(-1.-df_x[j]-df_x[i]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][3].push_back(zsmn_t);

		temp = dxy*dxy*(i_*(w_num1/w_Im1))*Diel_g_R(R_S,w_num1)*(df_y[j]*df_x[i]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[1][4].push_back(zsmn_t);
						
		//第三个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_x[i]-df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num1)*(df_y[i]-df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][3].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*((z[i]-z[j])-df_y[j]*(y[i]-y[j]))+df_y[i]*df_y[j]*(x[i]-x[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][4].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num1)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[2][5].push_back(zsmn_t);

		//第四个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(df_y[i]*df_y[j]+1.);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_y[i]*(y[i]-y[j])+df_x[j]*(x[i]-x[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][3].push_back(zsmn_t);

		//为零
		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_y[i]*(x[i]-x[j])+df_y[j]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][4].push_back(zsmn_t);

		temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(df_y[i]*(z[i]-z[j])+(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[3][5].push_back(zsmn_t);

		//第五个积分方程
		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-1.-df_x[i]*df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][0].push_back(zsmn_t);

		temp = dxy*dxy*(i_*w_num1*w_Im1)*Diel_g_R(R_S,w_num2)*(-df_x[i]*df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[j]*(y[i]-y[j])-df_x[i]*(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][3].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-(z[i]-z[j])+df_y[j]*(y[i]-y[j])+df_x[i]*(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][4].push_back(zsmn_t);

		temp = dxy*dxy*(1./permittivity)*Diel_G_R(R_S,w_num2)*(-df_x[i]*(z[i]-z[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[4][5].push_back(zsmn_t);

		//第六个积分方程
		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(-df_x[i]*df_x[j]*(y[i]-y[j])+df_y[i]*(df_x[j]*(x[i]-x[j])-(z[i]-z[j]))-(y[i]-y[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][0].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(z[i]-z[j])+df_y[j]*df_y[i]*(x[i]-x[j])-df_y[j]*df_x[i]*(y[i]-y[j])+(x[i]-x[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][1].push_back(zsmn_t);

		temp = dxy*dxy*Diel_G_R(R_S,w_num2)*(df_x[i]*(x[i]-x[j])+df_y[i]*(y[i]-y[j])-(z[i]-z[j]));
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][2].push_back(zsmn_t);
		
		temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_x[i]-df_x[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][3].push_back(zsmn_t);
		
		temp = dxy*dxy*(-i_*(w_num1*w_Im1*permittivity))*Diel_g_R(R_S,w_num2)*(df_y[i]-df_y[j]);
		//cout<<temp<<endl;
		zsmn_t.Zmn_z = temp;
		Zs[5][4].push_back(zsmn_t);		
	}
	//近场强相关矩阵填充完毕

    //////////////////////////////////////////////////////////////////////////
	//验证平面矩阵向量积
/*	Vec_Com t6[6];
	Vec_Zmn_Cell Zs_tran[6][6];
	for ( i=0; i<6; i++ )
	{
		for ( j=0; j<6; j++ )
		{
			for ( k=0; k<Zs[i][j].size(); k++ )
			{
				zsmn_t.Zmn_i = Zs[i][j][k].Zmn_j;
				zsmn_t.Zmn_j = Zs[i][j][k].Zmn_i;
				zsmn_t.Zmn_z = conj(Zs[i][j][k].Zmn_z);
				Zs_tran[j][i].push_back(zsmn_t);
			}
		}
	}
//	Two_D_Diel_Zs_Vector_Multiply('*',Zs,Vmn,t6);
	//	Two_D_Diel_Multi_Zw_Vec(0,0,x,y,z,df_x,df_y,Vmn,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,t6);
	Two_D_Diel_Multi_Zfs_Vec('*',x,y,z,df_x,df_y,Vmn,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,t6);
	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			outresult<<t6[k][i]<<"  "<<k<<"  "<<i<<endl;
		}
	}*/

    cout<<"近场强相关矩阵填充完毕；开始解矩阵方程；"<<endl;
//	cout<<Zs[0][0][0].Zmn_i<<"  "<<Zs[0][0][0].Zmn_j<<"  "<<Zs[0][0][0].Zmn_z<<endl;
	//迭代法解矩阵方程
//	Two_D_Diel_GPBI_CG(Zs,Vmn,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
//	Two_D_Diel_CGM(Zs,Vmn,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
//  Two_D_Diel_CGLS(Zs,Vmn,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
	Two_D_Diel_Bi_CGSTAB_Pre(Zs,Zs_Diag,Vmn,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);

	//右端弱相关矩阵向量积
	Vec_Com YM[6];
	Two_D_Diel_Multi_Zw_Vec(Para_a_Num,Para_b_Num,x,y,z,df_x,df_y,Diel_Uxy,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,YM);
/*  for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			outresult<<YM[k][i]<<"  "<<k<<"  "<<i<<endl;
		}
	}*/
	Vec_Com C[6],C0[6],r[6];
	for ( i=0; i<6; i++ )
	{
		C[i].assign(Num);
		C0[i].assign(Num);
		r[i].assign(Num);
	}
	
	for ( k=0; k<6; k++ )
	{
		for (i=0; i<Num; i++)
		{
			C[k][i] = Vmn[k][i] - YM[k][i];
			r[k][i] = -YM[k][i];
			//cout<<YM[i]<<endl;
		}
	}
	
	
	double Nom=L2_Norm(Vmn,6);
	double Crit=sqrt(L2_Norm(r,6)/Nom);
	cout<<Nom<<endl;
	while(Crit > criterion)
	{
		cout<<"criterion = "<<Crit<<endl;
		outstep<<Crit<<endl;
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				C0[k][i] = C[k][i];
			}
		}
	
		//解矩阵方程
//		Two_D_Diel_GPBI_CG(Zs,C0,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
//		Two_D_Diel_CGM(Zs,C0,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
//		Two_D_Diel_CGLS(Zs,C0,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);
		Two_D_Diel_Bi_CGSTAB_Pre(Zs,Zs_Diag,C0,x,y,z,df_x,df_y,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,eps,Diel_Uxy);

		//右端弱相关矩阵向量积
		Two_D_Diel_Multi_Zw_Vec(Para_a_Num,Para_b_Num,x,y,z,df_x,df_y,Diel_Uxy,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,YM);

//		Two_D_Cond_Vector_Bi_CGSTAB_Pre(Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zs_D_LU,Zs_D_RD,C0_U,C0_D,x,y,z,df_x,df_y,k_w,dxy,Nxy,RD_num,eps,Uxy);		
//		Two_D_Diel_Multi_Zw_Vec(Para_a_Num,Para_b_Num,x,y,z,df_x,df_y,Diel_Uxy,w_num1,w_num2,w_Im1,permittivity,dxy,Nxy,RD_num,YM);

		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				C[k][i] = Vmn[k][i] - YM[k][i];
				r[k][i] = C[k][i]-C0[k][i];
			}
		}
		Crit=sqrt(L2_Norm(r,6)/Nom);
	}

	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			out_Current<<Diel_Uxy[k][i]<<"  "<<k<<"  "<<i<<endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//
	for ( i=0; i<6; i++ )
	{
		Vmn[i].clear();
		Zs_Diag[i].clear();
		C[i].clear();
		C0[i].clear();
		r[i].clear();
		YM[i].clear();
		for ( j=0; j<6; j++ )
		{
			Zs[i][j].clear();
		}
	}
	
}
//////////////////////////////////////////////////////////////////////////
Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing(string LRUD,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk)
{
	Complex re;

	if (LRUD == "UL")
	{
		//re = G_R(rs,wave_num)*(-dfym*(ym-x[1])-(xm-x[0])*dfxn+zmn);
		re = 0.5 -dfym*Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('Y',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk)
			-dfxm*Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('X',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk)
			+Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('Z',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk);
	}
	else if (LRUD == "UR")
	{
		//re = G_R(rs,wave_num)*(dfym*(xm-x[0])-(xm-x[0])*dfyn);
		re = 0.0;
	}
	else if (LRUD == "DL")
	{
		//re = G_R(rs,wave_num)*(dfxm*(ym-x[1])-(ym-x[1])*dfxn);
		re = 0.0;
	}
	else
	{
		//re = G_R(rs,wave_num)*(-dfxm*(xm-x[0])-(ym-x[1])*dfyn+zmn);
		re = 0.5 -dfxm*Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('X',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk)
			-dfym*Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('Y',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk)
			+Two_D_Con_Vec_Wave_Integ_Sing_Unite_G('Z',xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk);
	}

	return re;
}
Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing_Unite_G(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk)
{
	Complex re;
    re = Two_D_Con_Vec_Wave_Integ_Sing_Left_No(xyz,xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,dxy,kk)
		+ Two_D_Con_Vec_Wave_Integ_Sing_Right(xyz,ym,dfxm,dfym,dxy,kk);
	return re;
}
void matrix::Two_D_Conductor_Vector_Integ_Up_Down(double dxy,int j,double ym,double y[2])
{
	switch(j)
	{
	case 0:{y[0]=ym-dxy/2.0;y[1]=ym+dxy/2.0;break;}
    default:{}
	}
}


Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing_Right(char xyz,double ym,double dfxm,double dfym,double dxy,Complex kk)
{
	int n=1;
    int js[1]={4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex *x,s,*a,*b;
	Complex p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Vector_Integ_Up_Down(dxy,j-1,ym,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Con_Vec_Wave_Integ_Sing_Right_Function(xyz,x,ym,dfxm,dfym,dxy,kk);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}

Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing_Right_Function(char xyz,Complex x[1],double ym,double dfxm,double dfym,double dxy,Complex kk)
{
	Complex re;
	double A;
	Complex B,C;
	A = 1. + pow(dfxm,2.);
	B = 2. * dfxm * dfym * (ym - x[0]);
	C = (1. + pow(dfym,2.)) * pow(ym-x[0],2.);
	Complex a2 = C/A - pow(B,2.)/4.*pow(A,2.);
	double U11,D12;
	Complex U21,D22;
	U11 = -dxy/2.;
	D12 = dxy/2.;
	U21 = -dxy/2. + B/(2.*A);
	D22 = dxy/2. + B/(2.*A);
	Complex IR1,IR3;
	switch(xyz)
	{
	case 'X':
		{
			IR1 = (sqrt(A*pow(U11,2.)+B*U11+C)/A-B/(2.*pow(A,1.5))*log(abs(2.*A*U11+B+2.*sqrt(A)*sqrt(A*pow(U11,2.)+B*U11+C))))
				-(sqrt(A*pow(D12,2.)+B*D12+C)/A-B/(2.*pow(A,1.5))*log(abs(2.*A*D12+B+2.*sqrt(A)*sqrt(A*pow(D12,2.)+B*D12+C))));
			IR3 = 1./pow(A,1.5) * ((-1./sqrt(pow(U21,2.)+a2)-B/(2.*A)*(U21/(a2*sqrt(pow(U21,2.)+a2)))) - (-1./sqrt(pow(D22,2.)+a2)-B/(2.*A)*(D22/(a2*sqrt(pow(D22,2.)+a2)))));
			re = -IR3/(4.*PI) - IR1*pow(kk,2.)/(8.*PI);
			break;
		}
	case 'Y':
		{
			IR1 = (log(abs(2.*A*U11+B+2.*sqrt(A)*sqrt(A*pow(U11,2.)+B*U11+C)))/sqrt(A)) - (log(abs(2.*A*D12+B+2.*sqrt(A)*sqrt(A*pow(D12,2.)+B*D12+C)))/sqrt(A));
			IR3 = 1./pow(A,1.5) * ((U21/(a2*sqrt(pow(U21,2.)+a2))) - (D22/(a2*sqrt(pow(D22,2.)+a2))));
			re = (ym-x[0]) * (-IR3/(4.*PI) - IR1*pow(kk,2.)/(8.*PI));
			break;
		}
	case 'Z':
		{
			IR1 = -dfxm * ((sqrt(A*pow(U11,2.)+B*U11+C)/A-B/(2.*pow(A,1.5))*log(abs(2.*A*U11+B+2.*sqrt(A)*sqrt(A*pow(U11,2.)+B*U11+C))))
				-(sqrt(A*pow(D12,2.)+B*D12+C)/A-B/(2.*pow(A,1.5))*log(abs(2.*A*D12+B+2.*sqrt(A)*sqrt(A*pow(D12,2.)+B*D12+C)))))
				- dfym*(ym-x[0]) * ((log(abs(2.*A*U11+B+2.*sqrt(A)*sqrt(A*pow(U11,2.)+B*U11+C)))/sqrt(A)) 
				- (log(abs(2.*A*D12+B+2.*sqrt(A)*sqrt(A*pow(D12,2.)+B*D12+C)))/sqrt(A)));
			IR3 = -dfxm/pow(A,1.5) * ((-1./sqrt(pow(U21,2.)+a2)) - (-1./sqrt(pow(D22,2.)+a2))) 
				+ ((B/(2.*A))*dfxm - dfym*(ym-x[0]))/pow(A,1.5) * ((U21/(a2*sqrt(pow(U21,2.)+a2))) - (D22/(a2*sqrt(pow(D22,2.)+a2))));
			re = -IR3/(4.*PI) - IR1*pow(kk,2.)/(8.*PI);
			break;
		}
	}

	return re;
}
Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing_Left_No(char xyz,double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,double dxy,Complex kk)
{
	int n=2;
    int js[2]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex *x,s,*a,*b;
	Complex p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Vector_Integ_Up_Down(dxy,j-1,xm,ym,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Con_Vec_Wave_Integ_Sing_Left_No_Function(xyz,x,xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,kk);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Two_D_Con_Vec_Wave_Integ_Sing_Left_No_Function(char a,Complex x[2],double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex kk)
{
	Complex re;
	Complex R,zxy;
    zxy = dfxm*(xm-x[0])+dfym*(ym-x[1])+ddfxxm*pow(xm-x[0],2.0)/2.0+ddfxym*(xm-x[0])*(ym-x[1])+ddfyym*pow(ym-x[1],2.0)/2.0;
	R = sqrt( pow(xm-x[0],2.0)+pow(ym-x[1],2.0)+pow(zxy,2.0) ) + 1.0e-18 ;
	Complex GR_Left = ((i_*kk*R-1.)*exp(i_*kk*R)+1.+pow(kk*R,2.0)/2.0)/(4.*PI*pow(R,3.0));
	switch( a )
	{
	case 'X':
		{
			re = GR_Left * (xm-x[0]);
			break;
		}
	case 'Y':
		{
			re = GR_Left * (ym-x[1]);
			break;
		}
	case 'Z':
		{
			re = GR_Left * (-zxy);
			break;
		}
	}
	
	return re;
}
Complex matrix::Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(string LRUD,double dxy,double xm,double xn,double ym,double yn,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num)
{
	int n=2;
    int js[2]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex *x,s,*a,*b;
	Complex p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Vector_Integ_Up_Down(dxy,j-1,xn,yn,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Conductor_Vector_Integ_Function_Zfs(LRUD,x,xm,ym,dfxm,dfxn,dfxn,dfyn,ddfxxm,ddfxym,ddfyym,wave_num);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Two_D_Conductor_Vector_Integ_Function_Zfs(string LRUD,Complex x[2],double xm,double ym,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num)
{
	Complex re;
	Complex zmn = - dfxm*(x[0]-xm) - dfym*(x[1]-ym) - ddfxxm*(x[0]-xm)*(x[0]-xm)/2.0 - ddfxym*(x[0]-xm)*(x[1]-ym) - ddfyym*(x[1]-ym)*(x[1]-ym)/2.0;
	Complex rs = sqrt(pow(x[0]-xm,2.0) + pow(x[1]-ym,2.0));
	
	if (LRUD == "UL")
	{
		re = G_R(rs,wave_num)*(-dfym*(ym-x[1])-(xm-x[0])*dfxn+zmn);
	}
	else if (LRUD == "UR")
	{
		re = G_R(rs,wave_num)*(dfym*(xm-x[0])-(xm-x[0])*dfyn);
	}
	else if (LRUD == "DL")
	{
		re = G_R(rs,wave_num)*(dfxm*(ym-x[1])-(ym-x[1])*dfxn);
	}
	else
	{
		re = G_R(rs,wave_num)*(-dfxm*(xm-x[0])-(ym-x[1])*dfyn+zmn);
	}	
	
	return re;
}
void matrix::Two_D_Conductor_Vector_SMFIA_CAG_Matrix_Integ(int Para_Num)
{//二维矢量波入射导体粗糙面的矩阵填充
	//当均方根高度h较小时，Bi_CGSTAB_Pre才收敛
	int i,j;
	int Num = Nxy*Nxy;//二维表面总的节点个数
	double dxy = RL/Nxy;
	cout<<"单位间隔为 dxy= "<<dxy/wave_len<<"*Wavelength ;"<<endl;
	
	cout<<"强弱相关距离 RD= "<<RD/wave_len<<"*Wavelength ;"<<endl;
	if ( RD > RL)
	{
		cout<<"RD must be smaller than RL;"<<endl;
		exit(0);//程序终止
	}
	
	//入射波:x,y,z三个分量
	Vec_Com Hix,Hiy,Hiz;
	Hix.assign(Num);
	Hiy.assign(Num);
	Hiz.assign(Num);
	char hix = 'x';
	char hiy = 'y';
	char hiz = 'z';

	ThreeD re;
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		re = Hi_xyz(x[i],y[i],z[i],Nxy,k_w,Inci_HL,Inci_AZ,g_p);
		Hix[i] = re.TD_x;
		Hiy[i] = re.TD_y;
		Hiz[i] = re.TD_z;
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<re.Module(re)<<endl;
	}

	//近似计算磁场：Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiy[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiz[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
	}
	*/
    //锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave(hix,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave(hiy,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave(hiz,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	//近似锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave_Approximation(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave_Approximation(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
//		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
/*
	//平面波
	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Plane_Wave_Hi(hix,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiy[i] = Incidence_Vector_Plane_Wave_Hi(hiy,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiz[i] = Incidence_Vector_Plane_Wave_Hi(hiz,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		outzm<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/

	Vec_Com Vmn,Vmn_U,Vmn_D;
	Vmn.assign(2*Num);
	Vmn_U.assign(Num);
	Vmn_D.assign(Num);
	for (i=0; i<Num; i++)
	{
		Vmn_U[i] = -df_y[i]*Hiz[i]-Hiy[i];
		Vmn_D[i] = df_x[i]*Hiz[i]+Hix[i];
		Vmn[i] = Vmn_U[i];
		Vmn[i+Num] = Vmn_D[i];
	}
	cout<<"入射波计算完毕; "<<endl;
	
	//强矩阵Zs
	int RD_num = int(RD/dxy);//整数倍
	cout<<"强弱相关距离相对单位间隔的整数倍 RD_num= "<<RD_num<<endl;
	//double R_S;
	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";

	//采用带状子矩阵存储
	//由于矩阵对称，则只存储一半
	Vec_Com Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zs_D_LU,Zs_D_RD;
	int Band_Ser,Band_Num;//带状矩阵的序号和子带矩阵的数目
	Complex temp;
	int Zs_bw;
	for( Band_Ser = 0; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			//cout<<"Band_Ser= "<<Band_Ser<<"; "<<" Zs_bw= "<<Zs_bw<<"; "<<" Band_Num= "<<Band_Num<<endl;
			for( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )  
			{
				for( j=Band_Num*Nxy+Band_Ser*Nxy; j < Nxy+Band_Num*Nxy+Band_Ser*Nxy; j++ )
				{
					if( abs((i-Band_Num*Nxy)-(j-Band_Num*Nxy-Band_Ser*Nxy)) <= Zs_bw )
					{
						if( i == j )
						{
							temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_LU.push_back(temp);
							Zs_D_LU.push_back(temp);
							temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_RU.push_back(temp);
							temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_LD.push_back(temp);
							temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_RD.push_back(temp);
							Zs_D_RD.push_back(temp);
						}
						else
						{
							//R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
								//dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
							temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							//outresult<<temp<<endl;
							Zs_LU.push_back(temp);
								//dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
							temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_RU.push_back(temp);
								//dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
							temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_LD.push_back(temp);
								//dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
							temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
							Zs_RD.push_back(temp);
						}
					}
				}
			}
		}
		//cout<<endl;
	}	
	//当Band_Ser=RD_num时；矩阵为一条线
	for( i=0; i < Num-RD_num*Nxy; i++ )
	{
		//R_S=sqrt(pow(x[i]-x[i+RD_num*Nxy],2.0)+pow(y[i]-y[i+RD_num*Nxy],2.0)+pow(z[i]-z[i+RD_num*Nxy],2.0));
			//dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[i+RD_num*Nxy])-(x[i]-x[i+RD_num*Nxy])*df_x[i+RD_num*Nxy]+z[i]-z[i+RD_num*Nxy]);
		temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i+RD_num*Nxy],y[i],y[i+RD_num*Nxy],df_x[i],df_x[i+RD_num*Nxy],df_y[i],df_y[i+RD_num*Nxy],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
		Zs_LU.push_back(temp);
			//dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[i+RD_num*Nxy])-(x[i]-x[i+RD_num*Nxy])*df_y[i+RD_num*Nxy]);
		temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i+RD_num*Nxy],y[i],y[i+RD_num*Nxy],df_x[i],df_x[i+RD_num*Nxy],df_y[i],df_y[i+RD_num*Nxy],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
		Zs_RU.push_back(temp);
			//dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[i+RD_num*Nxy])-(y[i]-y[i+RD_num*Nxy])*df_x[i+RD_num*Nxy]);
		temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i],x[i+RD_num*Nxy],y[i],y[i+RD_num*Nxy],df_x[i],df_x[i+RD_num*Nxy],df_y[i],df_y[i+RD_num*Nxy],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
		Zs_LD.push_back(temp);
			//dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[i+RD_num*Nxy])-(y[i]-y[i+RD_num*Nxy])*df_y[i+RD_num*Nxy]+z[i]-z[i+RD_num*Nxy]);
		temp = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i],x[i+RD_num*Nxy],y[i],y[i+RD_num*Nxy],df_x[i],df_x[i+RD_num*Nxy],df_y[i],df_y[i+RD_num*Nxy],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
		Zs_RD.push_back(temp);
	}
	//块状Toeplitz矩阵Zfs，由于对称，只存储一行
	Vec_Com Zfs_LU,Zfs_RU,Zfs_LD,Zfs_RD;
	Zfs_LU.assign(Num);
	Zfs_RU.assign(Num);
	Zfs_LD.assign(Num);
	Zfs_RD.assign(Num);
	double rou;
	for(i=0; i<Num; i++)
	{
		j=0;
		rou = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
		//cout<<rou<<endl;
		rou = rou / dxy;
		if(rou <= RD_num || fabs(rou-RD_num) < 1.0e-12){}
		else
		{
			//R_S = rou * dxy;
			Zfs_LU[i] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(ul,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				//dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
			Zfs_RU[i] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(ur,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				//dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
			Zfs_LD[i] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(dl,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				//dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
			Zfs_RD[i] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zfs(dr,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				//dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
			//outzm<<Zfs_LD[i]<<"     "<<i<<endl;
		}
	}
	cout<<"阻抗矩阵填充完毕；开始解矩阵方程；"<<endl;

	Two_D_Cond_Vector_Bi_CGSTAB_Pre(Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zfs_LU,Zfs_RU,Zfs_LD,Zfs_RD,Zs_D_LU,Zs_D_RD,Vmn_U,Vmn_D,Nxy,RD_num,eps,Uxy);
/*	for (i=0; i<2*Num; i++)
	{
		cout<<Uxy[i]<<"    "<<i<<endl;
	}*/

	//右端弱相关矩阵向量积
	Vec_Com YM;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(Para_Num,x,y,z,df_x,df_y,Uxy,k_w,dxy,Nxy,RD_num,YM);

	Vec_Com C,r,C0,C0_U,C0_D;
	C.assign(2*Num);
	r.assign(2*Num);
	C0.assign(2*Num);
	C0_U.assign(Num);
	C0_D.assign(Num);
	
	for (i=0; i<2*Num; i++)
	{
		C[i]=Vmn[i]-YM[i];
		r[i]=-YM[i];
		//cout<<YM[i]<<endl;
	}
	
	double Nom=L2_Norm(Vmn);
	double Crit=sqrt(L2_Norm(r)/Nom);
	//cout<<Nom<<endl;
	while(Crit > criterion)
	{
		cout<<"criterion = "<<Crit<<endl;
		outstep<<Crit<<endl;
		Assign_Vecor(C,C0);
		for (i=0; i<Num; i++)
		{
			C0_U[i] = C0[i];
			C0_D[i] = C0[i+Num];
		}
		
		//解矩阵方程
		Two_D_Cond_Vector_Bi_CGSTAB_Pre(Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zfs_LU,Zfs_RU,Zfs_LD,Zfs_RD,Zs_D_LU,Zs_D_RD,C0_U,C0_D,Nxy,RD_num,eps,Uxy);
		
		Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(Para_Num,x,y,z,df_x,df_y,Uxy,k_w,dxy,Nxy,RD_num,YM);

		for (i=0; i<2*Num; i++)
		{
			C[i]=Vmn[i]-YM[i];
			r[i]=C[i]-C0[i];
			//cout<<YM[i]<<"   "<<i<<endl;
		}
		Crit=sqrt(L2_Norm(r)/Nom);
	}

	for (i=0; i<2*Num; i++)
	{
		out_Current<<Uxy[i]<<endl;
	}
	
	Zs_LU.clear();
	Zs_RU.clear();
	Zs_LD.clear();
	Zs_RD.clear();
	Zs_D_LU.clear();
	Zs_D_RD.clear();
	Zfs_LU.clear();
	Zfs_RU.clear();
	Zfs_LD.clear();
	Zfs_RD.clear();
	Hix.clear();
	Hiy.clear();
	Hiz.clear();
	Vmn.clear();
	Vmn_U.clear();
	Vmn_D.clear();
	YM.clear();
	C.clear();
	r.clear();
	C0.clear();
	C0_U.clear();
	C0_D.clear();

}
void matrix::Two_D_Conductor_Vector_SMFIA_CAG_Matrix(int Para_Num)
{//二维矢量波入射导体粗糙面的矩阵填充
	//当均方根高度h较小时，Bi_CGSTAB_Pre才收敛
	int i,j;
	int Num = Nxy*Nxy;//二维表面总的节点个数
	double dxy = RL/Nxy;
	cout<<"单位间隔为 dxy= "<<dxy/wave_len<<"*Wavelength ;"<<endl;
	
	cout<<"强弱相关距离 RD= "<<RD/wave_len<<"*Wavelength ;"<<endl;
	if ( RD > RL)
	{
		cout<<"RD must be smaller than RL;"<<endl;
		exit(0);//程序终止
	}

	//入射波:x,y,z三个分量
	Vec_Com Hix,Hiy,Hiz;
	Hix.assign(Num);
	Hiy.assign(Num);
	Hiz.assign(Num);
	char hix = 'x';
	char hiy = 'y';
	char hiz = 'z';

	//近似计算入射磁场：Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiy[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiz[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
	}
	*/
    //锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave(hix,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave(hiy,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave(hiz,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	//近似锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave_Approximation(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave_Approximation(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
//		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
/*	Complex tt;
    for (double Value_UD = k_w; Value_UD < 30.0*k_w; Value_UD+=k_w/4.0 )
	{
		tt = Incidence_Vector_Wave_Test(hiz,x[1],y[1],z[1],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p,Value_UD);
		cout<<Value_UD/k_w<<"   "<<tt<<endl;
		outzm<<Value_UD/k_w<<"   "<<tt.real()<<endl;
		outvm<<Value_UD/k_w<<"   "<<tt.imag()<<endl;
	}
	//
	
	//平面波
	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Plane_Wave_Hi(hix,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiy[i] = Incidence_Vector_Plane_Wave_Hi(hiy,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiz[i] = Incidence_Vector_Plane_Wave_Hi(hiz,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		outzm<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/

	ThreeD re;
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		re = Hi_xyz(x[i],y[i],z[i],Nxy,k_w,Inci_HL,Inci_AZ,g_p);
		Hix[i] = re.TD_x;
		Hiy[i] = re.TD_y;
		Hiz[i] = re.TD_z;
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<re.Module(re)<<endl;
	}

	Vec_Com Vmn,Vmn_U,Vmn_D;
	Vmn.assign(2*Num);
	Vmn_U.assign(Num);
	Vmn_D.assign(Num);
	for (i=0; i<Num; i++)
	{
		Vmn_U[i] = -df_y[i]*Hiz[i]-Hiy[i];
		Vmn_D[i] = df_x[i]*Hiz[i]+Hix[i];
		Vmn[i] = Vmn_U[i];
		Vmn[i+Num] = Vmn_D[i];
	}
	cout<<"入射波计算完毕; "<<endl;

	//强矩阵Zs
	int RD_num = int(RD/dxy);//整数倍
	cout<<"强弱相关距离相对单位间隔的整数倍 RD_num= "<<RD_num<<endl;
	
	//////////////////////////////////////////////////////////////////////////
	//Zs和Zfs的完全矩阵填充
	//验证Zs与Vector相乘
/* 	int u;
	double rou1,R_S1;
    Vec_Com Zs_Full;
    Zs_Full.assign(2*Num*2*Num);
	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";
	//Zs_Full
	for(i=0; i<2*Num; i++)
	{
		for(j=0; j<2*Num; j++)
		{
			u=i*2*Num+j;
			if ( i<Num && j<Num )
			{
				rou1=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				rou1=rou1/dxy;
				if( rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12 )
				{
					if( i == j )
					{
						Zs_Full[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						//cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					else
					{
						R_S1=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						Zs_Full[u] = dxy*dxy*G_R(R_S1,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					//outresult<<Zs_Full[u]<<"   "<<u<<endl;
				}
			}
			else if ( i<Num && j>=Num )
			{
				rou1=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				rou1=rou1/dxy;
				if( rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12 )
				{
					if( i == j-Num )
					{//对角元素为零
						Zs_Full[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					else
					{//存在零元素
						R_S1=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
						Zs_Full[u] = dxy*dxy*G_R(R_S1,k_w)*(df_y[i]*(x[i]-x[j-Num])-(x[i]-x[j-Num])*df_y[j-Num]);
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					//outresult<<Zs_Full[u]<<"   "<<u<<endl;
				}
			}
			else if ( i>=Num && j<Num )
			{
				rou1=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				rou1=rou1/dxy;
				if( rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12 )
				{
					if ( i-Num == j )
					{//对角元素为零
						Zs_Full[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[j],x[j],y[j],y[j],df_x[j],df_x[j],df_y[j],df_y[j],ddf_x[j],ddf_xy[j],ddf_y[j],k_w);
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					else
					{//存在零元素
						R_S1 = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
						Zs_Full[u] = dxy*dxy*G_R(R_S1,k_w)*(df_x[i-Num]*(y[i-Num]-y[j])-(y[i-Num]-y[j])*df_x[j]);
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					//outresult<<Zs_Full[u]<<"   "<<u<<endl;
				}
			}
			else
			{
				rou1=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				rou1=rou1/dxy;
				if( rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12 )
				{
					if( i-Num == j-Num )
					{
						Zs_Full[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i-Num],x[i-Num],y[i-Num],y[i-Num],df_x[i-Num],df_x[i-Num],df_y[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],k_w);
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					else
					{
						R_S1 = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
						Zs_Full[u] = dxy*dxy*G_R(R_S1,k_w)*(-df_x[i-Num]*(x[i-Num]-x[j-Num])-(y[i-Num]-y[j-Num])*df_y[j-Num]+(z[i-Num]-z[j-Num]));
                        //cout<<Zs_Full[u]<<"    "<<u<<endl;  
					}
					//outresult<<Zs_Full[u]<<"   "<<u<<endl;
				}
			}			
			//display(Zs_Full[u]);
			//outresult<<"  ";
		}
		//outresult<<";"<<endl;
	}

	Vec_Com re1;
	Multiply_Full_Matrix_Vector(' ',Zs_Full,Vmn,re1);
	for ( i=0; i<2*Num; i++ )
	{
		outzm<<re1[i]<<"   "<<i<<endl;
	}*/
    //////////////////////////////////////////////////////////////////////////
    //平面矩阵的完全填充
	//
/*	double rou1;
	vector <Zmn_Cell> Zfs_full;
	Zmn_Cell Zfs_temp;
	for (i=0; i<2*Num; i++)
	{
		for (j=0; j<2*Num; j++)
		{
			if ( i<Num && j<Num )
			{
				rou1=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					Zfs_temp.Zmn_i = i;
					Zfs_temp.Zmn_j = j;
					Zfs_temp.Zmn_z = dxy*dxy*G_R(rou1,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
					Zfs_full.push_back(Zfs_temp);
				}
			}
			else if ( i<Num && j>=Num )
			{
				rou1=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					Zfs_temp.Zmn_i = i;
					Zfs_temp.Zmn_j = j;
					Zfs_temp.Zmn_z = dxy*dxy*G_R(rou1,k_w)*(df_y[i]*(x[i]-x[j-Num])-(x[i]-x[j-Num])*df_y[j-Num]);
					Zfs_full.push_back(Zfs_temp);
				}
			}
			else if ( i>=Num && j<Num )
			{
				rou1=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					Zfs_temp.Zmn_i = i;
					Zfs_temp.Zmn_j = j;
					Zfs_temp.Zmn_z = dxy*dxy*G_R(rou1,k_w)*(df_x[i-Num]*(y[i-Num]-y[j])-(y[i-Num]-y[j])*df_x[j]);
					Zfs_full.push_back(Zfs_temp);
				}
			}
			else
			{
				rou1=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					Zfs_temp.Zmn_i = i;
					Zfs_temp.Zmn_j = j;
					Zfs_temp.Zmn_z = dxy*dxy*G_R(rou1,k_w)*(-df_x[i-Num]*(x[i-Num]-x[j-Num])-(y[i-Num]-y[j-Num])*df_y[j-Num]+(z[i-Num]-z[j-Num]));
					Zfs_full.push_back(Zfs_temp);
				}
			}
		}
	}
	//检验弱相关矩阵向量积
	Vec_Com t3;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(0,x,y,z,df_x,df_y,Vmn,k_w,dxy,Nxy,RD_num,t3);
	//Two_D_Vector_NearMultiplyVector(' ',Zfs_full,Vmn,t3);
    for ( i=0; i<2*Num; i++ )
    {
		outresult<<t3[i]<<"   "<<i<<endl;
    }*/
	
	//远场弱相关矩阵的完全填充Zw
/*    double rou1,R_S1;
	vector <Zmn_Cell> Zw_full;
	Zmn_Cell Zw_temp;
	for (i=0; i<2*Num; i++)
	{
		for (j=0; j<2*Num; j++)
		{
			if ( i<Num && j<Num )
			{
				rou1=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					R_S1=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
					Zw_temp.Zmn_i = i;
					Zw_temp.Zmn_j = j;
					Zw_temp.Zmn_z = dxy*dxy*(G_R(R_S1,k_w)-G_R(rou1,k_w))*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
					Zw_full.push_back(Zw_temp);
				}
			}
			else if ( i<Num && j>=Num )
			{
				rou1=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					R_S1=sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
					Zw_temp.Zmn_i = i;
					Zw_temp.Zmn_j = j;
					Zw_temp.Zmn_z = dxy*dxy*(G_R(R_S1,k_w)-G_R(rou1,k_w))*(df_y[i]*(x[i]-x[j-Num])-(x[i]-x[j-Num])*df_y[j-Num]);
					Zw_full.push_back(Zw_temp);
				}
			}
			else if ( i>=Num && j<Num )
			{
				rou1=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					R_S1=sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
					Zw_temp.Zmn_i = i;
					Zw_temp.Zmn_j = j;
					Zw_temp.Zmn_z = dxy*dxy*(G_R(R_S1,k_w)-G_R(rou1,k_w))*(df_x[i-Num]*(y[i-Num]-y[j])-(y[i-Num]-y[j])*df_x[j]);
					Zw_full.push_back(Zw_temp);
				}
			}
			else
			{
				rou1=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0));
				//cout<<rou<<endl;
				rou1=rou1/dxy;
				if(rou1 <= RD_num || fabs(rou1-RD_num) < 1.0e-12){}
				else
				{
					rou1 = rou1 * dxy;
					R_S1=sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
					Zw_temp.Zmn_i = i;
					Zw_temp.Zmn_j = j;
					Zw_temp.Zmn_z = dxy*dxy*(G_R(R_S1,k_w)-G_R(rou1,k_w))*(-df_x[i-Num]*(x[i-Num]-x[j-Num])-(y[i-Num]-y[j-Num])*df_y[j-Num]+(z[i-Num]-z[j-Num]));
					Zw_full.push_back(Zw_temp);
				}
			}
		}
	}
	//检验弱相关矩阵向量积
	Vec_Com t1;
	//Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(Para_Num,x,y,z,df_x,df_y,Vmn,k_w,dxy,Nxy,RD_num,t1);
	Two_D_Vector_NearMultiplyVector(' ',Zw_full,Vmn,t1);
    for ( i=0; i<2*Num; i++ )
    {
		outresult<<t1[i]<<"   "<<i<<endl;
    }*/
	//////////////////////////////////////////////////////////////////////////
	//填充近场强相关矩阵
	//采用带状子矩阵存储
	//在矢量波入射的积分方程中，矩阵是不对称的，这要与标量波入射相区别，必须全部存储
	//子矩阵也为带状矩阵，所以子矩阵也进行了带状存储
	double R_S;
	Zmn_Cell zsmn_t;
	vector <Zmn_Cell> Zs_LU,Zs_RU,Zs_LD,Zs_RD;
	Vec_Com Zs_D_LU,Zs_D_RD;
	int Band_Ser,Band_Num;//带状矩阵的序号和子带矩阵的数目
	Complex temp;
	int Zs_bw;
	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";
	//第零带
	Band_Ser = 0;
	for( Band_Num = 0; Band_Num < Nxy; Band_Num++ )
	{
		Zs_bw = RD_num;
		for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
		{
			if (i-Band_Num*Nxy <= Zs_bw)
			{
				for ( j=Band_Num*Nxy; j <= Band_Num*Nxy+Zs_bw+(i-Band_Num*Nxy); j++ )
				{
					if ( i == j )
					{
						//temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = Two_D_Con_Vec_Wave_Integ_Sing(ul,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,k_w);
						//cout<<temp<<"   "<<0.5+Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w)<<endl;
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						Zs_D_LU.push_back(temp);
							
						//temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = 0.0;
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs_RU.push_back(zsmn_t);

						//temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = 0.0;
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);

						//temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = Two_D_Con_Vec_Wave_Integ_Sing(dr,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,k_w);
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);

						Zs_D_RD.push_back(temp);
					}
					else
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
			}
			else
			{
				for ( j=Band_Num*Nxy+(i-Band_Num*Nxy)-Zs_bw; j < Band_Num*Nxy+((i-Band_Num*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Nxy; j++ )
				{
					if ( i == j )
					{
						//temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = Two_D_Con_Vec_Wave_Integ_Sing(ul,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,k_w);
						//outresult<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						Zs_D_LU.push_back(temp);
							
						//temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = 0.0;
						zsmn_t.Zmn_z = temp;
						//cout<<temp<<endl;
						Zs_RU.push_back(zsmn_t);

						//temp=0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = 0.0;
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);

						//temp=0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
						temp = Two_D_Con_Vec_Wave_Integ_Sing(dr,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],dxy,k_w);
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);

						Zs_D_RD.push_back(temp);
					}
					else
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
			}
		}
	}			
	//从左到右的带
	for( Band_Ser = 1; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
			{
				if (i-Band_Num*Nxy <= Zs_bw)
				{
					for ( j=Band_Num*Nxy+Band_Ser*Nxy; j <= Band_Num*Nxy+Band_Ser*Nxy+Zs_bw+(i-Band_Num*Nxy); j++ )
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
				else
				{
					for ( j=Band_Num*Nxy+Band_Ser*Nxy+(i-Band_Num*Nxy)-Zs_bw; j < Band_Num*Nxy+Band_Ser*Nxy+((i-Band_Num*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Band_Ser*Nxy+Nxy; j++ )
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
			}				
		}
	}

	//当Band_Ser=RD_num时；矩阵为一条线
	for( i=0; i < Num-RD_num*Nxy; i++ )
	{
		R_S=sqrt(pow(x[i]-x[i+RD_num*Nxy],2.0)+pow(y[i]-y[i+RD_num*Nxy],2.0)+pow(z[i]-z[i+RD_num*Nxy],2.0));
		temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[i+RD_num*Nxy])-(x[i]-x[i+RD_num*Nxy])*df_x[i+RD_num*Nxy]+z[i]-z[i+RD_num*Nxy]);
		zsmn_t.Zmn_i = i;
		zsmn_t.Zmn_j = i+RD_num*Nxy;
		zsmn_t.Zmn_z = temp;
		Zs_LU.push_back(zsmn_t);
		//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
		//////////////////////////////////////////////////////////////////////////
		//Zs_LU.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[i+RD_num*Nxy])-(x[i]-x[i+RD_num*Nxy])*df_y[i+RD_num*Nxy]);
		zsmn_t.Zmn_z = temp;
		Zs_RU.push_back(zsmn_t);
		//Zs_RU.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[i+RD_num*Nxy])-(y[i]-y[i+RD_num*Nxy])*df_x[i+RD_num*Nxy]);
		zsmn_t.Zmn_z = temp;
		Zs_LD.push_back(zsmn_t);
		//cout<<temp<<endl;
		//Zs_LD.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[i+RD_num*Nxy])-(y[i]-y[i+RD_num*Nxy])*df_y[i+RD_num*Nxy]+z[i]-z[i+RD_num*Nxy]);
		zsmn_t.Zmn_z = temp;
		Zs_RD.push_back(zsmn_t);
		//Zs_RD.push_back(temp);
	}
	//从上到下的带
	for( Band_Ser = 1; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			for ( i=Band_Num*Nxy+Band_Ser*Nxy; i < Nxy+Band_Num*Nxy+Band_Ser*Nxy; i++ )
			{
				if (i-Band_Num*Nxy-Band_Ser*Nxy <= Zs_bw)
				{
					for ( j=Band_Num*Nxy; j <= Band_Num*Nxy+Zs_bw+(i-Band_Num*Nxy-Band_Ser*Nxy); j++ )
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
				else
				{
					for ( j=Band_Num*Nxy+(i-Band_Num*Nxy-Band_Ser*Nxy)-Zs_bw; j < Band_Num*Nxy+((i-Band_Num*Nxy-Band_Ser*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Nxy; j++ )
					{
						R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
						//cout<<temp<<endl;
						zsmn_t.Zmn_i = i;
						zsmn_t.Zmn_j = j;
						zsmn_t.Zmn_z = temp;
						Zs_LU.push_back(zsmn_t);
						//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
						//////////////////////////////////////////////////////////////////////////
						//Zs_LU.push_back(temp);
							
						temp=dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
						zsmn_t.Zmn_z = temp;
						Zs_RU.push_back(zsmn_t);
						//Zs_RU.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
						zsmn_t.Zmn_z = temp;
						Zs_LD.push_back(zsmn_t);
						//Zs_LD.push_back(temp);
						temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
						zsmn_t.Zmn_z = temp;
						Zs_RD.push_back(zsmn_t);
						//Zs_RD.push_back(temp);
					}
				}
			}				
		}
	}

	//当Band_Ser=RD_num时；矩阵为一条线
	for( j=0; j < Num-RD_num*Nxy; j++ )
	{
		R_S=sqrt(pow(x[j+RD_num*Nxy]-x[j],2.0)+pow(y[j+RD_num*Nxy]-y[j],2.0)+pow(z[j+RD_num*Nxy]-z[j],2.0));
		temp=dxy*dxy*G_R(R_S,k_w)*(-df_y[j+RD_num*Nxy]*(y[j+RD_num*Nxy]-y[j])-(x[j+RD_num*Nxy]-x[j])*df_x[j]+z[j+RD_num*Nxy]-z[j]);
		zsmn_t.Zmn_i = j+RD_num*Nxy;
		zsmn_t.Zmn_j = j;
		zsmn_t.Zmn_z = temp;
		Zs_LU.push_back(zsmn_t);
		//outresult<<Zs_LU[Zs_LU.size()-1].Zmn_z<<"  "<<i*2*Num+j<<"    "<<i<<"    "<<j<<endl;
		//////////////////////////////////////////////////////////////////////////
		//Zs_LU.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(df_y[j+RD_num*Nxy]*(x[j+RD_num*Nxy]-x[j])-(x[j+RD_num*Nxy]-x[j])*df_y[j]);
		zsmn_t.Zmn_z = temp;
		Zs_RU.push_back(zsmn_t);
		//Zs_RU.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(df_x[j+RD_num*Nxy]*(y[j+RD_num*Nxy]-y[j])-(y[j+RD_num*Nxy]-y[j])*df_x[j]);
		zsmn_t.Zmn_z = temp;
		Zs_LD.push_back(zsmn_t);
		//cout<<temp<<endl;
		//Zs_LD.push_back(temp);

		temp=dxy*dxy*G_R(R_S,k_w)*(-df_x[j+RD_num*Nxy]*(x[j+RD_num*Nxy]-x[j])-(y[j+RD_num*Nxy]-y[j])*df_y[j]+z[j+RD_num*Nxy]-z[j]);
		zsmn_t.Zmn_z = temp;
		Zs_RD.push_back(zsmn_t);
		//Zs_RD.push_back(temp);
	}
	//近场强相关矩阵填充完毕
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//验证强矩阵向量积：结果正确
/*	Vec_Com t1,t2,t3,t4,t5;
	t5.assign(2*Num);
	Two_D_Vector_NearMultiplyVector(' ',Zs_LU,Vmn_U,t1);
	Two_D_Vector_NearMultiplyVector(' ',Zs_RU,Vmn_D,t2);
	Two_D_Vector_NearMultiplyVector(' ',Zs_LD,Vmn_U,t3);
	Two_D_Vector_NearMultiplyVector(' ',Zs_RD,Vmn_D,t4);
		
	for (i=0; i<Num; i++)
	{
		t5[i] = t1[i] + t2[i];
		t5[i+Num] = t3[i] + t4[i];
	}
	for ( i=0; i<2*Num; i++ )
	{
		outresult<<t5[i]<<"  "<<i<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
	//块状Toeplitz矩阵Zfs，由于不对称，还需讨论
	//但可以利用2-D FFT快速计算，相当于按照零阶泰勒展开，然后运用2-D FFT
/*	Vec_Com Zfs_LU,Zfs_RU,Zfs_LD,Zfs_RD;
	Zfs_LU.assign(Num);
	Zfs_RU.assign(Num);
	Zfs_LD.assign(Num);
	Zfs_RD.assign(Num);
	double rou;
	for(i=0; i<Num; i++)
	{
		j=0;
		rou = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
		//cout<<rou<<endl;
		rou = rou/dxy;
		if(rou <= RD_num || fabs(rou-RD_num) < 1.0e-12){}
		else
		{
			R_S = rou * dxy;
			Zfs_LU[i] = dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
			Zfs_RU[i] = dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j])-(x[i]-x[j])*df_y[j]);
			Zfs_LD[i] = dxy*dxy*G_R(R_S,k_w)*(df_x[i]*(y[i]-y[j])-(y[i]-y[j])*df_x[j]);
			Zfs_RD[i] = dxy*dxy*G_R(R_S,k_w)*(-df_x[i]*(x[i]-x[j])-(y[i]-y[j])*df_y[j]+(z[i]-z[j]));
			//outzm<<Zfs_LD[i]<<"     "<<i<<endl;
		}
		
	}
	//////////////////////////////////////////////////////////////////////////
	//验证平面矩阵向量积
	Vec_Com t1,t2,t3,t4,t5;
	t5.assign(2*Num);
	Multi_BTTB_Vec(Zfs_LU,Vmn_U,Nxy,t1);
	Multi_BTTB_Vec(Zfs_RU,Vmn_D,Nxy,t2);
	Multi_BTTB_Vec(Zfs_LD,Vmn_U,Nxy,t3);
	Multi_BTTB_Vec(Zfs_RD,Vmn_D,Nxy,t4);
	for (i=0; i<Num; i++)
	{
		t5[i] = t1[i] + t2[i];
		t5[i+Num] = t3[i] + t4[i];
	}
	for (i=0; i<2*Num; i++)
	{
		outresult<<t5[i]<<"   "<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
	
    cout<<"阻抗矩阵填充完毕；开始解矩阵方程；"<<endl;
	//解矩阵方程
	Two_D_Cond_Vector_Bi_CGSTAB_Pre(Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zs_D_LU,Zs_D_RD,Vmn_U,Vmn_D,x,y,z,df_x,df_y,k_w,dxy,Nxy,RD_num,eps,Uxy);

	//右端弱相关矩阵向量积
	Vec_Com YM;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(Para_Num,x,y,z,df_x,df_y,Uxy,k_w,dxy,Nxy,RD_num,YM);
    
	Vec_Com C,r,C0,C0_U,C0_D;
	C.assign(2*Num);
	r.assign(2*Num);
	C0.assign(2*Num);
	C0_U.assign(Num);
	C0_D.assign(Num);
	
	for (i=0; i<2*Num; i++)
	{
		C[i]=Vmn[i]-YM[i];
		r[i]=-YM[i];
		//cout<<YM[i]<<endl;
	}
	
	double Nom=L2_Norm(Vmn);
	double Crit=sqrt(L2_Norm(r)/Nom);
	//cout<<Nom<<endl;
	while(Crit > criterion)
	{
		cout<<"criterion = "<<Crit<<endl;
		outstep<<Crit<<endl;
		Assign_Vecor(C,C0);
		for (i=0; i<Num; i++)
		{
			C0_U[i] = C0[i];
			C0_D[i] = C0[i+Num];
		}
		
		//解矩阵方程
		Two_D_Cond_Vector_Bi_CGSTAB_Pre(Zs_LU,Zs_RU,Zs_LD,Zs_RD,Zs_D_LU,Zs_D_RD,C0_U,C0_D,x,y,z,df_x,df_y,k_w,dxy,Nxy,RD_num,eps,Uxy);
		
		Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(Para_Num,x,y,z,df_x,df_y,Uxy,k_w,dxy,Nxy,RD_num,YM);

		for (i=0; i<2*Num; i++)
		{
			C[i]=Vmn[i]-YM[i];
			r[i]=C[i]-C0[i];
			//cout<<YM[i]<<"   "<<i<<endl;
		}
		Crit=sqrt(L2_Norm(r)/Nom);
	}

	for (i=0; i<2*Num; i++)
	{
		out_Current<<Uxy[i]<<"  "<<i<<endl;
	}
	for (i=0; i<Num; i++)
	{
		out_current_mag1<<x[i]<<"  "<<y[i]<<"  "<<abs(Uxy[i])<<endl;
		out_current_mag2<<x[i]<<"  "<<y[i]<<"  "<<abs(Uxy[i+Num])<<endl;
	}
	
	Zs_LU.clear();
	Zs_RU.clear();
	Zs_LD.clear();
	Zs_RD.clear();
	Zs_D_LU.clear();
	Zs_D_RD.clear();
	Hix.clear();
	Hiy.clear();
	Hiz.clear();
	Vmn.clear();
	Vmn_U.clear();
	Vmn_D.clear();
	YM.clear();
	C.clear();
	r.clear();
	C0.clear();
	C0_U.clear();
	C0_D.clear();

}

Complex matrix::Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(string LRUD,double dxy,double xm,double xn,double ym,double yn,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num)
{
	int n=2;
    int js[4]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Vector_Integ_Up_Down(dxy,j-1,xn,yn,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Conductor_Vector_Integ_Function_Zs(LRUD,x,xm,ym,dfxm,dfxn,dfym,dfyn,ddfxxm,ddfxym,ddfyym,wave_num);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Two_D_Conductor_Vector_Integ_Function_Zs(string LRUD,Complex x[2],double xm,double ym,double dfxm,double dfxn,double dfym,double dfyn,double ddfxxm,double ddfxym,double ddfyym,double wave_num)
{
	Complex re;
	Complex zmn = - dfxm*(x[0]-xm) - dfym*(x[1]-ym) - ddfxxm*(x[0]-xm)*(x[0]-xm)/2.0 - ddfxym*(x[0]-xm)*(x[1]-ym) - ddfyym*(x[1]-ym)*(x[1]-ym)/2.0;
//	Complex zmn = - dfxm*(x[0]-xm) - dfym*(x[1]-ym);
	Complex rs = sqrt(pow(x[0]-xm,2.0) + pow(x[1]-ym,2.0) + pow(zmn,2.0)) + 1.0e-18;
	
	if (LRUD == "UL")
	{
		re = G_R(rs,wave_num)*(-dfym*(ym-x[1])-(xm-x[0])*dfxn+zmn);
	}
	else if (LRUD == "UR")
	{
		re = G_R(rs,wave_num)*(dfym*(xm-x[0])-(xm-x[0])*dfyn);
	}
	else if (LRUD == "DL")
	{
		re = G_R(rs,wave_num)*(dfxm*(ym-x[1])-(ym-x[1])*dfxn);
	}
	else
	{
		re = G_R(rs,wave_num)*(-dfxm*(xm-x[0])-(ym-x[1])*dfyn+zmn);
	}
    
	/*
	switch(LRUD)
	{
	case "UL":{re = G_R(rs,wave_num)*(-dfyn*(yn-x[1])-(xn-x[0])*dfxm+zmn); break;}
	case "UR":{re = G_R(rs,wave_num)*(dfyn*(xn-x[0])-(xn-x[0])*dfym); break;}
	case "DL":{re = G_R(rs,wave_num)*(dfxn*(yn-x[1])-(yn-x[1])*dfxm); break;}
	case "DR":{re = G_R(rs,wave_num)*(-dfxn*(xn-x[0])-(yn-x[1])*dfym+zmn); break;}
	default:{}
	}*/
	
	return re;
}
void matrix::Two_D_Conductor_Vector_Integ_Up_Down(double dxy,int j,double xn,double yn,double y[2])
{
	switch(j)
	{
	case 0:{y[0] = xn-dxy/2.0;y[1] = xn+dxy/2.0;break;}
    case 1:{y[0] = yn-dxy/2.0;y[1] = yn+dxy/2.0;break;}
    case 2:{y[0]=-1.0;y[1]=1.0;break;}
    case 3:{y[0]=-1.0;y[1]=1.0;break;}
    default:{}
	}
}
void matrix::Two_D_Conductor_Vector_MoM_Matrix_Integ()
{
	int i,j,u;
	int Num = Nxy*Nxy;
	Vec_Com Zmn,Vmn;
	Zmn.assign(2*Num*2*Num);
	Vmn.assign(2*Num);
	double dxy = RL/Nxy;
    Complex R_S;

	//V矩阵的填充
	//入射波:x,y,z三个分量
	Vec_Com Hix,Hiy,Hiz;
	Hix.assign(Num);
	Hiy.assign(Num);
	Hiz.assign(Num);
	char hix = 'x';
	char hiy = 'y';
	char hiz = 'z';

	ThreeD re;
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		re = Hi_xyz(x[i],y[i],z[i],Nxy,k_w,Inci_HL,Inci_AZ,g_p);
		Hix[i] = re.TD_x;
		Hiy[i] = re.TD_y;
		Hiz[i] = re.TD_z;
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<re.Module(re)<<endl;
	}

	//近似计算磁场：Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiy[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiz[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
	}
	*/
    //锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave(hix,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave(hiy,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave(hiz,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	//近似锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave_Approximation(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave_Approximation(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
//		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
/*
	//平面波
	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Plane_Wave_Hi(hix,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiy[i] = Incidence_Vector_Plane_Wave_Hi(hiy,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiz[i] = Incidence_Vector_Plane_Wave_Hi(hiz,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		outzm<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	
	for (i=0; i<Num; i++)
	{
		Vmn[i] = -df_y[i]*Hiz[i]-Hiy[i];
		Vmn[i+Num] = df_x[i]*Hiz[i]+Hix[i];
	}

	cout<<"入射波计算完毕; "<<endl;

	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";
	//Z矩阵的填充
	for(i=0; i<2*Num; i++)
	{
		for(j=0; j<2*Num; j++)
		{
			u = i*2*Num+j;
			if(i < Num && j < Num)
			{
				if(i == j)
				{
					Zmn[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
		    	else
				{
					//R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
		        	Zmn[u] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[j],y[i],y[j],df_x[i],df_x[j],df_y[i],df_y[j],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
			}
			else if(i < Num && j >= Num)
			{
				if(i == (j-Num))
				{
    				Zmn[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
    			else
				{
     			    //R_S = sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
					Zmn[u] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[j-Num],y[i],y[j-Num],df_x[i],df_x[j-Num],df_y[i],df_y[j-Num],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
			}
			else if(i >= Num && j < Num)
			{
				if((i-Num) == j)
				{
	    			Zmn[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[j],x[j],y[j],y[j],df_x[j],df_x[j],df_y[j],df_y[j],ddf_x[j],ddf_xy[j],ddf_y[j],k_w);
				}
	    		else
				{
					//R_S = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
					Zmn[u] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[i-Num],x[j],y[i-Num],y[j],df_x[i-Num],df_x[j],df_y[i-Num],df_y[j],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],k_w);
				}
				
			}
			else
			{
				if((i-Num) == (j-Num))
				{
		    		Zmn[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i-Num],x[i-Num],y[i-Num],y[i-Num],df_x[i-Num],df_x[i-Num],df_y[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],k_w);
				}
	    		else
				{
		    		//R_S = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
					Zmn[u] = Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i-Num],x[j-Num],y[i-Num],y[j-Num],df_x[i-Num],df_x[j-Num],df_y[i-Num],df_y[j-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],k_w);
				}
			}
			//outzm<<Zmn[u]<<"   "<<u<<endl;
		}
		//outzm<<endl;
	}
	cout<<"阻抗矩阵填充完毕；开始解矩阵方程；"<<endl;
	/////////////////////////////////////////////////////////////////////
	//解方程组
	time_t ltime1,ltime2,ltime;//存储时间
	time( &ltime1 );
	cout<<"The time of starting resolving equation is "<<ctime( &ltime1 )<<endl;
	resolve_m(Zmn,Vmn,Uxy);
	time( &ltime2 );
	cout<<"The time of stoping resolving equation is "<<ctime( &ltime2 )<<endl;
	ltime=ltime2-ltime1;
	cout<<ltime<<" seconds"<<endl;
	////////////////////////////////////////////////////////////////////
    for (i=0; i<2*Num; i++)
	{
		out_Current<<Uxy[i]<<"   "<<i<<endl;
	}

	Zmn.clear();
	Vmn.clear();
	Hix.clear();
	Hiy.clear();
	Hiz.clear();
}
void matrix::Two_D_Conductor_Vector_MoM_Matrix()
{
	int i,j,u;
	int Num = Nxy*Nxy;
	Vec_Com Zmn,Vmn;
	Zmn.assign(2*Num*2*Num);
	Vmn.assign(2*Num);
	double dxy = RL/Nxy;
    double R_S;

	//V矩阵的填充
	//入射波:x,y,z三个分量
	Vec_Com Hix,Hiy,Hiz;
	Hix.assign(Num);
	Hiy.assign(Num);
	Hiz.assign(Num);
	char hix = 'x';
	char hiy = 'y';
	char hiz = 'z';

	ThreeD re;
	//徐鹏方法
	for (i=0; i<Num; i++)
	{
		re = Hi_xyz(x[i],y[i],z[i],Nxy,k_w,Inci_HL,Inci_AZ,g_p);
		Hix[i] = re.TD_x;
		Hiy[i] = re.TD_y;
		Hiz[i] = re.TD_z;
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
//		outresult<<x[i]<<"  "<<y[i]<<"  "<<re.Module(re)<<endl;
	}

	//近似计算磁场：Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiy[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
		Hiz[i] = Incidence_Vector_Wave_Approximation_Hi_All_Polarization(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p,Po_H,Po_V);
//		outresult<<Hix[i]<<"  "<<Hiy[i]<<"  "<<Hiz[i]<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Hix[i])+norm(Hiy[i])+norm(Hiz[i]))<<endl;
	}
	*/
    //锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave(hix,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave(hiy,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave(hiz,x[i],y[i],z[i],k_w,RL,Nxy,Inci_HL,Inci_AZ,g_p);
		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	//近似锥形波
/*	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Wave_Approximation(hix,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiy[i] = Incidence_Vector_Wave_Approximation(hiy,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
		Hiz[i] = Incidence_Vector_Wave_Approximation(hiz,k_w,x[i],y[i],z[i],Inci_HL,Inci_AZ,g_p);
//		cout<<Hix[i]<<"  "<<i<<endl;
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
/*
	//平面波
	for (i=0; i<Num; i++)
	{
		Hix[i] = Incidence_Vector_Plane_Wave_Hi(hix,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiy[i] = Incidence_Vector_Plane_Wave_Hi(hiy,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		Hiz[i] = Incidence_Vector_Plane_Wave_Hi(hiz,x[i],y[i],z[i],k_w,Inci_HL,Inci_AZ);
		outzm<<x[i]<<"  "<<y[i]<<"  "<<sqrt(pow(abs(Hix[i]),2.0)+pow(abs(Hiy[i]),2.0)+pow(abs(Hiz[i]),2.0))<<endl;
	}*/
	
	for (i=0; i<Num; i++)
	{
		Vmn[i] = -df_y[i]*Hiz[i]-Hiy[i];
		Vmn[i+Num] = df_x[i]*Hiz[i]+Hix[i];
	}

	cout<<"入射波计算完毕; "<<endl;

	string ul = "UL";
	string ur = "UR";
	string dl = "DL";
	string dr = "DR";
	//Z矩阵的填充
	for(i=0; i<2*Num; i++)
	{
		for(j=0; j<2*Num; j++)
		{
			u=i*2*Num+j;
			if(i < Num && j < Num)
			{
				if(i == j)
				{
					Zmn[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ul,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
		    	else
				{
					R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
		        	Zmn[u] = dxy*dxy*G_R(R_S,k_w)*(-df_y[i]*(y[i]-y[j])-(x[i]-x[j])*df_x[j]+(z[i]-z[j]));
				}
			}
			else if(i < Num && j >= Num)
			{
				if(i == (j-Num))
				{
    				Zmn[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(ur,dxy,x[i],x[i],y[i],y[i],df_x[i],df_x[i],df_y[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
				}
    			else
				{
     			    R_S = sqrt(pow(x[i]-x[j-Num],2.0)+pow(y[i]-y[j-Num],2.0)+pow(z[i]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*G_R(R_S,k_w)*(df_y[i]*(x[i]-x[j-Num])-(x[i]-x[j-Num])*df_y[j-Num]);
				}
			}
			else if(i >= Num && j < Num)
			{
				if((i-Num) == j)
				{
	    			Zmn[u] = 0.0 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dl,dxy,x[j],x[j],y[j],y[j],df_x[j],df_x[j],df_y[j],df_y[j],ddf_x[j],ddf_xy[j],ddf_y[j],k_w);
				}
	    		else
				{
					R_S = sqrt(pow(x[i-Num]-x[j],2.0)+pow(y[i-Num]-y[j],2.0)+pow(z[i-Num]-z[j],2.0));
					Zmn[u] = dxy*dxy*G_R(R_S,k_w)*(df_x[i-Num]*(y[i-Num]-y[j])-(y[i-Num]-y[j])*df_x[j]);
				}
				
			}
			else
			{
				if((i-Num) == (j-Num))
				{
		    		Zmn[u] = 0.5 + Two_D_Conductor_Vector_Integ_Matrix_Element_Zs(dr,dxy,x[i-Num],x[i-Num],y[i-Num],y[i-Num],df_x[i-Num],df_x[i-Num],df_y[i-Num],df_y[i-Num],ddf_x[i-Num],ddf_xy[i-Num],ddf_y[i-Num],k_w);
				}
	    		else
				{
		    		R_S = sqrt(pow(x[i-Num]-x[j-Num],2.0)+pow(y[i-Num]-y[j-Num],2.0)+pow(z[i-Num]-z[j-Num],2.0));
					Zmn[u] = dxy*dxy*G_R(R_S,k_w)*(-df_x[i-Num]*(x[i-Num]-x[j-Num])-(y[i-Num]-y[j-Num])*df_y[j-Num]+(z[i-Num]-z[j-Num]));
				}
			}
			//outzm<<Zmn[u]<<"  ";
		}
		//outzm<<endl;
	}
	cout<<"阻抗矩阵填充完毕，开始解矩阵方程；"<<endl;

	/////////////////////////////////////////////////////////////////////
	//解方程组
	time_t ltime1,ltime2,ltime;//存储时间
	time( &ltime1 );
	cout<<"The time of starting resolving equation is "<<ctime( &ltime1 )<<endl;
	resolve_m(Zmn,Vmn,Uxy);
	time( &ltime2 );
	cout<<"The time of stoping resolving equation is "<<ctime( &ltime2 )<<endl;
	ltime=ltime2-ltime1;
	cout<<ltime<<" seconds"<<endl;
	////////////////////////////////////////////////////////////////////
    for (i=0; i<2*Num; i++)
	{
		out_Current<<Uxy[i]<<"   "<<i<<endl;
	}

	Zmn.clear();
	Vmn.clear();
	Hix.clear();
	Hiy.clear();
	Hiz.clear();
}
Complex matrix::Incidence_Vector_Plane_Wave_Hi(char xyz,double xm,double ym,double zm,double wave_num,double In_HL,double In_AZ)
{
	Complex re;
	switch(xyz)
	{
	case 'x' :{re = cos(In_HL)*cos(In_AZ) * exp(-i_*wave_num*(xm*sin(In_HL)*cos(In_AZ)+ym*sin(In_HL)*sin(In_AZ)-zm*cos(In_HL))) / 377.0; break;}
    case 'y' :{re = cos(In_HL)*sin(In_AZ) * exp(-i_*wave_num*(xm*sin(In_HL)*cos(In_AZ)+ym*sin(In_HL)*sin(In_AZ)-zm*cos(In_HL))) / 377.0; break;}
    case 'z' :{re = sin(In_HL) * exp(-i_*wave_num*(xm*sin(In_HL)*cos(In_AZ)+ym*sin(In_HL)*sin(In_AZ)-zm*cos(In_HL))) / (120.0*PI); break;}
    default:{}
	}
	return re;
}
Complex matrix::Incidence_Vector_Plane_Wave_Ei(char xyz,double xm,double ym,double zm,double wave_num,double In_HL,double In_AZ)
{
	Complex re;
	switch(xyz)
	{
	case 'x' :{re = -sin(In_AZ) * exp(-i_*wave_num*(xm*sin(In_HL)*cos(In_AZ)+ym*sin(In_HL)*sin(In_AZ)-zm*cos(In_HL))); break;}
    case 'y' :{re =  cos(In_AZ) * exp(-i_*wave_num*(xm*sin(In_HL)*cos(In_AZ)+ym*sin(In_HL)*sin(In_AZ)-zm*cos(In_HL))); break;}
    case 'z' :{re = 0.0; break;}
    default:{}
	}
	return re;
}
Complex matrix::Incidence_Vector_Wave(char xyz,double xm,double ym,double zm,double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g)
{
	int n=4;
    int js[4]={2,2,2,2};
    int m,j,k,q,l,*is;
    double y[4];
	Complex s,*x,p,*a,*b;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Up_Down_Incidence_Vector_Wave(j-1,y,x,rl,nxy,wave_num,In_HL);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Incidence_Vector_Wave_Integral(xyz,x,wave_num,xm,ym,zm,In_HL,In_AZ,Value_g);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Incidence_Vector_Wave_Test(char xyz,double xm,double ym,double zm,double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g,double Value_UD)
{
	int n=4;
    int js[4]={2,2,2,2};
    int m,j,k,q,l,*is;
    double y[4];
	Complex s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Up_Down_Incidence_Vector_Wave_Test(j-1,y,rl,nxy,Value_UD);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Incidence_Vector_Wave_Integral(xyz,x,wave_num,xm,ym,zm,In_HL,In_AZ,Value_g);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
double matrix::Incidence_Vector_Wave_Pinc_h(double wave_num,double rl,int nxy,double In_HL,double In_AZ,double Value_g)
{
	int n=4;
    int js[4]={4,4,4,4};
    int m,j,k,q,l,*is;
    double y[4];
	double s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new double[n];
    a=new double[2*(n+1)];
    b=new double[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Up_Down_Incidence_Vector_Pinc_h(j-1,y,x,rl,nxy,wave_num);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Incidence_Vector_Integ_Pinc_h(x,wave_num,In_HL,In_AZ,Value_g);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Incidence_Vector_Wave_Integral(char xyz,Complex xup[4],double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g)
{
	Complex re;
	Complex Krou = sqrt(xup[0]*xup[0] + xup[1]*xup[1]);
	Complex Kz;
	if (Krou.real() <= wave_num)
	{
		Kz = sqrt(wave_num*wave_num - Krou*Krou);
	}
	else
	{
		Kz = - i_*sqrt(Krou*Krou - wave_num*wave_num);
	}
	
	Complex tx = pow(cos(In_HL)*cos(In_AZ)*xup[2] + cos(In_HL)*sin(In_AZ)*xup[3],2.0) / pow(Value_g*cos(In_HL),2.0);
	Complex ty = pow(-sin(In_AZ)*xup[2] + cos(In_AZ)*xup[3],2.0) / pow(Value_g,2.0);
	Complex Value_W = ((2.0*tx -1.0) / pow(Value_g*cos(In_HL),2.0) + (2.0*ty - 1.0) / pow(Value_g,2.0)) / pow(wave_num,2.0);
	Complex Value_t = (xup[2]*xup[2] + xup[3]*xup[3]) / pow(Value_g,2.0);
	Complex Inti_Function;
	Inti_Function = exp(i_*(xup[0]*xm + xup[1]*ym - Kz*zm)) * exp(-i_*(xup[0]*xup[2] + xup[1]*xup[3])) 
		* exp(i_*(wave_num*sin(In_HL)*cos(In_AZ)*xup[2] + wave_num*sin(In_HL)*sin(In_AZ)*xup[3]) * (1.0 + Value_W)) * exp(-Value_t);
	
	switch(xyz)
	{
	case 'x' :{re = Inti_Function * Kz * xup[0] / (wave_num * Krou); break;}
    case 'y' :{re = Inti_Function * Kz * xup[1] / (wave_num * Krou); break;}
    case 'z' :{re = Inti_Function * Krou / wave_num; break;}
    default:{}
	}
	return - re / (120.0 * PI * 4.0 * PI * PI);//377为空气的本征阻抗
}
Complex matrix::Incidence_Vector_Wave_Approximation(char xyz,double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g)
{
	Complex re;
	Complex Sr,Uz;
	double kix = wave_num*sin(In_HL)*cos(In_AZ);
	double kiy = wave_num*sin(In_HL)*sin(In_AZ);
	double kiz = wave_num*cos(In_HL);
	Sr = (1.0-2.0*i_*(pow(wave_num,2.0)-pow(kix,2.0))*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*pow(xm+kix*zm/kiz,2.0)
		+(1.0-2.0*i_*(pow(wave_num,2.0)-pow(kix,2.0))*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*pow(ym+kiy*zm/kiz,2.0)
		+4.0*i_*kix*kiy*zm*(xm+kix*zm/kiz)*(ym+kiy*zm/kiz)/(pow(Value_g,2.0)*pow(kiz,3.0));
	Uz = sqrt((1.0-2.0*i_*pow(wave_num,2.0)*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*(1.0-2.0*i_*zm/(pow(Value_g,2.0)*kiz)));
	Complex Inti_Function;
	Inti_Function = exp(i_*(kix*xm + kiy*ym - kiz*zm)) * exp(-Sr/pow(Value_g*Uz,2.0)) / Uz;
	
	switch(xyz)
	{//有问题
	case 'x' :{re = Inti_Function * cos(In_HL) * cos(In_AZ); break;}
    case 'y' :{re = Inti_Function * cos(In_HL) * sin(In_AZ); break;}
    case 'z' :{re = Inti_Function * sin(In_HL); break;}
    default:{}
	}
	return - re / (120.0 * PI);//377为空气的本征阻抗
}
Complex matrix::Incidence_Vector_Wave_Approximation_Hi_All_Polarization(char xyz,double wave_num,double xm,double ym,double zm,double In_HL,double In_AZ,double Value_g,int Pola_H,int Pola_V)
{//Pola_H=1，Pola_V=0为水平极化波；Pola_H＝0，Pola_V＝1为垂直极化波
	Complex re;
	int m,n,k,j,l;
	Complex Sr,Uz,U1z2,U2z2,Wyz;
	double kix = wave_num*sin(In_HL)*cos(In_AZ);
	double kiy = wave_num*sin(In_HL)*sin(In_AZ);
	double kiz = wave_num*cos(In_HL);
	Sr = (1.0-2.0*i_*(pow(wave_num,2.0)-pow(kix,2.0))*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*pow(xm+kix*zm/kiz,2.0)
		+(1.0-2.0*i_*(pow(wave_num,2.0)-pow(kiy,2.0))*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*pow(ym+kiy*zm/kiz,2.0)
		+4.0*i_*kix*kiy*zm*(xm+kix*zm/kiz)*(ym+kiy*zm/kiz)/(pow(Value_g,2.0)*pow(kiz,3.0));
	Uz = sqrt((1.0-2.0*i_*pow(wave_num,2.0)*zm/(pow(Value_g,2.0)*pow(kiz,3.0)))*(1.0-2.0*i_*zm/(pow(Value_g,2.0)*kiz)));
	U1z2 = 1.0 - 2.0*i_*(wave_num*wave_num-kix*kix)*zm / (pow(Value_g,2.0)*pow(kiz,3.0));
	U2z2 = pow(Uz,2.0) / U1z2;
	Wyz = 2.0*i_*kix*kiy*zm*(ym + kiy*zm/kiz) / (pow(Value_g,2.0)*U1z2*pow(kiz,3.0));

	//////////////////////////////////////////////////////////////////////////
	ThreeD a[3][3]; 
	ThreeD ei = ThreeD_hkp(wave_num,kix,kiy,In_HL,In_AZ) * Pola_H + ThreeD_vkp(wave_num,kix,kiy,In_HL,In_AZ) * Pola_V;
	ThreeD X(1.0,0.0,0.0),Y(0.0,1.0,0.0),Z(0.0,0.0,1.0);
	a[0][0] = ei;
	a[1][0] = ((X*(ei*X)-Z*(ei*Z))*(-2.0*kix*kiz)) / (pow(wave_num,2.0)*kiz) - ((Y*(ei*X)+X*(ei*Y))*(kiy*kiz)) / (pow(wave_num,2.0)*kiz) 
		+ ((Z*(ei*X)+X*(ei*Z))*(kiz*kiz-kix*kix)) / (pow(wave_num,2.0)*kiz) - ((Z*(ei*Y)+Y*(ei*Z))*(kix*kiy)) / (pow(wave_num,2.0)*kiz);

	a[0][1] = ((Y*(ei*Y)-Z*(ei*Z))*(-2.0*kiy*kiz)) / (pow(wave_num,2.0)*kiz) - ((Y*(ei*X)+X*(ei*Y))*(kix*kiz)) / (pow(wave_num,2.0)*kiz) 
		- ((Z*(ei*X)+X*(ei*Z))*(kix*kiy)) / (pow(wave_num,2.0)*kiz) + ((Z*(ei*Y)+Y*(ei*Z))*(kiz*kiz-kiy*kiy)) / (pow(wave_num,2.0)*kiz);
	
	a[1][1] = ((Y*(ei*X)+X*(ei*Y))*(-pow(kiz,3.0))) / (pow(wave_num,2.0)*pow(kiz,3.0)) - ((Z*(ei*X)+X*(ei*Z))*(kiy*(kix*kix+kiz*kiz))) / (pow(wave_num,2.0)*pow(kiz,3.0))
		- ((Z*(ei*Y)+Y*(ei*Z))*(kix*(kiy*kiy+kiz*kiz))) / (pow(wave_num,2.0)*pow(kiz,3.0));
	a[2][0] = ((X*(ei*X)-Z*(ei*Z))*(-2.0*pow(kiz,3.0))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0)) - ((Z*(ei*X)+X*(ei*Z))*(kix*(kix*kix+3.0*kiz*kiz))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0))
		- ((Z*(ei*Y)+Y*(ei*Z))*(kiy*(wave_num*wave_num-kiy*kiy))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0));
	a[0][2] = ((Y*(ei*Y)-Z*(ei*Z))*(-2.0*pow(kiz,3.0))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0)) - ((Z*(ei*X)+X*(ei*Z))*(kix*(wave_num*wave_num-kix*kix))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0))
		- ((Z*(ei*Y)+Y*(ei*Z))*(kiy*(kiy*kiy+3.0*kiz*kiz))) / (2.0*pow(wave_num,2.0)*pow(kiz,3.0));
/*	for (m=0; m<=2; m++)
	{
		for (n=0; n<=2-m; n++)
		{
			outzm<<a[m][n].TD_x<<"  "<<a[m][n].TD_y<<"  "<<a[m][n].TD_z<<endl;
		}
	}
	cout<<ei*Z<<endl;*/
	//////////////////////////////////////////////////////////////////////////
	
	ThreeD Pmn,Temp;
	Complex te1,te2,te3;
	//
	for (m=0; m<=2; m++)
	{
		for (n=0; n<=2-m; n++)
		{
			te1 = 0.0;
			for (k=0; k<=n/2; k++)
			{
				te2 = 0.0;
				for (j=0; j<= n-2*k; j++)
				{
					te3 = 0.0;
					for (l=0; l<=(m+j)/2; l++)
					{
						te3 = te3 + pow(-Value_g*Value_g*U2z2,l)/(pow(4.0,l)*factor(m+j-2*l)*factor(l))*pow(xm+kix*zm/kiz+Wyz,m+j-2*l);
					}
					te2 = te2 + factor(m+j)/(factor(j)*factor(n-2*k-j))*pow(ym+kiy*zm/kiz,n-2*k-j)
						*pow(kix*kiy*zm/pow(kiz,3.0),j)*pow(2.0*i_/(Value_g*Value_g*U2z2),m+j) * te3;
				}
				te1 = te1 + pow(-Value_g*Value_g*U1z2,k)/(pow(4.0,k)*factor(k)) * te2;
			}
			Pmn = Pmn + a[m][n] * (factor(n)*pow(2.0*i_/(Value_g*Value_g*U1z2),n) * te1);
		}
	}
	
	Complex Inti_Function;
	Inti_Function = exp(i_*(kix*xm + kiy*ym - kiz*zm)) * exp(-Sr/pow(Value_g*Uz,2.0)) / Uz;

	Temp = Pmn * Inti_Function;
	//磁场
	ThreeD kii(sin(In_HL)*cos(In_AZ),sin(In_HL)*sin(In_AZ),-cos(In_HL));
	Temp = kii.Xc(kii,Temp) / (120.0*PI);
//	cout<<Temp.TD_x<<"  "<<Temp.TD_y<<"  "<<Temp.TD_z<<endl;

	switch(xyz)
	{//有问题
	case 'x' :{re = Temp.TD_x; break;}
    case 'y' :{re = Temp.TD_y; break;}
    case 'z' :{re = Temp.TD_z; break;}
    default:{}
	}
	return re;//377为空气的本征阻抗
}
ThreeD matrix::ThreeD_hkp(double wave_num,double kx,double ky,double in_hl,double in_az)
{
	ThreeD re;
	double kp;
	kp = sqrt(kx*kx + ky*ky);
	if (kp == 0.0)
	{
		re.Assign(sin(in_az),-cos(in_az),0.0);
	}
	else
	{
		re.Assign(ky/kp,-kx/kp,0.0);
	}
	return re;
}
ThreeD matrix::ThreeD_vkp(double wave_num,double kx,double ky,double in_hl,double in_az)
{
	ThreeD re;
	double kp;
	kp = sqrt(kx*kx + ky*ky);
	Complex kz;
	if (kp <= wave_num)
	{
		kz = sqrt(wave_num*wave_num - kp*kp);
	}
	else
	{
		kz = - i_*sqrt(kp*kp - wave_num*wave_num);
	}
	if (kp == 0.0)
	{
		re.Assign(cos(in_az),sin(in_az),0.0);
	}
	else
	{
		re.Assign(kz*kx/(wave_num*kp),kz*ky/(wave_num*kp),kp/wave_num);
	}
	return re;
}
double matrix::Incidence_Vector_Integ_Pinc_h(double xup[4],double wave_num,double In_HL,double In_AZ,double Value_g)
{
	double re;
	double Krou = sqrt(xup[0]*xup[0] + xup[1]*xup[1]);
	double Kz = sqrt(wave_num*wave_num - Krou*Krou);
	double tx = pow(cos(In_HL)*cos(In_AZ)*xup[2] + cos(In_HL)*sin(In_AZ)*xup[3],2.0) / pow(Value_g*cos(In_HL),2.0);
	double ty = pow(-sin(In_AZ)*xup[2] + cos(In_AZ)*xup[3],2.0) / pow(Value_g,2.0);
	double Value_W = ((2.0*tx -1.0) / pow(Value_g*cos(In_HL),2.0) + (2.0*ty - 1.0) / pow(Value_g,2.0)) / pow(wave_num,2.0);
	double Value_t = (xup[2]*xup[2] + xup[3]*xup[3]) / pow(Value_g,2.0);
	Complex Inti_Function;
	Inti_Function = exp(-i_*(xup[0]*xup[2] + xup[1]*xup[3])) * exp(i_*(sin(In_HL)*cos(In_AZ)*xup[2] + sin(In_HL)*sin(In_AZ)*xup[3])
		* (1.0 + Value_W)) * exp(-Value_t);
	re = pow(abs(Inti_Function),2.0);
	re = re * Kz / wave_num;

	return re;
}
void matrix::Up_Down_Incidence_Vector_Wave(int j,double y[2],Complex xup[2],double rl,int nxy,double wave_num,double in_hl)
{
	Complex dky;
	double wave_num_sin;
	wave_num_sin = wave_num - wave_num * sin(in_hl);
//	wave_num_sin = wave_num;
	dky=sqrt(wave_num_sin*wave_num_sin-xup[0]*xup[0]);
	switch(j)
	{
//	case 0:{y[0] = -wave_num_sin; y[1] = wave_num_sin; break;}
//  case 1:{y[0] = -dky.real(); y[1] = dky.real(); break;}
	case 0:{y[0] = -PI*nxy/rl + 2.0*PI/rl; y[1] = PI*nxy/rl; break;}
    case 1:{y[0] = -PI*nxy/rl + 2.0*PI/rl; y[1] = PI*nxy/rl; break;}
//	case 0:{y[0] = -wave_num; y[1] = wave_num; break;}
//  case 1:{y[0] = -wave_num; y[1] = wave_num; break;}
	case 2:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    case 3:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    default:{}
	}
}
void matrix::Up_Down_Incidence_Vector_Wave_Test(int j,double y[2],double rl,int nxy,double Value_UD)
{
	switch(j)
	{
	case 0:{y[0] = -Value_UD; y[1] = Value_UD; break;}
    case 1:{y[0] = -Value_UD; y[1] = Value_UD; break;}
	case 2:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    case 3:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    default:{}
	}
}
void matrix::Up_Down_Incidence_Vector_Pinc_h(int j,double y[2],double xup[2],double rl,int nxy,double wave_num)
{
	double dky;
	dky=sqrt(wave_num*wave_num-xup[0]*xup[0]);
	switch(j)
	{
	case 0:{y[0] = -wave_num; y[1] = wave_num; break;}
    case 1:{y[0] = -dky; y[1] = dky; break;}
    case 2:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    case 3:{y[0] = -rl/2.0 + rl/nxy/2.0; y[1] = rl/2.0 - rl/nxy/2.0; break;}
    default:{}
	}
}
ThreeD matrix::Hi_xyz(double xi,double yi,double zi,int nxy,double wave_num,double In_HL,double In_AZ,double Value_g)
{
	Complex *Ete = new Complex[nxy*nxy];
	Complex *Ete_FFT = new Complex[nxy*nxy];
	int i,j,u;
	double tx,ty,Value_W,Value_t;
	for (i=0; i<nxy*nxy; i++)
	{
		tx = pow(cos(In_HL)*cos(In_AZ)*x[i] + cos(In_HL)*sin(In_AZ)*y[i],2.0) / pow(Value_g*cos(In_HL),2.0);
		ty = pow(-sin(In_AZ)*x[i] + cos(In_AZ)*y[i],2.0) / pow(Value_g,2.0);
		Value_W = ((2.0*tx - 1.0) / pow(Value_g*cos(In_HL),2.0) + (2.0*ty - 1.0) / pow(Value_g,2.0)) / pow(wave_num,2.0);
		Value_t = tx + ty;
		//Value_t = (x[i]*x[i] + y[i]*y[i]) / pow(Value_g,2.0);
		Ete[i] = exp(i_*wave_num*(sin(In_HL)*cos(In_AZ)*x[i] + sin(In_HL)*sin(In_AZ)*y[i]) * (1.0 + Value_W)) * exp(-Value_t);
	}

	//经过周期扩展后所得结果基本出现锥形，但结果仍有许多毛刺
	Complex *xt1=new Complex[nxy*nxy];
	Complex *xt2=new Complex[nxy*nxy];
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i <= nxy/2)
			{
				xt1[u]=Ete[(i+nxy/2-1)*nxy+j];
			}
			else
			{
				xt1[u]=Ete[(i-nxy/2-1)*nxy+j];
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

	fft2(nxy,xt2,xt2,1);
	//fft2(nxy,Ete,Ete_FFT,1);

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
				Ete_FFT[u]=xt1[i*nxy+j+nxy/2+1];
			}
	    	else
			{
				Ete_FFT[u]=xt1[i*nxy+j-nxy/2+1];
			}
		}
	}
	delete [] xt1;
	delete [] xt2;


    
/*	for (i=0; i<nxy*nxy; i++)
    {
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Ete_FFT[i]))<<endl;
		//display(Ete_FFT[i]);
		//outresult<<"   "<<i<<endl;
    }*/
	
	double kx,ky;
	double Krou;
	Complex Kz;
	double dkxy = 2.0*PI/RL;
	ThreeD re,hkz;
	for (i=-nxy/2+1; i<=nxy/2; i++)
	{
		for (j=-nxy/2+1; j<=nxy/2; j++)
		{
			u = (i+nxy/2-1)*nxy + (j+nxy/2-1);
			kx = wave_num*sin(In_HL)*cos(In_AZ) + 2.0*PI*i/RL;
			ky = wave_num*sin(In_HL)*sin(In_AZ) + 2.0*PI*j/RL;
			Krou = sqrt(kx*kx + ky*ky);
			Kz = Isqrt(wave_num*wave_num - Krou*Krou);//这个开方结果较好
			//Kz = sqrt(fabs(wave_num*wave_num - Krou*Krou));
			//cout<<wave_num<<"  "<<Krou<<"  "<<Kz<<endl;
/*			if (Krou <= wave_num)
			{
				Kz = sqrt(wave_num*wave_num - Krou*Krou);
			}
			else
			{
				Kz = - i_*sqrt(Krou*Krou - wave_num*wave_num);
			}*/
			hkz.Assign(Kz*kx/(wave_num*Krou),Kz*ky/(wave_num*Krou),Krou/wave_num);
			re = re - hkz * exp(i_*(kx*xi+ky*yi-Kz*zi)) * Ete_FFT[u] * dkxy * dkxy;
			//outresult<<wave_num<<"  ";
		}
		//outresult<<endl;
	}
	
	delete [] Ete;
	delete [] Ete_FFT;
	
	return re / (120.0*PI);
}
ThreeD matrix::Ei_xyz(double xi,double yi,double zi,int nxy,double wave_num,double In_HL,double In_AZ,double Value_g)
{
	Complex *Ete = new Complex[nxy*nxy];
	Complex *Ete_FFT = new Complex[nxy*nxy];
	int i,j,u;
	double tx,ty,Value_W,Value_t;
	for (i=0; i<nxy*nxy; i++)
	{
		tx = pow(cos(In_HL)*cos(In_AZ)*x[i] + cos(In_HL)*sin(In_AZ)*y[i],2.0) / pow(Value_g*cos(In_HL),2.0);
		ty = pow(-sin(In_AZ)*x[i] + cos(In_AZ)*y[i],2.0) / pow(Value_g,2.0);
		Value_W = ((2.0*tx - 1.0) / pow(Value_g*cos(In_HL),2.0) + (2.0*ty - 1.0) / pow(Value_g,2.0)) / pow(wave_num,2.0);
		Value_t = tx + ty;
		//Value_t = (x[i]*x[i] + y[i]*y[i]) / pow(Value_g,2.0);
		Ete[i] = exp(i_*wave_num*(sin(In_HL)*cos(In_AZ)*x[i] + sin(In_HL)*sin(In_AZ)*y[i]) * (1.0 + Value_W)) * exp(-Value_t);
	}

	//经过周期扩展后所得结果基本出现锥形，但结果仍有许多毛刺
	Complex *xt1=new Complex[nxy*nxy];
	Complex *xt2=new Complex[nxy*nxy];
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i <= nxy/2)
			{
				xt1[u]=Ete[(i+nxy/2-1)*nxy+j];
			}
			else
			{
				xt1[u]=Ete[(i-nxy/2-1)*nxy+j];
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

	fft2(nxy,xt2,xt2,1);
	//fft2(nxy,Ete,Ete_FFT,1);

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
				Ete_FFT[u]=xt1[i*nxy+j+nxy/2+1];
			}
	    	else
			{
				Ete_FFT[u]=xt1[i*nxy+j-nxy/2+1];
			}
		}
	}
	delete [] xt1;
	delete [] xt2;


    
/*	for (i=0; i<nxy*nxy; i++)
    {
		outresult<<x[i]<<"  "<<y[i]<<"  "<<sqrt(norm(Ete_FFT[i]))<<endl;
		//display(Ete_FFT[i]);
		//outresult<<"   "<<i<<endl;
    }*/
	
	double kx,ky;
	double Krou;
	Complex Kz;
	double dkxy = 2.0*PI/RL;
	ThreeD re,hkz;
	for (i=-nxy/2+1; i<=nxy/2; i++)
	{
		for (j=-nxy/2+1; j<=nxy/2; j++)
		{
			u = (i+nxy/2-1)*nxy + (j+nxy/2-1);
			kx = wave_num*sin(In_HL)*cos(In_AZ) + 2.0*PI*i/RL;
			ky = wave_num*sin(In_HL)*sin(In_AZ) + 2.0*PI*j/RL;
			Krou = sqrt(kx*kx + ky*ky);
			Kz = Isqrt(wave_num*wave_num - Krou*Krou);//这个开方结果较好
			//Kz = sqrt(fabs(wave_num*wave_num - Krou*Krou));
			//cout<<wave_num<<"  "<<Krou<<"  "<<Kz<<endl;
/*			if (Krou <= wave_num)
			{
				Kz = sqrt(wave_num*wave_num - Krou*Krou);
			}
			else
			{
				Kz = - i_*sqrt(Krou*Krou - wave_num*wave_num);
			}*/
			hkz.Assign(ky/Krou,-kx/Krou,0.);
			re = re - hkz * exp(i_*(kx*xi+ky*yi-Kz*zi)) * Ete_FFT[u] * dkxy * dkxy;
			//outresult<<wave_num<<"  ";
		}
		//outresult<<endl;
	}
	
	delete [] Ete;
	delete [] Ete_FFT;
	
	return re;
}
double matrix::Two_D_Conductor_Vector_Ph_inc(int nxy,double wave_num,double In_HL,double In_AZ,double Value_g)
{
	double Ph_inc;
	Complex *Ete = new Complex[nxy*nxy];
	int i,j,u;
	double tx,ty,Value_W,Value_t;
	for (i=0; i<nxy*nxy; i++)
	{
		tx = pow(cos(In_HL)*cos(In_AZ)*x[i] + cos(In_HL)*sin(In_AZ)*y[i],2.0) / pow(Value_g*cos(In_HL),2.0);
		ty = pow(-sin(In_AZ)*x[i] + cos(In_AZ)*y[i],2.0) / pow(Value_g,2.0);
		Value_W = ((2.0*tx - 1.0) / pow(Value_g*cos(In_HL),2.0) + (2.0*ty - 1.0) / pow(Value_g,2.0)) / pow(wave_num,2.0);
		Value_t = tx + ty;
		//Value_t = (x[i]*x[i] + y[i]*y[i]) / pow(Value_g,2.0);
		Ete[i] = exp(i_*wave_num*(sin(In_HL)*cos(In_AZ)*x[i] + sin(In_HL)*sin(In_AZ)*y[i]) * (1.0 + Value_W)) * exp(-Value_t);
	}

	//经过周期扩展后所得结果基本出现锥形，但结果仍有许多毛刺
	Complex *xt1=new Complex[nxy*nxy];
	Complex *xt2=new Complex[nxy*nxy];
	for(i=0; i<nxy; i++)
	{
		for(j=0; j<nxy; j++)
		{
			u=i*nxy+j;
			if(i <= nxy/2)
			{
				xt1[u]=Ete[(i+nxy/2-1)*nxy+j];
			}
			else
			{
				xt1[u]=Ete[(i-nxy/2-1)*nxy+j];
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

	fft2(nxy,xt2,xt2,1);
	//fft2(nxy,Ete,Ete_FFT,1);

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
				Ete[u]=xt1[i*nxy+j+nxy/2+1];
			}
	    	else
			{
				Ete[u]=xt1[i*nxy+j-nxy/2+1];
			}
		}
	}
	delete [] xt1;
	delete [] xt2;

	double kx,ky;
	double Krou;
	double Kz;
	double dkxy = 2.0*PI/RL;
	Ph_inc = 0.0;
	for (i=-nxy/2+1; i<=nxy/2; i++)
	{
		for (j=-nxy/2+1; j<=nxy/2; j++)
		{
			u = (i+nxy/2-1)*nxy + (j+nxy/2-1);
			kx = wave_num*sin(In_HL)*cos(In_AZ) + 2.0*PI*i/RL;
			ky = wave_num*sin(In_HL)*sin(In_AZ) + 2.0*PI*j/RL;
			Krou = sqrt(kx*kx + ky*ky);
			//Kz = Isqrt(wave_num*wave_num - Krou*Krou);
			if (Krou < wave_num)
			{
				Kz = sqrt(wave_num*wave_num - Krou*Krou);
				Ph_inc = Ph_inc + dkxy * dkxy * pow(abs(Ete[u]),2.0) * Kz / wave_num;
			}
		}
	}
	
	delete [] Ete;

	return Ph_inc * 2.0 * PI * PI / (120.0 * PI);
}
void matrix::Two_D_Conductor_Scalar_SMFIA_CAG_Matrix(int Para_Num)
{//二维标量波入射导体粗糙面的矩阵填充
	//Zs,Zfs,Zw矩阵的填充
	//Zs
	int i,j;
	int Num = Nxy*Nxy;//二维表面总的节点个数
	double dxy = RL/Nxy;
	cout<<"单位间隔为 dxy= "<<dxy/wave_len<<"*Wavelength ;"<<endl;
	
	cout<<"强相关距离 RD= "<<RD/wave_len<<"*Wavelength ;"<<endl;
	if ( RD > RL)
	{
		cout<<"RD must be smaller than RL;"<<endl;
		exit(0);//程序终止
	}
	int RD_num = int(RD/dxy);//整数倍
	cout<<"强相关距离 RD_num= "<<RD_num<<endl;
	double rou;
	double R_S;
    
	//v矩阵的填充
	Vec_Com Vmn;
	Vmn.assign(Num);
	//锥形波
	double xsz,ww,kg,tx,ty;
	for (i=0; i<Num; i++)
	{
		xsz=z[i]*cos(Inci_HL)-x[i]*sin(Inci_HL)*cos(Inci_AZ)-y[i]*sin(Inci_HL)*sin(Inci_AZ);
		kg=pow(g_p*cos(Inci_HL),2);
		tx=pow(x[i]*cos(Inci_HL)*cos(Inci_AZ)+y[i]*cos(Inci_HL)*sin(Inci_AZ)+z[i]*sin(Inci_HL),2)/kg;
		ty=pow((-x[i]*sin(Inci_AZ)+y[i]*cos(Inci_AZ))/g_p,2);        
    	ww=((2.0*tx-1.0)/kg+(2.0*ty-1.0)/pow(g_p,2))/pow(k_w,2);
		Vmn[i] = exp(-i_*k_w*xsz*(1.0+ww))*exp(-(tx+ty));
		//outresult<<x[i]/wave_len<<"  "<<y[i]/wave_len<<"  "<<abs(Vmn[i])<<endl;
		/*display(Vmn[i]);
		outresult<<";"<<endl;*/
	}
	//平面波：计算结果不正确，没有考虑因截断而引起的电流突变，故只能用锥形波作为入射波；所以计算结果与文献不一致可能是二维粗糙面的生成引起的。
/*  for (i=0; i<Num; i++)
    {
		Vmn[i] = exp(-i_*k_w*(x[i]*sin(Inci_HL)*cos(Inci_AZ)+y[i]*sin(Inci_HL)*sin(Inci_AZ)-z[i]*cos(Inci_HL)));
    }*/
	
	//////////////////////////////////////////////////////////////////////////
	//Zs和Zfs的完全矩阵填充
	//验证Zs与Vector相乘
/*	int u;
    Vec_Com Zs_Full;
    Zs_Full.assign(Num*Num);
	for(i=0; i<Num; i++)
	{
		for(j=0; j<Num; j++)
		{
			u=i*Num+j;
			rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
			rou=rou/dxy;
			if( rou <= RD_num || fabs(rou-RD_num) < 1.0e-12 )
			{
				if( i == j )
				{
					Zs_Full[u] = Two_D_Conductor_Scalar_Integ_Sing_Left(dxy,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w)
								+Two_D_Conductor_Scalar_Integ_Sing(dxy,df_x[i],df_y[i]);
				}
				else
				{
					R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
					Zs_Full[u] = dxy*dxy*G_0_fs_R(R_S,k_w);
				}
			}
			//display(Zs_Full[u]);
			//outresult<<"  ";
		}
		//outresult<<";"<<endl;
	}
	Vec_Com re1;
	Multiply_Full_Matrix_Vector(' ',Zs_Full,Vmn,re1);
	for ( i=0; i<Num; i++ )
	{
		outresult<<re1[i]<<"   "<<i<<endl;
	}

	//块状Toeplitz矩阵
	//显示块状 Toeplitz 矩阵
	Vec_Com Zfs_Full;
	Zfs_Full.assign(Num*Num);
	for(i=0;i<Num;i++)
	{
		for(j=0;j<Num;j++)
		{
			u=i*Num+j;
			rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
			rou=rou/dxy;
			if(rou <= RD_num || fabs(rou-RD_num) < 1.0e-12){}
			else
			{
				R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
				Zfs_Full[u]=dxy*dxy*G_0_fs_R(R_S,k_w);
			}
			display(Zfs_Full[u]);
			outresult<<"  ";
		}
        outresult<<";"<<endl;
	}
	//弱矩阵的完全填充Zw
	Vec_Com Zw;
	Zw.assign(Num*Num);

	for (i=0; i<Num; i++)
	{
		for (j=0; j<Num; j++)
		{
			rou=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
   	    	//cout<<rou<<endl;
   	    	rou=rou/dxy;
        	if(rou <= RD_num || fabs(rou-RD_num) < 1.0e-12){}
         	else
			{
	        	rou = rou * dxy;
	    		R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
	        	Zw[i*Num+j] = dxy*dxy*( G_0_fs_R(R_S,k_w) - G_0_fs_R(rou,k_w) );
			}
        	//display(Zw[i*Num+j]);
		}
		outzm<<Zw[i*Num+0]<<endl;
	}
	//////////////////////////////////////////////////////////////////////////
*/
	//////////////////////////////////////////////////////////////////////////
	//采用带状子矩阵存储
	//由于矩阵对称，则只存储一半
	//子矩阵也为带状矩阵，所以子矩阵也进行了带状存储
	Zmn_Cell zsmn_t;
	vector <Zmn_Cell> Zs;
	Vec_Com Zs_Diag;
	int Band_Ser,Band_Num;//带状矩阵的序号和各个带中子带矩阵的序号
	Complex temp;
	int Zs_bw;
	for( Band_Ser = 0; Band_Ser < RD_num; Band_Ser++ )
	{
		for( Band_Num = 0; Band_Num < Nxy-Band_Ser; Band_Num++ )
		{
			Zs_bw = int( sqrt(RD_num*RD_num-Band_Ser*Band_Ser) );
			if ( Band_Ser == 0 )
			{
				for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
				{
					for ( j=i; j <= i + Zs_bw && j < (Band_Num+1)*Nxy; j++ )
					{
						if ( i == j )
						{
							temp = Two_D_Conductor_Scalar_Integ_Sing_Left(dxy,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w)
								+Two_D_Conductor_Scalar_Integ_Sing(y[i],dxy,df_x[i],df_y[i]);
							//cout<<Two_D_Conductor_Scalar_Integ_Sing(y[i],dxy,df_x[i],df_y[i])<<"   "<<Two_D_Conductor_Scalar_Integ_Sing(dxy,df_x[i],df_y[i])<<endl;
							zsmn_t.Zmn_i = i;
							zsmn_t.Zmn_j = j;
							zsmn_t.Zmn_z = temp;
							Zs.push_back(zsmn_t);
							Zs_Diag.push_back(temp);
							//cout<<temp<<"   "<<i<<"  "<<j<<endl;
						}
						else
						{
							R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
							temp=dxy*dxy*G_0_fs_R(R_S,k_w);
							//temp = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(dxy,x[i],y[i],x[j],y[j],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);

							zsmn_t.Zmn_i = i;
							zsmn_t.Zmn_j = j;
							zsmn_t.Zmn_z = temp;
							Zs.push_back(zsmn_t);
							//cout<<temp<<"   "<<i<<"  "<<j<<endl;
						}
						//outresult<<temp<<"   "<<i<<"  "<<j<<endl;
					}
				}
			}
			else
			{
				for ( i=Band_Num*Nxy; i < Nxy+Band_Num*Nxy; i++ )
				{
					if (i-Band_Num*Nxy <= Zs_bw)
					{
						for ( j=Band_Num*Nxy+Band_Ser*Nxy; j <= Band_Num*Nxy+Band_Ser*Nxy+Zs_bw+(i-Band_Num*Nxy); j++ )
						{
							R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
							temp=dxy*dxy*G_0_fs_R(R_S,k_w);
							//temp = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(dxy,x[i],y[i],x[j],y[j],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);

							zsmn_t.Zmn_i = i;
							zsmn_t.Zmn_j = j;
							zsmn_t.Zmn_z = temp;
							Zs.push_back(zsmn_t);
							//cout<<temp<<"   "<<i<<"  "<<j<<endl;
						}
					}
					else
					{
						for ( j=Band_Num*Nxy+Band_Ser*Nxy+(i-Band_Num*Nxy)-Zs_bw; j < Band_Num*Nxy+Band_Ser*Nxy+((i-Band_Num*Nxy)-Zs_bw)+(2*Zs_bw+1) && j < Band_Num*Nxy+Band_Ser*Nxy+Nxy; j++ )
						{
							R_S=sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
							temp=dxy*dxy*G_0_fs_R(R_S,k_w);
							//temp = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(dxy,x[i],y[i],x[j],y[j],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);

							zsmn_t.Zmn_i = i;
							zsmn_t.Zmn_j = j;
							zsmn_t.Zmn_z = temp;
							Zs.push_back(zsmn_t);
							//cout<<temp<<"   "<<i<<"  "<<j<<endl;
						}
					}
				}				
			}
		}
	}

	//当Band_Ser=RD_num时；矩阵为一条线
	for( i=0; i < Num-RD_num*Nxy; i++ )
	{
		R_S=sqrt(pow(x[i]-x[i+RD_num*Nxy],2.0)+pow(y[i]-y[i+RD_num*Nxy],2.0)+pow(z[i]-z[i+RD_num*Nxy],2.0));
		temp=dxy*dxy*G_0_fs_R(R_S,k_w);
		//j = i+RD_num*Nxy;
		//temp = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(dxy,x[i],y[i],x[j],y[j],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);

		zsmn_t.Zmn_i = i;
		zsmn_t.Zmn_j = i+RD_num*Nxy;
		zsmn_t.Zmn_z = temp;
		Zs.push_back(zsmn_t);
		//cout<<temp<<"   "<<i<<"  "<<j<<endl;
		//Zs_Inver.push_back(conj(temp));
		//cout<<Zs[t]<<"   "<<t<<endl;
		//t++;
	}
	cout<<"Zs的维数= "<<Zs.size()<<endl;//Zs的维数
	out_Reflectivity<<"Zs的维数= "<<Zs.size()<<" ;"<<endl;

	//////////////////////////////////////////////////////////////////////////
    //检验矩阵向量积
/*	cout<<zsmn.size()<<endl;
	outresult<<zsmn.size()<<endl;
//	for (i=0; i<Zsmn.size(); i++)
//	{
//		outresult<<Zsmn[i].Zmn_z<<"   "<<Zsmn[i].Zmn_i<<"  "<<Zsmn[i].Zmn_j<<endl;
//	}

	Vec_Com tt;
	//Multi_Zs_Vec(Zs,Vmn,RD_num,Nxy,tt);
	NearMultiplyVector(' ',Zs,Vmn,tt);

	for (i=0; i<Num; i++)
	{
		outresult<<tt[i]<<"   "<<i<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
	//块状Toeplitz矩阵，由于对称，只存储一行
	Vec_Com Zfs,Zfs_Inver;
	Zfs.assign(Num);
	Zfs_Inver.assign(Num);
	for(i=0; i<Num; i++)
	{
		j=0;
		rou = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0));
		//cout<<rou<<endl;
		rou = rou/dxy;
		if(rou <= RD_num || fabs(rou-RD_num) < 1.0e-12){}
		else
		{
			R_S = rou * dxy;
			Zfs[i] = dxy*dxy*G_0_fs_R(R_S,k_w);
			//Zfs[i] = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zfs(dxy,x[i],y[i],x[j],y[j],k_w);
			//Zfs_Inver[i] = conj(Zfs[i]);
		}
		//cout<<Zfs[i]<<endl;
	}
	//////////////////////////////////////////////////////////////////////////
		
/*	Vec_Com result;
	Multi_Zs_Vec(Zs,Vmn,RD_num,Nxy,result);
	for(i=0;i<Num;i++)
		outvm<<result[i]<<"   "<<i<<endl;*/

	//////////////////////////////////////////////////////////////////////////
	cout<<"矩阵填充完毕，开始解矩阵方程;"<<endl;	
	//解矩阵方程
	//Two_D_Cond_CGS_Pre(Zs,Zfs,Zs_Diag,Vmn,nxy,RD_num,eps,Uxy);
	Two_D_Cond_Bi_CGSTAB_Pre(Zs,Zfs,Zs_Diag,Vmn,Nxy,RD_num,eps,Uxy);
	//Two_D_Cond_CGM(Zs,Zfs,Zs_Inver,Zfs_Inver,Zs_Diag,Vmn,nxy,RD_num,eps,Uxy);

	//右端弱相关矩阵向量积
	Vec_Com YM;
	Two_D_Conductor_Scalar_Multi_Zw_Vec(Para_Num,x,y,z,Uxy,k_w,dxy,Nxy,RD_num,YM);
	//////////////////////////////////////////////////////////////////////////
	Vec_Com C,r,C0;
	C.assign(Num);
	r.assign(Num);
	C0.assign(Num);
	
	for (i=0; i<Num; i++)
	{
		C[i]=Vmn[i]-YM[i];
		r[i]=-YM[i];
		//cout<<YM[i]<<endl;
	}
	
	double Nom=L2_Norm(Vmn);
	double Crit=sqrt(L2_Norm(r)/Nom);
	//cout<<Nom<<endl;
	while(Crit > criterion)
	{
		cout<<"criterion = "<<Crit<<endl;
		outstep<<Crit<<endl;
		Assign_Vecor(C,C0);
		
		//解矩阵方程
		//Two_D_Cond_CGS_Pre(Zs,Zfs,Zs_Diag,C0,nxy,RD_num,eps,Uxy);
		Two_D_Cond_Bi_CGSTAB_Pre(Zs,Zfs,Zs_Diag,C0,Nxy,RD_num,eps,Uxy);
		//Two_D_Cond_CGM(Zs,Zfs,Zs_Inver,Zfs_Inver,Zs_Diag,C0,nxy,RD_num,eps,Uxy);
		
		Two_D_Conductor_Scalar_Multi_Zw_Vec(Para_Num,x,y,z,Uxy,k_w,dxy,Nxy,RD_num,YM);
		
		for (i=0; i<Num; i++)
		{
			C[i]=Vmn[i]-YM[i];
			r[i]=C[i]-C0[i];
			//cout<<YM[i]<<"   "<<i<<endl;
		}
		Crit=sqrt(L2_Norm(r)/Nom);
	}

	for (i=0; i<Num; i++)
	{
		out_Current<<Uxy[i]<<endl;
	}

	Zs.clear();
	Zs_Diag.clear();
//	Zs_Inver.clear();
	Zfs.clear();
//	Zfs_Inver.clear();
	Vmn.clear();
	YM.clear();
	C.clear();
	r.clear();
	C0.clear();

}
Complex matrix::Diel_g_R(double R,Complex k)
{
	Complex re;
	re=exp(i_*k*R)/(4.0*PI*R);
	return re;
}
Complex matrix::G_R(double R,Complex k)
{
	Complex re,t0;
	t0 = i_*k*R;
	re = (t0 - 1.0)*exp(t0)/(4.0*PI*pow(R,3.0));
	return re;
}
Complex matrix::Diel_G_R(double R,Complex k)
{
	Complex re,t0;
	t0 = i_*k*R;
	re = (1. - t0)*exp(t0)/(4.*PI*pow(R,3.));
	return re;
}
Complex matrix::G_R(Complex R,Complex k)
{
	Complex re,t0;
	t0 = i_*k*R;
	re = (t0 - 1.0)*exp(t0)/(4.0*PI*pow(R,3.0));
	return re;
}
void matrix::Two_D_Conductor_Scalar_MoM_Matrix()
{
	int i,j,u;
	int Num = Nxy*Nxy;
	Vec_Com Zmn,Vmn;
	Zmn.assign(Num*Num);
	Vmn.assign(Num);
	double dxy = RL/Nxy;
    double R_S;
	//Z矩阵的填充
	for(i=0; i<Num; i++)
	{
		for(j=0; j<Num; j++)
		{
			u=i*Num+j;
			if( i == j )
			{
				Zmn[u] = Two_D_Conductor_Scalar_Integ_Sing_Left(dxy,x[i],y[i],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w)
					+Two_D_Conductor_Scalar_Integ_Sing(dxy,df_x[i],df_y[i]);
				//outzm<<Zmn[u]<<"  "<<u<<endl;
			}
			else
			{
				R_S = sqrt(pow(x[i]-x[j],2.0)+pow(y[i]-y[j],2.0)+pow(z[i]-z[j],2.0));
				Zmn[u] = dxy*dxy*G_0_fs_R(R_S,k_w);
				//Zmn[u] = Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(dxy,x[i],y[i],x[j],y[j],df_x[i],df_y[i],ddf_x[i],ddf_xy[i],ddf_y[i],k_w);
			}
		}
	}

	//V矩阵的填充
	double xsz,ww,kg,tx,ty;
	for (i=0; i<Num; i++)
	{
		xsz=z[i]*cos(Inci_HL)-x[i]*sin(Inci_HL)*cos(Inci_AZ)-y[i]*sin(Inci_HL)*sin(Inci_AZ);
		kg=pow(g_p*cos(Inci_HL),2);
		tx=pow(x[i]*cos(Inci_HL)*cos(Inci_AZ)+y[i]*cos(Inci_HL)*sin(Inci_AZ)+z[i]*sin(Inci_HL),2)/kg;
		ty=pow((-x[i]*sin(Inci_AZ)+y[i]*cos(Inci_AZ))/g_p,2);        
    	ww=((2.0*tx-1.0)/kg+(2.0*ty-1.0)/pow(g_p,2))/pow(k_w,2);
		Vmn[i] = exp(-i_*k_w*xsz*(1.0+ww))*exp(-(tx+ty));
		//cout<<Vmn[i]<<"  "<<i<<endl;
	}

	/////////////////////////////////////////////////////////////////////
	//解方程组
	time_t ltime1,ltime2,ltime;//存储时间
	time( &ltime1 );
	cout<<"The time of starting resolving equation is "<<ctime( &ltime1 )<<endl;
	resolve_m(Zmn,Vmn,Uxy);
	time( &ltime2 );
	cout<<"The time of stoping resolving equation is "<<ctime( &ltime2 )<<endl;
	ltime=ltime2-ltime1;
	cout<<ltime<<" seconds"<<endl;
	////////////////////////////////////////////////////////////////////
/*	for (i=0; i<Num; i++)
	{
		cout<<Uxy[i]<<endl;
	}*/

	Zmn.clear();
	Vmn.clear();
	
}
void matrix::One_D_Medium_PBTG_BMIA_CAG_Matrix()
{
	int i,j,u;
	k_w1=2.0*PI*fre*sqrt(((4.0*PI)*1.0e-7)*(1.0/(36.0*PI*1.0e9))*permittivity);//区域1的介电常数
	////////////////////////////////
	Ndg=N;//细格的剖分数
	double dx=RL/N;
	////////////////////介质的矩阵
	double d_rl;
	double Cons=1.5;//常数Cons
	d_rl=Cons/imag(k_w1);
	d_rl=d_rl/dx;
	int rl_bw=int(d_rl);//介质空间的子带状矩阵的半带宽
	int rl_mbw=2*rl_bw+1;//介质空间的子带状矩阵的带宽
	///////////////////////
	double d_rf=1.0*wave_len;
	d_rf=d_rf/dx;
	int rf_bw=int(d_rf);//空气空间的子带状矩阵的半带宽
	int rf_mbw=2*rf_bw+1;//空气空间的子带状矩阵的带宽
	//////////////////////
	if(rd_bw <= rf_bw)
	{
		cout<<"rd_bw must be larger than rf_bw;"<<endl;
		exit(0);//程序终止
		//break;//循环终止
	}
	
	//////////////////////
	Complex *a1mn=new Complex[Ndg*rl_mbw];
	Complex *b1mn=new Complex[Ndg*rl_mbw];
	Complex *asmn=new Complex[Ndg*rf_mbw];
	Complex *bsmn=new Complex[Ndg*rf_mbw];
	Complex tm5;
	double tm6,tm7;
	double p_e=2.7183;//自然对数的底数e
	//采用压缩矩阵存储带状矩阵:强矩阵
	//a1mn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<rl_mbw;j++)
		{
			u=i*rl_mbw+j;
			if(i <= rl_bw)
			{
				if(abs(i-j) <= rl_bw)
				{
					if(i == j)
					{
						tm5=k_w1*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
						tm5=Hankel(0,1,tm5);
						a1mn[u]=-dx*i_*tm5/4.0;
			        }
		        	else
					{
						tm5=k_w1*sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
						tm5=Hankel(0,1,tm5);
						a1mn[u]=-dx*i_*tm5/4.0;
		        	}
				}
			}
			else
			{
				if(abs(i-(j+i-rl_bw)) <= rl_bw && (j+i-rl_bw) < Ndg)
				{
					if(i == (j+i-rl_bw))
					{
						tm5=k_w1*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
						tm5=Hankel(0,1,tm5);
						a1mn[u]=-dx*i_*tm5/4.0;
					}
		        	else
					{
						tm5=k_w1*sqrt(pow(x[i]-x[j+i-rl_bw],2.0)+pow(z[i]-z[j+i-rl_bw],2.0));
						tm5=Hankel(0,1,tm5);
		        		a1mn[u]=-dx*i_*tm5/4.0;
					}
				}
			}
			a1mn[u]=a1mn[u]*r_em;
			//cout<<a1mn[u]<<"     ";
		}
		//cout<<endl;
	}
	/////////////////////////////
	//b1mn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<rl_mbw;j++)
		{
			u=i*rl_mbw+j;
			if(i <= rl_bw)
			{
				if(abs(i-j) <= rl_bw)
				{
					if(i == j)
					{
						b1mn[u]=0.5+ddf_x[i]*dx/(4*PI*(1+pow(df_x[i],2)));
					}
		        	else
					{
						tm6=sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
						tm7=df_x[j]*(x[j]-x[i])-(z[j]-z[i]);
						tm5=k_w1*tm6;
						tm5=Hankel(1,1,tm5);
						tm5=k_w1*tm5*tm7/tm6;
						b1mn[u]=dx*i_*tm5/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rl_bw)) <= rl_bw && (j+i-rl_bw) < Ndg)
				{
					if(i == (j+i-rl_bw))
					{
						b1mn[u]=0.5+ddf_x[i]*dx/(4*PI*(1+pow(df_x[i],2)));
					}
		        	else
					{
						tm6=sqrt(pow(x[i]-x[j+i-rl_bw],2)+pow(z[i]-z[j+i-rl_bw],2));
						tm7=df_x[j+i-rl_bw]*(x[j+i-rl_bw]-x[i])-(z[j+i-rl_bw]-z[i]);
						tm5=k_w1*tm6;
						tm5=Hankel(1,1,tm5);
						tm5=k_w1*tm5*tm7/tm6;
		        		b1mn[u]=dx*i_*tm5/4.0;
					}
				}
			}
			//cout<<b1mn[u]<<"     ";
		}
		//cout<<endl;
	}
	///////////空气的带状子矩阵
	//asmn
	
	
	//采用压缩矩阵存储带状矩阵:强矩阵
	//asmn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<rf_mbw;j++)
		{
			u=i*rf_mbw+j;
			if(i <= rf_bw)
			{
				if(abs(i-j) <= rf_bw)
				{
					if(i == j)
					{
						tm6=k_w*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
						asmn[u]=dx*i_*Hankel(0,1,tm6)/4.0;
					}
		        	else
					{
						tm6=k_w*sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
		        		asmn[u]=dx*i_*Hankel(0,1,tm6)/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rf_bw)) <= rf_bw && (j+i-rf_bw) < Ndg)
				{
					if(i == (j+i-rf_bw))
					{
						tm6=k_w*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
			         	asmn[u]=dx*i_*Hankel(0,1,tm6)/4.0;
					}
		        	else
					{
						tm6=k_w*sqrt(pow(x[i]-x[j+i-rf_bw],2.0)+pow(z[i]-z[j+i-rf_bw],2.0));
		        		asmn[u]=dx*i_*Hankel(0,1,tm6)/4.0;
					}
				}
			}
			//cout<<asmn[u]<<"     ";
		}
		//cout<<endl;
	}
	/////////////////////////////
	//bsmn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<rf_mbw;j++)
		{
			u=i*rf_mbw+j;
			if(i <= rf_bw)
			{
				if(abs(i-j) <= rf_bw)
				{
					if(i == j)
					{
						bsmn[u]=0.5-ddf_x[i]*dx/(4.0*PI*(1+pow(df_x[i],2.0)));
					}
		        	else
					{
						tm6=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
						tm7=df_x[j]*(x[j]-x[i])-(z[j]-z[i]);
						tm7=k_w*tm7/tm6;
						bsmn[u]=-dx*i_*tm7*Hankel(1,1,k_w*tm6)/4.0;
					}
				}
			}
			else
			{
				if(abs(i-(j+i-rf_bw)) <= rf_bw && (j+i-rf_bw) < Ndg)
				{
					if(i == (j+i-rf_bw))
					{
						bsmn[u]=0.5-ddf_x[i]*dx/(4.0*PI*(1+pow(df_x[i],2.0)));
					}
		        	else
					{
						tm6=sqrt(pow(x[i]-x[j+i-rf_bw],2.0)+pow(z[i]-z[j+i-rf_bw],2.0));
						tm7=df_x[j+i-rf_bw]*(x[j+i-rf_bw]-x[i])-(z[j+i-rf_bw]-z[i]);
						tm7=k_w*tm7/tm6;
		        		bsmn[u]=-dx*i_*tm7*Hankel(1,1,k_w*tm6)/4.0;
					}
				}
			}
			//cout<<bsmn[u]<<"     ";
		}
		//cout<<endl;
	}

	/////////////////
/*	//验证带状矩阵的解
	Complex *ta1mn=new Complex[Ndg*Ndg];
	Complex *tb1mn=new Complex[Ndg*Ndg];
	//采用压缩矩阵存储带状矩阵:强矩阵
	//a1mn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<Ndg;j++)
		{
			u=i*Ndg+j;
			if(abs(i-j) <= rl_bw)
			{
				if(i == j)
				{
					tm5=k_w1*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
					tm5=Hankel(0,1,tm5);
					ta1mn[u]=-dx*i_*tm5/4.0;
				}
		        else
				{
					tm5=k_w1*sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
					tm5=Hankel(0,1,tm5);
					ta1mn[u]=-dx*i_*tm5/4.0;
				}
			}
			else
			{
				ta1mn[u]=0.0;
			}
		}
		ta1mn[u]=ta1mn[u]*r_em;
	}
	/////////////////////////////
	//b1mn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<Ndg;j++)
		{
			u=i*Ndg+j;
			if(abs(i-j) <= rl_bw)
			{
				if(i == j)
				{
					tb1mn[u]=0.5+ddf_x[i]*dx/(4*PI*(1+pow(df_x[i],2)));
				}
		       	else
				{
					tm6=sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
					tm7=df_x[j]*(x[j]-x[i])-(z[j]-z[i]);
					tm5=k_w1*tm6;
					tm5=Hankel(1,1,tm5);
					tm5=k_w1*tm5*tm7/tm6;
					tb1mn[u]=dx*J*tm5/4.0;
				}
			}
			else
			{
				tb1mn[u]=0.0;
			}
		}
	}
	///////////空气的带状子矩阵
	//asmn
	Complex *tasmn=new Complex[Ndg*Ndg];
	Complex *tbsmn=new Complex[Ndg*Ndg];
	//采用压缩矩阵存储带状矩阵:强矩阵
	//asmn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<Ndg;j++)
		{
			u=i*Ndg+j;
			if(abs(i-j) <= rf_bw)
			{
				if(i == j)
				{
					tm6=k_w*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
					tasmn[u]=dx*J*Hankel(0,1,tm6)/4.0;
				}
		       	else
				{
					tm6=k_w*sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
		        	tasmn[u]=dx*J*Hankel(0,1,tm6)/4.0;
				}
			}
			else
			{
				tasmn[u]=0.0;
			}
		}
	}
	/////////////////////////////
	//bsmn
	for(i=0;i<Ndg;i++)
	{
		for(j=0;j<Ndg;j++)
		{
			u=i*Ndg+j;
			if(abs(i-j) <= rf_bw)
			{
				if(i == j)
				{
					tbsmn[u]=0.5-ddf_x[i]*dx/(4.0*PI*(1+pow(df_x[i],2.0)));
				}
		       	else
				{
					tm6=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
					tm7=df_x[j]*(x[j]-x[i])-(z[j]-z[i]);
					tm7=k_w*tm7/tm6;
					tbsmn[u]=-dx*J*tm7*Hankel(1,1,k_w*tm6)/4.0;
				}
			}
			else
			{
				tbsmn[u]=0.0;
			}
		}
	}
	Complex *allzm=new Complex[2*Ndg*2*Ndg];
	for(i=0;i<2*Ndg;i++)
	{
		for(j=0;j<2*Ndg;j++)
		{
			u=i*2*Ndg+j;
			if(i<Ndg && j<Ndg)
			{
				v=i*Ndg+j;
				allzm[u]=tasmn[v];
			}
			else if(i<Ndg && j>=Ndg)
			{
				v=i*Ndg+j-Ndg;
				allzm[u]=tbsmn[v];
			}
			else if(i>=Ndg && j<Ndg)
			{
				v=(i-Ndg)*Ndg+j;
				allzm[u]=ta1mn[v];
			}
			else
			{
				v=(i-Ndg)*Ndg+j-Ndg;
				allzm[u]=tb1mn[v];
			}
			outzm_R<<real(allzm[u])<<endl;
			outzm_I<<imag(allzm[u])<<endl;
		}
		//outzm<<";"<<endl;
	}*/
	/////////////////////////////
	//v矩阵的填充
	Complex *vm=new Complex[Ndg];
	double xsz,ww,kg,fac;
    for(i=0;i<Ndg;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		vm[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		//outvm<<vm[i]<<endl;
	}
	//L2_Norm(vm,Ndg)=93.9267;
	/////////////////////////////
	//采用CGM计算解矩阵，中间利用BMIA/CAG、插值等
	//粗格上的x,z坐标：在这里，n1只能是奇数，粗格上的点都在原细格上
	Ncg=Ndg/n1;
	double Cdx=n1*dx;
	Vec_Dou Cx,Cz,Cdf_x;
	Cx.assign(Ncg);
	Cz.assign(Ncg);
	Cdf_x.assign(Ncg);
	/*for(i=-Ncg/2+1;i<=Ncg/2;i++)
	{
		Cx[i+Ncg/2-1]=i*rl/Ncg;
	}
	double esl=1.0e-6;
	Cz=atk_Interpolate(z,x,Cx,Ndg,Ncg,esl);
	Cdf_x=atk_Interpolate(df_x,x,Cx,Ndg,Ncg,esl);*/
	
	
	Cx[0]=x[(n1-1)/2];
	Cz[0]=z[(n1-1)/2];
	Cdf_x[0]=df_x[(n1-1)/2];
	for(i=1;i<Ncg;i++)
	{
		Cx[i]=x[i*n1+(n1-1)/2];
		Cz[i]=z[i*n1+(n1-1)/2];
		Cdf_x[i]=df_x[i*n1+(n1-1)/2];
	}
	//for(i=0;i<Ncg;i++)
	//	outm<<Cx[i]<<"     "<<Cz[i]<<endl;
	//当r大于rf和小于rd时的矩阵元素
	//aabmn
	
	Complex *aabmn=new Complex[Ncg*Ncg];
	Complex *babmn=new Complex[Ncg*Ncg];
	for(i=0;i<Ncg;i++)
	{
		for(j=0;j<Ncg;j++)
		{
			u=i*Ncg+j;
			if(abs(i-j) < rd_bw/n1 && abs(i-j) > rf_bw/n1)
			{
				tm6=k_w*sqrt(pow(Cx[i]-Cx[j],2.0)+pow(Cz[i]-Cz[j],2.0));
				aabmn[u]=Cdx*i_*Hankel(0,1,tm6)/4.0;
			}
		}
	}
	//babmn
	for(i=0;i<Ncg;i++)
	{
		for(j=0;j<Ncg;j++)
		{
			u=i*Ncg+j;
			if(abs(i-j) < rd_bw/n1 && abs(i-j) > rf_bw/n1)
			{
				tm6=sqrt(pow(Cx[i]-Cx[j],2.0)+pow(Cz[i]-Cz[j],2.0));
				tm7=Cdf_x[j]*(Cx[j]-Cx[i])-(Cz[j]-Cz[i]);
		        tm7=k_w*tm7/tm6;
				babmn[u]=-Cdx*i_*tm7*Hankel(1,1,k_w*tm6)/4.0;
				//cout<<babmn[u]<<endl;
			}
		}
	}
	//FFT计算弱矩阵

	//矩阵的共轭转置	
/*	Complex *asmnt,*bsmnt,*a1mnt,*b1mnt;
	Complex *aabmnt=new Complex[Ncg*Ncg];
	Complex *babmnt=new Complex[Ncg*Ncg];
		
	asmnt=Transpose_BandMatrix(asmn,Ndg,rf_bw);
	bsmnt=Transpose_BandMatrix(bsmn,Ndg,rf_bw);
	a1mnt=Transpose_BandMatrix(a1mn,Ndg,rl_bw);
	b1mnt=Transpose_BandMatrix(b1mn,Ndg,rl_bw);
	for(i=0;i<Ncg;i++)
	{
		for(j=0;j<Ncg;j++)
		{
			if(abs(i-j) > rf_bw/n1 && abs(i-j) < rd_bw/n1)
			{
				u=i*Ncg+j;
				v=j*Ncg+i;
				aabmnt[u]=conj(aabmn[v]);
				babmnt[u]=conj(babmn[v]);
				//cout<<babmnt[u]<<"        "<<u<<endl;
			}
		}
	}*/
	////////矩阵转置结束/////////

	cout<<"矩阵填充完毕"<<endl;

	//解矩阵方程
	double eps_p=1.0e-2;//F_GMRES_K
	double eps_m=0.5;//F_GMRES_K

	//GPBI_CG_Pre(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	GPBI_CG(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//BI_CGSTAB_Pre(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//CGS(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//Bi_CG_Pre(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//F_GMRES_K(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw,eps_m,eps_p);
	//rex=GMRES_K(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw,eps_p);
	//CGLS(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//LSQR(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//CG(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//resolve_pbtg_cgm(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,vm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	//resolve_pbtg_HJG_cgm(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,Bm,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
	Complex *re0;
	//将2*Ndg个解由n1平均得到2*Ncg个解
	re0=Average(rex,n1,2*Ncg);
	
	Complex *ym=computer_ym(k_w,Cdx,re0,Cx,x,Cdf_x,Cz,Ncg,Ndg,rd_bw,r_em,n1);
	//for(i=0;i<2*Ndg;i++)
	//	cout<<ym[i]<<"    "<<i<<endl;
	    	
	Complex *C=new Complex[Ndg];
	Complex *r=new Complex[Ndg];
	Complex *C0=new Complex[Ndg];
	
	for(i=0;i<Ndg;i++)
	{
		C[i]=vm[i]-ym[i];
		r[i]=-ym[i];
		//outm<<r[i]<<"    "<<i<<endl;
	}
	
	//out<<L2_Norm(ym,n)<<endl;
	double Nom=L2_Norm(vm,Ndg);
	double te;
	te=sqrt(L2_Norm(r,Ndg)/Nom);
	cout<<te<<endl;
	while(sqrt(L2_Norm(r,Ndg)/Nom) > criterion)
	{
		te=sqrt(L2_Norm(r,Ndg)/Nom);
		cout<<"迭代准则为"<<te<<endl;
		outstep<<te<<endl;
		for(i=0;i<Ndg;i++)
		{
			C0[i]=C[i];
		}
		//解矩阵方程
		//GPBI_CG_Pre(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		GPBI_CG(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//BI_CGSTAB_Pre(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//CGS(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//Bi_CG_Pre(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//F_GMRES_K(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw,eps_m,eps_p);
		//rex=GMRES_K(a1mn,b1mn,asmn,bsmn,aabmn,babmn,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw,eps_p);
		//CGLS(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//LSQR(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//CG(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//resolve_pbtg_cgm(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);//
		//resolve_pbtg_HJG_cgm(a1mn,b1mn,asmn,bsmn,a1mnt,b1mnt,asmnt,bsmnt,aabmn,babmn,aabmnt,babmnt,Cx,x,C0,Ndg,Ncg,rl_bw,rf_bw,rd_bw);
		//将2*Ndg个解由n1平均得到2*Ncg个解
		re0=Average(rex,n1,2*Ncg);
		ym=computer_ym(k_w,Cdx,re0,Cx,x,Cdf_x,Cz,Ncg,Ndg,rd_bw,r_em,n1);
    	
		
		for(i=0;i<Ndg;i++)
		{
			C[i]=vm[i]-ym[i];
			r[i]=C[i]-C0[i];
		}
	}
	//////////////////////////////////////////////////
	Cx.clear();
	Cz.clear();
	Cdf_x.clear();
	delete [] re0;
	delete [] vm;
	delete [] ym;
	delete [] C;
	delete [] r;
	delete [] C0;
	delete [] a1mn;
	delete [] b1mn;
	delete [] asmn;
	delete [] bsmn;
	delete [] aabmn;
	delete [] babmn;

}
Complex matrix::Two_D_Conductor_Scalar_Integ_Sing_Left(double dxy,double xm,double ym,double dfx,double dfy,double ddfxx,double ddfxy,double ddfyy,Complex kk)
{
	int n=2;
    int js[2]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Scalar_Integ_Up_Down(dxy,j-1,xm,ym,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p=Two_D_Conductor_Scalar_Integ_Sing_Func_Left(x,xm,ym,dfx,dfy,ddfxx,ddfxy,ddfyy,kk);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }	
}
double matrix::Two_D_Conductor_Scalar_Integ_Sing(double dxy,double dfxm,double dfym)
{//dxy为x,y边上剖分的间隔
	int n=1;
    int js[1]={4};
    int m,j,k,q,l,*is;
    double y[2];
	double s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new double[n];
    a=new double[2*(n+1)];
    b=new double[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Scalar_Integ_Up_Down(dxy,j-1,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Conductor_Scalar_Integ_Sing_Function(x,dxy,dfxm,dfym);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}
Complex matrix::Two_D_Conductor_Scalar_Integ_Matrix_Element_Zs(double dxy,double xm,double ym,double xn,double yn,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex kk)
{
	int n=2;
    int js[2]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Scalar_Integ_Up_Down(dxy,j-1,xn,yn,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p=Two_D_Conductor_Scalar_Func_G_0_Zs(x,xm,ym,dfxm,dfym,ddfxxm,ddfxym,ddfyym,kk);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }	
}
Complex matrix::Two_D_Conductor_Scalar_Integ_Matrix_Element_Zfs(double dxy,double xm,double ym,double xn,double yn,Complex kk)
{
	int n=2;
    int js[2]={4,4};
    int m,j,k,q,l,*is;
    double y[2];
	Complex s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new Complex[n];
    a=new Complex[2*(n+1)];
    b=new Complex[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Scalar_Integ_Up_Down(dxy,j-1,xn,yn,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p=Two_D_Conductor_Scalar_Func_G_0_Zfs(x,xm,ym,kk);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }	
}
void matrix::Two_D_Conductor_Scalar_Integ_Up_Down(double dxy,int j,double y[2])
{
	switch(j)
	{
	case 0:{y[0]=-dxy/2.0;y[1]=dxy/2.0;break;}
    case 1:{y[0]=-dxy/2.0;y[1]=dxy/2.0;break;}
    case 2:{y[0]=-1.0;y[1]=1.0;break;}
    case 3:{y[0]=-1.0;y[1]=1.0;break;}
    default:{}
	}
}
void matrix::Two_D_Conductor_Scalar_Integ_Up_Down(double dxy,int j,double xn,double yn,double y[2])
{
	switch(j)
	{
	case 0:{y[0]=xn-dxy/2.0;y[1]=xn+dxy/2.0;break;}
    case 1:{y[0]=yn-dxy/2.0;y[1]=yn+dxy/2.0;break;}
    case 2:{y[0]=-1.0;y[1]=1.0;break;}
    case 3:{y[0]=-1.0;y[1]=1.0;break;}
    default:{}
	}
}
Complex matrix::Two_D_Conductor_Scalar_Integ_Sing_Func_Left(Complex x[2],double xm,double ym,double dfx,double dfy,double ddfxx,double ddfxy,double ddfyy,Complex k)
{
	Complex z;
	Complex r,temp;
	temp = dfx*(x[0]-xm) + dfy*(x[1]-ym) + ddfxx*(x[0]-xm)*(x[0]-xm)/2.0 + ddfxy*(x[0]-xm)*(x[1]-ym) + ddfyy*(x[1]-ym)*(x[1]-ym)/2.0;
	r=sqrt(pow(x[0]-xm,2) + pow(x[1]-ym,2) + pow(temp,2)) + 1.0e-18;
	z=(exp(i_*k*r)-1.0)/(4.0*PI*r);
	//z = (i_*k - k*k*r/2.0 - i_*k*k*k*r*r/factor(3)) / (4.0*PI);
	return z;
}
Complex matrix::Two_D_Conductor_Scalar_Func_G_0_Zs(Complex x[2],double xm,double ym,double dfxm,double dfym,double ddfxxm,double ddfxym,double ddfyym,Complex k)
{
	Complex z;
	Complex r,temp;
	temp = dfxm*(x[0]-xm) + dfym*(x[1]-ym) + ddfxxm*(x[0]-xm)*(x[0]-xm)/2.0+ddfxym*(x[0]-xm)*(x[1]-ym)+ddfyym*(x[1]-ym)*(x[1]-ym)/2.0;
	r=sqrt(pow(x[0]-xm,2) + pow(x[1]-ym,2) + pow(temp,2));
	z=exp(i_*k*r)/(4.0*PI*r);
	return z;
}
Complex matrix::Two_D_Conductor_Scalar_Func_G_0_Zfs(Complex x[2],double xm,double ym,Complex k)
{
	Complex z;
	Complex r;
	r=sqrt(pow(x[0]-xm,2) + pow(x[1]-ym,2));
	z=exp(i_*k*r)/(4.0*PI*r);
	return z;
}
double matrix::Two_D_Conductor_Scalar_Integ_Sing_Function(double x[1],double dxy,double dfxm,double dfym)
{
	double z;
	double a,b,c,t1,t2;
	a=1.0+dfxm*dfxm;
	b=2.0*dfxm*dfym*x[0];
	c=(1.0+dfym*dfym)*x[0]*x[0];
	t1=fabs(2.0*a*dxy/2.0+b+2.0*sqrt(a*(a*dxy*dxy/4.0+b*dxy/2.0+c)));
	t2=fabs(-2.0*a*dxy/2.0+b+2.0*sqrt(a*(a*dxy*dxy/4.0-b*dxy/2.0+c)));
	z=log(t1)-log(t2);
	
	return z/(4.0*PI*sqrt(a));
	
	/*double re;
	double t1,t2,a2;
	t1 = -dxy/2.0+dfxm*dfym*x[0]/(1.0+dfxm*dfxm);
	t2 =  dxy/2.0+dfxm*dfym*x[0]/(1.0+dfxm*dfxm);
	a2 = ((1.0+dfym*dfym)/(1.0+dfxm*dfxm)-dfxm*dfxm*dfym*dfym/pow(1.0+dfxm*dfxm,2.0))*pow(x[0],2.0);
	re = (log(t2+sqrt(t2*t2+a2))-log(t1+sqrt(t1*t1+a2))) / (4.0*PI);
	return re;*/	
}
double matrix::Two_D_Conductor_Scalar_Integ_Sing_Function(double x[1],double ym,double dxy,double dfxm,double dfym)
{
	double re;
	double A;
	double B,C;
	A = 1. + pow(dfxm,2.);
	B = 2. * dfxm * dfym * (x[0] - ym);
	C = (1. + pow(dfym,2.)) * pow(x[0]-ym,2.);
	double U11,D12;
	U11 = dxy/2.;
	D12 = -dxy/2.;
	double IR1;
	IR1 = (log(fabs(2.*A*U11+B+2.*sqrt(A)*sqrt(A*pow(U11,2.)+B*U11+C)))/sqrt(A)) - (log(fabs(2.*A*D12+B+2.*sqrt(A)*sqrt(A*pow(D12,2.)+B*D12+C)))/sqrt(A));
    re = IR1 / (4.*PI);

	return re;		
}
double matrix::Two_D_Conductor_Scalar_Integ_Sing(double ym,double dxy,double dfxm,double dfym)
{//dxy为x,y边上剖分的间隔
	int n=1;
    int js[1]={4};
    int m,j,k,q,l,*is;
    double y[2];
	double s,*x,*a,*b,p;
    static double t[5]={-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851};
    is=new int[2*(n+1)];
    x=new double[n];
    a=new double[2*(n+1)];
    b=new double[n+1];
    m=1; l=1;
    a[n]=1.0; a[2*n+1]=1.0;
    while (l==1)//非0时候执行以下语句
      { for (j=m;j<=n;j++)
          { Two_D_Conductor_Vector_Integ_Up_Down(dxy,j-1,ym,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];//都是乘以t[0];目的是先要计算出各个积分区间上的第一个子区间中第一个高斯型点;
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
          }
        j=n; q=1;
        while (q==1)
          { k=is[j-1];
            if (j==n) 
				p = Two_D_Conductor_Scalar_Integ_Sing_Function(x,ym,dxy,dfxm,dfym);
            else p=1.0;
             a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];//存放当前次和上一次插值求积公式之和//
            is[j-1]=is[j-1]+1;
            if (is[j-1]>5)
              if (is[n+j]>=js[j-1])
                { j=j-1; q=1;
                  if (j==0)
                    { s=a[n+1]*a[0]; 
				      delete [] is; 
				      delete [] x;
                      delete [] a; 
					  delete [] b; 
					  return s;
                    }
                }
              else
                { is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
                }
            else
              { k=is[j-1];
                x[j-1]=a[j-1]*t[k-1]+b[j-1];
                if (j==n) q=1;
                else q=0;
              }
          }
        m=j+1;
      }
}			
void matrix::One_D_Conductor_MoM_Matrix()
{
	int i,j,u;
	Vec_Com Zmn,Vmn;
	Zmn.assign(N*N);
	Vmn.assign(N);
	double dx = RL/N;
	double euler=1.78107;
	double dl;
	double abx;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			u=i*N+j;
			if(i==j)
			{
				dl=dx*sqrt(1.0+pow(df_x[i],2));
			    Zmn[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
			}
			else
			{
				abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        Zmn[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
			}
		}
	}
	//v矩阵的填充
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Vmn[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}
	resolve_m(Zmn,Vmn,N);
	for (i=0;i<N;i++)
	{
		cout<<rex[i]<<endl;
	}
	
	Zmn.clear();
	Vmn.clear();
}
void matrix::One_D_Conductor_BMIA_CAG_Matrix()
{
	int i,j,u;
	int mbw=2*rd_bw+1;
	Complex *zm,*vm;
	zm=new Complex[N*mbw];
	vm=new Complex[N];
	k_w=2*PI/wave_len;
	double dx=RL/N;
	double euler=1.78107;
	double dl;
	double abx;
	//采用压缩矩阵存储带状矩阵:强矩阵
	for(i=0;i<N;i++)
	{
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= rd_bw)
			{
				if(abs(i-j) <= rd_bw)
				{
					if(i == j)
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j],2)+pow(z[i]-z[j],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
				/*else
				{
					zm[u]=0.0;
				}*/
			}
			else
			{
				if(abs(i-(j+i-rd_bw)) <= rd_bw && (j+i-rd_bw) < N)
				{
					if(i == (j+i-rd_bw))
					{
						dl=dx*sqrt(1.0+pow(df_x[i],2));
			         	zm[u]=-dx*(log(euler*k_w*dl/4.0)-1.0)/(2.0*PI)+i_*dx/4.0;
					}
		        	else
					{
		        		abx=k_w*sqrt(pow(x[i]-x[j+i-rd_bw],2)+pow(z[i]-z[j+i-rd_bw],2));
		        		zm[u]=-dx*_y0(abx)/4.0+i_*dx*_j0(abx)/4.0;
					}
				}
			}
			//cout<<zm[u]<<"   "<<u<<endl;
		}
	}
	//FFT计算弱矩阵
	
	//v矩阵的填充
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		vm[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		//vm[i].display();
		//outresult<<";"<<endl;
		//outvm<<vm[i]<<"  "<<i<<endl;
	}
	
	//迭代法解矩阵方程
	//resolve_m(zm,n,vm);
	//resolve_bm(zm,vm,n,rd_bw);
	resolve_cgm(zm,vm,N,rd_bw);
	//for(i=0;i<n;i++)
	//	outzm<<rex[i]<<endl;
	/////////////////////////////////////////////
	Complex *ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
    //Complex *ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);

	Complex *C,*r,*C0;
	C=new Complex[N];
	r=new Complex[N];
	C0=new Complex[N];
	
	for(i=0;i<N;i++)
	{
		C[i]=vm[i]-ym[i];
		r[i]=-ym[i];
		
	}
	//cout<<Norm2(ym,n)<<endl;
	double Nom=L2_Norm(vm,N);
	double te;
	//te=sqrt(Norm2(r,N)/Nom);
	//cout<<te<<endl;
	while(sqrt(L2_Norm(r,N)/Nom) > 0.01)
	{
		te=sqrt(L2_Norm(r,N)/Nom);
		cout<<"BMIA迭代标准＝"<<te<<endl;
		for(i=0;i<N;i++)
		{
			C0[i]=C[i];
			//C0[i].display();
		}
		//C0=C;//??
		//resolve_bm(zm,C,N,rd_bw);
		resolve_cgm(zm,C,N,rd_bw);
		//ym=computer_ym(temp0,temp1,temp2,temp3,rex,z,N);
		ym=computer_ym_te(k_w,dx,rex,z,N,rd_bw,1);
		for(i=0;i<N;i++)
		{
			C[i]=vm[i]-ym[i];
			r[i]=C[i]-C0[i];
			//C[i].display();
		}
		te=sqrt(L2_Norm(r,N)/Nom);
	}
	//////////////////////////////////////////////////
	delete [] zm;
	delete [] vm;
	delete [] ym;
	for(i=0;i<N;i++)
	{
		rex[i]=rex[i]*i_/(k_w*120.*PI)/sqrt(1+pow(df_x[i],2));
		outzm<<x[i]/wave_len<<" "<<real(rex[i])<<" "<<imag(rex[i])<<"  "<<abs(rex[i])<<endl;
	}

}
void matrix::One_D_Medium_MoM_Matrix()
{
	int i,j,u;
	k_w1=2.0*PI*fre*sqrt(((4.0*PI)*1.0e-7)*(1.0/(36.0*PI*1.0e9))*permittivity);//区域1的介电常数
	//cout<<k_w1<<endl;
	Vec_Com Zmn,Vmn;
	Zmn.assign(2*N*2*N);
	Vmn.assign(2*N);
	//Complex *Zmn=new Complex[2*N*2*N];
	//Complex *Vmn=new Complex[2*N];
	double dx=RL/N;
	double p_e=2.7183;
	Complex t3;
	double t1,t2;
	
	//填充阻抗矩阵
	for(i=0;i<2*N;i++)
	{
		for(j=0;j<2*N;j++)
		{
			u=i*2*N+j;
			if(i < N && j < N)
			{
				if(i == j)
				{
			    	t1=k_w*dx*sqrt(1.0+pow(df_x[i],2.0))/(2.0*p_e);
					Zmn[u]=dx*i_*Hankel(0,1,t1)/4.0;
				}
		    	else
				{
			    	t1=k_w*sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
		        	Zmn[u]=dx*i_*Hankel(0,1,t1)/4.0;
				}
			}
			else if(i < N && j >= N)
			{
				if(i == (j-N))
				{
    				Zmn[u]=0.5-ddf_x[i]*dx/(4.0*PI*(1+pow(df_x[i],2.0)));
				}
    			else
				{
     			    t1=sqrt(pow(x[i]-x[j-N],2.0)+pow(z[i]-z[j-N],2.0));
					t2=df_x[j-N]*(x[j-N]-x[i])-(z[j-N]-z[i]);
					t2=k_w*t2/t1;
					Zmn[u]=-dx*i_*t2*Hankel(1,1,k_w*t1)/4.0;
				}
			}
			else if(i >= N && j < N)
			{
				if((i-N) == j)
				{
	    			t3=k_w1*dx*sqrt(1.0+pow(df_x[i-N],2.0))/(2.0*p_e);
					t3=Hankel(0,1,t3);
					Zmn[u]=-dx*i_*t3/4.0;
					Zmn[u]=Zmn[u]*r_em;
				}
	    		else
				{
					t3=k_w1*sqrt(pow(x[i-N]-x[j],2.0)+pow(z[i-N]-z[j],2.0));
					t3=Hankel(0,1,t3);
					Zmn[u]=-dx*i_*t3/4.0;
					Zmn[u]=Zmn[u]*r_em;
				}
				
			}
			else
			{
				if((i-N) == (j-N))
				{
		    		Zmn[u]=0.5+ddf_x[i-N]*dx/(4*PI*(1+pow(df_x[i-N],2)));
				}
	    		else
				{
		    		t1=sqrt(pow(x[i-N]-x[j-N],2)+pow(z[i-N]-z[j-N],2));
					t2=df_x[j-N]*(x[j-N]-x[i-N])-(z[j-N]-z[i-N]);
					t3=k_w1*t1;
					t3=Hankel(1,1,t3);
					t3=k_w1*t3*t2/t1;
					Zmn[u]=dx*i_*t3/4.0;
				}
			}
			//cout<<Zmn[u]<<"      "<<u<<endl;
		}
	}
	///////////////////////////////////////////////////////////////////
	//v矩阵的填充
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Vmn[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}
	
	/////////////////////////////////////////////////////////////////////
	//解方程组
	cout<<"开始解方程"<<endl;
	resolve_m(Zmn,Vmn,2*N);
	Zmn.clear();
	Vmn.clear();
}

//自由空间左边弱子矩阵的泰勒展开的系数,及其和列向量的乘积
Complex matrix::para_a_te(int Bool,Complex k,double x)
{
	Complex re;
	Complex t0,t1;
	Complex t2=0.25*i_;
	switch(Bool)//系数a0,a1,a2,a3
	{
	case 0:
		{
			re=Hankel(0,1,k*x);
			break;
		}
	case 1: 
		{
			t1=Hankel(1,1,k*x);
			re=-k*x*t1/2.0;
			break;
		}
	case 2: 
		{
			t0=Hankel(0,1,k*x);
			t1=Hankel(1,1,k*x);
			re=-pow(k*x,2.0)*t0/8.0+k*x*t1/4.0;
			break;
		}
	case 3: 
		{
			t0=Hankel(0,1,k*x);
			t1=Hankel(1,1,k*x);
			re=-pow(k*x,2.0)*t0/12.0+pow(k*x,3.0)*t1/48.0-k*x*t1/6.0;
			break;
		}
	}
	return re*t2;
}

Complex *matrix::computer_ym_te(Complex k_w,double dx,Complex *x,Vec_Dou &fmn,int nn,int bw,int n1)
{
	int i;
	Complex *ym,*fu1,*fu2,*fu3,*fu4,*fu5,*fu6;
	Complex *temp0=new Complex[nn];
	Complex *temp1=new Complex[nn];
	Complex *temp2=new Complex[nn];
	Complex *temp3=new Complex[nn];
	double temp;
	for(i=0;i<nn;i++)
	{
		if(i > bw/n1)
		{
			temp=i*dx;
			temp0[i]=para_a_te(0,k_w,temp);
			temp1[i]=para_a_te(1,k_w,temp)/pow(temp,2.0);
			temp2[i]=para_a_te(2,k_w,temp)/pow(temp,4.0);
			temp3[i]=para_a_te(3,k_w,temp)/pow(temp,6.0);
		}
	}
	fu1=new Complex[nn];
	fu2=new Complex[nn];
	fu3=new Complex[nn];
	fu4=new Complex[nn];
	fu5=new Complex[nn];
	fu6=new Complex[nn];
	 ym=new Complex[nn];
	
	Complex *t0,*t1,*t2,*t3,*t4,*t5,*t6,*t7,*t8,*t9,*t10,*t11,*t12,*t13,*t14,*t15;
	
	for(i=0;i<nn;i++)
	{
		fu1[i]=fmn[i]*x[i];
		fu2[i]=fu1[i]*fmn[i];
		fu3[i]=fu2[i]*fmn[i];
		fu4[i]=fu3[i]*fmn[i];
		fu5[i]=fu4[i]*fmn[i];
		fu6[i]=fu5[i]*fmn[i];
	}
	t0=fft_toematrix_vector(temp0,x,nn);//0
	t1=fft_toematrix_vector(temp1,x,nn);//1
	t2=fft_toematrix_vector(temp1,fu1,nn);//2
	t3=fft_toematrix_vector(temp1,fu2,nn);//3
	t4=fft_toematrix_vector(temp2,x,nn);//4
	t5=fft_toematrix_vector(temp2,fu1,nn);//5
	t6=fft_toematrix_vector(temp2,fu2,nn);//6
	t7=fft_toematrix_vector(temp2,fu3,nn);//7
	t8=fft_toematrix_vector(temp2,fu4,nn);//8
	t9=fft_toematrix_vector(temp3,x,nn);//9
	t10=fft_toematrix_vector(temp3,fu1,nn);//10
	t11=fft_toematrix_vector(temp3,fu2,nn);//11
	t12=fft_toematrix_vector(temp3,fu3,nn);//12
	t13=fft_toematrix_vector(temp3,fu4,nn);//13
	t14=fft_toematrix_vector(temp3,fu5,nn);//14
	t15=fft_toematrix_vector(temp3,fu6,nn);//15

	for(i=0;i<nn;i++)
	{
		ym[i]=t0[i]+pow(fmn[i],2.0)*t1[i]-2.0*fmn[i]*t2[i]+t3[i]+pow(fmn[i],4.0)*t4[i]
			-4.0*pow(fmn[i],3.0)*t5[i]+6.0*pow(fmn[i],2.0)*t6[i]-4.0*fmn[i]*t7[i]+t8[i]
			+pow(fmn[i],6.0)*t9[i]-6.0*pow(fmn[i],5.0)*t10[i]+15.0*pow(fmn[i],4.0)*t11[i]
			-20.0*pow(fmn[i],3.0)*t12[i]+15.0*pow(fmn[i],2.0)*t13[i]-6.0*fmn[i]*t14[i]+t15[i];
		ym[i] = dx * ym[i];
		//ym[i].display();
	}

	delete [] temp0;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] fu1;
	delete [] fu2;
	delete [] fu3;
	delete [] fu4;
	delete [] fu5;
	delete [] fu6;
	delete [] t1;
	delete [] t2;
	delete [] t3;
	delete [] t4;
	delete [] t5;
	delete [] t6;
	delete [] t7;
	delete [] t8;
	delete [] t9;
	delete [] t10;
	delete [] t11;
	delete [] t12;
	delete [] t13;
	delete [] t14;
	delete [] t15;

	return ym;
}
//自由空间右边弱子矩阵的泰勒展开的系数,及其和列向量的乘积
Complex matrix::para_a_tm(int Bool,Complex k,double x)
{
	Complex re;
	Complex t0,t1;
	Complex t2=-i_*k/4.0;
	switch(Bool)//系数a0,a1,a2,a3
	{
	case 0:
		{
			t1=Hankel(1,1,k*x);
			re=t1/x;
			break;
		}
	case 1: 
		{
			t0=Hankel(0,1,k*x);
			t1=Hankel(1,1,k*x);
			re=(0.5*k*x*t0-t1)/pow(x,3.0);
			break;
		}
	case 2: 
		{
			t0=Hankel(0,1,k*x);
			t1=Hankel(1,1,k*x);
			re=(8.0*t1-4.0*k*x*t0-k*k*x*x*t1)/(8.0*pow(x,5.0));
			break;
		}
	case 3: 
		{
			t0=Hankel(0,1,k*x);
			t1=Hankel(1,1,k*x);
			re=(-48.0*t1+24.0*k*x*t0+8.0*k*k*x*x*t1-pow(k*x,3.0)*t0)/(48.0*pow(x,7.0));
			break;
		}
	}
	return t2*re;
}

Complex *matrix::computer_ym_tm(Complex k_w,double dx,Complex *rx,Vec_Dou &xmn,Vec_Dou &dfx,Vec_Dou &fmn,int nn,int bw,int n1)
{
	int i;
	Complex *ym,*ym1,*ym2,*ym3,*ym4,*fu0,*fu1,*fu2,*fu3,*fu4,*fu5,*fu6;
	Complex *temp0=new Complex[nn];
	Complex *temp1=new Complex[nn];
	Complex *temp2=new Complex[nn];
	Complex *temp3=new Complex[nn];
	double temp;
	for(i=0;i<nn;i++)
	{
		if(i > bw/n1)
		{
			temp=i*dx;
			temp0[i]=para_a_tm(0,k_w,temp);
			temp1[i]=para_a_tm(1,k_w,temp);
			temp2[i]=para_a_tm(2,k_w,temp);
			temp3[i]=para_a_tm(3,k_w,temp);
			//cout<<temp0[i]<<endl;
		}
	}
	fu0=new Complex[nn];
	fu1=new Complex[nn];
	fu2=new Complex[nn];
	fu3=new Complex[nn];
	fu4=new Complex[nn];
	fu5=new Complex[nn];
	fu6=new Complex[nn];
	 ym=new Complex[nn];
	ym1=new Complex[nn];
	ym2=new Complex[nn];
	ym3=new Complex[nn];
	ym4=new Complex[nn];
	
	Complex *t0,*t1,*t2,*t3,*t4,*t5,*t6,*t7,*t8,*t9,*t10,*t11,*t12,*t13,*t14,*t15;
	
	//ym4问题在这
	for(i=0;i<nn;i++)
	{
		fu0[i]=rx[i];
		fu1[i]=fu0[i]*fmn[i];
		fu2[i]=fu1[i]*fmn[i];
		fu3[i]=fu2[i]*fmn[i];
		fu4[i]=fu3[i]*fmn[i];
		fu5[i]=fu4[i]*fmn[i];
		fu6[i]=fu5[i]*fmn[i];
	}
	t0=fft_toematrix_vector(temp0,fu0,nn);//0
	t1=fft_toematrix_vector(temp1,fu0,nn);//1
	t2=fft_toematrix_vector(temp1,fu1,nn);//2
	t3=fft_toematrix_vector(temp1,fu2,nn);//3
	t4=fft_toematrix_vector(temp2,fu0,nn);//4
	t5=fft_toematrix_vector(temp2,fu1,nn);//5
	t6=fft_toematrix_vector(temp2,fu2,nn);//6
	t7=fft_toematrix_vector(temp2,fu3,nn);//7
	t8=fft_toematrix_vector(temp2,fu4,nn);//8
	t9=fft_toematrix_vector(temp3,fu0,nn);//9
	t10=fft_toematrix_vector(temp3,fu1,nn);//10
	t11=fft_toematrix_vector(temp3,fu2,nn);//11
	t12=fft_toematrix_vector(temp3,fu3,nn);//12
	t13=fft_toematrix_vector(temp3,fu4,nn);//13
	t14=fft_toematrix_vector(temp3,fu5,nn);//14
	t15=fft_toematrix_vector(temp3,fu6,nn);//15

	for(i=0;i<nn;i++)//如何处理前面的项????????????
	{
		ym4[i]=t0[i]+pow(fmn[i],2.0)*t1[i]-2.0*fmn[i]*t2[i]+t3[i]+pow(fmn[i],4.0)*t4[i]
			-4.0*pow(fmn[i],3.0)*t5[i]+6.0*pow(fmn[i],2.0)*t6[i]-4.0*fmn[i]*t7[i]+t8[i]
			+pow(fmn[i],6.0)*t9[i]-6.0*pow(fmn[i],5.0)*t10[i]+15.0*pow(fmn[i],4.0)*t11[i]
			-20.0*pow(fmn[i],3.0)*t12[i]+15.0*pow(fmn[i],2.0)*t13[i]-6.0*fmn[i]*t14[i]+t15[i];
		ym4[i]=fmn[i]*ym4[i];
		//cout<<ym4[i]<<endl;
	}
	//ym1
	for(i=0;i<nn;i++)
	{
		fu0[i]=dfx[i]*xmn[i]*rx[i];
		fu1[i]=fu0[i]*fmn[i];
		fu2[i]=fu1[i]*fmn[i];
		fu3[i]=fu2[i]*fmn[i];
		fu4[i]=fu3[i]*fmn[i];
		fu5[i]=fu4[i]*fmn[i];
		fu6[i]=fu5[i]*fmn[i];
	}
	t0=fft_toematrix_vector(temp0,fu0,nn);//0
	t1=fft_toematrix_vector(temp1,fu0,nn);//1
	t2=fft_toematrix_vector(temp1,fu1,nn);//2
	t3=fft_toematrix_vector(temp1,fu2,nn);//3
	t4=fft_toematrix_vector(temp2,fu0,nn);//4
	t5=fft_toematrix_vector(temp2,fu1,nn);//5
	t6=fft_toematrix_vector(temp2,fu2,nn);//6
	t7=fft_toematrix_vector(temp2,fu3,nn);//7
	t8=fft_toematrix_vector(temp2,fu4,nn);//8
	t9=fft_toematrix_vector(temp3,fu0,nn);//9
	t10=fft_toematrix_vector(temp3,fu1,nn);//10
	t11=fft_toematrix_vector(temp3,fu2,nn);//11
	t12=fft_toematrix_vector(temp3,fu3,nn);//12
	t13=fft_toematrix_vector(temp3,fu4,nn);//13
	t14=fft_toematrix_vector(temp3,fu5,nn);//14
	t15=fft_toematrix_vector(temp3,fu6,nn);//15

	for(i=0;i<nn;i++)//如何处理前面的项????????????
	{
		ym1[i]=t0[i]+pow(fmn[i],2.0)*t1[i]-2.0*fmn[i]*t2[i]+t3[i]+pow(fmn[i],4.0)*t4[i]
			-4.0*pow(fmn[i],3.0)*t5[i]+6.0*pow(fmn[i],2.0)*t6[i]-4.0*fmn[i]*t7[i]+t8[i]
			+pow(fmn[i],6.0)*t9[i]-6.0*pow(fmn[i],5.0)*t10[i]+15.0*pow(fmn[i],4.0)*t11[i]
			-20.0*pow(fmn[i],3.0)*t12[i]+15.0*pow(fmn[i],2.0)*t13[i]-6.0*fmn[i]*t14[i]+t15[i];
		//cout<<ym1[i]<<endl;
	}
	//ym2
	for(i=0;i<nn;i++)
	{
		fu0[i]=dfx[i]*rx[i];
		fu1[i]=fu0[i]*fmn[i];
		fu2[i]=fu1[i]*fmn[i];
		fu3[i]=fu2[i]*fmn[i];
		fu4[i]=fu3[i]*fmn[i];
		fu5[i]=fu4[i]*fmn[i];
		fu6[i]=fu5[i]*fmn[i];
	}
	t0=fft_toematrix_vector(temp0,fu0,nn);//0
	t1=fft_toematrix_vector(temp1,fu0,nn);//1
	t2=fft_toematrix_vector(temp1,fu1,nn);//2
	t3=fft_toematrix_vector(temp1,fu2,nn);//3
	t4=fft_toematrix_vector(temp2,fu0,nn);//4
	t5=fft_toematrix_vector(temp2,fu1,nn);//5
	t6=fft_toematrix_vector(temp2,fu2,nn);//6
	t7=fft_toematrix_vector(temp2,fu3,nn);//7
	t8=fft_toematrix_vector(temp2,fu4,nn);//8
	t9=fft_toematrix_vector(temp3,fu0,nn);//9
	t10=fft_toematrix_vector(temp3,fu1,nn);//10
	t11=fft_toematrix_vector(temp3,fu2,nn);//11
	t12=fft_toematrix_vector(temp3,fu3,nn);//12
	t13=fft_toematrix_vector(temp3,fu4,nn);//13
	t14=fft_toematrix_vector(temp3,fu5,nn);//14
	t15=fft_toematrix_vector(temp3,fu6,nn);//15

	for(i=0;i<nn;i++)//如何处理前面的项????????????
	{
		ym2[i]=t0[i]+pow(fmn[i],2.0)*t1[i]-2.0*fmn[i]*t2[i]+t3[i]+pow(fmn[i],4.0)*t4[i]
			-4.0*pow(fmn[i],3.0)*t5[i]+6.0*pow(fmn[i],2.0)*t6[i]-4.0*fmn[i]*t7[i]+t8[i]
			+pow(fmn[i],6.0)*t9[i]-6.0*pow(fmn[i],5.0)*t10[i]+15.0*pow(fmn[i],4.0)*t11[i]
			-20.0*pow(fmn[i],3.0)*t12[i]+15.0*pow(fmn[i],2.0)*t13[i]-6.0*fmn[i]*t14[i]+t15[i];
		ym2[i]=xmn[i]*ym2[i];
		//cout<<ym2[i]<<endl;
	}
	//ym3
	for(i=0;i<nn;i++)
	{
		fu0[i]=fmn[i]*rx[i];
		fu1[i]=fu0[i]*fmn[i];
		fu2[i]=fu1[i]*fmn[i];
		fu3[i]=fu2[i]*fmn[i];
		fu4[i]=fu3[i]*fmn[i];
		fu5[i]=fu4[i]*fmn[i];
		fu6[i]=fu5[i]*fmn[i];
	}
	t0=fft_toematrix_vector(temp0,fu0,nn);//0
	t1=fft_toematrix_vector(temp1,fu0,nn);//1
	t2=fft_toematrix_vector(temp1,fu1,nn);//2
	t3=fft_toematrix_vector(temp1,fu2,nn);//3
	t4=fft_toematrix_vector(temp2,fu0,nn);//4
	t5=fft_toematrix_vector(temp2,fu1,nn);//5
	t6=fft_toematrix_vector(temp2,fu2,nn);//6
	t7=fft_toematrix_vector(temp2,fu3,nn);//7
	t8=fft_toematrix_vector(temp2,fu4,nn);//8
	t9=fft_toematrix_vector(temp3,fu0,nn);//9
	t10=fft_toematrix_vector(temp3,fu1,nn);//10
	t11=fft_toematrix_vector(temp3,fu2,nn);//11
	t12=fft_toematrix_vector(temp3,fu3,nn);//12
	t13=fft_toematrix_vector(temp3,fu4,nn);//13
	t14=fft_toematrix_vector(temp3,fu5,nn);//14
	t15=fft_toematrix_vector(temp3,fu6,nn);//15

	for(i=0;i<nn;i++)//如何处理前面的项????????????
	{
		ym3[i]=t0[i]+pow(fmn[i],2.0)*t1[i]-2.0*fmn[i]*t2[i]+t3[i]+pow(fmn[i],4.0)*t4[i]
			-4.0*pow(fmn[i],3.0)*t5[i]+6.0*pow(fmn[i],2.0)*t6[i]-4.0*fmn[i]*t7[i]+t8[i]
			+pow(fmn[i],6.0)*t9[i]-6.0*pow(fmn[i],5.0)*t10[i]+15.0*pow(fmn[i],4.0)*t11[i]
			-20.0*pow(fmn[i],3.0)*t12[i]+15.0*pow(fmn[i],2.0)*t13[i]-6.0*fmn[i]*t14[i]+t15[i];
	}

	for(i=0;i<nn;i++)
	{
		ym[i]=ym1[i]-ym2[i]-ym3[i]+ym4[i];
		ym[i] = dx * ym[i];
	}

	delete [] temp0;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] fu0;
	delete [] fu1;
	delete [] fu2;
	delete [] fu3;
	delete [] fu4;
	delete [] fu5;
	delete [] fu6;
	delete t0;
	delete t1;
	delete t2;
	delete t3;
	delete t4;
	delete t5;
	delete t6;
	delete t7;
	delete t8;
	delete t9;
	delete t10;
	delete t11;
	delete t12;
	delete t13;
	delete t14;
	delete t15;
	delete [] ym1;
	delete [] ym2;
	delete [] ym3;
	delete [] ym4;

	return ym;
}
Complex *matrix::computer_ym(Complex k_w,double dx,Complex *rx,Vec_Dou &Cx,Vec_Dou &x,Vec_Dou &dfx,vector <double> &fmn,int ncg,int ndg,int rd_bw,Complex r_em,int n1)
{//ym输出为2*ndg维
	int i;
	Complex *re1=new Complex[ncg];
	Complex *re2=new Complex[ncg];
	for(i=0;i<ncg;i++)
	{
	   	re1[i]=rx[i];
	   	re2[i]=rx[i+ncg];
	}
	Complex *yml_Up=computer_ym_te(k_w,dx,re1,fmn,ncg,rd_bw,n1);
	Complex *ymr_Up=computer_ym_tm(k_w,dx,re2,Cx,dfx,fmn,ncg,rd_bw,n1);
	
	Complex *ym_Up=Plus_Matrix(yml_Up,ymr_Up,ncg);//列向量的上半部分
	//ym_Up=Linear_Interpolate(ym_Up,Cx,x,ncg,ndg);
	ym_Up=Lg3_Interpolate(ym_Up,Cx,x,ncg,ndg);
	
	

	delete [] re1;
	delete [] re2;
	delete yml_Up;
	delete ymr_Up;

	return ym_Up;
}