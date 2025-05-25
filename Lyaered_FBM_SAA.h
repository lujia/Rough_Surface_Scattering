

#include "surface.h "
#include "Target.h"

class Lyaered_FBM_SAA:public post_cp, public CylinderQuestion
{
public:
	     Lyaered_FBM_SAA();
		 ~Lyaered_FBM_SAA();

		 void Gauss_Layered_Rough_surface_FBM();
		 void Gauss_Layered_Rough_surface_FBM_SAA();

         
		 void Gauss_Layered_Rough_surface_FBM_Bried_Target();
		 void Gauss_Layered_Rough_surface_FBM_Above_Target();
		 void Gauss_Layered_Rough_surface_FBM_Bried_conduct_Target();



};

Lyaered_FBM_SAA::Lyaered_FBM_SAA()
{

}

Lyaered_FBM_SAA::~Lyaered_FBM_SAA()
{

}


void Lyaered_FBM_SAA::Gauss_Layered_Rough_surface_FBM()
{
//////////////////2008.11.13///////////////////////
//////////////////////////////////////////////////////////////////////////
/////        如果媒质介电常数虚数部分过大则TM波收敛性更好           //////
//////////////////////////////////////////////////////////////////////////


////////////////////////////////姬伟杰////////////////////////////////////
	
     int NUMBER_REALIZATION =1; // The number of realization.
	//一维导体粗糙面
	 wave_len=3.0e8/(3.0e9);      //wave length
	 h=0.05*wave_len;            //rms height 
	 lc=0.05*wave_len;            //correlation length 
	 i_an=0.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

    double g_p =5.*wave_len;

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	N=(rl/wave_len)*20;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;


	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828183;

	double d=5.0*wave_len;

	double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta0;
	double Delta1;

	Vec_Dou Delta;
	        Delta.assign(15);


	Complex NR(0.,0.);
	
	Complex er2(6.0,0.6);
	Complex er3(25,2.5);

//	Complex rou_2 = er2;             //The TM wave
//	Complex rou_3 = er3/er2;           

	double rou_2 = 1.;               //The TE wave
	double rou_3 = 1.;

	double K1;
	Complex K2;
	Complex K3;

	Complex zero(0.,0.);

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);
	Complex Sum_E0(0.,0.);
	Complex Sum_E1(0.,0.);                                                                                                           

	int i,j,k,l,m,s;

	

	Vec_Com A;
	        A.assign(4*4);
	Vec_Com X;
	        X.assign(4);
	Vec_Com V;
	        V.assign(4);
	
	Vec_Com Ei;
	        Ei.assign(N);

	Vec_Dou Z1;
	        Z1.assign(N);
	Vec_dou Z2;
	        Z2.assign(N);

	Vec_Dou df_1;
	        df_1.assign(N);
	Vec_dou df_2;
	        df_2.assign(N);
			
	Vec_Dou ddf_1;
	        ddf_1.assign(N);
	Vec_dou ddf_2;
	        ddf_2.assign(N);

	Vec_Com U1_f;
	        U1_f.assign(N);
	Vec_Com U1_b;
	        U1_b.assign(N);
	Vec_Com U2_f;
	        U2_f.assign(N);
	Vec_Com U2_b;
	        U2_b.assign(N);

	Vec_Com W1_f;
	        W1_f.assign(N);
	Vec_Com W1_b;
	        W1_b.assign(N);
	Vec_Com W2_f;
	        W2_f.assign(N);
	Vec_Com W2_b;
	        W2_b.assign(N);

	Vec_Com Jf1;
	        Jf1.assign(4*N);
	Vec_Com Jf2;
	        Jf2.assign(4*N);
	Vec_Com Jf3;
	        Jf3.assign(4*N);

	Vec_Com Resi_V;
            Resi_V.assign(4*N);




	Vec_Com A_1_1_1;
	        A_1_1_1.assign(N*N);
	Vec_Com B_1_1_1;
	        B_1_1_1.assign(N*N);

			
	Vec_Com A_2_1_1;
	        A_2_1_1.assign(N*N);
	Vec_Com B_2_1_1;
	        B_2_1_1.assign(N*N);

	Vec_Com A_2_1_2;
	        A_2_1_2.assign(N*N);
	Vec_Com B_2_1_2;
	        B_2_1_2.assign(N*N);

	Vec_Com A_2_2_1;
	        A_2_2_1.assign(N*N);
	Vec_Com B_2_2_1;
	        B_2_2_1.assign(N*N);

	Vec_Com A_2_2_2;
	        A_2_2_2.assign(N*N);
	Vec_Com B_2_2_2;
	        B_2_2_2.assign(N*N);

	Vec_Com A_3_2_2;
	        A_3_2_2.assign(N*N);
	Vec_Com B_3_2_2;
	        B_3_2_2.assign(N*N);

	Vec_Com ZT;
	        ZT.assign(4*N*4*N);
	Vec_Com EE;
	        EE.assign(4*N);
	Vec_Com Jf;
	        Jf.assign(4*N);


	
	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);

	 int N0=180;
     double R_D = 10000.*wave_len;

	 Vec_Com Es;
		 Es.assign(N0);
	 Vec_dou D_Es;
		  D_Es.assign(N0);
     Vec_dou DB_Es;
		  DB_Es.assign(N0);
				 
	 Vec_Com Es1;
		  Es1.assign(N0);	
	 Vec_dou D_Es1;
		  D_Es1.assign(N0);
     Vec_dou DB_Es1;
		  DB_Es1.assign(N0);			
		
			
	 Vec_dou Sita_D;
		  Sita_D.assign(N0);
			
	 Vec_Dou x_D;
		  x_D.assign(N0);
	 Vec_Dou z_D;
		  z_D.assign(N0);




	 int N1 = N;
		Vec_Dou ZZ;
		        ZZ.assign(N1);
	 Vec_Dou XX;
		        XX.assign(N1);

	 Vec_Dou EE_Total;
		        EE_Total.assign(N1*N1);

	Vec_Com U_in;    U_in.assign(N);			

					 
	 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
		    x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);
		 }

	double i_an_1 = i_an*PI/180.;

	mean=0.;              //均值
	vari=1.;              //方差
	seed=123456.;            //seed for random number generator 
	RL=rl;             //rough surface length
	k_w=2.*PI/wave_len;  //wave number

	 
	
for (s=0;s<NUMBER_REALIZATION;s++)
{

	cout<<"第"<<s+1<<"次生成粗糙面"<<endl;
	
//	h=0.2*wave_len/(2.0*PI);            //rms height 
//	lc=1.0*wave_len/(2.0*PI);            //correlation length 
	h=0.05*wave_len;            //rms height 
	lc=1.0*wave_len;            //correlation length 

	gauss_surface();	
//	One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	dx_Suf = RL/N;




/*    for (i=0;i<N;i++)
    {
		z[i] = 0.;
    }*/

	for (i=0;i<N;i++)
	{
		Z1[i] = z[i];
//		Z2[i] = z[i] - d;df_2[i]=ddf_2[i] =
	}

	for (i=0;i<N;i++)
	{
		df_1[i] = df_x[i];
		
		ddf_1[i] = ddf_x[i];
		 
	}

		 
//	h=0.1*wave_len/(2.0*PI);            //rms height 
//	lc=1.0*wave_len/(2.0*PI);            //correlation length 

	h=0.05*wave_len;            //rms height 
	lc=1.0*wave_len;            //correlation length 

	gauss_surface();
//	One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	for (i=0;i<N;i++)
	{
		Z2[i] = z[i] - d;
	}

	for (i=0;i<N;i++)
	{
		df_2[i]= df_x[i];		
		ddf_2[i] = ddf_x[i];		 
	}



	double Yita = 120.*PI;

	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an_1)-Z1[i]*cos(i_an_1);
		kg=pow(K1*g_p*cos(i_an_1),2.0);
		fac=pow((x[i]+Z1[i]*tan(i_an_1))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);
		U_in[i] = Ei[i]*(J*K1*(1.0+ww)*(-df_x[i]*sin(i_an_1)-cos(i_an_1))+4.0*J*K1*xsz*(x[i]+z[i]*tan(i_an_1))/pow(g_p*K1*g_p*cos(i_an_1),2.0)*(tan(i_an_1)-df_x[i])+2.0*(x[i]+z[i]*tan(i_an_1))/pow(g_p,2.0)*(df_x[i]-tan(i_an_1)));
	}

/////////////////////////////////////////////////////////////////////////////////////////


	cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}




   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (-df_1[j]*(x[i]-x[j])+(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}

//////////////////////////////////////2_1_1///////////////////////////////////////////////
   	
   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}



   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (-df_1[j]*(x[i]-x[j])+(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////2_1_2///////////////////////////////////////////////

   
   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
//			dl = sqrt(1.+pow(df_2[j],2.0));

		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
		    A_2_1_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
				NR = (-df_2[j]*(x[i]-x[j])+(Z1[i]-Z2[j]))/rm;
				B_2_1_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_1///////////////////////////////////////////////


   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
//			dl = sqrt(1.+pow(df_1[i],2.0));

		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
		    A_2_2_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

		      rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
			  NR = (-df_1[j]*(x[i]-x[j])+(Z2[i]-Z1[j]))/rm;
			  B_2_2_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_2///////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (-df_2[j]*(x[i]-x[j])+(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////3_2_2///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 + ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (-df_2[j]*(x[i]-x[j])+(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}


		

/////////////////////////////The impedance matrix for the MOM/////////////////////////////////////////////

 /*   for (i=0;i<N;i++)
    {
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = A_1_1_1[i*N+j];
		}
        for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_1_1_1[i*N+j-N];
		}
	}
	for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = rou_2*A_2_1_1[(i-N)*N+j];

		}
		for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_2_1_1[(i-N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = A_2_1_2[(i-N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_2_1_2[(i-N)*N+j-3*N];
		}
	}

	for (i=2*N;i<3*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = rou_2*A_2_2_1[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_2_2_1[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = A_2_2_2[(i-2*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_2_2_2[(i-2*N)*N+j-3*N];
		}
	}

	for (i=3*N;i<4*N;i++)
	{
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = rou_3*A_3_2_2[(i-3*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_3_2_2[(i-3*N)*N+j-3*N];
		}
	}

	for (i=0;i<N;i++)
	{
		EE[i] = Ei[i];
	}

*/


	cout<<"The matrix has been filled."<<endl;


//	Bi_CG(ZT,EE,Jf);                     //The MOM method


//CLinearEquations linear;



  for (l=0;l<N;l++)
  {

	cout<<"Begin the "<<l+1<<"th step of the FBM......------>>"<<endl;

/////////////////////////The forward sweep/////////////////////////////////////////////////


    for (i=0;i<N;i++)
    {
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = Ei[i] - Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

		


			GAUSS(A,X,V,4);

		

			U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

    }


/////////////////////////The backward sweep/////////////////////////////////////////////////
	
	
   for (i=N-1;i>=0;i--)
    {
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;

			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

			
				
			GAUSS(A,X,V,4);


			U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];

    }

//////////////////////////The convergence condition////////////////////////////////////////////////

	
	for (i=0;i<N;i++)
	{
		Jf1[i]   = U1_b[i];
		Jf1[i+N] = W1_b[i];
		Jf1[i+2*N] = U2_b[i];
		Jf1[i+3*N] = W2_b[i];
	}

	for (i=0;i<4*N;i++)
	{
		Jf3[i] = Jf1[i] - Jf2[i];
	}




	Delta0 = sqrt(L2_Norm(Jf3))/sqrt(L2_Norm(Jf1)); 
	cout<<"The error of the FBM is: Delta ="<<Delta0<<endl;
	

	  if (Delta0<0.01)
	  {
		  	
		  for (i=0;i<4*N;i++)
		  {
		    for (j=0;j<4*N;j++)
			{
             Sum1  = Sum1 + ZT[i*4*N+j]*(Jf1[j] - Jf2[j]);
			}
		     Resi_V[i] = Sum1;
		     Sum1  =  zero;
		  }
         
  
		 Delta1 = sqrt(L2_Norm(Resi_V))/sqrt(L2_Norm(EE));
//		 Delta[l] = Delta1;

	    if(Delta0<=0.01)
		{
	       cout<<"The total number of the FBM step is :"<<l+1<<"."<<"and the final error is :"<<Delta1<<endl;
		    break;
		}
				
	     else  
		     continue;
	  
	  
	  }
	  
	
		Assign_Vecor(Jf1,Jf2); 	  

}
       
 



  
			
/*	   double max_z;
			for (i=0;i<N;i++)
			{
				if (max_z<z[i])
				{
					max_z = z[i];
				}
				else 
					continue;
			}

		   for (i=0;i<N1;i++)
		   {
			XX[i] = x[i];
			ZZ[i] = max_z+(i+1)*dx_Suf;
		   }


/*		   for (i=0;i<N1;i++)
		   {
			   for (j=0;j<N1;j++)
			   {
				    xsz=XX[i]*sin(i_an)-ZZ[j]*cos(i_an);
		            kg=pow(K1*g_p*cos(i_an),2.0);
		            fac=pow((XX[i]+ZZ[j]*tan(i_an))/g_p,2.0);
		            ww=(2.0*fac-1.0)/kg;
				   
				   for (k=0;k<N;k++)
				   {
				   	  rm = sqrt(pow(XX[i]-x[k],2.0)+pow(ZZ[j]-z[k],2.0));
					  NR = (-1.*df_x[j]*(XX[i]-x[k])+(ZZ[j]-z[k]))/(rm);
					 /* if (rm==0.)
					  {
						 dl = sqrt(1.+pow(df_1[k],2.0));

                         Sum0 = Sum0 + dx_Suf*((0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0))))*(W1_f[k]+W1_b[k])-J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)))*(U1_f[k]+U1_b[k]));
					  }

					  else*/
/*					  {	}				 						  
						  Sum0 = Sum0 + dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[k]+W1_b[k])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[k]+U1_b[k]);
					  

				   }                   
					EE_Total[i*N1+j] = abs(exp(J*K1*xsz*(1.0+ww))*exp(-fac) + Sum0);
					Sum0 = zero;
			     double r = sqrt(pow(XX[i],2.0)+pow(ZZ[j],2.0));
			     EE_Total[i*N1+j]= r*pow(EE_Total[i*N1+j],2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 EE_Total[i*N1+j]=10.*log10( EE_Total[i*N1+j]);
			   }

		   }*/


			
	double KG = 1-(1.+2.*pow(tan(i_an_1),2.0))/(2.*pow(K1*g_p*cos(i_an_1),2.0));

	
	
		 
/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {					 
			   Es[i] = Es[i] + dx_Suf*((J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i])))*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			     
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

           
		 }*/

		double abx   =  8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an*PI/180)*(1.0-((1.0+2.0*pow(tan(i_an*PI/180),2))/(2.0*pow(k_w*g_p*cos(i_an*PI/180),2))));

/*	     for (i=0;i<N0;i++)
		 {
			 			
			 Complex integ=0.0;

			 for (j=0;j<N;j++)
			 {
				double ksc=-k_w*(x[j]*sin(Sita_D[i] )+z[j]*cos(Sita_D[i] ));
				integ=integ+JJ[j]*exp(i_*ksc);
			 }

			Es1[i] = integ;

			double temp = norm(integ)/abx;
			D_Es[i] = (temp+nn*D_Es[i])/(nn+1);
			
			temp = 10.*log10(D_Es[i]);
			DB_Es[i] = (temp + nn * DB_Es[i]) / (nn + 1);*/


		    for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-Z1[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
//					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					NR = k_w*(df_1[j]*sin(Sita_D[i])-cos(Sita_D[i]))/sqrt(1+pow(df_1[j],2.0));


				    double ksc=-k_w*(x[j]*sin(Sita_D[i] )+Z1[j]*cos(Sita_D[i] ));
				    Sum_E0=Sum_E0+dx_Suf*(NR*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(i_*ksc);
//					Sum_E0 = Sum_E0+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

                 Es[i] = Sum_E0;				 
				 D_Es[i]  = (s*D_Es[i]+norm(Sum_E0)/abx)/(s+1);
                 DB_Es[i] = (s*DB_Es[i]+10.*log10(D_Es[i]))/(s+1);
		         Sum_E0   = zero;
			}
         
			
/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Sum_E0 = Sum_E0+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(Jf[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(Jf[j]);
				}


                   Es1[i] = Sum_E0;
				   Sum_E0 = zero;
				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an_1)*KG);
                 DB_Es1[i]=DB_Es1[i]+10.*log10(D_Es1[i]);
		
			}*/
		 
}  


/*			for (i=0;i<N0;i++)
			{
				                
				D_Es[i] = D_Es[i]/(s);
                D_Es1[i] = D_Es1[i]/(s);

                DB_Es[i] = DB_Es[i]/(s);
                DB_Es1[i] = DB_Es1[i]/(s);
			}*/


	Complex Tou1, Fan1,Tou2,Fan2;
	double R1,T1,R2,T2;
	double KG;
	         KG = (k_w*g_p*sqrt(PI/2.0)*cos(i_an*PI/180)*(1.0-((1.0+2.0*pow(tan(i_an*PI/180),2))/(2.0*pow(k_w*g_p*cos(i_an*PI/180),2)))));
	         Tou1 = 0.0;  Fan1 = 0.0; Tou2=0.0; Fan2=0.0;
	         for (j=0;j<N;j++)
	         {
				 Tou1 = Tou1 + dx_Suf*(W1_f[j]+W1_b[j])*conj(U1_f[j]+U1_b[j]);
				 Fan1 = Fan1 + dx_Suf*(W1_f[j]+W1_b[j]-Ei[j])*conj(U1_f[j]+U1_b[j]-U_in[j]);

				 Tou2 = Tou2 + dx_Suf*(W2_f[j]+W2_b[j])*conj(U2_f[j]+U2_b[j]);
	         }

			 T1 = Tou1.imag()/KG;
			 T2 = Tou2.imag()/KG;

			 R1 = Fan1.imag()/KG;

			 cout<<"T1 = "<<T1<<endl;
			 cout<<"T2 = "<<T2<<endl;
			 cout<<"R1 = "<<R1<<endl;
			 cout<<"R1+T1 = "<<T1-R1<<endl;

	ofstream out_U1;out_U1.open("Layered_FBM_U1.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<abs(Jf[i])<<endl;
			 }
   
    ofstream out_W1;out_W1.open("Layered_FBM_W1.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<abs(Jf[i+N])<<endl;
			 }

	ofstream out_U2;out_U2.open("Layered_FBM_U2.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<abs(Jf[i+2*N])<<endl;
			 }
	ofstream out_W2;out_W2.open("Layered_FBM_W2.txt");	        
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<abs(Jf[i+3*N])<<endl;
			 }
	ofstream out_Es;out_Es.open("Layered_FBM_Es.txt");	         
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	ofstream out_Es1;out_Es1.open("Layered_FBM_Es1.txt");	         
			 for (i=0;i<N0;i++)
			 {
				 out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<D_Es1[i]<<"   "<<DB_Es1[i]<<endl;
			 }
	ofstream out_z;out_z.open("Layered_FBM_Z.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_z<<x[i]/wave_len<<"   "<<Z1[i]/wave_len<<"   "<<Z2[i]/wave_len<<endl;
			 }
	ofstream out_Delta;out_Delta.open("Lyaered_FBM_Delta.txt");	         
			 for (i=0;i<15;i++)
			 {
				 out_Delta<<i<<"   "<<Delta[i]<<endl;
			 }
}

void Lyaered_FBM_SAA::Gauss_Layered_Rough_surface_FBM_SAA()
{
//////////////////////////////2008.11.20//////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////        如果媒质介电常数虚数部分过大则TM波收敛性更好           //////
//////////////////////////////////////////////////////////////////////////
/////          只是对场点和源点在同一粗糙面的情况进行了加速         ////// 
/////	     而场点和源点不在同一粗糙面的情况还用的是常规的FBM法    //////
////////////////////////////////姬伟杰////////////////////////////////////
	
	//一维导体粗糙面
	double wave_len=3.0e8/3e9;         //wave length
	double h=0.5*wave_len;             //rms height 
	double lc=1.0*wave_len;             //correlation length 
	double i_an=30.0;                  //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

    double g_p =10.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	int N=(rl/wave_len)*10;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                                      //total number of surface realizations
	int n_sc=180;                                   //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));          //rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;


	
	double eps=1.0e-2;                                 //CGM的收敛准则
	double criterion=1.0e-2;                           //BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828183;

	double d=5.*wave_len;

	double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta;


	Complex NR(0.,0.);
	
	Complex er2(2.0,0.05);
	Complex er3(25,0);

	Complex rou_2 = er2;                                     //The TM wave
	Complex rou_3 = er3/er2;           

//	double rou_2 = 1.;                                       //The TE wave
//	double rou_3 = 1.;

	double K1;
	Complex K2;
	Complex K3;

	Complex zero(0.,0.);

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);
	

	int i,j,k,l;

	Vec_Com A;
	        A.assign(4*4);
	Vec_Com X;
	        X.assign(4);
	Vec_Com V;
	        V.assign(4);
	
	Vec_Com Ei;
	        Ei.assign(N);

	Vec_Dou Z1;
	        Z1.assign(N);
	Vec_dou Z2;
	        Z2.assign(N);

	Vec_Dou df_1;
	        df_1.assign(N);
	Vec_dou df_2;
	        df_2.assign(N);
			
	Vec_Dou ddf_1;
	        ddf_1.assign(N);
	Vec_dou ddf_2;
	        ddf_2.assign(N);

	Vec_Com U1_f;
	        U1_f.assign(N);
	Vec_Com U1_b;
	        U1_b.assign(N);
	Vec_Com U2_f;
	        U2_f.assign(N);
	Vec_Com U2_b;
	        U2_b.assign(N);

	Vec_Com W1_f;
	        W1_f.assign(N);
	Vec_Com W1_b;
	        W1_b.assign(N);
	Vec_Com W2_f;
	        W2_f.assign(N);
	Vec_Com W2_b;
	        W2_b.assign(N);

	Vec_Com Jf1;
	        Jf1.assign(4*N);
	Vec_Com Jf2;
	        Jf2.assign(4*N);
	Vec_Com Jf3;
	        Jf3.assign(4*N);





	Vec_Com A_1_1_1;
	        A_1_1_1.assign(N*N);
	Vec_Com B_1_1_1;
	        B_1_1_1.assign(N*N);

			
	Vec_Com A_2_1_1;
	        A_2_1_1.assign(N*N);
	Vec_Com B_2_1_1;
	        B_2_1_1.assign(N*N);

	Vec_Com A_2_1_2;
	        A_2_1_2.assign(N*N);
	Vec_Com B_2_1_2;
	        B_2_1_2.assign(N*N);

	Vec_Com A_2_2_1;
	        A_2_2_1.assign(N*N);
	Vec_Com B_2_2_1;
	        B_2_2_1.assign(N*N);

	Vec_Com A_2_2_2;
	        A_2_2_2.assign(N*N);
	Vec_Com B_2_2_2;
	        B_2_2_2.assign(N*N);

	Vec_Com A_3_2_2;
	        A_3_2_2.assign(N*N);
	Vec_Com B_3_2_2;
	        B_3_2_2.assign(N*N);

	Vec_Com ZT;
	        ZT.assign(4*N*4*N);
	Vec_Com EE;
	        EE.assign(4*N);
	Vec_Com Jf;
	        Jf.assign(4*N);
	          		
	One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	dx_Suf = RL/N;

	i_an = i_an*PI/180.;

	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);

//////////////////////////////////////////////////////////////////////////

	double max_z = 0.;
	double min_z = 0.;
	double L_s_1,L_s_2;
	double Bita_s_1,Bita_s_2,Bita_s_3;
	double R_s_1,R_s_2;
	double Fai_S_max_1,Fai_S_max_2;
	double Delta_0_1,Delta_0_2,Delta_0_3;
	double Delta_Fai_1,Delta_Fai_2,Delta_Fai_3;
	Complex Fai_p;
	
	int Q1,Q2,Q3;
	int N_s_1,N_s_2;
		
	for (i=0;i<N;i++)
	{
		if (z[i]>max_z)
		{
			max_z=z[i];
		}
	   if (z[i]<=min_z)
	    {
			min_z=z[i];
	    }
	}

//////////////////////The parameter of the SAA////////////////////////////////////////////////////


	L_s_1  = 10.*wave_len ;                             //The strong interaction zone.
	L_s_2 =  L_s_1;

	R_s_1 = sqrt(pow(L_s_1,2.0)+pow(max_z-min_z,2.0)); 		
	R_s_2 = sqrt(pow(L_s_2,2.0)+pow(max_z-min_z,2.0));
	
//	cout<<"R_S_1 = "<<R_s_1<<endl;
//	cout<<"R_S_2 = "<<R_s_2<<endl;


	
	N_s_1 = int(L_s_1/dx_Suf);                                      //The number of the strong interaction point.                          
    N_s_2 = int(L_s_2/dx_Suf);

//	cout<<"N_s_1 = "<<N_s_1<<endl;
//	cout<<"N_s_2 = "<<N_s_2<<endl;


//  cout<<"L_s_1="<<L_s_1/wave_len<<"wave_len"<<endl; 	
//	cout<<"L_s_2="<<L_s_2/wave_len<<"wave_len"<<endl;

	Bita_s_1 = sqrt(10./(K1*L_s_1));                                //The geometric optics lit region.
	Bita_s_2 = sqrt(10./(abs(K2)*L_s_2));
	Bita_s_3 = sqrt(10./(abs(K3)*L_s_1));
	

	
	Fai_S_max_1 = atan((max_z-min_z)/L_s_1);                       //The outermost possible saddle point.
	Fai_S_max_2 = atan((max_z-min_z)/L_s_2);

//	cout<<"Fai_S_max_1 = "<<Fai_S_max_1<<endl;
//	cout<<"Fai_S_max_2 = "<<Fai_S_max_2<<endl;
	
	

	if (sqrt(K1*R_s_1/20.)*Fai_S_max_1-1>1)
	{
		Delta_0_1 = atan(1./(sqrt(K1*R_s_1/20.)*Fai_S_max_1-1.));     
	}
	else
	{
		Delta_0_1 = PI/5.;
	}

    	
	if (sqrt(abs(K2)*R_s_2/20.)*Fai_S_max_2-1>1)
	{
		Delta_0_2 = atan(1./(sqrt(abs(K2)*R_s_2/20.)*Fai_S_max_2-1.));     
	}
	else
	{
		Delta_0_2 = PI/5.;
	}

		
	if (sqrt(abs(K3)*R_s_1/20.)*Fai_S_max_1-1>1)
	{
		Delta_0_3 = atan(1./(sqrt(abs(K3)*R_s_1/20.)*Fai_S_max_1-1.));     
	}
	else
	{
		Delta_0_3 = PI/5.;
	}
	

	cout<<"Delta_0_1 = "<<Delta_0_1<<endl;
    cout<<"Delta_0_2 = "<<Delta_0_2<<endl;
    cout<<"Delta_0_3 = "<<Delta_0_3<<endl;

	
	Delta_Fai_1 = sqrt(5./(K1*R_s_1))/22.; 
	Delta_Fai_2 = sqrt(5./(abs(K2)*R_s_2))/22.; 
	Delta_Fai_3 = sqrt(5./(abs(K3)*R_s_1))/22.; 

	cout<<"Delta_Fai_1 = "<<Delta_Fai_1<<endl;
    cout<<"Delta_Fai_2 = "<<Delta_Fai_2<<endl;
    cout<<"Delta_Fai_3 = "<<Delta_Fai_3<<endl;


	Q1 = int(Bita_s_1/Delta_Fai_1)+5;
	Q2 = int(Bita_s_2/Delta_Fai_2)+5;
	Q3 = int(Bita_s_3/Delta_Fai_3)+5;

	cout<<"Q1 = "<<Q1<<endl;
    cout<<"Q2 = "<<Q2<<endl;
    cout<<"Q3 = "<<Q3<<endl;

    Vec_Com FN_f_1;
           FN_f_1.assign(2*Q1+1);
    
    Vec_Com FN_f_2_1;
           FN_f_2_1.assign(2*Q2+1);
	Vec_Com FN_f_3_1;
           FN_f_3_1.assign(2*Q2+1);
		   
	Vec_Com FN_f_2_2;
           FN_f_2_2.assign(2*Q2+1);
	Vec_Com FN_f_3_2;
           FN_f_3_2.assign(2*Q2+1);

	Vec_Com FN_f_4;
           FN_f_4.assign(2*Q3+1);
	
	Vec_Com FN_b_1;
           FN_b_1.assign(2*Q1+1);
    
	Vec_Com FN_b_2_1;
           FN_b_2_1.assign(2*Q2+1);
	Vec_Com FN_b_3_1;
           FN_b_3_1.assign(2*Q2+1);

	    
	Vec_Com FN_b_2_2;
           FN_b_2_2.assign(2*Q2+1);
	Vec_Com FN_b_3_2;
           FN_b_3_2.assign(2*Q2+1);
	   
	Vec_Com FN_b_4;
           FN_b_4.assign(2*Q3+1);	   
		   
 
		   
    Vec_dou w_Sita_1;
	       w_Sita_1.assign(2*Q1+1);
    Vec_dou w_Sita_2;
	       w_Sita_2.assign(2*Q2+1);
    Vec_dou w_Sita_3;
	       w_Sita_3.assign(2*Q3+1);

	Vec_Com Fai;
	        Fai.assign(2*Q1+1);
  

	for (k=-Q1;k<=Q1;k++)
	{   
	     Fai_p = k*Delta_Fai_1*exp(-1.*J*Delta_0_1);                //The window function.
         

		 if (abs(Fai_p)<=Bita_s_1)				  
			w_Sita_1[Q1+k] = 1.;				  
		 else				  
		    w_Sita_1[Q1+k] = 0.;
	}

		
	for (k=-Q2;k<=Q2;k++)
	{   
	     Fai_p = k*Delta_Fai_2*exp(-1.*J*Delta_0_2);                //The window function.
         		 
		 
		 if (abs(Fai_p)<=Bita_s_2)				  
			w_Sita_2[Q2+k] = 1.;				  
		 else				  
		    w_Sita_2[Q2+k] = 0.;
	}

		
	for (k=-Q3;k<=Q3;k++)
	{   
	     Fai_p = k*Delta_Fai_3*exp(-1.*J*Delta_0_3);                //The window function.
         Fai[k+Q1] = Fai_p;
		 
		 if (abs(Fai_p)<=Bita_s_3)				  
			w_Sita_3[Q3+k] = 1.;				  
		 else				  
		    w_Sita_3[Q3+k] = 0.;
		 
	}




//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		Z1[i] = z[i];
		Z2[i] = z[i] - d;
	}

	for (i=0;i<N;i++)
	{
		df_1[i] =df_2[i]= df_x[i];
		
		ddf_1[i] =ddf_2[i] = ddf_x[i];
		 
	}

	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(K1*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);	
	}

//////////////////////////////////////////////////////////////////////////


	cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////


	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}

	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}

	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}
	







	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}

	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}


	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}








//////////////////////////////////////2_1_1///////////////////////////////////////////////
   	
	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}


	
	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}
  

	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}









//////////////////////////////////////2_1_2///////////////////////////////////////////////

   
   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
		    A_2_1_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z1[i]-Z2[j]))/rm;
				B_2_1_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_1///////////////////////////////////////////////


   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
		    A_2_2_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

		      rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
			  NR = (df_1[j]*(x[i]-x[j])-(Z2[i]-Z1[j]))/rm;
			  B_2_2_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_2///////////////////////////////////////////////

	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}








	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

		
	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

		
	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}


//////////////////////////////////////3_2_2///////////////////////////////////////////////


	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}

		
	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}
		

	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}





	for (i=0;i<N_s_1+1;i++)
	{
		for (j=0;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 - ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}

		
	for (i=N_s_1+1;i<N-N_s_1;i++)
	{
		for (j=i-N_s_1;j<=i+N_s_1;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 - ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}
	
	
	
	for (i=N-N_s_1;i<N;i++)
	{
		for (j=i-N_s_1;j<=i;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 - ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}


/////////////////////////////The impedance matrix for the MOM/////////////////////////////////////////////
/*
    for (i=0;i<N;i++)
    {
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = A_1_1_1[i*N+j];
		}
        for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_1_1_1[i*N+j-N];
		}
	}
	for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = rou_2*A_2_1_1[(i-N)*N+j];

		}
		for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_2_1_1[(i-N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = A_2_1_2[(i-N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_2_1_2[(i-N)*N+j-3*N];
		}
	}

	for (i=2*N;i<3*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT[i*4*N+j] = rou_2*A_2_2_1[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT[i*4*N+j] = B_2_2_1[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = A_2_2_2[(i-2*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_2_2_2[(i-2*N)*N+j-3*N];
		}
	}

	for (i=3*N;i<4*N;i++)
	{
		for (j=2*N;j<3*N;j++)
		{
			ZT[i*4*N+j] = rou_3*A_3_2_2[(i-3*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT[i*4*N+j] = B_3_2_2[(i-3*N)*N+j-3*N];
		}
	}

	for (i=0;i<N;i++)
	{
		EE[i] = Ei[i];
	}*/


	cout<<"The matrix has been filled."<<endl;

//	Bi_CG(ZT,EE,Jf);                     //The MOM method




//	CLinearEquations linear;

for (l=0;l<N;l++)
{

	cout<<"Begin the "<<l+1<<"th step of the FBM_SAA......------>>"<<endl;

/////////////////////////The forward sweep/////////////////////////////////////////////////


    for (i=0;i<=N_s_1;i++)
    {
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = Ei[i] - Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

		


			GAUSS(A,X,V,4);

		

			U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

    }

//////////////////////////////////////////////////////////////////////////

	    i = N_s_1+1;
	
	    for (j=i-N_s_1;j<i;j++)
		{
			
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum2 = Sum2 + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
            Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);		

		}
	    
		for (j=0;j<i;j++)
		{					
			Sum1 = Sum1 + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]);
		
		}
		
            V[0] = Ei[i] - Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

				
			for (k=-Q1;k<=Q1;k++)
			{
			   Fai_p = k*Delta_Fai_1*exp(-1.*J*Delta_0_1);

			   FN_f_1[k+Q1] = (-J*K1*(-sin(Fai_p)+df_1[i-N_s_1-1]*cos(Fai_p))*(W1_f[i-N_s_1-1]+W1_b[i-N_s_1-1])+(U1_f[i-N_s_1-1]+U1_b[i-N_s_1-1]))*exp(J*K1*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K1*Z1[i-N_s_1-1]*sin(Fai_p));

			   Sum0 = Sum0 + J*dx_Suf*exp(-J*Delta_0_1)/(4.*PI)*w_Sita_1[k+Q1]*FN_f_1[k+Q1]*exp(J*K1*Z1[i]*sin(Fai_p))*Delta_Fai_1;

			}

					   
			for (k=-Q2;k<=Q2;k++)
			{
				Fai_p = k*Delta_Fai_2*exp(-1.*J*Delta_0_2);

				FN_f_2_1[k+Q2] = (-rou_2*(U1_f[i-N_s_2-1]+U1_b[i-N_s_2-1])+J*K2*(-sin(Fai_p)+df_1[i-N_s_2-1]*cos(Fai_p))*(W1_f[i-N_s_2-1]+W1_b[i-N_s_2-1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i-N_s_2-1]*sin(Fai_p));

			    FN_f_3_1[k+Q2] =  ((U2_b[i-N_s_2-1]+U2_f[i-N_s_2-1])-J*K2*(sin(Fai_p)+df_2[i-N_s_2-1]*cos(Fai_p))*(W2_b[i-N_s_2-1]+W2_f[i-N_s_2-1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i-N_s_2-1]*sin(Fai_p));

				
				Sum1 = Sum1 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_f_2_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;
			   
				Sum2 = Sum2 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_f_3_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;

			}

					   
		    for (k=-Q3;k<=Q3;k++)
			{
   	            Fai_p = k*Delta_Fai_3*exp(-1.*J*Delta_0_3);
				
				FN_f_4[k+Q3] = (-rou_3*(U2_f[i-N_s_1-1]+U2_b[i-N_s_1-1])+J*K3*(-sin(Fai_p)+df_2[i-N_s_1-1]*cos(Fai_p))*(W2_f[i-N_s_1-1]+W2_b[i-N_s_1-1]))
					           *exp(J*K3*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K3*Z1[i-N_s_1-1]*sin(Fai_p));

				Sum3 = Sum3 + J*dx_Suf*exp(-J*Delta_0_3)/(4.*PI)*w_Sita_3[k+Q3]*FN_f_4[k+Q3]*exp(J*K3*Z1[i]*sin(Fai_p))*Delta_Fai_3;
			}

			
            V[0] = V[0] - Sum0;
			V[1] = V[1] - Sum1;
			V[2] = V[2] - Sum2;
			V[3] = V[3] - Sum3;

 		    Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			GAUSS(A,X,V,4);

	        U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

//////////////////////////////////////////////////////////////////////////


	for (i=N_s_1+2;i<N;i++)
	{
		for (j=i-N_s_1;j<i;j++)
		{							
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);	
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) ;
			Sum2 = Sum2 + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
	        Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);		
		}				
		for (j=0;j<i;j++)
		{			
			Sum1 = Sum1 + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]) ;
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) ;				
		}

            V[0] = Ei[i] - Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

						
			for (k=-Q1;k<=Q1;k++)
			{
			   Fai_p = k*Delta_Fai_1*exp(-1.*J*Delta_0_1);

			   FN_f_1[k+Q1] = FN_f_1[k+Q1]*exp(J*K1*dx_Suf*cos(Fai_p)) + (-J*K1*(-sin(Fai_p)+df_1[i-N_s_1-1]*cos(Fai_p))*(W1_f[i-N_s_1-1]+W1_b[i-N_s_1-1])+(U1_f[i-N_s_1-1]+U1_b[i-N_s_1-1]))*exp(J*K1*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K1*Z1[i-N_s_1-1]*sin(Fai_p));

			   Sum0 = Sum0 + J*dx_Suf*exp(-J*Delta_0_1)/(4.*PI)*w_Sita_1[k+Q1]*FN_f_1[k+Q1]*exp(J*K1*Z1[i]*sin(Fai_p))*Delta_Fai_1;

			}

			for (k=-Q2;k<=Q2;k++)
			{
				Fai_p = k*Delta_Fai_2*exp(-1.*J*Delta_0_2);

				FN_f_2_1[k+Q2] =  FN_f_2_1[k+Q2]*exp(J*K2*dx_Suf*cos(Fai_p)) + (-rou_2*(U1_f[i-N_s_2-1]+U1_b[i-N_s_2-1])+J*K2*(-sin(Fai_p)+df_1[i-N_s_2-1]*cos(Fai_p))*(W1_f[i-N_s_2-1]+W1_b[i-N_s_2-1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i-N_s_2-1]*sin(Fai_p));

			    FN_f_3_1[k+Q2] =  FN_f_3_1[k+Q2]*exp(J*K2*dx_Suf*cos(Fai_p)) + ((U2_b[i-N_s_2-1]+U2_f[i-N_s_2-1])-J*K2*(sin(Fai_p)+df_2[i-N_s_2-1]*cos(Fai_p))*(W2_b[i-N_s_2-1]+W2_f[i-N_s_2-1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i-N_s_2-1]*sin(Fai_p));

				
				Sum1 = Sum1 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_f_2_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;
			   
				Sum2 = Sum2 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_f_3_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;

			}


			for (k=-Q3;k<=Q3;k++)
			{
   	            Fai_p = k*Delta_Fai_3*exp(-1.*J*Delta_0_3);
				
				FN_f_4[k+Q3] = (-rou_3*(U2_f[i-N_s_1-1]+U2_b[i-N_s_1-1])+J*K3*(-sin(Fai_p)+df_2[i-N_s_1-1]*cos(Fai_p))*(W2_f[i-N_s_1-1]+W2_b[i-N_s_1-1]))
					           *exp(J*K3*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K3*Z1[i-N_s_1-1]*sin(Fai_p));

				Sum3 = Sum3 + J*dx_Suf*exp(-J*Delta_0_3)/(4.*PI)*w_Sita_3[k+Q3]*FN_f_4[k+Q3]*exp(J*K3*Z1[i]*sin(Fai_p))*Delta_Fai_3;
			}

	
			

			
            V[0] = V[0] - Sum0;
			V[1] = V[1] - Sum1;
			V[2] = V[2] - Sum2;
			V[3] = V[3] - Sum3;

 		    Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			GAUSS(A,X,V,4);

	        U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

		}

//////////////////////////////////////////////////////////////////////////
	    
	for (i=N-1;i>=N-N_s_1-1;i--)
    {
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

		


			GAUSS(A,X,V,4);

		

			U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];

    }

//////////////////////////////////////////////////////////////////////////

	    i = N-N_s_1-2;
	
	    for (j=i+1;j<i+N_s_1+1;j++)
		{
			
						
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
		    Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);	
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum2 = Sum2 + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);

			
		}

		for (j=i+1;j<N;j++)
		{	 		
			Sum1 = Sum1 + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]) ;
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) ;				
		}


            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

						
		   for (k=-Q1;k<=Q1;k++)
			{
			   Fai_p = k*Delta_Fai_1*exp(-1.*J*Delta_0_1);

			   FN_b_1[k+Q1] = (-J*K1*(-sin(Fai_p)+df_1[i+N_s_1+1]*cos(Fai_p))*(W1_f[i+N_s_1+1]+W1_b[i+N_s_1+1])+(U1_f[i+N_s_1+1]+U1_b[i+N_s_1+1]))*exp(J*K1*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K1*Z1[i+N_s_1+1]*sin(Fai_p));

			   Sum0 = Sum0 + J*dx_Suf*exp(-J*Delta_0_1)/(4.*PI)*w_Sita_1[k+Q1]*FN_b_1[k+Q1]*exp(J*K1*Z1[i]*sin(Fai_p))*Delta_Fai_1;

			}

					
		
		   for (k=-Q2;k<=Q2;k++)
			{
				Fai_p = k*Delta_Fai_2*exp(-1.*J*Delta_0_2);

				FN_b_2_1[k+Q2] =  (-rou_2*(U1_f[i+N_s_2+1]+U1_b[i+N_s_2+1])+J*K2*(-sin(Fai_p)+df_1[i+N_s_2+1]*cos(Fai_p))*(W1_f[i+N_s_2+1]+W1_b[i+N_s_2+1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i+N_s_2+1]*sin(Fai_p));
			    FN_b_3_1[k+Q2] =  ((U2_b[i+N_s_2+1]+U2_f[i+N_s_2+1])-J*K2*(sin(Fai_p)+df_2[i+N_s_2+1]*cos(Fai_p))*(W2_b[i+N_s_2+1]+W2_f[i+N_s_2+1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i+N_s_2+1]*sin(Fai_p));
			   
				Sum1 = Sum1 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_b_2_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;
			   
				Sum2 = Sum2 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_b_3_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;

			}

		   			
		   for (k=-Q3;k<=Q3;k++)
			{
   	            Fai_p = k*Delta_Fai_3*exp(-1.*J*Delta_0_3);
				
				FN_b_4[k+Q3] = (-rou_3*(U2_f[i+N_s_1+1]+U2_b[i+N_s_1+1])+J*K3*(-sin(Fai_p)+df_2[i+N_s_1+1]*cos(Fai_p))*(W2_f[i+N_s_1+1]+W2_b[i+N_s_1+1]))
					           *exp(J*K3*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K3*Z1[i+N_s_1+1]*sin(Fai_p));

				Sum3 = Sum3 + J*dx_Suf*exp(-J*Delta_0_3)/(4.*PI)*w_Sita_3[k+Q3]*FN_b_4[k+Q3]*exp(J*K3*Z1[i]*sin(Fai_p))*Delta_Fai_3;

			}	

			
            V[0] = V[0] - Sum0;
			V[1] = V[1] - Sum1;
			V[2] = V[2] - Sum2;
			V[3] = V[3] - Sum3;

 		    Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;
			
			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

			GAUSS(A,X,V,4);

	        U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];


//////////////////////////////////////////////////////////////////////////


	for (i=N-N_s_1-3;i>=0;i--)
	{
		for (j=i+1;j<i+N_s_1+1;j++)
		{
			
					
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
		    Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);	
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum2 = Sum2 + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);

			
		}

		for (j=i+1;j<N;j++)
		{				
			Sum1 = Sum1 + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]);
		
		}

            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;
			
			
			for (k=-Q1;k<=Q1;k++)
			{
			   Fai_p = k*Delta_Fai_1*exp(-1.*J*Delta_0_1);

			   FN_b_1[k+Q1] = FN_b_1[k+Q1]*exp(J*K1*dx_Suf*cos(Fai_p)) + (-J*K1*(-sin(Fai_p)+df_1[i+N_s_1+1]*cos(Fai_p))*(W1_f[i+N_s_1+1]+W1_b[i+N_s_1+1])+(U1_f[i+N_s_1+1]+U1_b[i+N_s_1+1]))
				              *exp(J*K1*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K1*Z1[i+N_s_1+1]*sin(Fai_p));

			   Sum0 = Sum0 + J*dx_Suf*exp(-J*Delta_0_1)/(4.*PI)*w_Sita_1[k+Q1]*FN_b_1[k+Q1]*exp(J*K1*Z1[i]*sin(Fai_p))*Delta_Fai_1;

			}


			for (k=-Q2;k<=Q2;k++)
			{
				Fai_p = k*Delta_Fai_2*exp(-1.*J*Delta_0_2);

				FN_b_2_1[k+Q2] = FN_b_2_1[k+Q2]*exp(J*K2*dx_Suf*cos(Fai_p)) + (-rou_2*(U1_f[i+N_s_2+1]+U1_b[i+N_s_2+1])+J*K2*(-sin(Fai_p)+df_1[i+N_s_2+1]*cos(Fai_p))*(W1_f[i+N_s_2+1]+W1_b[i+N_s_2+1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i+N_s_2+1]*sin(Fai_p));
			    FN_b_3_1[k+Q2] = FN_b_3_1[k+Q2]*exp(J*K2*dx_Suf*cos(Fai_p)) + ((U2_b[i+N_s_2+1]+U2_f[i+N_s_2+1])-J*K2*(sin(Fai_p)+df_2[i+N_s_2+1]*cos(Fai_p))*(W2_b[i+N_s_2+1]+W2_f[i+N_s_2+1]))*exp(J*K2*(N_s_2+1.)*dx_Suf*cos(Fai_p))*exp(-J*K2*Z1[i+N_s_2+1]*sin(Fai_p));			   
				
				Sum1 = Sum1 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_b_2_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;			   
				Sum2 = Sum2 + J*dx_Suf*exp(-J*Delta_0_2)/(4.*PI)*w_Sita_2[k+Q2]*FN_b_3_1[k+Q2]*exp(J*K2*Z1[i]*sin(Fai_p))*Delta_Fai_2;

			}

			for (k=-Q3;k<=Q3;k++)
			{
   	            Fai_p = k*Delta_Fai_3*exp(-1.*J*Delta_0_3);
				
				FN_b_4[k+Q3] = FN_b_4[k+Q3]*exp(J*K3*dx_Suf*cos(Fai_p)) + (-rou_3*(U2_f[i+N_s_1+1]+U2_b[i+N_s_1+1])+J*K3*(-sin(Fai_p)+df_2[i+N_s_1+1]*cos(Fai_p))*(W2_f[i+N_s_1+1]+W2_b[i+N_s_1+1]))
					           *exp(J*K3*(N_s_1+1.)*dx_Suf*cos(Fai_p))*exp(-J*K3*Z1[i+N_s_1+1]*sin(Fai_p));

				Sum3 = Sum3 + J*dx_Suf*exp(-J*Delta_0_3)/(4.*PI)*w_Sita_3[k+Q3]*FN_b_4[k+Q3]*exp(J*K3*Z1[i]*sin(Fai_p))*Delta_Fai_3;

			}

							
            V[0] = V[0] - Sum0;
			V[1] = V[1] - Sum1;
			V[2] = V[2] - Sum2;
			V[3] = V[3] - Sum3;

 		    Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

						
			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];


			GAUSS(A,X,V,4);

	        U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];

		}

//////////////////////////The convergence condition////////////////////////////////////////////////

	
	for (i=0;i<N;i++)
	{
		Jf1[i]   = U1_b[i];
		Jf1[i+N] = W1_b[i];
		Jf1[i+2*N] = U2_b[i];
		Jf1[i+3*N] = W2_b[i];
	}

	for (i=0;i<4*N;i++)
	{
		Jf3[i] = Jf1[i] - Jf2[i];
	}

	  Assign_Vecor(Jf1,Jf2);  
	  Delta = sqrt(L2_Norm(Jf3))/sqrt(L2_Norm(Jf1));      
	  
	  cout<<"The error is: Delta ="<<Delta<<endl;
			  
	  if(Delta<=0.01)
	  {
	    cout<<"The total number of the step is :"<<l+1<<"."<<"and the final error is :"<<Delta<<endl;
		  break;
	  }
				
	     else  
		     continue;

}


    	 int N0=180;
		 double R_D = 10000.*wave_len;

		 Vec_Com Es;
		    Es.assign(N0);
		 Vec_dou D_Es;
		    D_Es.assign(N0);
    	 Vec_dou DB_Es;
		    DB_Es.assign(N0);
				 
		 Vec_Com Es1;
		    Es1.assign(N0);	
		 Vec_dou D_Es1;
		    D_Es1.assign(N0);
    	 Vec_dou DB_Es1;
		    DB_Es1.assign(N0);			
		
			
		 Vec_dou Sita_D;
		    Sita_D.assign(N0);
			
		 Vec_Dou x_D;
		         x_D.assign(N0);
		 Vec_Dou z_D;
		         z_D.assign(N0);


		int N1 = N;
		Vec_Dou ZZ;
		        ZZ.assign(N1);
		Vec_Dou XX;
		        XX.assign(N1);

		Vec_Dou EE_Total;
		        EE_Total.assign(N1*N1);

	  double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));
			

		   for (i=0;i<N1;i++)
		   {
			XX[i] = x[i];
			ZZ[i] = max_z+(i+1)*dx_Suf;
		   }


/*		   for (i=0;i<N1;i++)
		   {
			   for (j=0;j<N1;j++)
			   {
				    xsz=XX[i]*sin(i_an)-ZZ[j]*cos(i_an);
		            kg=pow(K1*g_p*cos(i_an),2.0);
		            fac=pow((XX[i]+ZZ[j]*tan(i_an))/g_p,2.0);
		            ww=(2.0*fac-1.0)/kg;
				   
				   for (k=0;k<N;k++)
				   {
				   	  rm = sqrt(pow(XX[i]-x[k],2.0)+pow(ZZ[j]-z[k],2.0));
					  NR = (-1.*df_x[j]*(XX[i]-x[k])+(ZZ[j]-z[k]))/(rm);
					 /* if (rm==0.)
					  {
						 dl = sqrt(1.+pow(df_1[k],2.0));

                         Sum0 = Sum0 + dx_Suf*((0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0))))*(W1_f[k]+W1_b[k])-J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)))*(U1_f[k]+U1_b[k]));
					  }

					  else*/
/*					  {	}				 						  
						  Sum0 = Sum0 + dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[k]+W1_b[k])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[k]+U1_b[k]);
					  

				   }                   
					EE_Total[i*N1+j] = abs(exp(J*K1*xsz*(1.0+ww))*exp(-fac) + Sum0);
					Sum0 = zero;
			     double r = sqrt(pow(XX[i],2.0)+pow(ZZ[j],2.0));
			     EE_Total[i*N1+j]= r*pow(EE_Total[i*N1+j],2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 EE_Total[i*N1+j]=10.*log10( EE_Total[i*N1+j]);
			   }

		   }*/


			
	
		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
		    x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);
		 }

	
		 
/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {					 
			   Es[i] = Es[i] + dx_Suf*((J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i])))*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			     
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

           
		 }*/

		   /* for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es[i] = Es[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

				 D_Es[i]=R_D*pow(abs(Es[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es[i]=10.*log10(D_Es[i]);
		
			}*/

			double abx   =  8.0*PI*K1*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-((1.0+2.0*pow(tan(i_an),2))/(2.0*pow(K1*g_p*cos(i_an),2))));

	Vec_Dou Sita_S; Sita_S.assign(N0);
	Complex Sum_E0(0.0,0.0);
	     for (i=0;i<N0;i++)
		 {
		    Sita_S[i] = i*PI/180 - PI/2;
		 }
		    for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					NR = K1*(df_x[j]*sin(Sita_S[i])-cos(Sita_S[i]));
				    double ksc=-K1*(x[j]*sin(Sita_S[i])+z[j]*cos(Sita_S[i]));
				    Sum_E0=Sum_E0+dx_Suf*(i_*NR*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(i_*ksc);
//				    Sum_E0=Sum_E0+dx_Suf*(NR*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(i_*ksc);
				}

                 Es[i] = Sum_E0;

//				 D_ES[i]  = (s*D_ES[i]+2*PI*norm(Sum_E0*exp(i_*K1*R_S))/abx)/(s+1);
				 D_Es[i]  = norm(Sum_E0)/abx;

                 DB_Es[i] = 10.*log10(D_Es[i]);
		         Sum_E0   = zero;
			}
         
			
/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es1[i] = Es1[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(Jf[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(Jf[j]);
				}

				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es1[i]=10.*log10(D_Es1[i]);
		
			}*/
	//	abs(Jf[i])<<<<abs(Jf[i+N])	<<abs(Jf[i+2*N])<<abs(Jf[i+3*N])

	ofstream out_U1;
	         out_U1.open("Layered_FBM_SAA_U1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<endl;
			 }
   
    ofstream out_W1;
	         out_W1.open("Layered_FBM_SAA_W1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<endl;
			 }

	ofstream out_U2;
	         out_U2.open("Layered_FBM_SAA_U2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<endl;
			 }
	ofstream out_W2;
	         out_W2.open("Layered_FBM_SAA_W2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<endl;
			 }



	ofstream out_Es;
	         out_Es.open("Layered_FBM_SAA_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }

	ofstream out_FN_1;
	         out_FN_1.open("lyaered_FN_f_1.txt");
			 for (i=-Q1;i<=Q1;i++)
			 {
				 out_FN_1<<i<<"   "<<abs(FN_f_1[i+Q1])<<endl;
			 }
	ofstream out_FN_2;
	         out_FN_2.open("lyaered_FN_f_2.txt");
			 for (i=-Q2;i<=Q2;i++)
			 {
				 out_FN_2<<i<<"   "<<abs(FN_f_2_1[i+Q2])<<endl;
			 }
	ofstream out_FN_3;
	         out_FN_3.open("lyaered_FN_f_3.txt");
			 for (i=-Q2;i<=Q2;i++)
			 {
				 out_FN_3<<i<<"   "<<abs(FN_f_3_1[i+Q2])<<endl;
			 }
			 	
	ofstream out_FN_4;
	         out_FN_4.open("lyaered_FN_f_4.txt");
			 for (i=-Q3;i<=Q3;i++)
			 {
				 out_FN_4<<i<<"   "<<abs(FN_f_4[i+Q3])<<endl;
			 }
			 
	ofstream out_Fai;
	         out_Fai.open("lyaered_Fai.txt");
			 for (i=-Q1;i<=Q1;i++)
			 {
				 out_Fai<<i<<"   "<<real(Fai[i+Q3])<<"   "<<imag(Fai[i+Q3])<<"   "<<abs(Fai[i+Q3])<<endl;
			 }



/*	ofstream out_Es;
	         out_Es.open("Layered_FBM_SAA_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }*/

}





void Lyaered_FBM_SAA::Gauss_Layered_Rough_surface_FBM_Bried_Target()
{
//////////////////2008.11.13///////////////////////
//////////////////////////////////////////////////////////////////////////
/////        如果媒质介电常数虚数部分过大则TM波收敛性更好           //////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/////                  2008.11.20                                   //////
/////                 加入了地下目标                               //////
//////////////////////////////////////////////////////////////////////////

///  若需要计算TM波入射情况。则需要将Z_Target Z_S_T_1 Z_S_T_2 改变。再将rou_2 rou_3改变///
////////////////////////////////姬伟杰////////////////////////////////////
	
	int NUMBER_REALIZATION=1;   //The number of realization.

	//一维导体粗糙面
	wave_len=3.0e8/(3e8);      //wave length
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

    g_p =10.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	N=(rl/wave_len)*10;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;


	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828183;

	double d=10.0*wave_len;
	
	double Sita_Target = 0.;
	double R_Target = 0.;
    R_Target=3.0*wave_len;
	H=-5*wave_len;

  
    double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta;


	Complex NR(0.,0.);
	
	Complex er2(2.0,0.05);         //如果应用迭代法求解，该介电常数虚部不能小于0.05,而矩量法则无此限制。
	Complex er3(25.0,0.0);
	Complex er_Target(2.25,0);

//	Complex rou_2 = er2;             //The TM wave
//	Complex rou_3 = er3/er2;  
//  Complex roru_Target = er_Target/er2;	

	double rou_2 = 1.;               //The TE wave
	double rou_3 = 1.;
	double rou_Target = 1.;

	double K1;
	Complex K2;
	Complex K3;
	Complex K_Target;


	Complex zero(0.,0.);

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);
	Complex Sum_E0(0.,0.);
	Complex Sum_E1(0.,0.);

	int i,j,k,l,m;

	

	Vec_Com A;A.assign(4*4);
	        
	Vec_Com X;X.assign(4);
	        
	Vec_Com V;V.assign(4);
	        
	
	Vec_Com Ei;Ei.assign(N);
	        

	Vec_Dou Z1;Z1.assign(N);
	        
	Vec_dou Z2;Z2.assign(N);
	        

	Vec_Dou df_1;df_1.assign(N);
	        
	Vec_dou df_2;df_2.assign(N);
	        
			
	Vec_Dou ddf_1;ddf_1.assign(N);
	        
	Vec_dou ddf_2;ddf_2.assign(N);
	        

	Vec_Com U1_f;U1_f.assign(N);
	        
	Vec_Com U1_b; U1_b.assign(N);
	       
	Vec_Com U2_f;U2_f.assign(N);
	        
	Vec_Com U2_b;U2_b.assign(N);
	        

	Vec_Com W1_f;W1_f.assign(N);
	        
	Vec_Com W1_b;W1_b.assign(N);
	        
	Vec_Com W2_f;W2_f.assign(N);
	        
	Vec_Com W2_b;W2_b.assign(N);
	        

	Vec_Com Jf1;Jf1.assign(4*N);
	        
	Vec_Com Jf2;Jf2.assign(4*N);
	        
	Vec_Com Jf3;Jf3.assign(4*N);
	        
	Vec_Com A_1_1_1;A_1_1_1.assign(N*N);
	        
	Vec_Com B_1_1_1;B_1_1_1.assign(N*N);
	        

			
	Vec_Com A_2_1_1; A_2_1_1.assign(N*N);
	       
	Vec_Com B_2_1_1;B_2_1_1.assign(N*N);
	        

	Vec_Com A_2_1_2;A_2_1_2.assign(N*N);
	        
	Vec_Com B_2_1_2;B_2_1_2.assign(N*N);
	        

	Vec_Com A_2_2_1; A_2_2_1.assign(N*N);
	       
	Vec_Com B_2_2_1;B_2_2_1.assign(N*N);
	        

	Vec_Com A_2_2_2;A_2_2_2.assign(N*N);
	        
	Vec_Com B_2_2_2; B_2_2_2.assign(N*N);
	       

	Vec_Com A_3_2_2;A_3_2_2.assign(N*N);
	        
	Vec_Com B_3_2_2;B_3_2_2.assign(N*N);
	        

	Vec_Com ZT_MOM; ZT_MOM.assign((4*N+2*M)*(4*N+2*M));
	       
	Vec_Com EE_MOM; EE_MOM.assign(4*N+2*M);
	       
	Vec_Com Jf_MOM; Jf_MOM.assign(4*N+2*M);
	       

	Vec_Com Z_Target; Z_Target.assign(M*M);
	       
	Vec_Com Z_T_S_A_1; Z_T_S_A_1.assign(M*N);
	       
	Vec_Com Z_T_S_B_1; Z_T_S_B_1.assign(M*N);
	       

			
	Vec_Com Z_T_S_A_2;Z_T_S_A_2.assign(M*N);
	        
	Vec_Com Z_T_S_B_2; Z_T_S_B_2.assign(M*N);
	       

	Vec_Com E_Target_Surf;E_Target_Surf.assign(M);
	        
	Vec_Com E_Surf1_Target;E_Surf1_Target.assign(N);
	        
	Vec_Com E_Surf2_Target; E_Surf2_Target.assign(N);
	       
	
	Vec_Com J_Target;J_Target.assign(M);
	       

	Vec_Com Z_S_T_1;Z_S_T_1.assign(N*M);
	        
	Vec_Com Z_S_T_2;Z_S_T_2.assign(N*M);
	        

    Vec_Com J_Target_1;J_Target_1.assign(M);
		        
    Vec_Com J_Target_2;J_Target_2.assign(M);

	Vec_Com Delta__J_Target; Delta__J_Target.assign(M);

	Vec_Dou Delta_Matrix; Delta_Matrix.assign(15);
		        
    double Delta_Target_Surface = 0.;
	
	int N0=180;
	double R_D = 10000.*wave_len;
	Vec_Com Es;		    Es.assign(N0);
	Vec_dou D_Es;		D_Es.assign(N0);
	Vec_dou DB_Es;		DB_Es.assign(N0);
	Vec_Com Es1;		Es1.assign(N0);	
	Vec_dou D_Es1;		D_Es1.assign(N0);
	Vec_dou DB_Es1;		DB_Es1.assign(N0);
	Vec_dou Sita_D;		Sita_D.assign(N0);
	Vec_Dou x_D;		x_D.assign(N0);
	Vec_Dou z_D;		z_D.assign(N0);

	Vec_Com ZT_I; ZT_I.assign(M*M);
	Vec_Com ZT_J; ZT_J.assign(M*M);
	Vec_Com ZT_K; ZT_K.assign(M*M);
	Vec_Com ZT_L; ZT_L.assign(M*M);

	
	int N1 = N;
	Vec_Dou ZZ;		ZZ.assign(N1);
    Vec_Dou XX;	    XX.assign(N1);
	Vec_Dou EE_Total; EE_Total.assign(N1*N1);
	double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));

	mean = 0.;
	vari = 1.0;
	RL = rl;
	k_w = 2*PI/wave_len;
	i_an = i_an*PI/180.;

	seed=123456;

	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);
	K_Target = K1*sqrt(er_Target);

			 
	for (i=0;i<N0;i++)
	{
       Sita_D[i] = i*PI/180-PI/2.;
	   x_D[i] = R_D*sin(Sita_D[i]);
	   z_D[i] = R_D*cos(Sita_D[i]);
	}

	
	for (int mm=0;mm<NUMBER_REALIZATION;mm++)
	{

		cout<<"Begin the "<<mm+1<<"th realization."<<endl;

//	 h=0.2*wave_len/(2.*PI);            //rms height 
//	 lc=1.0*wave_len/(2.*PI);           //correlation length 

	 h=0.5*wave_len;            //rms height 
	 lc=1.0*wave_len;           //correlation length 


	 gauss_surface();
	
//	One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	dx_Suf = RL/N;


/*    for (i=0;i<N;i++)
    {
		z[i] = 0.;
    }*/

	for (i=0;i<N;i++)
	{
		Z1[i] = z[i];
	}

	for (i=0;i<N;i++)
	{
		df_1[i] = df_x[i];		
		ddf_1[i]  = ddf_x[i];
		 
	}

//	 h=0.1*wave_len/(2.*PI);            //rms height 
//	 lc=1.0*wave_len/(2.*PI);           //correlation length 

	 h=0.5*wave_len;            //rms height 
	 lc=1.0*wave_len;           //correlation length 

	 gauss_surface();

	for (i=0;i<N;i++)
	{
		Z2[i] = z[i] - d;
	}

	for (i=0;i<N;i++)
	{
		df_2[i]= df_x[i];		
	    ddf_2[i] = ddf_x[i];		 
	}

/*	for (i=0;i<N;i++)
	{
		Z2[i] = -d;
		df_2[i] = 0.;
		ddf_2[i] = 0.;
	}*/

//	double Yita = 120.*PI;


	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-Z1[i]*cos(i_an);
		kg=pow(K1*g_p*cos(i_an),2.0);
		fac=pow((x[i]+Z1[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);	
	}

/////////////////////////////////////////////////////////////////////////////////////////


	cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}




   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}

//////////////////////////////////////2_1_1///////////////////////////////////////////////
   	
   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}



   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////2_1_2///////////////////////////////////////////////

   
   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
//			dl = sqrt(1.+pow(df_2[j],2.0));

		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
		    A_2_1_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z1[i]-Z2[j]))/rm;
				B_2_1_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_1///////////////////////////////////////////////


   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
//			dl = sqrt(1.+pow(df_1[i],2.0));

		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
		    A_2_2_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{

		      rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
			  NR = (df_1[j]*(x[i]-x[j])-(Z2[i]-Z1[j]))/rm;
			  B_2_2_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_2///////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////3_2_2///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 - ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}

	
//////////////////////////////////// The impedance matrix of the Target ///////////////////////////////////////////////
	

		
	cout<<"进行目标的剖分......------>>"<<endl;
	
	Lamda=wave_len; 
  
		

	ds=2.*PI*R_Target/M;
	Sita0=i_an;

	dx_t = 4*R/M;


	for(i=0;i<M;i++)
	{
		Sita_Target=i*2.*PI/M ;
		xtt[i]=R_Target*cos(Sita_Target) ;
		ztt[i]=R_Target*sin(Sita_Target)+H;		
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


	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
			if (i==j)
			{
				Z_Target[i*M+j] = ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K2/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));				
				Z_Target[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);
			}

		}
	}

/*	for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		            ZT_I[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);
					ZT_J[i*M+j] = J/4.*ds*Hankel(0,1,K2*rm);
				  }
		          else
				  {

					ZT_I[i*M+j] = 0.5-ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_J[i*M+j] = J/4.*ds*(1.+J*2./PI*log(exp(GAMA)*ds*K2/(4.*ee)));
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			       double  Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_K[i*M+j] = ds*NR*J*K_Target/4.*Hankel(1,1,K_Target*Rm);
					ZT_L[i*M+j] = -rou_Target*ds*J/4.*Hankel(0,1,K_Target*Rm);
				  }
		          else
				  {

					ZT_K[i*M+j] = 0.5+ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_L[i*M+j] = -rou_Target*ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K_Target/(4.*ee)));
                    
				  }
		   }
	   }*/




	/*     for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射情况
		            Z_Target[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);

				  }
		          else
				  {
					Z_Target[i*M+j] = -0.5;                    
				  }
		   }
	   }*/


//////////////////////////////////////////////////////////////////////////
	Vec_Com Z_S_T_3;  Z_S_T_3.assign(N*M);
	Vec_Com Z_S_T_4;  Z_S_T_4.assign(N*M);

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));		   		   
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_3[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);		   		   
		   Z_S_T_1[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);

		}
	}

		
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));		   		   
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_4[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);
		   Z_S_T_2[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);

		}
	}
	


/*	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));		   			        
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_3[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K2*rm);

		}
	}

		
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));		   			        
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_4[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K2*rm);
		}
	}*/

	


/*		for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射时情况
		   Z_S_T_1[i*M+j] = -1.*NR*J*ds*K2/4.*Hankel(1,1,K2*rm);
		}
	}

		for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射时情况
		   Z_S_T_2[i*M+j] = -1.*NR*J*ds*K2/4.*Hankel(1,1,K2*rm);
		}
	}*/



//////////////////////////////////////////////////////////////////////////


	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{		
				rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z1[j],2.0));
				NR = (df_1[j]*(xtt[i]-x[j])-(ztt[i]-Z1[j]))/rm;
				
				
				Z_T_S_A_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
				Z_T_S_B_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				
		}
	}

	
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
			    rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z2[j],2.0));
				NR = (df_2[j]*(xtt[i]-x[j])-(ztt[i]-Z2[j]))/rm;
				
				Z_T_S_A_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);				
				Z_T_S_B_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				

		}
	}




		

/////////////////////////////The impedance matrix for the MOM/////////////////////////////////////////////
/*   int N_M = 4*N+2*M;
    
   for (i=0;i<N;i++)
    {
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = A_1_1_1[i*N+j];
		}
        for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_1_1_1[i*N+j-N];
		}
	}

    for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_1_1[(i-N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_1[(i-N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_1_2[(i-N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_2[(i-N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_3[(i-N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_1[(i-N)*M+j-4*N-M];			
		}
	}

	for (i=2*N;i<3*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_2_1[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_1[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_2_2[(i-2*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_2[(i-2*N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_4[(i-2*N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_2[(i-2*N)*M+j-4*N-M];
		}

	}

	for (i=3*N;i<4*N;i++)
	{
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_3*A_3_2_2[(i-3*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_3_2_2[(i-3*N)*N+j-3*N];
		}
	}

	for (i=4*N;i<4*N+M;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*Z_T_S_A_1[(i-4*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_B_1[(i-4*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_A_2[(i-4*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_B_2[(i-4*N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_I[(i-4*N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_J[(i-4*N)*M+j-4*N-M];
		}
	}

	for (i=4*N+M;i<N_M;i++)
	{
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_K[(i-4*N-M)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_L[(i-4*N-M)*M+j-4*N-M];
		}
	}

	

	for (i=0;i<N;i++)
	{
		EE_MOM[i] = Ei[i];
	}

	cout<<"The matrix has been filled."<<endl;
	
   GAUSS(ZT_MOM,Jf_MOM,EE_MOM,N_M);*/






//	Bi_CG(ZT_MOM,EE_MOM,Jf_MOM);                     //The MOM method


//CLinearEquations linear;

for (m=0;m<15;m++)
{

cout<<"Begin the "<<m<<"th step of the iteration......------>>"<<endl;
 
 for (l=0;l<N;l++)
  {

	cout<<"Begin the "<<l+1<<"th step of the FBM......------>>"<<endl;

/////////////////////////The forward sweep/////////////////////////////////////////////////


    for (i=0;i<N;i++)
    {
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = Ei[i] - Sum0;
			V[1] = -E_Surf1_Target[i]-Sum1;
			V[2] = -E_Surf2_Target[i]-Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

		


			GAUSS(A,X,V,4);

		

			U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

    }


/////////////////////////The backward sweep/////////////////////////////////////////////////
	
	
   for (i=N-1;i>=0;i--)
    {
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;

			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

			
				
			GAUSS(A,X,V,4);


			U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];

    }

//////////////////////////The convergence condition////////////////////////////////////////////////

	
	for (i=0;i<N;i++)
	{
		Jf1[i]   = U1_f[i];
		Jf1[i+N] = W1_f[i];
		Jf1[i+2*N] = U2_f[i];
		Jf1[i+3*N] = W2_f[i];
	}

	for (i=0;i<4*N;i++)
	{
		Jf3[i] = Jf1[i] - Jf2[i];
	}

	  Assign_Vecor(Jf1,Jf2);  
	  Delta = sqrt(L2_Norm(Jf3))/sqrt(L2_Norm(Jf1));      
	  
	  cout<<"The error of the FBM is: Delta ="<<Delta<<endl;
			  
	  if(Delta<=0.01)
	  {
	    cout<<"The total number of the FBM step is :"<<l+1<<"."<<"and the final error is :"<<Delta<<endl;
		  break;
	  }
				
	     else  
		     continue;

}



	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
            Sum_E0  = Sum_E0 + rou_2*Z_T_S_A_1[i*N+j]*(U1_f[j] + U1_b[j]) + Z_T_S_B_1[i*N+j]*(W1_f[j]+W1_b[j]) + Z_T_S_A_2[i*N+j]*(U2_f[j] + U2_b[j]) + Z_T_S_B_2[i*N+j]*(W2_f[j] + W2_b[j]);
		}
		
           E_Target_Surf[i] = -Sum_E0;
		   Sum_E0 = zero;
	}


//     GAUSS(Z_Target,J_Target,E_Target_Surf,M);
	
//   Bi_CG(Z_Target,E_Target_Surf,J_Target); 


   for (i=0;i<N;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   Sum_E0 = Sum_E0 + Z_S_T_1[i*M+j]*J_Target[j];
		   Sum_E1 = Sum_E1 + Z_S_T_2[i*M+j]*J_Target[j];
	   }
	       E_Surf1_Target[i] = Sum_E0;
		   E_Surf2_Target[i] = Sum_E1;

		   Sum_E0 = zero;
		   Sum_E1 = zero;

   }



		for (i=0;i<M;i++)
		{
			J_Target_1[i] = J_Target[i] - J_Target_2[i]; 
		}

		 Assign_Vecor(J_Target,J_Target_2);

		 for (i=0;i<M;i++)
		 {
			 for (j=0;j<M;j++)
			 {
				 Sum_E0 = Sum_E0 + Z_Target[i*M+j]*J_Target_1[j];
			 }
			 Delta__J_Target[i] = Sum_E0;
			 Sum_E0 = zero;
             
		 }

		 Delta_Target_Surface = sqrt(L2_Norm(Delta__J_Target))/sqrt(L2_Norm(E_Target_Surf));  
		 
//		 Delta_Matrix[m] = Delta_Target_Surface;
	  
	     cout<<"The error of the alternate is: Delta_Target_Surface ="<<Delta_Target_Surface<<endl;
			  
	    if(Delta_Target_Surface<=0.01)
		{
	     cout<<"The total number of the alternate step is :"<<m+1<<"."<<"and the final error is :"<<Delta_Target_Surface<<endl;
		  break;
		}
				
	     else  
		     continue;


}
			
/*	   double max_z;
			for (i=0;i<N;i++)
			{
				if (max_z<z[i])
				{
					max_z = z[i];
				}
				else 
					continue;
			}

		   for (i=0;i<N1;i++)
		   {
			XX[i] = x[i];
			ZZ[i] = max_z+(i+1)*dx_Suf;
		   }


/*		   for (i=0;i<N1;i++)
		   {
			   for (j=0;j<N1;j++)
			   {
				    xsz=XX[i]*sin(i_an)-ZZ[j]*cos(i_an);
		            kg=pow(K1*g_p*cos(i_an),2.0);
		            fac=pow((XX[i]+ZZ[j]*tan(i_an))/g_p,2.0);
		            ww=(2.0*fac-1.0)/kg;
				   
				   for (k=0;k<N;k++)
				   {
				   	  rm = sqrt(pow(XX[i]-x[k],2.0)+pow(ZZ[j]-z[k],2.0));
					  NR = (-1.*df_x[j]*(XX[i]-x[k])+(ZZ[j]-z[k]))/(rm);
					 /* if (rm==0.)
					  {
						 dl = sqrt(1.+pow(df_1[k],2.0));

                         Sum0 = Sum0 + dx_Suf*((0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0))))*(W1_f[k]+W1_b[k])-J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)))*(U1_f[k]+U1_b[k]));
					  }

					  else*/
/*					  {	}				 						  
						  Sum0 = Sum0 + dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[k]+W1_b[k])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[k]+U1_b[k]);
					  

				   }                   
					EE_Total[i*N1+j] = abs(exp(J*K1*xsz*(1.0+ww))*exp(-fac) + Sum0);
					Sum0 = zero;
			     double r = sqrt(pow(XX[i],2.0)+pow(ZZ[j],2.0));
			     EE_Total[i*N1+j]= r*pow(EE_Total[i*N1+j],2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 EE_Total[i*N1+j]=10.*log10( EE_Total[i*N1+j]);
			   }

		   }*/


			
	

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
		    x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);
		 }

	
		 
/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {					 
			   Es[i] = Es[i] + dx_Suf*((J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i])))*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			     
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

           
		 }*/

/*		    for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es[i] = Es[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

		
			}*/

		double abx   =  8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-((1.0+2.0*pow(tan(i_an),2))/(2.0*pow(k_w*g_p*cos(i_an),2))));
		   
		for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-Z1[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
//					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					NR = k_w*(df_1[j]*sin(Sita_D[i])-cos(Sita_D[i]))/sqrt(1+pow(df_1[j],2.0));


				    double ksc=-k_w*(x[j]*sin(Sita_D[i] )+Z1[j]*cos(Sita_D[i] ));
				    Sum_E0=Sum_E0+dx_Suf*(NR*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(i_*ksc);
//				    Sum_E0=Sum_E0+dx_Suf*(NR*(Jf_MOM[j+N])-(Jf_MOM[j]))*exp(i_*ksc);

//					Sum_E0 = Sum_E0+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

                 Es[i] = Sum_E0;
				 
				 D_Es[i]  = (mm*D_Es[i]+norm(Sum_E0)/abx)/(mm+1);
                 DB_Es[i] = (mm*DB_Es[i]+10.*log10(D_Es[i]))/(mm+1);
		         Sum_E0   = zero;
			}

   	}      

/*		    for (i=0;i<N0;i++)
			{
				Es[i] = Es[i]/double(NUMBER_REALIZATION);
				D_Es[i]=R_D*pow(abs(Es[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
				DB_Es[i]=10.*log10(D_Es[i]);
			}

	        
			
			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es1[i] = Es1[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(Jf_MOM[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(Jf_MOM[j]);
				}

				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es1[i]=10.*log10(D_Es1[i]);
		
			}*/


		 
    

     

     


/*    ofstream out_A_1_1_1;
	          out_A_1_1_1.open("B_3_2_2.txt");
			  for (i=0;i<N;i++)
			  {
				  for (j=0;j<N;j++)
				  {  
				    out_A_1_1_1<<B_3_2_2[i*N+j]<<"   ";
					if (j%N==0)
					{
						out_A_1_1_1<<endl;
					}
				  }

			  }*/

	ofstream out_U1;
	         out_U1.open("Layered_Bried_Target_U1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<abs(Jf_MOM[i])<<endl;
			 }
   
    ofstream out_W1;
	         out_W1.open("Layered_Bried_Target_W1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<abs(Jf_MOM[i+N])<<endl;
			 }

	ofstream out_U2;
	         out_U2.open("Layered_Bried_Target_U2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<abs(Jf_MOM[i+2*N])<<endl;
			 }
	ofstream out_W2;
	         out_W2.open("Layered_Bried_Target_W2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<abs(Jf_MOM[i+3*N])<<endl;
			 }
	ofstream out_Es;
	         out_Es.open("Layered_Bried_Target_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	ofstream out_Es1;
	         out_Es1.open("Layered_Bried_Target_Es1.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<D_Es1[i]<<"   "<<DB_Es1[i]<<endl;
			 }

     ofstream out_J_Target;
	         out_J_Target.open("Layered_Bried_Target_J.txt");
		     for (i=0;i<M;i++)
		     {
				 out_J_Target<<i<<"   "<<real(J_Target[i])<<"   "<<imag(J_Target[i])<<"   "<<abs(J_Target[i])<<"  "<<abs(Jf_MOM[i+4*N])<<endl;
		     }
     ofstream out_J_Target_Z;
	         out_J_Target_Z.open("Layered_Bried_Target_Z.txt");
		     for (i=0;i<N;i++)
		     {
				 out_J_Target_Z<<i<<"   "<<Z1[i]<<"   "<<Z2[i]<<endl;
		     }

/*     ofstream out_Delta;
	         out_Delta.open("Layered_Bried_Target_Delta.txt");
		     for (i=0;i<15;i++)
		     {
				 out_Delta<<i<<"   "<<Delta_Matrix[i]<<endl;
		     }


/*	ofstream out_EE_Total;
	         out_EE_Total.open("EE_ToTal.txt");
			 for (i=0;i<N1;i++)
			 {
				 for (j=0;j<N1;j++)
				 {
					 out_EE_Total<<i<<"   "<<j<<"    "<<EE_Total[i*N1+j]<<endl;
				 }
			 }*/
}




void Lyaered_FBM_SAA::Gauss_Layered_Rough_surface_FBM_Above_Target()
{

//////////////////////////////////////////////////////////////////////////
/////                   2008.11.21                                  //////
/////                 加入了地上目标   ----TE                       //////
//////////////////////////////////////////////////////////////////////////

//   IF you want compute the TM wave,you must change the T_target    ///
///	       and Z_S_T_1.Then change the rou_2 and rou_3               ///
////////////////////////////////姬伟杰////////////////////////////////////
	

	//一维导体粗糙面
	double wave_len=3.0e8/3.0e8;      //wave length
	double h=0.2*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=20.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=102.4*wave_len;
	double R=1.*wave_len;
	       H=2.*wave_len; 
//	double g_p = 4.*R/cos(i_an*PI/180.);
           g_p = 10.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	int N=(rl/wave_len)*10;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;


	
	double eps=1.0e-2;        //CGM的收敛准则
	double criterion=1.0e-2;  //BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828183;

	double d=5.*wave_len;

	double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta;


	Complex NR(0.,0.);
	
	Complex er2(2.5,0.18);
	Complex er3(9.8,1.2);

//	Complex rou_2 = er2;             //The TM wave
//	Complex rou_3 = er3/er2;           

	double rou_2 = 1.;               //The TE wave
	double rou_3 = 1.;

	double K1;
	Complex K2;
	Complex K3;

	Complex zero(0.,0.);

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);
	Complex Sum_E0(0.,0.);
	Complex Sum_E1(0.,0.);

	int i,j,k,l,m;

	

	Vec_Com A;A.assign(4*4);	        
	Vec_Com X;X.assign(4);	        
	Vec_Com V; V.assign(4);
	       
	
	Vec_Com Ei; Ei.assign(N);
	       

	Vec_Dou Z1;Z1.assign(N);	        
	Vec_dou Z2;Z2.assign(N);
	        

	Vec_Dou df_1;df_1.assign(N);	        
	Vec_dou df_2; df_2.assign(N);
	       
			
	Vec_Dou ddf_1;ddf_1.assign(N);	        
	Vec_dou ddf_2;ddf_2.assign(N);
	        

	Vec_Com U1_f;U1_f.assign(N);	        
	Vec_Com U1_b;U1_b.assign(N);	        
	Vec_Com U2_f;U2_f.assign(N);	        
	Vec_Com U2_b;U2_b.assign(N);
	        

	Vec_Com W1_f;W1_f.assign(N);	        
	Vec_Com W1_b;W1_b.assign(N);	        
	Vec_Com W2_f;W2_f.assign(N);	        
	Vec_Com W2_b; W2_b.assign(N);
	       

	Vec_Com Jf1;Jf1.assign(4*N);	        
	Vec_Com Jf2;Jf2.assign(4*N);	        
	Vec_Com Jf3;Jf3.assign(4*N);
	        

	Vec_Com U1_f_1;U1_f_1.assign(N);	        
	Vec_Com U1_b_1;U1_b_1.assign(N);
	        
	Vec_Com W1_f_1;W1_f_1.assign(N);
	Vec_Com W1_b_1;W1_b_1.assign(N);
	        




	Vec_Com A_1_1_1;A_1_1_1.assign(N*N);	        
	Vec_Com B_1_1_1;B_1_1_1.assign(N*N);
	        

			
	Vec_Com A_2_1_1; A_2_1_1.assign(N*N);
	Vec_Com B_2_1_1;B_2_1_1.assign(N*N);
	        

	Vec_Com A_2_1_2;A_2_1_2.assign(N*N);	        
	Vec_Com B_2_1_2; B_2_1_2.assign(N*N);
	       

	Vec_Com A_2_2_1;A_2_2_1.assign(N*N);	        
	Vec_Com B_2_2_1;B_2_2_1.assign(N*N);
	        

	Vec_Com A_2_2_2; A_2_2_2.assign(N*N);	       
	Vec_Com B_2_2_2;B_2_2_2.assign(N*N);
	        

	Vec_Com A_3_2_2; A_3_2_2.assign(N*N);	       
	Vec_Com B_3_2_2;B_3_2_2.assign(N*N);
	        

	Vec_Com ZT; ZT.assign(4*N*4*N);
	       
	Vec_Com ZT_MOM; ZT_MOM.assign((4*N+M)*(4*N+M));	       
	Vec_Com EE;EE.assign(4*N);	        
	Vec_Com EE_M0M; EE_M0M.assign(4*N+M);	       
	Vec_Com Jf;Jf.assign(4*N);	        
	Vec_Com JJ_MOM; JJ_MOM.assign(4*N+M);
	       
	        

	Vec_Com Z_Target;Z_Target.assign(M*M);	        
	Vec_Com Z_T_S_A_1;Z_T_S_A_1.assign(M*N);	        
	Vec_Com Z_T_S_B_1; Z_T_S_B_1.assign(M*N);
	       			
	Vec_Com Z_T_S_A_2; Z_T_S_A_2.assign(M*N);	       
	Vec_Com Z_T_S_B_2; Z_T_S_B_2.assign(M*N);
	       

	Vec_Com E_Target_Surf;E_Target_Surf.assign(M);	        
	Vec_Com E_Surf1_Target; E_Surf1_Target.assign(N);	       
	Vec_Com E_Surf2_Target;  E_Surf2_Target.assign(N);
	      
	
	Vec_Com J_Target;J_Target.assign(M);
	       
	Vec_Com Z_S_T_1;Z_S_T_1.assign(N*M);        
	Vec_Com Z_S_T_2; Z_S_T_2.assign(N*M);	       
	Vec_Com Ei_T; Ei_T.assign(M);
	       

    Vec_Com J_Target_1;J_Target_1.assign(M);		        
    Vec_Com J_Target_2;J_Target_2.assign(M);
		        
    double Delta_Target_Surface = 0.;
	Vec_Dou Delta_Matrix;  Delta_Matrix.assign(15);
	Vec_Com Delta__J_Target; Delta__J_Target.assign(M);

	   	 int N0=181;
		 double R_D = 1000000.*wave_len;

		 Vec_Com Es;
		    Es.assign(N0);
		 Vec_dou D_Es;
		    D_Es.assign(N0);
    	 Vec_dou DB_Es;
		    DB_Es.assign(N0);
				 
		 Vec_Com Es1;
		    Es1.assign(N0);	
		 Vec_dou D_Es1;
		    D_Es1.assign(N0);
    	 Vec_dou DB_Es1;
		    DB_Es1.assign(N0);			
		
			
		 Vec_dou Sita_D;
		    Sita_D.assign(N0);
			
		 Vec_Dou x_D;
		         x_D.assign(N0);
		 Vec_Dou z_D;
		         z_D.assign(N0);


		int N1 = N;
		Vec_Dou ZZ;
		        ZZ.assign(N1);
		Vec_Dou XX;
		        XX.assign(N1);

		Vec_Dou EE_Total;
		        EE_Total.assign(N1*N1);


	    double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));
	i_an = i_an*PI/180.;

	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);


	int NUMBER_REALIZATION=1;

	Vec_dou BI_RCS; BI_RCS.assign(NUMBER_REALIZATION);

	for (int nn=0;nn<NUMBER_REALIZATION;nn++)
	{	
	   cout<<"Begin the "<<nn+1<<"th realization."<<endl;
	   
	   One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);
	   dx_Suf = RL/N;
/*	      {
			  z[i] = 0.;
			  df_x[i]=0.;
			  ddf_x[i]=0.;
		  }*/
		  for (i=0;i<N;i++)
		  {
			  Z1[i] = z[i];
			  Z2[i] = z[i] - d;
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  df_1[i] =df_2[i]= df_x[i];
			  ddf_1[i] =ddf_2[i] = ddf_x[i];
		  }
		  
/*		  for (i=0;i<N;i++)
		  {
			  Z2[i] = -d;
			  df_2[i]= 0;
			  ddf_2[i] = 0;
		  }*/
		  
		  //	double Yita = 120.*PI;+H*tan(i_an)H+*tan(i_an)
		  
		  double xsz,ww,kg,fac;
		  for(i=0;i<N;i++)
		  {
			  xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
			  kg=pow(k_w*g_p*cos(i_an),2.0);
			  fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
			  ww=(2.0*fac-1.0)/kg;
			  Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);
		  }
	
/////////////////////////////////////////////////////////////////////////////////////////
		  
		  cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  dl = sqrt(1.+pow(df_1[i],2.0));
					  A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
					  A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
				  }
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
					  NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
					  B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
				  }
			  }
		  }

//////////////////////////////////////2_1_1///////////////////////////////////////////////
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  dl = sqrt(1.+pow(df_1[i],2.0));
					  A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
					  A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
				  }
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
					  NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
					  B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				  }
			  }
		  }

//////////////////////////////////////2_1_2///////////////////////////////////////////////
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  
				  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
				  A_2_1_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
				  NR = (df_2[j]*(x[i]-x[j])-(Z1[i]-Z2[j]))/rm;
				  B_2_1_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			  }
		  }

//////////////////////////////////////2_2_1///////////////////////////////////////////////
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
				  A_2_2_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
				  NR = (df_1[j]*(x[i]-x[j])-(Z2[i]-Z1[j]))/rm;
				  B_2_2_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			  }
		  }

//////////////////////////////////////2_2_2///////////////////////////////////////////////
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  dl = sqrt(1.+pow(df_2[i],2.0));
					  A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
					  A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
				  }
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
					  NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
					  B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				  }
			  }
		  }

//////////////////////////////////////3_2_2///////////////////////////////////////////////
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  dl = sqrt(1.+pow(df_2[i],2.0));
					  A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
					  A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
				  }
			  }
		  }
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<N;j++)
			  {
				  if (i==j)
				  {
					  B_3_2_2[i*N+j] = 0.5 - ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
				  }
				  else
				  {
					  rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
					  NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
					  B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
				  }
			  }
		  }
	
//////////////////////////////////// The impedance matrix of the Target ///////////////////////////////////////////////
		  double Sita_Target = 0.;
		  double R_Target = 0.;
		  
		  cout<<"进行目标的剖分......------>>"<<endl;
		  
		  Lamda=wave_len; 
		  R_Target=R;
		  Yita=120.*PI;	
		  Gama=1.781072418;
		  ds=2.*PI*R_Target/M;
		  Sita0=i_an;
		  dx_t = 4*R/M;
		  
		  for(i=0;i<M;i++)
		  {
			  Sita_Target=i*2.*PI/M ;
			  xtt[i]=R_Target*cos(Sita_Target) ;
			  ztt[i]=R_Target*sin(Sita_Target)+H;
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
		  
		  for(i=0;i<M;i++)
		  {
			  xsz=xtt[i]*sin(i_an)-ztt[i]*cos(i_an);
			  Ei_T[i]=exp(J*K1*xsz);		  
		  }
	
	
/*	for( i=0;i<=M-1;i++)
	{
		for( j=0;j<=M-1;j++)
		{
			if(i!=j)
			{
				double kk=K1*sqrt((xtt[i]-xtt[j])*(xtt[i]-xtt[j])+(ztt[i]-ztt[j])*(ztt[i]-ztt[j]));
				Z_Target[i*M+j]=(K1*Yita/4.)*ds*Hankel(0,1,kk);
			}
			else
			{
                Z_Target[i*M+j]=(K1*Yita/4.)*ds*(1.+2.*J/PI*(log(Gama*K1*ds/4.)-1.));	
			}
		}
	}*/

/*	for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射情况
		            Z_Target[i*M+j] = -1.*NR*ds*J*K1/4.*Hankel(1,1,K1*rm);

				  }
		          else
				  {
					Z_Target[i*M+j] = -0.5;
                    
				  }
		   }
	   }*/
		  for (i=0;i<M;i++)
		  {
			  for (j=0;j<M;j++)
			  {
				  if (i==j)
				  {
					  Z_Target[i*M+j] = ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K1/(4.*ee)));    //TE波入射情况
				  }
				  else
				  {
					  rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
					  Z_Target[i*M+j] = ds*J/4.*Hankel(0,1,K1*rm);
				  }
			  }
		  }

//////////////////////////////////////////////////////////////////////////
		  
		  for (i=0;i<N;i++)
		  {
			  for (j=0;j<M;j++)
			  {
				  rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));                      //TE波入射情况
				  Z_S_T_1[i*M+j] = ds*J/4.*Hankel(0,1,K1*rm);
			  }
		  }

/*	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射时情况
		   Z_S_T_1[i*M+j] = -1.*NR*ds*J*K1/4.*Hankel(1,1,K1*rm);
		}
	}*/


		
/*	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));
		   Z_S_T_2[i*M+j] = J/4.*Hankel(0,1,K2*rm);

		}
	}*/

//////////////////////////////////////////////////////////////////////////
		  for (i=0;i<M;i++)
		  {
			  for (j=0;j<N;j++)
			  {	
				  rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z1[j],2.0));
				  NR = (df_1[j]*(xtt[i]-x[j])-(ztt[i]-Z1[j]))/rm;
				  
				  Z_T_S_A_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
				  Z_T_S_B_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			  }
		  }

	
/*	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
			    rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z2[j],2.0));
				NR = (df_2[j]*(xtt[i]-x[j])-(ztt[i]-Z2[j]))/rm;
				
				Z_T_S_A_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);				
				Z_T_S_B_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				

		}
	}*/

/////////////////////////////The impedance matrix for the MOM/////////////////////////////////////////////

/*	int N_M = 4*N+M;
	   
	for (i=0;i<N;i++)
    {
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = A_1_1_1[i*N+j];
		}
        for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_1_1_1[i*N+j-N];
		}
		
	}


	for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_1_1[(i-N)*N+j];

		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_1[(i-N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_1_2[(i-N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_2[(i-N)*N+j-3*N];
		}
	}

	for (i=2*N;i<3*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_2_1[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_1[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_2_2[(i-2*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_2[(i-2*N)*N+j-3*N];
		}
	}

	for (i=3*N;i<4*N;i++)
	{
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_3*A_3_2_2[(i-3*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_3_2_2[(i-3*N)*N+j-3*N];
		}
	}


		
	for (i=0;i<N;i++)
	{
		for (j=4*N;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_1[i*M+j-4*N];
		}
	}
	
	for (i=4*N;i<N_M;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_A_1[(i-4*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_B_1[(i-4*N)*N+j-N];
		}
		for (j=4*N;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_Target[(i-4*N)*M+j-4*N];
		}
		
	}

		
	for (i=0;i<N;i++)
	{
		EE_M0M[i] = Ei[i];
        
	}

	for (i=4*N;i<N_M;i++)
	{
		EE_M0M[i] = Ei_T[i-4*N];
	}



  GAUSS(ZT_MOM,JJ_MOM,EE_M0M,N_M);  */

	cout<<"The matrix has been filled."<<endl;
	for (m=0;m<15;m++)
	{
		cout<<"Begin the "<<m<<"th step of the iteration......------>>"<<endl;
		for (l=0;l<N;l++)
		{
			
			cout<<"Begin the "<<l+1<<"th step of the FBM......------>>"<<endl;

/////////////////////////The forward sweep/////////////////////////////////////////////////
			
			for (i=0;i<N;i++)
			{
				for (j=0;j<i;j++)
				{
					Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
					Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			        Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			        Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
				}
				
				V[0] = Ei[i] - Sum0 - E_Surf1_Target[i];
				V[1] = -Sum1;
				V[2] = -Sum2;
				V[3] = -Sum3;
				
				Sum0 = zero;
				Sum1 = zero;
				Sum2 = zero;
				Sum3 = zero;
				
				A[0] = A_1_1_1[i*N+i];
				A[1] = B_1_1_1[i*N+i];
				A[4] = rou_2*A_2_1_1[i*N+i];
				A[5] = B_2_1_1[i*N+i];
				A[6] = A_2_1_2[i*N+i];
				A[7] = B_2_1_2[i*N+i];
				A[8] = rou_2*A_2_2_1[i*N+i];
				A[9] = B_2_2_1[i*N+i];
				A[10] = A_2_2_2[i*N+i];
				A[11] = B_2_2_2[i*N+i];
				A[14] = rou_3*A_3_2_2[i*N+i];
				A[15] = B_3_2_2[i*N+i];
				
				GAUSS(A,X,V,4);
				
				U1_f[i] = X[0];
				W1_f[i] = X[1];
				U2_f[i] = X[2];
				W2_f[i] = X[3];
			}


/////////////////////////The backward sweep/////////////////////////////////////////////////
			
			for (i=N-1;i>=0;i--)
			{
				for (j=i+1;j<N;j++)
				{
					Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			        Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			        Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		   	        Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
				}
				
				V[0] = -Sum0;
				V[1] = -Sum1;
				V[2] = -Sum2;
				V[3] = -Sum3;
				
				Sum0 = zero;
				Sum1 = zero;
				Sum2 = zero;
				Sum3 = zero;
				
				A[0] = A_1_1_1[i*N+i];
				A[1] = B_1_1_1[i*N+i];
				A[4] = rou_2*A_2_1_1[i*N+i];
				A[5] = B_2_1_1[i*N+i];
				A[6] = A_2_1_2[i*N+i];
				A[7] = B_2_1_2[i*N+i];
				A[8] = rou_2*A_2_2_1[i*N+i];
				A[9] = B_2_2_1[i*N+i];
				A[10] = A_2_2_2[i*N+i];
				A[11] = B_2_2_2[i*N+i];
				A[14] = rou_3*A_3_2_2[i*N+i];
				A[15] = B_3_2_2[i*N+i];
				
				GAUSS(A,X,V,4);
				
				U1_b[i] = X[0];
				W1_b[i] = X[1];
				U2_b[i] = X[2];
				W2_b[i] = X[3];
			
			}

//////////////////////////The convergence condition////////////////////////////////////////////////
			
			for (i=0;i<N;i++)
			{
				Jf1[i]   = U1_f[i];
				Jf1[i+N] = W1_f[i];
				Jf1[i+2*N] = U2_f[i];
				Jf1[i+3*N] = W2_f[i];
			}
			
			for (i=0;i<4*N;i++)
			{
				Jf3[i] = Jf1[i] - Jf2[i];
			}
			
			Assign_Vecor(Jf1,Jf2);
			Delta = sqrt(L2_Norm(Jf3))/sqrt(L2_Norm(Jf1));
			cout<<"The error of the FBM is: Delta ="<<Delta<<endl;
			
			if(Delta<=0.01)
			{
				cout<<"The total number of the FBM step is :"<<l+1<<"."<<"and the final error is :"<<Delta<<endl;
				break;
			}
			else  
				continue;
}

  if (m==0)
  {
	  for (i=0;i<N;i++)
	  {
		  U1_f_1[i] = U1_f[i];
          U1_b_1[i] = U1_b[i];
		  W1_f_1[i] = W1_f[i];
		  W1_b_1[i] = W1_b[i];
	  }
  }

   
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
            Sum_E0  = Sum_E0 + Z_T_S_A_1[i*N+j]*(U1_f[j] + U1_b[j]) + Z_T_S_B_1[i*N+j]*(W1_f[j]+W1_b[j]);
		}
		
           E_Target_Surf[i] = Ei_T[i]-Sum_E0;
		   Sum_E0 = zero;
	}

    GAUSS(Z_Target,J_Target,E_Target_Surf,M);

	
//   Bi_CG(Z_Target,E_Target_Surf,J_Target); 


   for (i=0;i<N;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   Sum_E0 = Sum_E0 + Z_S_T_1[i*M+j]*J_Target[j];

	   }
	       E_Surf1_Target[i] = Sum_E0;

		   Sum_E0 = zero;
   }

		for (i=0;i<M;i++)
		{
			J_Target_1[i] = J_Target[i] - J_Target_2[i]; 
		}

		 Assign_Vecor(J_Target,J_Target_2);
		 for (i=0;i<M;i++)
		 {
			 for (j=0;j<M;j++)
			 {
				 Sum_E0 = Sum_E0 + Z_Target[i*M+j]*J_Target_1[j];
			 }
			 Delta__J_Target[i] = Sum_E0;
			 Sum_E0 = zero;             
		 }

		 Delta_Target_Surface = sqrt(L2_Norm(Delta__J_Target))/sqrt(L2_Norm(E_Target_Surf));
		 Delta_Matrix[m] = Delta_Target_Surface;
	  
	     cout<<"The error of the alternate is: Delta_Target_Surface ="<<Delta_Target_Surface<<endl;
			  
	     if(Delta_Target_Surface<=0.01)
		 {
	       cout<<"The total number of the alternate step is :"<<m+1<<"."<<"and the final error is :"<<Delta_Target_Surface<<endl;
		   break;
		 }				
	     else  
		     continue;
}


 



			
/*	   double max_z;
			for (i=0;i<N;i++)
			{
				if (max_z<z[i])
				{
					max_z = z[i];
				}
				else 
					continue;
			}

		   for (i=0;i<N1;i++)
		   {
			XX[i] = x[i];
			ZZ[i] = max_z+(i+1)*dx_Suf;
		   }


/*		   for (i=0;i<N1;i++)
		   {
			   for (j=0;j<N1;j++)
			   {
				    xsz=XX[i]*sin(i_an)-ZZ[j]*cos(i_an);
		            kg=pow(K1*g_p*cos(i_an),2.0);
		            fac=pow((XX[i]+ZZ[j]*tan(i_an))/g_p,2.0);
		            ww=(2.0*fac-1.0)/kg;
				   
				   for (k=0;k<N;k++)
				   {
				   	  rm = sqrt(pow(XX[i]-x[k],2.0)+pow(ZZ[j]-z[k],2.0));
					  NR = (-1.*df_x[j]*(XX[i]-x[k])+(ZZ[j]-z[k]))/(rm);
					 /* if (rm==0.)
					  {
						 dl = sqrt(1.+pow(df_1[k],2.0));

                         Sum0 = Sum0 + dx_Suf*((0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0))))*(W1_f[k]+W1_b[k])-J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)))*(U1_f[k]+U1_b[k]));
					  }

					  else*/
/*					  {	}				 						  
						  Sum0 = Sum0 + dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[k]+W1_b[k])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[k]+U1_b[k]);
					  

				   }                   
					EE_Total[i*N1+j] = abs(exp(J*K1*xsz*(1.0+ww))*exp(-fac) + Sum0);
					Sum0 = zero;
			     double r = sqrt(pow(XX[i],2.0)+pow(ZZ[j],2.0));
			     EE_Total[i*N1+j]= r*pow(EE_Total[i*N1+j],2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 EE_Total[i*N1+j]=10.*log10( EE_Total[i*N1+j]);
			   }

		   }*/


			
	

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
		    x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);
		 }

	
		 
/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {					 
			   Es[i] = Es[i] + dx_Suf*((J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i])))*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			     
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

           
		 }*/


		    for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es[i] = Es[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

				for (j=0;j<M;j++)
				{
					rm = sqrt(pow(x_D[i]-xtt[j],2.0)+pow(z_D[i]-ztt[j],2.0));
                    NR = ((x_D[i]-xtt[j])*(-1.*df_z[j])+(z_D[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
					Es[i] = Es[i] - ds*J/4.*Hankel(0,1,K1*rm)*J_Target[j];
				}

			}
    	}
		
		    for (i=0;i<N0;i++)
			{
								    		

				Es[i] = Es[i]/double(NUMBER_REALIZATION);
                D_Es[i] = R_D*pow(abs(Es[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                DB_Es[i]=10.*log10(D_Es[i]);
			}
				
/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es1[i] = Es1[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(JJ_MOM[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(JJ_MOM[j]);
				}
								
				for (j=0;j<M;j++)
				{
					rm = sqrt(pow(x_D[i]-xtt[j],2.0)+pow(z_D[i]-ztt[j],2.0));
					
					                    
					NR = ((x_D[i]-xtt[j])*(-1.*df_z[j])+(z_D[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
//					Es[i] = Es[i] - ds*J/4.*Hankel(0,1,K1*rm)*J_Target[j];
//					Es[i] = Es[i] - (K1*Yita/4.)*ds*Hankel(0,1,K1*rm)*J_Target[j];/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
//					Es1[i] = Es1[i] + NR*ds*J*K1/4.*Hankel(1,1,K1*rm)*JJ_MOM[j+4*N];

					Es1[i] = Es1[i] - ds*J/4.*Hankel(0,1,K1*rm)*JJ_MOM[j+4*N];
				}

				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es1[i]=10.*log10(D_Es1[i]);
		
			}*/

	ofstream out_U1;out_U1.open("Layered_Above_Target_U1.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<abs(JJ_MOM[i])<<endl;
			 }
   
    ofstream out_W1;out_W1.open("Layered_Above_Target_W1.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<abs(JJ_MOM[i+N])<<endl;
			 }

	ofstream out_U2;out_U2.open("Layered_Above_Target_U2.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<abs(JJ_MOM[i+2*N])<<endl;
			 }
	ofstream out_W2;out_W2.open("Layered_Above_Target_W2.txt");	         
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<abs(JJ_MOM[i+3*N])<<endl;
			 }
	ofstream out_Es;out_Es.open("Layered_Above_Target_Es.txt");	         
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	ofstream out_Es1;out_Es1.open("Layered_Above_Target_Es1.txt");	         
			 for (i=0;i<N0;i++)
			 {
				 out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<D_Es1[i]<<"   "<<DB_Es1[i]<<endl;
			 }

     ofstream out_J_Target;out_J_Target.open("Layered_Above_Target_J.txt");	         
		     for (i=0;i<M;i++)
		     {
 				 out_J_Target<<i<<"   "<<real(J_Target[i])<<"   "<<imag(J_Target[i])<<"   "<<abs(J_Target[i])<<endl;
		     }
	 ofstream out_JJ_Target_MOM;out_JJ_Target_MOM.open("Layered_Above_Target_J_MOM.txt");	         
		     for (i=0;i<M;i++)
		     {
                 out_JJ_Target_MOM<<i<<"   "<<real(JJ_MOM[i+4*N])<<"   "<<imag(JJ_MOM[i+4*N])<<"   "<<abs(JJ_MOM[i+4*N])<<endl;
		     }

/*	 ofstream out_Delta;out_Delta.open("Layered_Above_Target_Delta.txt");	         
		     for (i=0;i<15;i++)
		     {
                 out_Delta<<i<<"   "<<Delta_Matrix[i]<<endl;
		     }*/
		 
/*	ofstream out_EE_Total;
	         out_EE_Total.open("EE_ToTal.txt");
			 for (i=0;i<N1;i++)
			 {
				 for (j=0;j<N1;j++)
				 {
					 out_EE_Total<<i<<"   "<<j<<"    "<<EE_Total[i*N1+j]<<endl;
				 }
			 }*/

}

void Lyaered_FBM_SAA::Gauss_Layered_Rough_surface_FBM_Bried_conduct_Target()
{
//////////////////2008.11.13///////////////////////
//////////////////////////////////////////////////////////////////////////
/////        如果媒质介电常数虚数部分过大则TM波收敛性更好           //////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/////                  2008.11.20                                   //////
/////                 加入了地下目标                               //////
//////////////////////////////////////////////////////////////////////////

///  若需要计算TM波入射情况。则需要将Z_Target Z_S_T_1 Z_S_T_2 改变。再将rou_2 rou_3改变///
////////////////////////////////姬伟杰////////////////////////////////////
	
	int NUMBER_REALIZATION=1;   //The number of realization.

	//一维导体粗糙面
	wave_len=3.0e8/(3e8);      //wave length
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

    double g_p =10.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	N=(rl/wave_len)*10;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;


	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828183;

	double d=5.0*wave_len;
	
	double Sita_Target = 0.;
	double R_Target = 0.;
    R_Target=1.0*wave_len;
	H=-1.5*wave_len;

  
    double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta;


	Complex NR(0.,0.);
	
	Complex er2(2.5,0.18);         //如果应用迭代法求解，该介电常数虚部不能小于0.05,而矩量法则无此限制。
	Complex er3(9.8,1.2);
	Complex er_Target(2.25,0);

//	Complex rou_2 = er2;             //The TM wave
//	Complex rou_3 = er3/er2;  
//  Complex roru_Target = er_Target/er2;	

	double rou_2 = 1.;               //The TE wave
	double rou_3 = 1.;
	double rou_Target = 1.;

	double K1;
	Complex K2;
	Complex K3;
	Complex K_Target;


	Complex zero(0.,0.);

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);
	Complex Sum_E0(0.,0.);
	Complex Sum_E1(0.,0.);

	int i,j,k,l,m;

	

	Vec_Com A;A.assign(4*4);
	        
	Vec_Com X;X.assign(4);
	        
	Vec_Com V;V.assign(4);
	        
	
	Vec_Com Ei;Ei.assign(N);
	        

	Vec_Dou Z1;Z1.assign(N);
	        
	Vec_dou Z2;Z2.assign(N);
	        

	Vec_Dou df_1;df_1.assign(N);
	        
	Vec_dou df_2;df_2.assign(N);
	        
			
	Vec_Dou ddf_1;ddf_1.assign(N);
	        
	Vec_dou ddf_2;ddf_2.assign(N);
	        

	Vec_Com U1_f;U1_f.assign(N);
	        
	Vec_Com U1_b; U1_b.assign(N);
	       
	Vec_Com U2_f;U2_f.assign(N);
	        
	Vec_Com U2_b;U2_b.assign(N);
	        

	Vec_Com W1_f;W1_f.assign(N);
	        
	Vec_Com W1_b;W1_b.assign(N);
	        
	Vec_Com W2_f;W2_f.assign(N);
	        
	Vec_Com W2_b;W2_b.assign(N);
	        

	Vec_Com Jf1;Jf1.assign(4*N);
	        
	Vec_Com Jf2;Jf2.assign(4*N);
	        
	Vec_Com Jf3;Jf3.assign(4*N);
	        
	Vec_Com A_1_1_1;A_1_1_1.assign(N*N);
	        
	Vec_Com B_1_1_1;B_1_1_1.assign(N*N);
	        

			
	Vec_Com A_2_1_1; A_2_1_1.assign(N*N);
	       
	Vec_Com B_2_1_1;B_2_1_1.assign(N*N);
	        

	Vec_Com A_2_1_2;A_2_1_2.assign(N*N);
	        
	Vec_Com B_2_1_2;B_2_1_2.assign(N*N);
	        

	Vec_Com A_2_2_1; A_2_2_1.assign(N*N);
	       
	Vec_Com B_2_2_1;B_2_2_1.assign(N*N);
	        

	Vec_Com A_2_2_2;A_2_2_2.assign(N*N);
	        
	Vec_Com B_2_2_2; B_2_2_2.assign(N*N);
	       

	Vec_Com A_3_2_2;A_3_2_2.assign(N*N);
	        
	Vec_Com B_3_2_2;B_3_2_2.assign(N*N);
	        

	Vec_Com ZT_MOM; ZT_MOM.assign((4*N+2*M)*(4*N+2*M));
	       
	Vec_Com EE_MOM; EE_MOM.assign(4*N+2*M);
	       
	Vec_Com Jf_MOM; Jf_MOM.assign(4*N+2*M);
	       

	Vec_Com Z_Target; Z_Target.assign(M*M);
	       
	Vec_Com Z_T_S_A_1; Z_T_S_A_1.assign(M*N);
	       
	Vec_Com Z_T_S_B_1; Z_T_S_B_1.assign(M*N);
	       

			
	Vec_Com Z_T_S_A_2;Z_T_S_A_2.assign(M*N);
	        
	Vec_Com Z_T_S_B_2; Z_T_S_B_2.assign(M*N);
	       

	Vec_Com E_Target_Surf;E_Target_Surf.assign(M);
	        
	Vec_Com E_Surf1_Target;E_Surf1_Target.assign(N);
	        
	Vec_Com E_Surf2_Target; E_Surf2_Target.assign(N);
	       
	
	Vec_Com J_Target;J_Target.assign(M);
	       

	Vec_Com Z_S_T_1;Z_S_T_1.assign(N*M);
	        
	Vec_Com Z_S_T_2;Z_S_T_2.assign(N*M);
	        

    Vec_Com J_Target_1;J_Target_1.assign(M);
		        
    Vec_Com J_Target_2;J_Target_2.assign(M);

	Vec_Com Delta__J_Target; Delta__J_Target.assign(M);

	Vec_Dou Delta_Matrix; Delta_Matrix.assign(15);
		        
    double Delta_Target_Surface = 0.;
	
	int N0=180;
	double R_D = 10000.*wave_len;
	Vec_Com Es;		    Es.assign(N0);
	Vec_dou D_Es;		D_Es.assign(N0);
	Vec_dou DB_Es;		DB_Es.assign(N0);
	Vec_Com Es1;		Es1.assign(N0);	
	Vec_dou D_Es1;		D_Es1.assign(N0);
	Vec_dou DB_Es1;		DB_Es1.assign(N0);
	Vec_dou Sita_D;		Sita_D.assign(N0);
	Vec_Dou x_D;		x_D.assign(N0);
	Vec_Dou z_D;		z_D.assign(N0);

	Vec_Com ZT_I; ZT_I.assign(M*M);
	Vec_Com ZT_J; ZT_J.assign(M*M);
	Vec_Com ZT_K; ZT_K.assign(M*M);
	Vec_Com ZT_L; ZT_L.assign(M*M);

	
	int N1 = N;
	Vec_Dou ZZ;		ZZ.assign(N1);
    Vec_Dou XX;	    XX.assign(N1);
	Vec_Dou EE_Total; EE_Total.assign(N1*N1);
	double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));

	mean = 0.;
	vari = 1.0;
	RL = rl;
	seed=123456;
	k_w = 2*PI/wave_len;
	i_an = i_an*PI/180.;

	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);
	K_Target = K1*sqrt(er_Target);

			 
	for (i=0;i<N0;i++)
	{
        Sita_D[i] = i*PI/180-PI/2.;
		x_D[i] = R_D*sin(Sita_D[i]);
	    z_D[i] = R_D*cos(Sita_D[i]);
	}

	
	for (int mm=0;mm<NUMBER_REALIZATION;mm++)
	{

		cout<<"Begin the "<<mm+1<<"th realization."<<endl;

//	 h=0.2/K1;            //rms height 
//	 lc=1.0/K1;           //correlation length 

	 h=0.02*wave_len;            //rms height 
	 lc=1.0*wave_len;           //correlation length 


	 gauss_surface();
	
//	One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	dx_Suf = RL/N;


/*    for (i=0;i<N;i++)
    {
		z[i] = 0.;
    }*/

	for (i=0;i<N;i++)
	{
		Z1[i] = z[i];
	}

	for (i=0;i<N;i++)
	{
		df_1[i] = df_x[i];
		
		ddf_1[i]  = ddf_x[i];
		 
	}

//	 h=0.1/K1;            //rms height 
//	 lc=1.0/K1;           //correlation length 


	 h=0.02*wave_len;            //rms height 
	 lc=1.0*wave_len;           //correlation length 


	 gauss_surface();

	for (i=0;i<N;i++)
	{
		Z2[i] = z[i] - d;
	}

	for (i=0;i<N;i++)
	{
		df_2[i]= df_x[i];		
	    ddf_2[i] = ddf_x[i];		 
	}

	for (i=0;i<N;i++)
	{
		Z2[i] = -d;
		df_2[i] = 0.;
		ddf_2[i] = 0.;
	}

	double Yita = 120.*PI;

	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-Z1[i]*cos(i_an);
		kg=pow(K1*g_p*cos(i_an),2.0);
		fac=pow((x[i]+Z1[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);	
	}

/////////////////////////////////////////////////////////////////////////////////////////


	cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K1*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_1_1_1[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K1*rm);
			}
		}
	}




   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_1_1_1[i*N+j] = 0.5 - ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_1_1_1[i*N+j] = -J*(K1*dx_Suf/4.)*NR*Hankel(1,1,K1*rm);
			}
		}
	}

//////////////////////////////////////2_1_1///////////////////////////////////////////////
   	
   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_1[i],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				A_2_1_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}



   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_1_1[i*N+j] = 0.5 + ddf_1[j]*dx_Suf/(4.*PI*(1.+pow(df_1[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z1[j],2.0));
				NR = (df_1[j]*(x[i]-x[j])-(Z1[i]-Z1[j]))/rm;
				B_2_1_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////2_1_2///////////////////////////////////////////////

   
   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
		    A_2_1_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z1[i]-Z2[j],2.0));
			NR = (df_2[j]*(x[i]-x[j])-(Z1[i]-Z2[j]))/rm;
			B_2_1_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_1///////////////////////////////////////////////


   	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		    rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
		    A_2_2_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);			
		}
	}

   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		      rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z1[j],2.0));
			  NR = (df_1[j]*(x[i]-x[j])-(Z2[i]-Z1[j]))/rm;
			  B_2_2_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			
		}
	}

//////////////////////////////////////2_2_2///////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K2*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_2_2_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_2_2_2[i*N+j] = 0.5 - ddf_2[j]*dx_Suf/(4.*PI*(1.+pow(df_2[j],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_2_2_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
			}
		}
	}

//////////////////////////////////////3_2_2///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			
			if (i==j)
			{
				dl = sqrt(1.+pow(df_2[i],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*(1.+J*2./PI*log(exp(GAMA)*K3*dx_Suf*dl/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				A_3_2_2[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K3*rm);
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (i==j)
			{
				B_3_2_2[i*N+j] = 0.5 + ddf_2[i]*dx_Suf/(4.*PI*(1.+pow(df_2[i],2.0)));
			}
			else
			{
				rm = sqrt(pow(x[i]-x[j],2.0)+pow(Z2[i]-Z2[j],2.0));
				NR = (df_2[j]*(x[i]-x[j])-(Z2[i]-Z2[j]))/rm;
				B_3_2_2[i*N+j] = J*(K3*dx_Suf/4.)*NR*Hankel(1,1,K3*rm);
			}
		}
	}
	
//////////////////////////////////// The impedance matrix of the Target ///////////////////////////////////////////////
	

		
	cout<<"进行目标的剖分......------>>"<<endl;
	
	Lamda=wave_len; 
  
		

	ds=2.*PI*R_Target/M;
	Sita0=i_an;

	dx_t = 4*R/M;


	for(i=0;i<M;i++)
	{
		Sita_Target=i*2.*PI/M ;
		xtt[i]=R_Target*cos(Sita_Target) ;
		ztt[i]=R_Target*sin(Sita_Target)+H;		
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


	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
			if (i==j)
			{
				Z_Target[i*M+j] = ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K2/(4.*ee)));
			}
			else
			{
				rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));				
				Z_Target[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);
			}

		}
	}

/*	for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		            ZT_I[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);
					ZT_J[i*M+j] = J/4.*ds*Hankel(0,1,K2*rm);
				  }
		          else
				  {

					ZT_I[i*M+j] = 0.5-ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_J[i*M+j] = J/4.*ds*(1.+J*2./PI*log(exp(GAMA)*ds*K2/(4.*ee)));
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			       double  Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_K[i*M+j] = ds*NR*J*K_Target/4.*Hankel(1,1,K_Target*Rm);
					ZT_L[i*M+j] = -rou_Target*ds*J/4.*Hankel(0,1,K_Target*Rm);
				  }
		          else
				  {

					ZT_K[i*M+j] = 0.5+ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_L[i*M+j] = -rou_Target*ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K_Target/(4.*ee)));
                    
				  }
		   }
	   }*/




	/*     for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射情况
		            Z_Target[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);

				  }
		          else
				  {
					Z_Target[i*M+j] = -0.5;                    
				  }
		   }
	   }*/


//////////////////////////////////////////////////////////////////////////
	Vec_Com Z_S_T_3;  Z_S_T_3.assign(N*M);
	Vec_Com Z_S_T_4;  Z_S_T_4.assign(N*M);

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));		   		   
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_3[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);		   		   
		   Z_S_T_1[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);

		}
	}

		
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));		   		   
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_4[i*M+j] = -1.*NR*ds*J*K2/4.*Hankel(1,1,K2*rm);
		   Z_S_T_2[i*M+j] = ds*J/4.*Hankel(0,1,K2*rm);

		}
	}
	


/*	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));		   			        
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_3[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K2*rm);

		}
	}

		
	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));		   			        
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));
		   Z_S_T_4[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K2*rm);
		}
	}*/

	


/*		for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z1[i]-ztt[j],2.0));
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z1[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射时情况
		   Z_S_T_1[i*M+j] = -1.*NR*J*ds*K2/4.*Hankel(1,1,K2*rm);
		}
	}

		for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		   rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(Z2[i]-ztt[j],2.0));
		   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(Z2[i]-ztt[j]))/(rm*(1+pow(df_z[j],2.0)));  //TM波入射时情况
		   Z_S_T_2[i*M+j] = -1.*NR*J*ds*K2/4.*Hankel(1,1,K2*rm);
		}
	}*/



//////////////////////////////////////////////////////////////////////////


	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{		
				rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z1[j],2.0));
				NR = (df_1[j]*(xtt[i]-x[j])-(ztt[i]-Z1[j]))/rm;
				
				
				Z_T_S_A_1[i*N+j] = -J*dx_Suf/4.*Hankel(0,1,K2*rm);
				Z_T_S_B_1[i*N+j] = J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				
		}
	}

	
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
			    rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-Z2[j],2.0));
				NR = (df_2[j]*(xtt[i]-x[j])-(ztt[i]-Z2[j]))/rm;
				
				Z_T_S_A_2[i*N+j] = J*dx_Suf/4.*Hankel(0,1,K2*rm);				
				Z_T_S_B_2[i*N+j] = -J*(K2*dx_Suf/4.)*NR*Hankel(1,1,K2*rm);
				

		}
	}




		

/////////////////////////////The impedance matrix for the MOM/////////////////////////////////////////////
/*   int N_M = 4*N+2*M;
    
   for (i=0;i<N;i++)
    {
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = A_1_1_1[i*N+j];
		}
        for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_1_1_1[i*N+j-N];
		}
	}

    for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_1_1[(i-N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_1[(i-N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_1_2[(i-N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_1_2[(i-N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_3[(i-N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_1[(i-N)*M+j-4*N-M];			
		}
	}

	for (i=2*N;i<3*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*A_2_2_1[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_1[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = A_2_2_2[(i-2*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_2_2_2[(i-2*N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_4[(i-2*N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = Z_S_T_2[(i-2*N)*M+j-4*N-M];
		}

	}

	for (i=3*N;i<4*N;i++)
	{
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_3*A_3_2_2[(i-3*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = B_3_2_2[(i-3*N)*N+j-3*N];
		}
	}

	for (i=4*N;i<4*N+M;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_MOM[i*N_M+j] = rou_2*Z_T_S_A_1[(i-4*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_B_1[(i-4*N)*N+j-N];
		}
		for (j=2*N;j<3*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_A_2[(i-4*N)*N+j-2*N];
		}
		for (j=3*N;j<4*N;j++)
		{
			ZT_MOM[i*N_M+j] = Z_T_S_B_2[(i-4*N)*N+j-3*N];
		}
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_I[(i-4*N)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_J[(i-4*N)*M+j-4*N-M];
		}
	}

	for (i=4*N+M;i<N_M;i++)
	{
		for (j=4*N;j<4*N+M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_K[(i-4*N-M)*M+j-4*N];
		}
		for (j=4*N+M;j<N_M;j++)
		{
			ZT_MOM[i*N_M+j] = ZT_L[(i-4*N-M)*M+j-4*N-M];
		}
	}

	

	for (i=0;i<N;i++)
	{
		EE_MOM[i] = Ei[i];
	}

	cout<<"The matrix has been filled."<<endl;
	
   GAUSS(ZT_MOM,Jf_MOM,EE_MOM,N_M);*/






//	Bi_CG(ZT_MOM,EE_MOM,Jf_MOM);                     //The MOM method


//CLinearEquations linear;


for (m=0;m<15;m++)
{

cout<<"Begin the "<<m<<"th step of the iteration......------>>"<<endl;
 
for (l=0;l<N;l++)
{

	cout<<"Begin the "<<l+1<<"th step of the FBM......------>>"<<endl;

/////////////////////////The forward sweep/////////////////////////////////////////////////


    for (i=0;i<N;i++)
    {
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = Ei[i] - Sum0;
			V[1] = -E_Surf1_Target[i]-Sum1;
			V[2] = -E_Surf2_Target[i]-Sum2;
			V[3] = -Sum3;



			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

		


			GAUSS(A,X,V,4);

		

			U1_f[i] = X[0];
			W1_f[i] = X[1];
			U2_f[i] = X[2];
			W2_f[i] = X[3];

    }


/////////////////////////The backward sweep/////////////////////////////////////////////////
	
	
   for (i=N-1;i>=0;i--)
    {
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + A_1_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_1_1_1[i*N+j]*(W1_f[j]+W1_b[j]);
			Sum1 = Sum1 + rou_2*A_2_1_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_1_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_1_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_1_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum2 = Sum2 + rou_2*A_2_2_1[i*N+j]*(U1_f[j]+U1_b[j]) + B_2_2_1[i*N+j]*(W1_f[j]+W1_b[j]) + A_2_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_2_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
			Sum3 = Sum3 + rou_3*A_3_2_2[i*N+j]*(U2_f[j]+U2_b[j]) + B_3_2_2[i*N+j]*(W2_f[j]+W2_b[j]);
		}
            V[0] = -Sum0;
			V[1] = -Sum1;
			V[2] = -Sum2;
			V[3] = -Sum3;

			Sum0 = zero;
			Sum1 = zero;
			Sum2 = zero;
			Sum3 = zero;

			A[0] = A_1_1_1[i*N+i];
			A[1] = B_1_1_1[i*N+i];
			A[4] = rou_2*A_2_1_1[i*N+i];
			A[5] = B_2_1_1[i*N+i];
			A[6] = A_2_1_2[i*N+i];
			A[7] = B_2_1_2[i*N+i];
			A[8] = rou_2*A_2_2_1[i*N+i];
			A[9] = B_2_2_1[i*N+i];
			A[10] = A_2_2_2[i*N+i];
			A[11] = B_2_2_2[i*N+i];
			A[14] = rou_3*A_3_2_2[i*N+i];
			A[15] = B_3_2_2[i*N+i];

			
				
			GAUSS(A,X,V,4);


			U1_b[i] = X[0];
			W1_b[i] = X[1];
			U2_b[i] = X[2];
			W2_b[i] = X[3];

    }

//////////////////////////The convergence condition////////////////////////////////////////////////

	
	for (i=0;i<N;i++)
	{
		Jf1[i]   = U1_f[i];
		Jf1[i+N] = W1_f[i];
		Jf1[i+2*N] = U2_f[i];
		Jf1[i+3*N] = W2_f[i];
	}

	for (i=0;i<4*N;i++)
	{
		Jf3[i] = Jf1[i] - Jf2[i];
	}

	  Assign_Vecor(Jf1,Jf2);  
	  Delta = sqrt(L2_Norm(Jf3))/sqrt(L2_Norm(Jf1));      
	  
	  cout<<"The error of the FBM is: Delta ="<<Delta<<endl;
			  
	  if(Delta<=0.01)
	  {
	    cout<<"The total number of the FBM step is :"<<l+1<<"."<<"and the final error is :"<<Delta<<endl;
		  break;
	  }
				
	     else  
		     continue;

}



	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
            Sum_E0  = Sum_E0 + rou_2*Z_T_S_A_1[i*N+j]*(U1_f[j] + U1_b[j]) + Z_T_S_B_1[i*N+j]*(W1_f[j]+W1_b[j]) + Z_T_S_A_2[i*N+j]*(U2_f[j] + U2_b[j]) + Z_T_S_B_2[i*N+j]*(W2_f[j] + W2_b[j]);
		}
		
           E_Target_Surf[i] = -Sum_E0;
		   Sum_E0 = zero;
	}


  GAUSS(Z_Target,J_Target,E_Target_Surf,M);
	
//   Bi_CG(Z_Target,E_Target_Surf,J_Target); 


   for (i=0;i<N;i++)
   {

	   for (j=0;j<M;j++)
	   {
		   Sum_E0 = Sum_E0 + Z_S_T_1[i*M+j]*J_Target[j];
		   Sum_E1 = Sum_E1 + Z_S_T_2[i*M+j]*J_Target[j];
	   }
	       E_Surf1_Target[i] = Sum_E0;
		   E_Surf2_Target[i] = Sum_E1;
		   Sum_E0 = zero;
		   Sum_E1 = zero;

   }



		for (i=0;i<M;i++)
		{
			J_Target_1[i] = J_Target[i] - J_Target_2[i]; 
		}

		 Assign_Vecor(J_Target,J_Target_2);

		 for (i=0;i<M;i++)
		 {
			 for (j=0;j<M;j++)
			 {
				 Sum_E0 = Sum_E0 + Z_Target[i*M+j]*J_Target_1[j];
			 }
			 Delta__J_Target[i] = Sum_E0;
			 Sum_E0 = zero;
             
		 }

		 Delta_Target_Surface = sqrt(L2_Norm(Delta__J_Target))/sqrt(L2_Norm(E_Target_Surf));  
		 
		 Delta_Matrix[m] = Delta_Target_Surface;
	  
	     cout<<"The error of the alternate is: Delta_Target_Surface ="<<Delta_Target_Surface<<endl;
			  
/**/	    if(Delta_Target_Surface<=0.01)
		{
	     cout<<"The total number of the alternate step is :"<<m+1<<"."<<"and the final error is :"<<Delta_Target_Surface<<endl;
		  break;
		}
				
	     else  
		     continue;


}
			

		double abx   =  8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an)*(1.0-((1.0+2.0*pow(tan(i_an),2))/(2.0*pow(k_w*g_p*cos(i_an),2))));
		   
		for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-Z1[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
//					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					NR = k_w*(df_1[j]*sin(Sita_D[i])-cos(Sita_D[i]))/sqrt(1+pow(df_1[j],2.0));


				    double ksc=-k_w*(x[j]*sin(Sita_D[i])+Z1[j]*cos(Sita_D[i]));
				    Sum_E0=Sum_E0+dx_Suf*(NR*(W1_f[j]+W1_b[j])-(U1_f[j]+U1_b[j]))*exp(i_*ksc);
//				    Sum_E0=Sum_E0+dx_Suf*(NR*(Jf_MOM[j+N])-(Jf_MOM[j]))*exp(i_*ksc);

//					Sum_E0 = Sum_E0+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}
                 Es[i] = Sum_E0;				 
				 D_Es[i]  = (mm*D_Es[i]+norm(Sum_E0)/abx)/(mm+1);
                 DB_Es[i] = (mm*DB_Es[i]+10.*log10(D_Es[i]))/(mm+1);
		         Sum_E0   = zero;
			}
   	}      

/*		    for (i=0;i<N0;i++)
			{
				Es[i] = Es[i]/double(NUMBER_REALIZATION);
				D_Es[i]=R_D*pow(abs(Es[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
				DB_Es[i]=10.*log10(D_Es[i]);
			}

	        
			
			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es1[i] = Es1[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(Jf_MOM[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(Jf_MOM[j]);
				}

				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es1[i]=10.*log10(D_Es1[i]);
		
			}*/
		 

	ofstream out_U1;
	         out_U1.open("Layered_Bried_Target_U1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<abs(Jf_MOM[i])<<endl;
			 }
   
    ofstream out_W1;
	         out_W1.open("Layered_Bried_Target_W1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<abs(Jf_MOM[i+N])<<endl;
			 }

	ofstream out_U2;
	         out_U2.open("Layered_Bried_Target_U2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<abs(Jf_MOM[i+2*N])<<endl;
			 }
	ofstream out_W2;
	         out_W2.open("Layered_Bried_Target_W2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<abs(Jf_MOM[i+3*N])<<endl;
			 }
	ofstream out_Es;
	         out_Es.open("Layered_Bried_Target_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	ofstream out_Es1;
	         out_Es1.open("Layered_Bried_Target_Es1.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<D_Es1[i]<<"   "<<DB_Es1[i]<<endl;
			 }

     ofstream out_J_Target;
	         out_J_Target.open("Layered_Bried_Target_J.txt");
		     for (i=0;i<M;i++)
		     {
				 out_J_Target<<i<<"   "<<real(J_Target[i])<<"   "<<imag(J_Target[i])<<"   "<<abs(J_Target[i])<<"  "<<abs(Jf_MOM[i+4*N])<<endl;
		     }
     ofstream out_J_Target_Z;
	         out_J_Target_Z.open("Layered_Bried_Target_Z.txt");
		     for (i=0;i<N;i++)
		     {
				 out_J_Target_Z<<i<<"   "<<Z1[i]<<"   "<<Z2[i]<<endl;
		     }

    ofstream out_Delta;
	         out_Delta.open("Layered_Bried_Target_Delta.txt");
		     for (i=0;i<15;i++)
		     {
				 out_Delta<<i<<"   "<<Delta_Matrix[i]<<endl;
		     }


/*	ofstream out_EE_Total;
	         out_EE_Total.open("EE_ToTal.txt");
			 for (i=0;i<N1;i++)
			 {
				 for (j=0;j<N1;j++)
				 {
					 out_EE_Total<<i<<"   "<<j<<"    "<<EE_Total[i*N1+j]<<endl;
				 }
			 }*/
}
