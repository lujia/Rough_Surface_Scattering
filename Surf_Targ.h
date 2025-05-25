#include"surface.h"
#include "Target.h"

class Surf_Targ_Question: public post_cp, public CylinderQuestion
{

public:
	Surf_Targ_Question();
	~Surf_Targ_Question();
	void Gauss_Surface_Cylinder_Target_J();
	void Gauss_Surface_Cylinder_Target_DRCS(Vec_Com &Jsd_in,Vec_Com &Jto_in);
	void Gauss_Surface_Cylinder_Target_Jsd1_Jsd2();
	void Gauss_Surface_Cylinder_Target_Muti_Target();
    void Gauss_Surface_medium_Cylinder_Target_J();
	void Gauss_Del_Surface_Del_Cylinder__Above_Target();
	void Gauss_Del_Surface_Con_Cylinder__Above_Target();
    void Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Es1();
    void Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Esd();
	void Gauss_Surface_Below_Muti_Target();
    

	
};

Surf_Targ_Question::Surf_Targ_Question()
{

}

Surf_Targ_Question::~Surf_Targ_Question()
{

}

//////////////////////////////////////////////////////////////////////////


void Surf_Targ_Question::Gauss_Surface_Cylinder_Target_J()
{
	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=102.4*wave_len;
	double R=5.*wave_len;
	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double dl=0.;
	double abx1=0.;
	double abx2=0.;
	double abx3=0.;
	double ET01_M=0.;
	double ET0_M=0.;
	double Delta=0.;
	double xpp=0.;

	int i,j;

	Vec_Com EST;         // 粗糙面在目标表面的散射场；
	   EST.assign(M);
	Vec_Com EFT;         //差值电流在目标表面的散射场；
	   EFT.assign(M);
	Vec_Com ET0;         // 目标总的入射场；
	   ET0.assign(M); 
	Vec_Com ETF;
	   ETF.assign(N);
	Vec_Com Jsd1;
	   Jsd1.assign(N);
    Vec_Com Jsd2;
	   Jsd2.assign(N);
    Vec_Com J01;
       J01.assign(M);
	Vec_Com ET01;
	   ET01.assign(M);
	Vec_Com JT0;
	   JT0.assign(N);
	Vec_Com JT01;
	   JT01.assign(N);
   Vec_Com Jsd;
	   Jsd.assign(N);

	Complex zero(0.,0.);
	Complex max =zero;
    Complex max1 =zero;
	Complex Sum_EST=zero;
	Complex Sum_ETF=zero;
	Complex Sum_EFT=zero;

    cout<<"进行目标的剖分......------>>"<<endl;
    Cy_PoufenMatix();
   
    cout<<"进行目标阻抗矩阵的填充......------>>"<<endl;
	Cy_ImplenceMatix();
	 
	cout<<"计算Js0......------>>"<<endl;

	One_D_Conductor_Gauss_Js0(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);	

    double dx_Suf=RL/N;

	for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
		   dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
		   abx1=k_w*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
		   Sum_EST=Sum_EST-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Js0[j];
		}
			EST[i]=Sum_EST;
            Sum_EST=zero;
	}


	cout<<"计算粗糙面差值感应电流Jsd......------>>"<<endl;
    for (int k=0;k<20;k++)
	{
		cout<<"进行第"<<k+1<<"次目标与粗糙面的耦合迭代计算......------>>"<<endl;	

//////////////////////////////计算ES0/////////////////////////////////////

		  

//////////////////////////////计算差值电流在目标表面产生的电场/////////////////////////////////////


		  for(i=0;i<M;i++)
		  {
		       for (j=0;j<N;j++)
			   {
			     dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			     abx3=k_w*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
			     Sum_EFT=Sum_EFT-(k_w*Yita/4.)*Hankel(0,1,abx3)*dl*Jsd[j];
			   }
			   EFT[i]=Sum_EFT;
               Sum_EFT=zero;
		  }




//////////////////////////////计算目标表面总激励场/////////////////////////////////////
		 

	     for(i=0;i<M;i++)
		 {
           ET0[i]= exp(J*K*(xtt[i]*sin(Sita0)-ztt[i]*cos(Sita0)))+EST[i]+EFT[i];
		 }

//////////////////////////////计算有粗糙面的情况下目标表面的总感应电流/////////////////////////////////////

            cout<<"用CG法求目标表面电流分布......------>>"<<endl;;
	        CLinearEquations linear;
	           linear.CG(ZT,JT0,ET0,M,M,0.01);


//////////////////////////////计算目标表面的总感应电流在粗糙表面产生的电场/////////////////////////////////////
		    
	       for(i=0;i<N;i++)
		   {
		       for (j=0;j<M;j++)
			   {
				   abx2=k_w*sqrt(pow(xtt[j]-x[i],2)+pow(ztt[j]-z[i],2));
			       Sum_ETF = Sum_ETF-(k_w*Yita/4.)*Hankel(0,1,abx2)*ds*JT0[j];						
			   }
			   ETF[i] = Sum_ETF;
			   Sum_ETF = zero;

		   }

//////////////////////////////计算差值电流/////////////////////////////////////

		   		             	       
		   One_D_Conductor_BMIA_CAG_Matrix_Jsd(ETF,Jsd);


//////////////////////////////收敛判断条件/////////////////////////////////////	         
	       
		   max =zero;
           max1 =zero;

	       for (i=0;i<N;i++)
		   {
		        if (abs(max)<abs(JT0[i]-JT01[i]))
			               max=JT0[i]-JT01[i];
				else 
					continue;
			
		   }

		    for (i=0;i<N;i++)
			{
		       JT01[i] = JT0[i];
			}

		   
		   for (i=0;i<N;i++)
		   {
		        if (abs(max1)<abs(JT0[i]))
			               max1=JT0[i];
				else 
					continue;
			
		   }



		      Delta = abs(max)/abs(max1);
              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.001)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;

	  
	 
	}


//////////////////////////////输出计算结果/////////////////////////////////////
	        
             ofstream outEST;
		         outEST.open("EST.txt");
			 for (i=0;i<M;i++)
			  {
				  outEST<<i<<" "<<real(EST[i])<<"  "<<imag(EST[i])<<"  "<<abs(EST[i])<<endl;

			  }	   



	          ofstream outJT0;
		         outJT0.open("JT0.txt");
			  for (i=0;i<M;i++)
			  {
				  outJT0<<i<<" "<<real(JT0[i])<<"  "<<imag(JT0[i])<<"  "<<abs(JT0[i])<<endl;

			  }

	


			  ofstream outETF;
		         outETF.open("ETF.txt");
			  for (i=0;i<N;i++)
			  {
				outETF<<x[i]/Lamda<<"   "<<real(ETF[i])<<"   "<<imag(ETF[i])<<"   "<<ETF[i]<<endl;

			  }

              ofstream outJsd_1;
		         outJsd_1.open("Jsd_1.txt");
		      
			  for(i=0;i<N;i++)
			  {
				 outJsd_1<<x[i]/Lamda<<"     "<<real(Jsd[i])<<"     "<<imag(Jsd[i])<<"     "<<abs(Jsd[i])<<endl;
			  }

			  ofstream outJsd_2;
			     outJsd_2.open("Jsd_2.txt");
			  
			  for (i=0;i<N;i++)
			  {
			    outJsd_2<<Jsd[i]<<endl;
			  }


//	           Gauss_Surface_Cylinder_Target_DRCS(Jsd,JT0);


////////////////////////////////计算差场散射截面/////////////////////////////

   Vec_dou Jsd_IN_1;
	  Jsd_IN_1.assign(4);
   Vec_dou JT0_IN_1;
	  JT0_IN_1.assign(4);
   Vec_Com Jsd_IN;
	  Jsd_IN.assign(N);
   Vec_Com JT0_IN;
	  JT0_IN.assign(N);  
	ifstream in_Jsd;
	      in_Jsd.open("Jsd1.txt");
	    for(i=0;i<N;i++)
		{
			for (j=0;j<4;j++)
			{
				in_Jsd>>Jsd_IN_1[j];

			}
		Jsd_IN[i] = Jsd_IN_1[1]+J*Jsd_IN_1[2];

		}

	ifstream in_JT0;
	      in_JT0.open("JT0.txt");

       for (i=0;i<M;i++)
       {
		   for (j=0;j<4;j++)
			{
				in_JT0>>JT0_IN_1[j];

			}
		JT0_IN[i] = JT0_IN_1[1]+J*JT0_IN_1[2];
		   
       } 

/*
	   for (i=0;i<N;i++)
	   {
		   cout<<x[i]<<endl;
	   }
	   
*/
	   


	int N0=180;
	double R_Dsc = 100000.*Lamda;
	double Sita_Dsc = 0.;
	Complex Sum_Esd;
	double dl_Esd=0.;
	double abx_Esd_S=0.;
    double abx_Esd_T=0.;
	
	Vec_Com Esd;
	   Esd.assign(N0);
	Vec_dou x_Dsc;
	  x_Dsc.assign(N0);
	Vec_dou z_Dsc;
	  z_Dsc.assign(N0);
	Vec_Com Ei_Dsc;
	  Ei_Dsc.assign(N0);
	Vec_dou Dsc;
	  Dsc.assign(N0);


	for (i=N0-1;i>=0;i--)
	{  
		Sita_Dsc = i*PI/N0;
		
		x_Dsc[N0-i-1] = R_Dsc*cos(Sita_Dsc);
		z_Dsc[N0-i-1] = R_Dsc*sin(Sita_Dsc);

	}
/*    
	i_an = i_an*PI/180;
//	cout<<"i_an="<<i_an<<endl;
 
	double xsz,ww,kg,fac;
    for(i=0;i<N0;i++)
	{
	//	xsz=x_Dsc[i]*sin(i_an)-z_Dsc[i]*cos(i_an);
	//	kg=pow(k_w*g_p*cos(i_an),2.0);
	//	fac=pow((x_Dsc[i]+H*tan(i_an)+z_Dsc[i]*tan(i_an))/g_p,2.0);
	//	ww=(2.0*fac-1.0)/kg;
	//	Ei_Dsc[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		Ei_Dsc[i]=exp(i_*k_w*xsz);

		
	}

	ofstream outEi;
	  outEi.open("Ei.txt");
	for (i=0;i<N0;i++)
	{
		outEi<<i-90<<"   "<<real(Ei_Dsc[i])<<"   "<<imag(Ei_Dsc[i])<<"   "<<abs(Ei_Dsc[i])<<endl;
	}
*/
	dx_Suf=0.01;


	for (k=0;k<N0;k++)
	{

        for (i=0;i<M;i++)
        {
			abx_Esd_T = k_w*sqrt(pow(x_Dsc[k]-xtt[i],2)+pow(z_Dsc[k]-ztt[i],2));
			Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_T)*ds*JT0_IN[i];						
		}

		for (j=0;j<N;j++)
		{
			 dl_Esd=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			 abx_Esd_S=k_w*sqrt(pow(x_Dsc[k]-x[j],2)+pow(z_Dsc[k]-z[j],2));
			 Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_S)*dl_Esd*Jsd_IN[j];
		}
		
	}

	for (i=0;i<N0;i++)
	{  
//		Dsc[i] = PI*R_Dsc*pow(abs(Esd[i]),2.0)/pow(abs(Ei_Dsc[i]),2.0);
      Dsc[i] = 2.*PI*R_Dsc*pow(abs(Esd[i]),2.0);
	}

    cout<<"输出差场散射计算结果......------>>"<<endl;
	ofstream outEsd;
	 outEsd.open("Esd.txt");
	
	ofstream outDsc;
	 outDsc.open("Dsc.txt");
	 
	ofstream outDscDB;
	 outDscDB.open("DscDB.txt");

		 for (i=0;i<N0;i++)
		 {
			 outDsc<<i-90<<" "<<Dsc[i]<<endl;
			 outDscDB<<i-90<<" "<<10.*log10(Dsc[i])<<endl;
			 outEsd<<i-90<<" "<<abs(Esd[i])<<endl;
		 }	
		 


/////////////////////////////////////////////////////////////////////////////////


}


void Surf_Targ_Question::Gauss_Surface_Cylinder_Target_DRCS(Vec_Com &Jsd_in,Vec_Com &Jt0_in)
{
	int i,j,k;
	int N0=180;
	double R_Dsc = 10000.*Lamda;
	double Sita_Dsc = 0.;
	Complex Sum_Esd;
	double dl_Esd=0.;
	double abx_Esd_S=0.;
    double abx_Esd_T=0.;
	
	Vec_Com Esd;
	   Esd.assign(N0);
	Vec_dou x_Dsc;
	  x_Dsc.assign(N0);
	Vec_dou z_Dsc;
	  z_Dsc.assign(N0);
	Vec_Com Ei_Dsc;
	  Ei_Dsc.assign(N0);
	Vec_dou Dsc;
	  Dsc.assign(N0);

	for (i=N0-1;i>=0;i--)
	{  
		Sita_Dsc = i*PI/N0;
		
		x_Dsc[N0-i-1] = R_Dsc*cos(Sita_Dsc);
		z_Dsc[N0-i-1] = R_Dsc*sin(Sita_Dsc);

	}

	double xsz,ww,kg,fac;
    for(i=0;i<N0;i++)
	{
		xsz=x_Dsc[i]*sin(i_an)-z_Dsc[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x_Dsc[i]+H*tan(i_an)+z_Dsc[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei_Dsc[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		
	}

	dx_Suf=RL/N;


	for (k=0;k<N0;k++)
	{

//////////////////////////////////目标表面电流的散射场////////////////////////

        for (i=0;i<M;i++)
        {
			abx_Esd_T = k_w*sqrt(pow(x_Dsc[k]-xtt[i],2)+pow(z_Dsc[k]-ztt[i],2));
			Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_T)*ds*Jt0_in[i];						
		}

/////////////////////////////////粗糙面差场电流的散射场////////////////////////

		for (j=0;j<N;j++)
		{
			 dl_Esd=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			 abx_Esd_S=k_w*sqrt(pow(x_Dsc[k]-x[j],2)+pow(z_Dsc[k]-z[j],2));
			 Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_S)*dl_Esd*Jsd_in[j];
		}
		
	}

////////////////////////////////差场散射截面////////////////////////////////////
	
	for (i=0;i<N0;i++)
	{  
//		Dsc[i] = PI*R_Dsc*pow(abs(Esd[i]),2.0)/pow(abs(Ei_Dsc[i]),2.0);
        Dsc[i] = 2.*PI*R_Dsc*pow(abs(Esd[i]),2.0);
	}

/////////////////////////////////////////////////////////////////////////////////////////////////

    cout<<"输出差场散射......------>>"<<endl;
	ofstream outEsd;
	 outEsd.open("Esd.txt");
	
	ofstream outDsc;
	 outDsc.open("Dsc.txt");
	 
	ofstream outDscDB;
	 outDscDB.open("DscDB.txt");

		 for (i=0;i<N0;i++)
		 {
			 outDsc<<i-90<<" "<<Dsc[i]<<endl;
			 outDscDB<<i-90<<" "<<10.*log10(Dsc[i])<<endl;
			 outEsd<<i-90<<" "<<abs(Esd[i])<<endl;
		 }
}

void Surf_Targ_Question::Gauss_Surface_Cylinder_Target_Jsd1_Jsd2()
{
	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=102.4*wave_len;
	double R=5.*wave_len;
	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double dl=0.;
	double abx1=0.;
	double abx2=0.;
	double abx3=0.;
	double ET01_M=0.;
	double ET0_M=0.;
	double Delta=0.;

	int i,j,k;

	Vec_Com EFT_Js0;         // 粗糙面在目标表面的散射场；
	   EFT_Js0.assign(M);
	Vec_Com EFT_Jsd2;         //差值电流在目标表面的散射场；
	   EFT_Jsd2.assign(M);
	Vec_Com ET0;         // 目标总的入射场；
	   ET0.assign(M);
	Vec_Com Jsd1;
	   Jsd1.assign(M);
	Vec_Com EF0;
	   EF0.assign(N);

	Vec_Com ETF_J0;
	   ETF_J0.assign(N);
    Vec_Com ETF_Jsd1;
	   ETF_Jsd1.assign(N);

    Vec_Com Jsd2;
	   Jsd2.assign(N);
    Vec_Com J01;
       J01.assign(M);
	Vec_Com ET01;
	   ET01.assign(M);
	Vec_Com JT0;
	   JT0.assign(N);
	Vec_Com JT01;
	   JT01.assign(N);
	Vec_Com Jsd;   
	   Jsd.assign(M);

	Complex zero(0.,0.);
	Complex max =zero;
    Complex max1 =zero;
	Complex Sum_EST=zero;
	Complex Sum_ETF=zero;
	Complex Sum_EFT=zero;


    cout<<"计算J0......------>>"<<endl;
	
	Cy_ResultofJ();
	 
	cout<<"计算Js0.....------>>"<<endl;

	One_D_Conductor_Gauss_Js0(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);	
	
	double dx_Suf=RL/N;

	cout<<"计算目标表面差值感应电流Jsd1与粗糙面差值感应电流Jsd2......------>>"<<endl;

	 for(i=0;i<M;i++)
	{ 
	    for (j=0;j<N;j++)
		{
		   dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
		   abx1=k_w*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
		   Sum_EST=Sum_EST-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Js0[j];
		}
			 EFT_Js0[i]=Sum_EST;
             Sum_EST=zero;
	}



	for(i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			abx2=k_w*sqrt(pow(xtt[j]-x[i],2)+pow(ztt[j]-z[i],2));
			Sum_ETF = Sum_ETF-(k_w*Yita/4.)*Hankel(0,1,abx2)*ds*J0[j];						
		}
		    ETF_J0[i] = Sum_ETF;
			Sum_ETF = zero;

	}


	for (k=0;k<20;k++)
	{

	    cout<<"进行第"<<k+1<<"次目标与粗糙面的耦合迭代计算......------>>"<<endl;	


	   

	  for(i=0;i<M;i++)
	  {
		  for (j=0;j<N;j++)
		  {
			  dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			  abx1=k_w*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
			  Sum_EFT=Sum_EFT-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Jsd2[j];
		  }
			   EFT_Jsd2[i]=Sum_EFT;
               Sum_EFT=zero;
	  }

	  for(i=0;i<M;i++)
	  {
           ET0[i]= EFT_Js0[i]+EFT_Jsd2[i];
	  }


	   cout<<"用CG法求目标表面电流分布......------>>"<<endl;;
	        CLinearEquations linear;
	           linear.CG(ZT,Jsd1,ET0,M,M,0.01);


	  

	    
	   for(i=0;i<N;i++)
		   {
		       for (j=0;j<M;j++)
			   {
				   abx2=k_w*sqrt(pow(xtt[j]-x[i],2)+pow(ztt[j]-z[i],2));
			       Sum_ETF = Sum_ETF-(k_w*Yita/4.)*Hankel(0,1,abx2)*ds*Jsd1[j];						
			   }
			   ETF_Jsd1[i] = Sum_ETF;
			   Sum_ETF = zero;

		   }

	   for (i=0;i<N;i++)
	   {
		   EF0[i] = ETF_J0[i] + ETF_Jsd1[i] ;
	   }

	   One_D_Conductor_BMIA_CAG_Matrix_Jsd(EF0,Jsd2);


//////////////////////////////收敛判断条件/////////////////////////////////////	         

		   max =zero;
           max1 =zero;

	       for (i=0;i<M;i++)
		   {
		        if (abs(max)<abs(Jsd1[i]-Jsd[i]))
			               max=Jsd1[i]-Jsd[i];
				else 
					continue;
			
		   }

		    for (i=0;i<M;i++)
			{
		       Jsd[i] = Jsd1[i];
			}

		   
		   for (i=0;i<M;i++)
		   {
		        if (abs(max1)<abs(Jsd1[i]))
			               max1=Jsd1[i];
				else 
					continue;
			
		   }



		      Delta = abs(max)/abs(max1);
              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.001)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;
	 
	}

	          ofstream outJsd1;
		         outJsd1.open("Jsd1.txt");
		      
			  for(i=0;i<M;i++)
			  {
				 outJsd1<<i<<"     "<<real(Jsd1[i])<<"     "<<imag(Jsd1[i])<<"     "<<abs(Jsd1[i])<<endl;
			  }

			  ofstream outJsd2;
			     outJsd2.open("Jsd2.txt");
			  
			  for (i=0;i<N;i++)
			  {
			   	 outJsd2<<x[i]/Lamda<<"     "<<real(Jsd2[i])<<"     "<<imag(Jsd2[i])<<"     "<<abs(Jsd2[i])<<endl;

			  }

			  ofstream outJsd1_J0;
			     outJsd1_J0.open("Jsd1_J0.txt");
			  
			  for (i=0;i<M;i++)
			  {
			   	 outJsd1_J0<<i<<"     "<<real(Jsd1[i]+J0[i])<<"     "<<imag(Jsd1[i]+J0[i])<<"     "<<abs(Jsd1[i]+J0[i])<<endl;

			  }

}

void Surf_Targ_Question::Gauss_Surface_Cylinder_Target_Muti_Target()
{
	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=102.4*wave_len;
	double R=5.*wave_len;
	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double dl=0.;
	double abx1=0.;
	double abx2=0.;
	double abx3=0.;
	double ET01_M=0.;
	double ET0_M=0.;
    double Delta=0.;
	double Delta1=0.;
	double Delta2=0.;

	Lamda = 3e8/3e9;
	i_an = i_an*PI/180;

	double H = 10.*Lamda;

	double a = 5.*Lamda;
	double b =3.*Lamda/2.;

	int i,j,k;
	
	double xp1 = -25.*Lamda;
	double xp2 = 25*Lamda;

	double xpp1=-25.*Lamda-H*tan(i_an);
	double xpp2=25*Lamda-H*tan(i_an);

	Vector_double xtt_1;
	   xtt_1.assign(M);
	Vector_double ztt_1;
	   ztt_1.assign(M);
	Vector_double xtt_2;
	   xtt_2.assign(M);
	Vector_double ztt_2;
	   ztt_2.assign(M);
	Vector_double ds_1;
	   ds_1.assign(M);
	Vector_double ds_2;
	   ds_2.assign(M);


	Vector_complex ZT_1;
	   ZT_1.assign(M*M);
	Vector_complex ZT_2;
	   ZT_2.assign(M*M);

	Vector_complex Jsd_T1;
	   Jsd_T1.assign(N);
	Vector_complex Jsd_T2;
	   Jsd_T2.assign(N);

	Vector_complex J0_T1;
	   J0_T1.assign(M);
	Vector_complex J0_T1_1;
	   J0_T1_1.assign(M);
	Vector_complex J0_T2;
	   J0_T2.assign(M);
	Vector_complex J0_T2_1;
	   J0_T2_1.assign(M);

	Vector_complex Js0_1;
	   Js0_1.assign(N);
	Vector_complex Js0_2;
	   Js0_2.assign(N);

	Vector_complex E_Js0_T1;
	    E_Js0_T1.assign(M);
	Vector_complex E_Js0_T2;
	    E_Js0_T2.assign(M);
	Vector_complex E_J0_T1;
	    E_J0_T1.assign(M);
	Vector_complex E_J0_T2;
	    E_J0_T2.assign(M);
   	
	Vector_complex E_J0_T1_S;
	    E_J0_T1_S.assign(N);
	Vector_complex E_J0_T2_S;
	    E_J0_T2_S.assign(N);


	Vector_complex E_Jsd_T1_T1;
	    E_Jsd_T1_T1.assign(M);
    Vector_complex E_Jsd_T2_T1;
	    E_Jsd_T2_T1.assign(M);

	Vector_complex E_Jsd_T1_T2;
	    E_Jsd_T1_T2.assign(M);
    Vector_complex E_Jsd_T2_T2;
	    E_Jsd_T2_T2.assign(M);


	Vector_complex Ei_T1;
		   Ei_T1.assign(M);
    Vector_complex Ei_T2;
		   Ei_T2.assign(M);
	
    Vector_complex Ei_s1;
		   Ei_s1.assign(N);
    Vector_complex Ei_s2;
		   Ei_s2.assign(N);

	Vector_complex E0_T1;
		   E0_T1.assign(M);
    Vector_complex E0_T2;
		   E0_T2.assign(M);
	Vector_complex ET0_T1;
		   ET0_T1.assign(M);
    Vector_complex ET0_T2;
		   ET0_T2.assign(M);


	
		   
	for (i=0;i<M;i++)
	{
		Ei_T1[i] = exp(J*k_w*(xtt_1[i]*sin(Sita0)-ztt_1[i]*cos(Sita0)));
		Ei_T1[i] = exp(J*k_w*(xtt_2[i]*sin(Sita0)-ztt_2[i]*cos(Sita0)));	    
	}



	Complex zero(0.,0.);

	Complex Sum_Js0 = zero;
	Complex Sum_J0 = zero;
	Complex Sum_Jsd = zero;
	Complex max_J0_T1=zero;
	Complex max1_J0_T1=zero;

	Complex max_J0_T2=zero;
	Complex max1_J0_T2=zero;
		

  cout<<"对目标1进行剖分......------>>"<<endl;	   

    Muti_target_poufen(xp1,xtt_1,ztt_1,ds_1, a, b,M);

//	for (i=0;i<N;i++)
//	{
//		cout<<xtt_1[i]<<endl;
			
//	}
	
  cout<<"对目标2进行剖分......------>>"<<endl;	   
	
	Muti_target_poufen(xp2,xtt_2,ztt_2,ds_2,b,a,M);

  cout<<"填充目标1的阻抗矩阵......------>>"<<endl;	   

    Muti_target_ImplenceMatix(xtt_1,ztt_1,ZT_1,ds_1,M);

  cout<<"填充目标2的阻抗矩阵......------>>"<<endl;	

    Muti_target_ImplenceMatix(xtt_2,ztt_2,ZT_2,ds_2,M);
/*
     cout<<"Sita0 = "<<Sita0<<endl;

	 cout<<"K="<<K<<endl;

    for (i=0;i<M;i++)
    {
		ET0_T1[i] = exp(J*K*(xtt_1[i]*sin(Sita0)-ztt_1[i]*cos(Sita0)));
    }
	



	CLinearEquations linear;
	           linear.CG(ZT_1,J0_T1,ET0_T1,M,M,0.01);
//    resolve_FBM(ZT_1,ET0_T1,J0_T1,M);
	for (i=0;i<M;i++)
	{
		cout<<J0_T1[i]<<endl;
			
	}


	ofstream outztt_1;
	   outztt_1.open("ztt_1.txt" );
	ofstream outztt_2;
	   outztt_2.open("ztt_2.txt" );   

	for (i=0;i<M;i++)
	{
		outztt_1<<xtt_1[i]/Lamda<<"   "<<ztt_1[i]/Lamda<<"   "<<ds_1[i]/Lamda<<endl;
	}

		for (i=0;i<M;i++)
	{
		outztt_2<<xtt_2[i]/Lamda<<"   "<<ztt_2[i]/Lamda<<"   "<<ds_2[i]/Lamda<<endl;
	}


  
    ofstream outJ0_T1;
			  outJ0_T1.open("J0_T1.txt");
			  
			  for (i=0;i<M;i++)
			  {
			   	 outJ0_T1<<i<<"     "<<real(J0_T1[i])<<"     "<<imag(J0_T1[i])<<"     "<<abs(J0_T1[i])<<endl;

			  }




 */ 

//	One_D_Conductor_Gauss_Js0(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);
 	
   cout<<"生成粗糙面......------>>"<<endl;	
  	
   One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

	

   double dx_Suf=RL/N;

////////////////////设置两个粗糙面的锥形入射波////////////////////////////////

	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]-xpp1*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei_s1[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		
	}



   for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]-xpp2*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei_s2[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		
	}

///////////////////////////////////////////////////////////////////////////////

  cout<<"计算Js0_1......------>>"<<endl;	

   One_D_Conductor_BMIA_CAG_Matrix_Jsd(Ei_s1,Js0_1);

  cout<<"计算Js0_2......------>>"<<endl;
  
   One_D_Conductor_BMIA_CAG_Matrix_Jsd(Ei_s2,Js0_2);



  cout<<"计算 E_Js0_T1 和 E_Js0_T2......------>>"<<endl;
  

/////////////////////////目标1的总入射波/////////////////////////////////

   for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_1[i]-x[j],2)+pow(ztt_1[i]-z[j],2));
			Sum_Js0=Sum_Js0-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Js0_1[j];
		}
			E_Js0_T1[i]=Sum_Js0;
            Sum_Js0=zero;
	}


/////////////////////////目标2的总入射波/////////////////////////////////


    for (i=0;i<M;i++)
	{
	   E0_T1[i] = exp(J*k_w*(xtt_1[i]*sin(Sita0)-ztt_1[i]*cos(Sita0)))+E_Js0_T1[i];
	}

	
   
  
    for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_2[i]-x[j],2)+pow(ztt_2[i]-z[j],2));
			Sum_Js0=Sum_Js0-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Js0_2[j];
		}
			E_Js0_T2[i]=Sum_Js0;
            Sum_Js0=zero;
	}
	


	for (i=0;i<M;i++)
	{
	   E0_T2[i] = exp(J*k_w*(xtt_2[i]*sin(Sita0)-ztt_2[i]*cos(Sita0)))+E_Js0_T2[i];
	}


////////////////////////////////////////////////////////////////////////////////////////

cout<<"开始耦合迭代运算......------>>"<<endl;

for(k=0;k<20;k++)
{

    cout<<"进行第"<<k+1<<"次迭代"<<"......------>>"<<endl;

/////////////////////////目标2上的感应电流在目标1表面的入射场/////////////////////////////////


	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
			abx1=k_w*sqrt(pow(xtt_1[i]-xtt_2[j],2.0)+pow(ztt_1[i]-ztt_2[j],2.0));
			Sum_J0 = Sum_J0-(k_w*Yita/4.)*Hankel(0,1,abx1)*ds_2[j]*J0_T2[j];  
		}
	     	E_J0_T2[i]=Sum_J0;
            Sum_J0=zero;

	}

/////////////////////////目标1引起的感应差值电流在目标1表面的入射场/////////////////////////////////


	for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_1[i]-x[j],2)+pow(ztt_1[i]-z[j],2));
			Sum_Jsd=Sum_Jsd-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Jsd_T1[j];
		}
			E_Jsd_T1_T1[i]=Sum_Jsd;
            Sum_Jsd=zero;
	}

/////////////////////////目标2引起的感应差值电流在目标1表面的入射场/////////////////////////////////


	for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_1[i]-x[j],2)+pow(ztt_1[i]-z[j],2));
			Sum_Jsd=Sum_Jsd-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Jsd_T2[j];
		}
			E_Jsd_T2_T1[i]=Sum_Jsd;
            Sum_Jsd=zero;
	}

	for (i=0;i<M;i++)
	{
		ET0_T1[i] = E_Jsd_T1_T1[i]+E_Jsd_T2_T1[i]+E_J0_T2[i]-E0_T1[i];
	}

//	for (i=0;i<N;i++)
//	{
//		cout<<ZT_1[i]<<endl;
			
//	}

cout<<"计算J0_T1......------>>"<<endl;
	           

	CLinearEquations linear;
	           linear.CG(ZT_1,J0_T1,ET0_T1,M,M,0.01);

/////////////////////////目标1上的感应电流在目标2表面的入射场/////////////////////////////////

	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
			abx1=k_w*sqrt(pow(xtt_2[i]-xtt_1[j],2.0)+pow(ztt_2[i]-ztt_1[j],2.0));
			Sum_J0 = Sum_J0-(k_w*Yita/4.)*Hankel(0,1,abx1)*ds_1[j]*J0_T1[j];  
		}
		    E_J0_T1[i]=Sum_J0;
            Sum_J0=zero;

	}
	
/////////////////////////目标1引起的感应差值电流在目标2表面的入射场/////////////////////////////////


	for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_2[i]-x[j],2)+pow(ztt_2[i]-z[j],2));
			Sum_Jsd=Sum_Jsd-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Jsd_T1[j];
		}
			E_Jsd_T1_T2[i]=Sum_Jsd;
            Sum_Jsd=zero;
	}

/////////////////////////目标2引起的感应差值电流在目标2表面的入射场/////////////////////////////////
	
	for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx1=k_w*sqrt(pow(xtt_2[i]-x[j],2)+pow(ztt_2[i]-z[j],2));
			Sum_Jsd=Sum_Jsd-(k_w*Yita/4.)*Hankel(0,1,abx1)*dl*Jsd_T2[j];
		}
			E_Jsd_T2_T2[i]=Sum_Jsd;
            Sum_Jsd=zero;
	}

	for (i=0;i<M;i++)
	{
		ET0_T2[i] = E_Jsd_T1_T2[i]+E_Jsd_T2_T2[i]+E_J0_T1[i]-E0_T2[i];
	}
	
cout<<"计算J0_T2......------>>"<<endl;

	           linear.CG(ZT_2,J0_T2,ET0_T2,M,M,0.01);

//////////////////////计算目标1引起的差值感应电流//////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			abx1=k_w*sqrt(pow(x[i]-xtt_1[j],2.0)+pow(z[i]-ztt_1[j],2.0));
			Sum_J0 = Sum_J0-(k_w*Yita/4.)*Hankel(0,1,abx1)*ds_1[j]*J0_T1[j];  
		}
		   E_J0_T1_S[i]=Sum_J0;
           Sum_J0=zero;
	}

cout<<"计算Jsd_T1......------>>"<<endl;


     One_D_Conductor_BMIA_CAG_Matrix_Jsd(E_J0_T1_S,Jsd_T1);

//////////////////////计算目标2引起的差值感应电流//////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			abx1=k_w*sqrt(pow(x[i]-xtt_2[j],2.0)+pow(z[i]-ztt_2[j],2.0));
			Sum_J0 = Sum_J0-(k_w*Yita/4.)*Hankel(0,1,abx1)*ds_2[j]*J0_T2[j];  
		}
		    E_J0_T2_S[i]=Sum_J0;
            Sum_J0=zero;
	}

cout<<"计算Jsd_T2......------>>"<<endl;


	 One_D_Conductor_BMIA_CAG_Matrix_Jsd(E_J0_T2_S,Jsd_T2);
	
///////////////////////判断收敛条件///////////////////////////////
    

	       max_J0_T1 =zero;
           max1_J0_T1 =zero;

	       for (i=0;i<M;i++)
		   {
		        if (abs( max_J0_T1)<abs(J0_T1_1[i]-J0_T1[i]))
			                max_J0_T1=J0_T1_1[i]-J0_T1[i];
				else 
					continue;
			
		   }

		    for (i=0;i<M;i++)
			{
		       J0_T1_1[i] = J0_T1[i];
			}

		   
		   for (i=0;i<M;i++)
		   {
		        if (abs(max1_J0_T1)<abs(J0_T1[i]))
			               max1_J0_T1 =J0_T1[i];
				else 
					continue;
			
		   }
		   
		   max_J0_T2 =zero;
           max1_J0_T2 =zero;

	       for (i=0;i<M;i++)
		   {
		        if (abs( max_J0_T2)<abs(J0_T2_1[i]-J0_T2[i]))
			                max_J0_T2=J0_T2_1[i]-J0_T2[i];
				else 
					continue;
			
		   }

		    for (i=0;i<M;i++)
			{
		       J0_T2_1[i] = J0_T2[i];
			}

		   
		   for (i=0;i<M;i++)
		   {
		        if (abs(max1_J0_T2)<abs(J0_T2[i]))
			               max1_J0_T2=J0_T2[i];
				else 
					continue;
			
		   }





		      Delta1 = abs(max_J0_T1)/abs(max1_J0_T1);
			  Delta2 = abs(max_J0_T2)/abs(max1_J0_T2);
              
			  
///////////////////////////////两步收敛条件///////////////////////////////////////			  
			  
			  if (Delta1>Delta2)
				    Delta = Delta1;
			  else
				    Delta = Delta2;
			                			  			  			  			  			  			  			  			  
			  cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;
	          			  			  			  
			  if (Delta<0.001)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;

	  
	 
	}
     
////////////////////////////////////输出计算结果////////////////////////////////////////////
	

    ofstream outJsd_T1;
		      outJsd_T1.open("Jsd_T1.txt");
		      
			  for(i=0;i<N;i++)
			  {
				 outJsd_T1<<x[i]/Lamda<<"     "<<real(Jsd_T1[i])<<"     "<<imag(Jsd_T1[i])<<"     "<<abs(Jsd_T1[i])<<endl;
			  }

	ofstream outJsd_T2;
			  outJsd_T2.open("Jsd_T2.txt");
			  
			  for (i=0;i<N;i++)
			  {
			   	 outJsd_T2<<x[i]/Lamda<<"     "<<real(Jsd_T2[i])<<"     "<<imag(Jsd_T2[i])<<"     "<<abs(Jsd_T2[i])<<endl;

			  }

    ofstream outJ0_T1;
			  outJ0_T1.open("J0_T1.txt");
			  
			  for (i=0;i<M;i++)
			  {
			   	 outJ0_T1<<i<<"     "<<real(J0_T1[i])<<"     "<<imag(J0_T1[i])<<"     "<<abs(J0_T1[i])<<endl;

			  }

	ofstream outJ0_T2;
			  outJ0_T2.open("J0_T2.txt");


			  for (i=0;i<M;i++)
			  {
				  if (i<270)
				  {
			   	     outJ0_T2<<i<<"     "<<real(J0_T2[i+90])<<"     "<<imag(J0_T2[i+90])<<"     "<<abs(J0_T2[i+90])<<endl;
				  }
				  else
   			   	     outJ0_T2<<i<<"     "<<real(J0_T2[i-270])<<"     "<<imag(J0_T2[i-270])<<"     "<<abs(J0_T2[i-270])<<endl;
			  }
			  
			  
	ofstream outJs0_1;
			  outJs0_1.open("Js0_1.txt");
			  
			  for (i=0;i<N;i++)
			  {
			   	 outJs0_1<<x[i]/Lamda<<"     "<<real(Js0_1[i])<<"     "<<imag(Js0_1[i])<<"     "<<abs(Js0_1[i])<<endl;

			  }

	ofstream outJs0_2;
			  outJs0_2.open("Js0_2.txt");
			  
			  for (i=0;i<N;i++)
			  {
			   	 outJs0_2<<x[i]/Lamda<<"     "<<real(Js0_2[i])<<"     "<<imag(Js0_2[i])<<"     "<<abs(Js0_2[i])<<endl;

			  }
	
}

void Surf_Targ_Question::Gauss_Surface_medium_Cylinder_Target_J()
{

	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=0.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=102.4*wave_len;
	double R=5.*wave_len;
	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double dl;
	double abx;
	double Delta;
	
	double Rm;
	double sn;
	double cn;

	Complex zero(0.,0.);


	Complex Epcle(2.,50.);   //介电常数
	Complex Sum_E_Jsd;
	Complex Sum_E_JT;
	Complex Sum_E_MT;
	Complex SE;
    Complex kk;

	int i,j,k;

    Lamda=3e8/3e9;
	K=2*PI/Lamda;
	Yita=120.*PI;
    Complex Kd;
	Complex Yita_d;
      Kd=K*sqrt(Epcle);
	  Yita_d=Yita/sqrt(Epcle);
	Complex max;
	Complex max1;
 
	
	Vec_Com ET0;
	   ET0.assign(2*M);
	Vec_Com JT;	
	   JT.assign(2*M);
	Vec_Com JT1;	
	   JT1.assign(2*M);
	Vec_Com JT2;	
	   JT2.assign(2*M);
	Vec_Com ET;
	   ET.assign(2*M);   
	Vec_Com ET_Js0;
	   ET_Js0.assign(M);
	Vec_Com ET_Jsd;
	   ET_Jsd.assign(M);

   	Vec_Com ES0;
	   ES0.assign(N);
	Vec_Com Es_JT;
	   Es_JT.assign(N);
	Vec_Com Es_MT;
	   Es_MT.assign(N);

	Vec_Com Jsd;
	   Jsd.assign(N);
	Vec_Com Jsd1;
	   Jsd1.assign(N);
	Vec_Com Vector_temp;


	Vec_Com Z_T_S;
	   Z_T_S.assign(2*M*2*M);

	cout<<"进行目标的剖分......------>>"<<endl;
    Cy_PoufenMatix();
   
	 
	cout<<"计算Js0......------>>"<<endl;

	One_D_Conductor_Gauss_Js0(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);	


    dx_Suf = RL/N;

    for(i=0;i<M;i++)
	{
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx=K*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
		    ET_Js0[i]=ET_Js0[i]-(K*Yita/4.)*Hankel(0,1,abx)*dl*Js0[j];
		}
	}

	
     cout<<"开始填充矩阵......------>>"<<endl;

////////////////////////////////////////ZJtt_out//////////////////////////////////////////////
	   
	  for (i=0;i<M;i++)
	  {
		for (j=0;j<M;j++)
		{
		
	    	if(i!=j)
			{
				kk=K*sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
				Z_T_S[i*2*M+j]=(K*Yita/4.)*ds*Hankel(0,1,kk);
			}
			else
			{
                Z_T_S[i*2*M+j]=(K*Yita/4.)*ds*(1.+2.*J/PI*(log(Gama*K*ds/4.)-1.));	
			}
		 }

	  }

/////////////////////////////////////////ZKtt_out/////////////////////////////////////////////


    	for (i=0;i<M;i++)
		{
		  for(j=M;j<2*M;j++)
		  {
			if (j-i==M)
			{
                Z_T_S[i*2*M+j]=0.5;				
			}
			else
			{
				Rm=sqrt(pow(xtt[i]-xtt[j-M],2.0)+pow(ztt[i]-ztt[j-M],2.0));
				if (j<2*M-1)
				{
			       cn=(xtt[j+1-M]-xtt[j-M])/sqrt(pow(xtt[j+1-M]-xtt[j-M],2.0)+pow(ztt[j+1-M]-ztt[j-M],2.0)); 
			       sn=(ztt[j+1-M]-ztt[j-M])/sqrt(pow(xtt[j+1-M]-xtt[j-M],2.0)+pow(ztt[j+1-M]-ztt[j-M],2.0)); 
				}
				else
				{
                   cn=(xtt[0]-xtt[j-M])/sqrt(pow(xtt[0]-xtt[j-M],2.0)+pow(ztt[0]-ztt[j-M],2.0)); 
			       sn=(ztt[0]-ztt[j-M])/sqrt(pow(xtt[0]-xtt[j-M],2.0)+pow(ztt[0]-ztt[j-M],2.0)); 
				}
				
			   
				Z_T_S[i*2*M+j]=(J*K/4.)*ds*((ztt[i]-ztt[j-M])/Rm*cn-(xtt[i]-xtt[j-M])/Rm*sn)*Hankel(1,1,K*Rm);
			}
		  }
		
		} 

/////////////////////////////////////////ZJtt_ins/////////////////////////////////////////////

	for (i=M;i<2*M;i++)
	{
		for (j=0;j<M;j++)
		{
					
	    	if(i-j==M)
			{
			    Z_T_S[i*2*M+j]=(Kd*Yita_d/4.)*ds*(1.+2.*J/PI*(log(Gama*Kd*ds/4.)-1.));
			}
								
			else 
			{
               kk=Kd*sqrt(pow(xtt[i-M]-xtt[j],2.0)+pow(ztt[i-M]-ztt[j],2.0));
			   Z_T_S[i*2*M+j]=(Kd*Yita_d/4.)*ds*Hankel(0,1,kk);
			}

		}
	}

////////////////////////////////////////////////////ZKtt_ins///////////////////////////////////////////////////////


	for (i=M;i<2*M;i++)
	{
		for (j=M;j<2*M;j++)
		{
			if (i==j)
			{ 
				Z_T_S[i*2*M+j]=-0.5;
		
			}
			              			
			else
			{
				Rm=sqrt(pow(xtt[i-M]-xtt[j-M],2.0)+pow(ztt[i-M]-ztt[j-M],2.0));
				if (j<2*M-1)
				{				
			    cn=(xtt[j+1-M]-xtt[j-M])/sqrt(pow(xtt[j+1-M]-xtt[j-M],2.0)+pow(ztt[j+1-M]-ztt[j-M],2.0));
				sn=(ztt[j+1-M]-ztt[j-M])/sqrt(pow(xtt[j+1-M]-xtt[j-M],2.0)+pow(ztt[j+1-M]-ztt[j-M],2.0)); 
				}
				else
				{
				cn=(xtt[0]-xtt[j-M])/sqrt(pow(xtt[0]-xtt[j-M],2.0)+pow(ztt[0]-ztt[j-M],2.0)); 
			    sn=(ztt[0]-ztt[j-M])/sqrt(pow(xtt[0]-xtt[j-M],2.0)+pow(ztt[0]-ztt[j-M],2.0)); 

				}

				Z_T_S[i*2*M+j]=	(J*Kd/4.)*ds*((ztt[i-M]-ztt[j-M])/Rm*cn-(xtt[i-M]-xtt[j-M])/Rm*sn)*Hankel(1,1,Kd*Rm);
			}
		}
	}

	cout<<"矩阵填充完毕。"<<endl;


   
	
    for (k=0;k<50;k++)
	{

//////////////////////////开始进行粗糙面与目标的耦合迭代////////////////////////////////////
      
	cout<<"进行粗糙面与目标的第"<<k+1<<"次耦合迭代......------>>"<<endl;

      for(i=0;i<M;i++)
	  {
	    for (j=0;j<N;j++)
		{
			dl=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			abx=K*sqrt(pow(xtt[i]-x[j],2)+pow(ztt[i]-z[j],2));
		    Sum_E_Jsd=Sum_E_Jsd-(K*Yita/4.)*Hankel(0,1,abx)*dl*Jsd[j];
		}
		    ET_Jsd[i]=Sum_E_Jsd;
			Sum_E_Jsd=zero;
		
	  }

	   for (i=0;i<M;i++)
	   {

		ET0[i]=exp(J*K*(xtt[i]*sin(Sita0)-ztt[i]*cos(Sita0)))+ET_Js0[i]+ET_Jsd[i];
			
	
	   }




    Bi_CG(Z_T_S,ET0,JT);

//////////////////////////////////////////////Esf_out//////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
		  	kk=K*sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
			Sum_E_JT = Sum_E_JT-(K*Yita/4.)*ds*JT[j]*Hankel(0,1,kk);
		}
		    Es_JT[i] = Sum_E_JT;
			Sum_E_JT = zero;
	}

    for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{

	        	Rm=sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));

				if (j<M-1)
				{
				  cn=(xtt[j+1]-xtt[j])/sqrt(pow(xtt[j+1]-xtt[j],2.0)+pow(ztt[j+1]-ztt[j],2.0)); 
			      sn=(ztt[j+1]-ztt[j])/sqrt(pow(xtt[j+1]-xtt[j],2.0)+pow(ztt[j+1]-ztt[j],2.0)); 
				}
				else
				{
				   cn=(xtt[0]-xtt[j])/sqrt(pow(xtt[0]-xtt[j],2.0)+pow(ztt[0]-ztt[j],2.0)); 
			       sn=(ztt[0]-ztt[j])/sqrt(pow(xtt[0]-xtt[j],2.0)+pow(ztt[0]-ztt[j],2.0)); 
				}

			Sum_E_MT = Sum_E_MT-(J*K/4.)*JT[j+M]*ds*((z[i]-ztt[j])/Rm*cn-(x[i]-xtt[j])/Rm*sn)*Hankel(1,1,K*Rm);
		}
		    Es_MT[i] = Sum_E_MT;
			Sum_E_MT = zero;
	}

	

	for (i=0;i<N;i++)
	{  
         ES0[i]= Es_JT[i]+Es_MT[i];
	}

	

	One_D_Conductor_BMIA_CAG_Matrix_Jsd(ES0,Jsd);

////////////////////////////////判//////////////////断收敛条件/////////////////////////////////////////////////
	
/*
      max =zero;
      max1 =zero;

	       for (i=0;i<2*M;i++)
		   {
		        if (abs(max)<abs(JT1[i]-JT[i]))
			               max=JT1[i]-JT[i];
				else 
					continue;
			
		   }

		    for (i=0;i<2*M;i++)
			{
		       JT1[i] = JT[i];
			}

		   
		   for (i=0;i<2*M;i++)
		   {
		        if (abs(max1)<abs(JT[i]))
			               max1=JT[i];
				else 
					continue;
			
		   }

		     


		      Delta = abs(max)/abs(max1);

 */             

           CLinearEquations num0;

		    for (i=0;i<2*M;i++)
		    {
				JT2[i] = JT[i]-JT1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',Z_T_S,JT2,Vector_temp);
              Assign_Vecor(JT,JT1);
              Delta = num0.VectorNorm2(Vector_temp,2*M)/num0.VectorNorm2(ET0,2*M);

              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.01)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;

			  
	     


    }

///////////////////////////////////////输出计算结果//////////////////////////////////////////////

	ofstream outJs0;
	         outJs0.open("medium_Js0.txt");
			 for (i=0;i<N;i++)
			 {
				 outJs0<<x[i]/Lamda<<"   "<<real(Js0[i])<<"   "<<imag(Js0[i])<<"   "<<abs(Js0[i])<<endl;
			 }

	ofstream outJsd;
	         outJsd.open("medium_Jsd.txt");
			 for (i=0;i<N;i++)
			 {
				 outJsd<<x[i]/Lamda<<"   "<<real(Jsd[i])<<"   "<<imag(Jsd[i])<<"   "<<abs(Jsd[i])<<endl;
			 }

	ofstream outJt;
	         outJt.open("medium_Jt.txt");
			 for (i=0;i<M;i++)
			 {
				 outJt<<i<<"   "<<real(JT[i])<<"   "<<imag(JT[i])<<"   "<<abs(JT[i])<<endl;
			 }

 	ofstream outKt;
	         outKt.open("medium_Kt.txt");
			 for (i=M;i<2*M;i++)
			 {
				 outKt<<i-M<<"   "<<real(JT[i])<<"   "<<imag(JT[i])<<"   "<<abs(JT[i])<<endl;
			 }
	
////////////////////////////////计算差场散射截面/////////////////////////////

   Vec_dou Jsd_IN_1;
	  Jsd_IN_1.assign(4);
   Vec_dou JT0_IN_1;
	  JT0_IN_1.assign(4);
   Vec_Com Jsd_IN;
	  Jsd_IN.assign(N);
   Vec_Com JT0_IN;
	  JT0_IN.assign(N); 
   Vec_Com kt_IN;
	  kt_IN.assign(N); 
   Vec_Com kt_IN_1;
	  kt_IN_1.assign(N);  

//////////////////////////////////读取数据//////////////////////////////////////////	  
	  
   ifstream in_Jsd;
	      in_Jsd.open("medium_Jsd.txt");
	    for(i=0;i<N;i++)
		{
			for (j=0;j<4;j++)
			{
				in_Jsd>>Jsd_IN_1[j];

			}
		Jsd_IN[i] = Jsd_IN_1[1]+J*Jsd_IN_1[2];

		}

	ifstream in_JT0;
	      in_JT0.open("medium_JT.txt");

       for (i=0;i<M;i++)
       {
		   for (j=0;j<4;j++)
			{
				in_JT0>>JT0_IN_1[j];

			}
		JT0_IN[i] = JT0_IN_1[1]+J*JT0_IN_1[2];
		   
       } 
   ifstream in_medium_kT;
	      in_medium_kT.open("medium_kT.txt");

       for (i=0;i<M;i++)
       {
		   for (j=0;j<4;j++)
			{
				in_JT0>>kt_IN_1[j];

			}
		kt_IN[i] = kt_IN_1[1]+J*kt_IN_1[2];
		   
       } 



//////////////////////////////////////////////////////////////////////

	int N0=180;
	double R_Dsc = 100000.*Lamda;
	double Sita_Dsc = 0.;
	Complex Sum_Esd;
	double dl_Esd=0.;
	double abx_Esd_S=0.;
    double abx_Esd_T=0.;
	
	Vec_Com Esd;
	   Esd.assign(N0);
	Vec_dou x_Dsc;
	  x_Dsc.assign(N0);
	Vec_dou z_Dsc;
	  z_Dsc.assign(N0);
	Vec_Com Ei_Dsc;
	  Ei_Dsc.assign(N0);
	Vec_dou Dsc;
	  Dsc.assign(N0);


	for (i=N0-1;i>=0;i--)
	{  
		Sita_Dsc = i*PI/N0;
		
		x_Dsc[N0-i-1] = R_Dsc*cos(Sita_Dsc);
		z_Dsc[N0-i-1] = R_Dsc*sin(Sita_Dsc);

	}
   
	i_an = i_an*PI/180;
//	cout<<"i_an="<<i_an<<endl;
 
/*	double xsz,ww,kg,fac;
    for(i=0;i<N0;i++)
	{
	//	xsz=x_Dsc[i]*sin(i_an)-z_Dsc[i]*cos(i_an);
	//	kg=pow(k_w*g_p*cos(i_an),2.0);
	//	fac=pow((x_Dsc[i]+H*tan(i_an)+z_Dsc[i]*tan(i_an))/g_p,2.0);
	//	ww=(2.0*fac-1.0)/kg;
	//	Ei_Dsc[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
		Ei_Dsc[i]=exp(i_*k_w*xsz);

		
	}

	ofstream outEi;
	  outEi.open("Ei.txt");
	for (i=0;i<N0;i++)
	{
		outEi<<i-90<<"   "<<real(Ei_Dsc[i])<<"   "<<imag(Ei_Dsc[i])<<"   "<<abs(Ei_Dsc[i])<<endl;
	}*/

	dx_Suf=0.01;


	for (k=0;k<N0;k++)
	{

        for (i=0;i<M;i++)
        {
			abx_Esd_T = k_w*sqrt(pow(x_Dsc[k]-xtt[i],2)+pow(z_Dsc[k]-ztt[i],2));
			Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_T)*ds*JT0_IN[i];						
		}

		for (j=0;j<N;j++)
		{
			 dl_Esd=dx_Suf*sqrt(1.0+pow(df_x[j],2));
			 abx_Esd_S=k_w*sqrt(pow(x_Dsc[k]-x[j],2)+pow(z_Dsc[k]-z[j],2));
			 Esd[k] = Esd[k] - (k_w*Yita/4.)*Hankel(0,1,abx_Esd_S)*dl_Esd*Jsd_IN[j];
		}
		
	}
	

	for (i=0;i<N0;i++)
	{
		for (j=0;j<M;j++)
		{

	        	Rm=sqrt(pow(x_Dsc[i]-xtt[j],2.0)+pow(x_Dsc[i]-ztt[j],2.0));

				if (j<M-1)
				{
				  cn=(xtt[j+1]-xtt[j])/sqrt(pow(xtt[j+1]-xtt[j],2.0)+pow(ztt[j+1]-ztt[j],2.0)); 
			      sn=(ztt[j+1]-ztt[j])/sqrt(pow(xtt[j+1]-xtt[j],2.0)+pow(ztt[j+1]-ztt[j],2.0)); 
				}
				else
				{
				   cn=(xtt[0]-xtt[j])/sqrt(pow(xtt[0]-xtt[j],2.0)+pow(ztt[0]-ztt[j],2.0)); 
			       sn=(ztt[0]-ztt[j])/sqrt(pow(xtt[0]-xtt[j],2.0)+pow(ztt[0]-ztt[j],2.0)); 
				}

			Esd[i] = Esd[i]-(J*K/4.)*kt_IN[j+M]*ds*((x_Dsc[i]-ztt[j])/Rm*cn-(x_Dsc[i]-xtt[j])/Rm*sn)*Hankel(1,1,K*Rm);
		}
		    
	}




	for (i=0;i<N0;i++)
	{  
//		Dsc[i] = PI*R_Dsc*pow(abs(Esd[i]),2.0)/pow(abs(Ei_Dsc[i]),2.0);
      Dsc[i] = 2.*PI*R_Dsc*pow(abs(Esd[i]),2.0);
	}

    cout<<"输出差场散射计算结果......------>>"<<endl;
	ofstream outEsd;
	 outEsd.open("medium_Esd.txt");
	
	ofstream outDsc;
	 outDsc.open("medium_Dsc.txt");
	 
	ofstream outDscDB;
	 outDscDB.open("medium_DscDB.txt");

		 for (i=0;i<N0;i++)
		 {
			 outDsc<<i-90<<" "<<Dsc[i]<<endl;
			 outDscDB<<i-90<<" "<<10.*log10(Dsc[i])<<endl;
			 outEsd<<i-90<<" "<<abs(Esd[i])<<endl;
		 }	
		 



}




void Surf_Targ_Question::Gauss_Del_Surface_Del_Cylinder__Above_Target()
{

	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e8;      //wave length
	double h=0.5*wave_len;            //rms height 
	double lc=1.0*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=5.*wave_len;
	double R=1.*wave_len;
		   H=3.*wave_len;

		   
	double g_p = 4.*R/cos(i_an*PI/180.);
//		   g_p = 20.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double e=2.71828138;
	double GAMA = 0.577216;
    
	double Epsle_0 = 8.85e-7;
	Complex Epsle(9.8,1.0);
	Complex Epsle_Target(3,0.03);
//	Complex rou_0= Epsle;
//	Complex rou_Target = Epsle_Target;
	Complex rou_0= 1;
    Complex rou_Target = 1;

//	cout<<"rou="<<rou_0<<endl;
	

	Complex Sum(0.,0);
	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum_V3(0.,0.);
	Complex Sum_V4(0.,0.);

	double abx;
	double NR;
	double Rm;
	double dl;
	double Delta;
	double Delta0;
	double Delta1;

	

 	int N0 = 181;
	int i,j,k,m,l;
	Vec_Com V1; V1.assign(N);	  
	Vec_Com V2; V2.assign(N);	  
	Vec_Com V3; V3.assign(2*M);	  
	Vec_Com HT_S;	HT_S.assign(M);	
	Vec_Com V4; V4.assign(M);
	Vec_Com HV1_T; HV1_T.assign(M);	   
	Vec_Com HV2_T;  HV2_T.assign(M);	  
	Vec_Com HT0;  HT0.assign(2*M);	  
	Vec_Com H0;  H0.assign(M);	  
	Vec_Com H_in; H_in.assign(N);	  
    Vec_Com HS0; HS0.assign(2*N);	  
	Vec_Com V12; V12.assign(2*N);	  
	Vec_Com V12_1;V12_1.assign(2*N);	   
	Vec_Com V12_2;V12_2.assign(2*N);	   
	Vec_Com Vector_temp;Vector_temp.assign(2*N);	   
	Vec_Com  Hs; Hs.assign(N0);	  
	Vec_Com Hs_T; Hs_T.assign(N);	  
	Vec_Com Hs_S; Hs_S.assign(N);
	Vec_dou Sita_Dsc; Sita_Dsc.assign(N0);	   
	Vec_Com ZT_S; ZT_S.assign(2*N*2*N);
	Vec_Com ZT_T; ZT_T.assign(2*M*2*M);
	   

   Vec_Com ZT_MOM;  ZT_MOM.assign((2*N+2*M)*(2*N+2*M));
   Vec_Com JJ_MOM;  JJ_MOM.assign(2*N+2*M);
   Vec_Com EE_MOM;  EE_MOM.assign(2*N+2*M);

   Vec_Com ZT_A; ZT_A.assign(N*N);           
   Vec_Com ZT_B; ZT_B.assign(N*N);          
   Vec_Com ZT_C; ZT_C.assign(N*M);           
   Vec_Com ZT_D; ZT_D.assign(N*M);           
   Vec_Com ZT_E; ZT_E.assign(N*N);          
   Vec_Com ZT_F; ZT_F.assign(N*N);           
   Vec_Com ZT_G; ZT_G.assign(N*M);           
   Vec_Com ZT_H; ZT_H.assign(N*M);           
   Vec_Com ZT_I; ZT_I.assign(M*M);         
   Vec_Com ZT_J; ZT_J.assign(M*M);           		   
   Vec_Com ZT_K; ZT_K.assign(M*M);
   Vec_Com ZT_L; ZT_L.assign(M*M);

    Vec_Com Jf_1; Jf_1.assign(N);
	Vec_Com Jf_2; Jf_2.assign(N);
	Vec_Com Jb_1; Jb_1.assign(N);
    Vec_Com Jb_2; Jb_2.assign(N);

	Vec_Com Jf1;  Jf1.assign(2*N);
	Vec_Com Jf2;  Jf2.assign(2*N);
	Vec_Com Jf3;  Jf3.assign(2*N);

	Vec_Com J_Target_1; J_Target_1.assign(2*M);
	Vec_Com J_Target_2; J_Target_2.assign(2*M);
	Vec_Com Delta_V3;   Delta_V3.assign(2*M);

	Vec_dou Delta_Matrix; Delta_Matrix.assign(15);


	cout<<"进行目标的剖分......------>>"<<endl;

	Lamda=wave_len; 
 
		

	ds=2.*PI*R/M;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;
    Yita=120.*PI;

	for( i=0;i<M;i++)
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
   
//   Cy_PoufenMatix_above_target();

	Complex K1 = K*sqrt(Epsle);
	Complex K_Target = K*sqrt(Epsle_Target);
	
    cout<<"生成粗糙面......------>>"<<endl;	
  	
    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    dx_Suf=RL/N;
	
/*	for (i=0;i<N;i++)
	{
		z[i] = 0.;
		df_x[i] = 0.;
		ddf_x[i] = 0.;
	}*/



	double xsz,ww,kg,fac;
//	double d=4.*Lamda;

	i_an = i_an*PI/180;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+H*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		H_in[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}


	
	for (i=0;i<M;i++)
	{
		xsz=xtt[i]*sin(i_an)-ztt[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((xtt[i]+H*tan(i_an)+ztt[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
 //		H0[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
        H0[i] = exp(i_*k_w*xsz);
		 
	}






	cout<<" 开始填充阻抗矩阵......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

/*	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = (x[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = -1.*dl*(J*K/4.)*NR*Hankel(1,1,K*Rm);
//			   cout<<ZT_S[i*2*N+j]<<endl;
			}
			else
			{
			   ZT_S[i*2*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
		}
	}*/


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_A[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);
//			   cout<<ZT_S[i*2*N+j]<<endl;
			}
			else
			{
			   ZT_A[i*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
			
		}
	}

//////////////////////////////////////////////////////////////////////////

/*	for (i=0;i<N;i++)
	{
		for (j=N;j<2*N;j++)
		{
			if (j-i==N)
			{
               ZT_S[i*2*N+j] = dx_Suf*sqrt(1+pow(df_x[i],2.0))*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm            = sqrt(pow(x[i]-x[j-N],2.0)+pow(z[i]-z[j-N],2.0));
			   dl           = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_S[i*2*N+j] = dl*J/4.*Hankel(0,1,K*Rm);

			}			
		}
		
	}*/


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (j==i)
			{
               ZT_B[i*N+j] = dx_Suf*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_B[i*N+j]   = dx_Suf*J/4.*Hankel(0,1,K*Rm);

			}			

		}
	}
//////////////////////////////////////////////////////////////////////////

for (i=0;i<N;i++)
{
	for (j=0;j<M;j++)
	{
					       
		Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
//		NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(z[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
//		ZT_I[i*M+j] = -1.*ds*NR*J*K/4.*Hankel(1,1,K*Rm);
		ZT_D[i*M+j] = ds*J/4.*Hankel(0,1,K*Rm);
			
	}
}

for (i=0;i<N;i++)
{
	for (j=0;j<M;j++)
	{
		Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
		NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(z[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		ZT_C[i*M+j] = -1.*NR*J*K/4.*ds*Hankel(1,1,K*Rm);
	}
}

//////////////////////////////////////////////////////////////////////////

/*   for (i=N;i<2*N;i++)
   {
	   for(j=0;j<N;j++)
	   {
		   if (i-j==N)
		   {
			   ZT_S[i*2*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i-N]-x[j],2.0)+pow(z[i-N]-z[j],2.0));
			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = (x[i-N]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i-N]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = dl*Rm*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }
	   }
   }*/

   for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   		   
		   if (i==j)
		   {
			   ZT_E[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_E[i*N+j]   = dx_Suf*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }

	   }
   }
   

//////////////////////////////////////////////////////////////////////////

/*   for (i=N;i<2*N;i++)
   {
	   for (j=N;j<2*N;j++)
	   {
		   if (i==j)
		   {
              ZT_S[i*2*N+j] = -1.*dx_Suf*rou_0*sqrt(1+pow(df_x[j-N],2.0))*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j-N],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i-N]-x[j-N],2.0)+pow(z[i-N]-z[j-N],2.0));
			   dl            = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_S[i*2*N+j] = -1.*rou_0*dl*J/4.*Hankel(0,1,K1*Rm);

	   
		   }
	   }
   }*/

   for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   		   
		   if (i==j)
		   {
              ZT_F[i*N+j] = -1.*dx_Suf*rou_0*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_F[i*N+j]   = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);

	   
		   }

	   }
   }
//////////////////////////////////////////////////////////////////////////
for (i=0;i<M;i++)
{
	for (j=0;j<N;j++)
	{
		 			  
		       Rm          = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = ((xtt[i]-x[j])*(-1.*df_x[j])+(ztt[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_G[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);			   
		
	}
}

for (i=0;i<M;i++)
{
	for (j=0;j<N;j++)
	{			       
		       Rm            = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_H[i*N+j]   = dx_Suf*J/4.*Hankel(0,1,K*Rm);	 			   				       		
	}
}



//////////////////////////////////////////////////////////////////////////


/*       for (i=0;i<M;i++)
		  {
	            for(j=0;j<M;j++)
				{
		          if (i!=j)
				  {
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = (xtt[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(ztt[i]-ztt[j])/Rm/(1+pow(df_z[j],2.0));
		            ZT_T[i*M+j] = -1.*ds*NR*J*K/4.*Hankel(1,1,K*Rm); 
				  }
		          else
				  {
			        ZT_T[i*M+j] = 0.5-ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));

				  }
		   
		   
				}
		  }*/

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_I[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K*Rm);
					ZT_J[i*M+j] = J/4.*ds*Hankel(0,1,K*Rm);
				  }
		          else
				  {

					ZT_I[i*M+j] = 0.5-ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_J[i*M+j] = J/4.*ds*(1.+J*2./PI*log(exp(GAMA)*ds*K/(4.*e)));
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_K[i*M+j] = ds*NR*J*K_Target/4.*Hankel(1,1,K_Target*Rm);
					ZT_L[i*M+j] = -rou_Target*ds*J/4.*Hankel(0,1,K_Target*Rm);
				  }
		          else
				  {

					ZT_K[i*M+j] = 0.5+ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_L[i*M+j] = -rou_Target*ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K_Target/(4.*e)));
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_I[i*M+j];
		   }
		   for (j=M;j<2*M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_J[i*M+j-M];
		   }
	   }

	   for (i=M;i<2*M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
               ZT_T[i*2*M+j] = ZT_K[(i-M)*M+j];
		   }
		   for (j=M;j<2*M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_L[(i-M)*M+j-M];
		   }
	   }
	   
/*	   for (i=0;i<M;i++)
	   {
		   HT0[i] = H0[i];
	   }*/

//	   Bi_CG(ZT_T,H0,V3);
//	   Bi_CG(ZT_L,HT0,V4);



//////////////////////////////////////////////////////////////////////////
  int N_M = 2*N+2*M;
	   
  for (i=0;i<N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_A[i*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_B[i*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_C[i*M+j-2*N];
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_D[i*M+j-2*N-M];
	  }

  }

  for (i=N;i<2*N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_E[(i-N)*N+j]; 
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_F[(i-N)*N+j-N];
	  }

  }

  for (i=2*N;i<2*N+M;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_G[(i-2*N)*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_H[(i-2*N)*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
           ZT_MOM[i*N_M+j] = ZT_I[(i-2*N)*M+j-2*N];
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_J[(i-2*N)*M+j-2*N-M];
	  }
  }
  for (i=2*N+M;i<N_M;i++)
  {
	  for (j=2*N;j<2*N+M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_K[(i-2*N-M)*M+j-2*N];		  
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_L[(i-2*N-M)*M+j-2*N-M];
	  }
  }

  for (i=0;i<N;i++)
  {
	  EE_MOM[i] = H_in[i];
  }
  for (i=0;i<M;i++)
  {
	  EE_MOM[i+2*N] = H0[i];
  }
  
   cout<<"阻抗矩阵填充完毕。"<<endl; 

//     Bi_CG(ZT_MOM,EE_MOM,JJ_MOM);
//   Bi_CGSTAB(ZT_MOM,EE_MOM,JJ_MOM);
	 GAUSS(ZT_MOM,JJ_MOM,EE_MOM,N_M);

//////////////////////////////////////////////////////////////////////////

        


for (l=0;l<15;l++)
{

cout<<"Begin the "<<l+1<<"step of the iteration......------>>"<<endl;
  
   for (m=0;m<N;m++)
   {


	for (i=0;i<N;i++)
	{
		for (j=0;j<i;j++)
		{
						
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_E[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_F[i*N+j];
			
		}
	
		Jf_1[i] = ((H_in[i]-Sum0-Hs_T[i])*ZT_F[i*N+i]-(-1.*Sum1)*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Jf_2[i] = ((-1.*Sum1)*ZT_A[i*N+i]-(H_in[i]-Sum0-Hs_T[i])*ZT_E[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
        

		Sum0 = zero;
		Sum1 = zero;
	}

		 
	 for (i=N-1;i>=0;i--)
	 {
		for (j=i+1;j<N;j++)
		{
		    Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_E[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_F[i*N+j];	
		}


     	Jb_1[i] = ((-1.*Sum0)*ZT_F[i*N+i]-(-1.*Sum1)*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*ZT_A[i*N+i]-(-1.*Sum0)*ZT_E[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;

	}

    		 CLinearEquations num;
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta0 = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta0 ="<<Delta0<<endl;
			  
			 if(Delta0<=0.01)
			 {
			    cout<<"The total number of the FBM step is :"<<m+1<<"."<<"and the final error is :"<<Delta0<<endl;
		        break;
			 }
				
	          else  
		           continue;



   }
   

  
  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   Sum0 = Sum0 + ZT_G[i*N+j]*(Jf_1[j]+Jb_1[j])+ZT_H[i*N+j]*(Jf_2[j]+Jb_2[j]);
	   }
	   HT_S[i] = Sum0;
	   Sum0 = zero;
   }
   for (i=0;i<M;i++)
   {
	   HT0[i] = H0[i]-HT_S[i];
   }
 //    Bi_CG(ZT_T,HT0,V3);

    GAUSS(ZT_T,V3,HT0,2*M);
//   Bi_CG(ZT_J,HT0,V3);
	 
//    GAUSS(ZT_J,V3,H0,M);
    
   for (i=0;i<N;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   Sum0 = Sum0 + ZT_C[i*M+j]*V3[j]+ZT_D[i*M+j]*V3[j+M];
	   }
	   Hs_T[i] = Sum0;
	   Sum0 = zero;
   }
   
		  

 		     for (i=0;i<2*M;i++)
			 {
				 J_Target_2[i] = J_Target_1[i] - V3[i];
			 }
			 for (i=0;i<M;i++)
			 {
				 for (j=0;j<M;j++)
				 {
					  Sum_V3 = Sum_V3+ZT_I[i*M+j]*J_Target_2[j]+ZT_J[i*M+j]*J_Target_2[j+M];
					  Sum_V4 = Sum_V4+ZT_K[i*M+j]*J_Target_2[j]+ZT_L[i*M+j]*J_Target_2[j+M];
				 }
                  Delta_V3[i] = Sum_V3;
				  Delta_V3[i+M] = Sum_V4;
				  Sum_V3 = zero;
                  Sum_V4 = zero;
			 }
			  CLinearEquations num;

             Assign_Vecor(V3,J_Target_1);
             Delta1 = num.VectorNorm2(Delta_V3,2*M)/num.VectorNorm2(HT_S,M);

			 Delta_Matrix[l] = Delta1;

             cout<<"The error of the iteration is: Delta ="<<Delta1<<endl;
			  
			 if(Delta1<=0.01)
			 {
			    cout<<"The total number of the iteration step is :"<<l+1<<"."<<"and the final error is :"<<Delta1<<endl;
		        break;
			 }
				
	          else  
		           continue;

}
/*
   for (k=0;k<1;k++)
   {

        cout<<"开始粗糙面与目标的第"<<k+1<<"次耦合迭代......------>>"<<endl;

        for (i=0;i<N;i++)
		  {
	         for (j=0;j<M;j++)
			 {
		        Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
		        NR = (x[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(z[i]-ztt[j])/Rm/sqrt(1+pow(df_z[j],2.0));
		        Sum = Sum - ds*(J*K/4.)*NR*V3[j]*Hankel(1,1,K*Rm);

			 }
   
                HT_S[i] = Sum;
	            Sum = zero;
		  }     
 
          for (i=0;i<N;i++)
		  {
	          HS0[i] = H_in[i] - HT_S[i];
				 
		  }

 	Bi_CG(ZT_S,HS0,V12);
 //       resolve_FBM(ZT_S,HS0,V12,2*N);

  
         for (i=0;i<M;i++)
		 {
	         for (j=0;j<N;j++)
			 {
                Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		        dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        NR = (xtt[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(ztt[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum - dl*NR*J*K/4.*Hankel(1,1,K*Rm)*V12[j];
			 }
	            HV1_T[i] = Sum;
	            Sum = zero;

		 }

        for (i=0;i<M;i++)
		{
	         for (j=0;j<N;j++)
			 {
		         Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		         dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum+dl*J/4.*Hankel(0,1,K*Rm)*V12[j+N];
			 }
	            HV2_T[i] = Sum;
	            Sum = zero;
		}

        for (i=0;i<M;i++)
		{
	        HT0[i] = H0[i] - HV1_T[i] -HV2_T[i];
		}
 
		Bi_CG(ZT_T,HT0,V3);
   	   
	   CLinearEquations num0;

		    for (i=0;i<2*N;i++)
		    {
				V12_2[i] = V12[i]-V12_1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',ZT_S,V12_2,Vector_temp);
              Assign_Vecor(V12,V12_1);
              Delta = num0.VectorNorm2(Vector_temp,2*N)/num0.VectorNorm2(HS0,2*N);

              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.01)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;   
   
 }
*/
/////////////////////////////////////////////////////////////////////////////////
         
          

          double R_Dsc = 1000000*Lamda;
		 
		  Vec_dou Dsc;
		    Dsc.assign(N0);
		  Vec_dou Dsc_DB;
		    Dsc_DB.assign(N0);
		  Vec_dou x_Dsc;
		    x_Dsc.assign(N0);
		  Vec_dou z_Dsc;
		    z_Dsc.assign(N0);
		  Vec_Com Hs_MOM;  Hs_MOM.assign(N0);
		  Vec_dou Dsc_DB_MOM;  Dsc_DB_MOM.assign(N0);
		  Vec_dou Dsc_MOM;     Dsc_MOM.assign(N0);

		    


			for (i=0;i<N0;i++)
			{
			  Sita_Dsc[i] = (i-90)*PI/N0;
		
		      x_Dsc[i] = R_Dsc*sin(Sita_Dsc[i]);
		      z_Dsc[i] = R_Dsc*cos(Sita_Dsc[i]);
				
			}

/*			for (i=0;i<N0;i++)
			{
				 KSR = K*(x_Dsc[i]*sin(Sita_Dsc[i])+z_Dsc[i]*cos(Sita_Dsc[i]));
				for (j=0;j<N;j++)
				{
					NKS = K*((-1.*df_x[j]/sqrt(1+pow(df_x[j],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_x[j],2.0))));                   
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					Hs[i] = Hs[i] + J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*((-1.*J)*NKS*V12[j]-V12[j+N])*exp(-1.*J*KSR)*dl;
				}

				for (k=0;k<M;k++)
				{
				    NKS = K*((-1.*df_z[k]/sqrt(1+pow(df_z[k],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_z[k],2.0))));
                    Hs[i] = Hs[i] - J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*J*NKS*V3[k]*exp(-1.*J*KSR)*ds;
                     
				}
				    Hs[i] =   Hs[i]*exp(J*K*R_Dsc)/sqrt(R_Dsc);
			}*/

			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs[i] = Hs[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*Rm)*(Jf_2[j]+Jb_2[j]));
					Hs_MOM[i] = Hs_MOM[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(JJ_MOM[j])-J/4.*Hankel(0,1,K*Rm)*(JJ_MOM[j+N]));
				}
				for (k=0;k<M;k++)
				{
					Rm = sqrt(pow(x_Dsc[i]-xtt[k],2.0)+pow(z_Dsc[i]-ztt[k],2.0));
					NR = (-1.*df_z[k]*(x_Dsc[i]-xtt[k])+(z_Dsc[i]-ztt[k]))/(Rm*sqrt(1+pow(df_z[k],2.0)));
 //                   Hs[i] = Hs[i] + J*K/4.*NR*Hankel(1,1,K*Rm)*V3[k]*ds; 
                    Hs[i] = Hs[i] - ds*J/4.*Hankel(0,1,K*Rm)*V3[k]+NR*ds*J*K/4.*Hankel(1,1,K*Rm)*V3[k+M]; 
                    Hs_MOM[i] = Hs_MOM[i]-ds*J/4.*Hankel(0,1,K*Rm)*JJ_MOM[k+2*N]+NR*ds*J*K/4.*Hankel(1,1,K*Rm)*JJ_MOM[k+2*N+M];
				}		
			}

            double KG = 1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0)); 
			for (i=0;i<N0;i++)
			{
				
				Dsc[i] = R_Dsc*pow(abs(Hs[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB[i] = 10.*log10(Dsc[i]);

				Dsc_MOM[i] = R_Dsc*pow(abs(Hs_MOM[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB_MOM[i] = 10.*log10(Dsc_MOM[i]);
			}
//*2*PI*R_Dsc
			  ofstream outHs;
			         outHs.open("DEL_HS.txt");
					   for (i=0;i<N0;i++)
					  {
						  outHs<<i-90<<"   "<<real(Hs[i])<<"   "<<imag(Hs[i])<<"   "<<abs(Hs[i])<<endl;
					  }


			  ofstream outDsc;
			          outDsc.open("Del_S_Above_Target_DSC.txt");
					  for (i=0;i<N0;i++)
					  {
						  outDsc<<i-90<<"   "<<Dsc[i]<<"   "<<Dsc_DB[i]<<"   "<<Dsc_DB_MOM[i]<<endl;
					  }
			 

              ofstream outV1;
	                   outV1.open("Del_S_Above_Target_V1.txt");
			        for (i=0;i<N;i++)
					{
				       outV1<<x[i]/Lamda<<"   "<<real(V12[i])<<"   "<<imag(V12[i])<<"   "<<abs(V12[i])<<"  "<<abs(Jf_1[i]+Jb_1[i])<<"   "<<abs(JJ_MOM[i])<<endl;
					}


               ofstream outV2;
	                    outV2.open("Del_S_Above_Target_V2.txt");
			            for (i=0;i<N;i++)
						{
				            outV2<<x[i]/Lamda<<"   "<<real(V12[i+N])<<"   "<<imag(V12[i+N])<<"   "<<abs(V12[i+N])<<"  "<<abs(Jf_2[i]+Jb_2[i])<<"   "<<abs(JJ_MOM[i+N])<<endl;
						}

	  
                ofstream outV3;
	                     outV3.open("Del_S_Above_Target_V3.txt");
			             for (i=0;i<M;i++)
						 {
				            outV3<<i<<"   "<<real(V3[i])<<"   "<<imag(V3[i])<<"   "<<abs(V3[i])<<"   "<<abs(JJ_MOM[i+2*N])<<"   "<<abs(V3[i+M])<<"   "<<abs(JJ_MOM[i+2*N+M])<<endl;
						 }

               ofstream outDeDelta;
	                     outDeDelta.open("Del_S_Above_Target_Delta.txt");
			             for (i=0;i<15;i++)
						 {
				            outDeDelta<<i<<"   "<<Delta_Matrix[i]<<endl;
						 }

	
}

void Surf_Targ_Question::Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Es1()
{
/////////////////////////////////////////////////////////////////////////
//// 2008.12.3 修改，增加了FBM计算粗糙面方程，用Bi_CG计算目标方程   /////
////	            并与完全用Bi-CG得到的结果相比较                 /////
//////////////////////////////////////////////////////////////////////////

	//一维导体粗糙面
	double wave_len=3.0e8/3.0e9;      //wave length
	double h=0.5*wave_len;            //rms height 
	double lc=1.*wave_len;           //correlation length 
	double i_an=30.;                 //incident angle in degree 

  	double Epsle_0 = 8.85e-7;
//	Complex Epsle(3.7,0.13);
    Complex Epsle(2,0.05);
//	Complex rou_0= Epsle;
    Complex rou_0= 1;

//	cout<<"rou="<<rou_0<<endl;

	double T_an = asin(abs(sin(i_an*PI/180)/sqrt(Epsle)));


	
	if (T_an == atan(sqrt(abs(Epsle))))
	{
		cout<<"发生全反射。"<<endl;
	}

	double R=2.*wave_len;              //The R of the Target
		   H=-3.*wave_len;             //The z station of the Target
	double xp = 0*Lamda;              //The x station of the target
//	double g_p = 4.*R/cos(T_an);

//	cout<<"g_p = "<<g_p<<endl;


	
//	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
	double g_p=10.*wave_len;
//	double R=5.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double e=2.71828138;
    



	Complex Sum(0,0);
	Complex Sum_V3(0,0);

	double abx;
	double NR;
	double Rm;
	double dl;
	double Delta0;
	double Delta1;

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);

	

 	int N0 = 180;

	int i,j,k,m,l;

	Vec_dou Delda_Matrix;
	       Delda_Matrix.assign(15);
	Vec_Com Delta_V3;
	        Delta_V3.assign(M);

	Vec_Com V1;
	   V1.assign(N);
	Vec_Com V2;
	   V2.assign(N);
	Vec_Com V3;
	   V3.assign(M);
	Vec_Com HT_S;
		HT_S.assign(M);
	Vec_Com HS_T;
		HS_T.assign(N);

	Vec_Com HV1_T;
	    HV1_T.assign(M);
	Vec_Com HV2_T;
	    HV2_T.assign(M);
	Vec_Com HT0;
	    HT0.assign(M);
	Vec_Com H0;
	    H0.assign(M);

	Vec_Com H_in;
	   H_in.assign(N);
    Vec_Com H_in_Target;
	   H_in_Target.assign(M);
		
    Vec_Com HS0;
	   HS0.assign(2*N);
	Vec_Com V12;
	   V12.assign(2*N);
	Vec_Com V12_1;
	   V12_1.assign(2*N);

    Vec_Com Jf_1;
	       Jf_1.assign(N);
	Vec_Com Jf_2;
	       Jf_2.assign(N);
        
    Vec_Com Jb_1;
	       Jb_1.assign(N);
	Vec_Com Jb_2;
	       Jb_2.assign(N);

    Vec_Com Jf1;
	        Jf1.assign(2*N);	
    Vec_Com Jf2;
	        Jf2.assign(2*N);
	Vec_Com Jf3;
	        Jf3.assign(2*N);

	
	   
	Vec_Com V12_2;
	   V12_2.assign(2*N);
	Vec_Com Vector_temp;
	   Vector_temp.assign(2*N);

	Vec_Com  Hs;
	   Hs.assign(N0);
	Vec_Com  Hs_1;
	   Hs_1.assign(N0);

    Vec_Com  Hs_0;
	   Hs_0.assign(N);
	Vec_Com Hs_T;
	   Hs_T.assign(M);
	Vec_Com Hs_S;
	   Hs_S.assign(N);


	Vec_Com ZT_A;
	    ZT_A.assign(N*N);
	Vec_Com ZT_B;
	    ZT_B.assign(N*N);
	Vec_Com ZT_C;
	    ZT_C.assign(N*N);
	Vec_Com ZT_D;
	    ZT_D.assign(N*N);
	Vec_Com ZT_E;
	    ZT_E.assign(N*M);
	Vec_Com ZT_F;
	    ZT_F.assign(N*M);
	Vec_Com ZT_G;
	    ZT_G.assign(N*M);
	Vec_Com ZT_H;
	    ZT_H.assign(M*M);

	Vec_Com ZT_E_TE;
	    ZT_E_TE.assign(N*M);
	Vec_Com ZT_H_TE;
	    ZT_H_TE.assign(M*M);

	Vec_Com JJ_MOM;
	    JJ_MOM.assign(2*N+M);
	Vec_Com EE_MOM;
	    EE_MOM.assign(2*N+M);

	Vec_Com JJ;
	    JJ.assign(2*N+M);
	Vec_Com JJ_1;
	    JJ_1.assign(2*N+M);
	Vec_Com JJ_2;
	    JJ_2.assign(2*N+M);	

	Vec_Com J_Target_1;
	        J_Target_1.assign(M);
	Vec_Com J_Target_2;
	        J_Target_2.assign(M);
	   
	   
	   
   Vec_dou Sita_Dsc;
	    Sita_Dsc.assign(N0);

	Vec_Com ZZ_MOM;
	    ZZ_MOM.assign((2*N+M)*(2*N+M));

	Vec_Com ZT_S;
	   ZT_S.assign(2*N*2*N);
	Vec_Com ZT_T;
	   ZT_T.assign(M*M);
	Vec_Com Z_S_T;
	     Z_S_T.assign(M*N);

	Lamda=wave_len; 
		

	ds=2.*PI*R/M;
	Sita0=30.*PI/180.;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;




	for(i=0;i<=M-1;i++)
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

//	cout<<"进行目标的剖分......------>>"<<endl;
   
//    Cy_PoufenMatix_buried_target();

	Complex K1 = K*sqrt(Epsle);
	
    cout<<"生成粗糙面......------>>"<<endl;	
  	
    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

/*	for (i=0;i<N;i++)
	{
		z[i] = 0.;
		df_x[i] = 0.;
		ddf_x[i] = 0.;
	}*/

    dx_Suf=RL/N;
//+abs(H)*tan(T_an)
	double xsz,ww,kg,fac;
	i_an = i_an*PI/180;

    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(K*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		H_in[i]=exp(J*k_w*xsz*(1.0+ww))*exp(-fac);
	}

/*	for (i=0;i<M;i++)
	{
		xsz=xtt[i]*sin(i_an)-ztt[i]*cos(i_an);
		H_in_Target[i]=exp(J*K1*xsz);
	}*/



	cout<<" 开始填充阻抗矩阵......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

/*	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = (x[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);
			}
			else
			{
			   ZT_S[i*2*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
		}
	}*/


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_A[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);
			   
			}
			else
			{
			   ZT_A[i*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}

		}
	}

//////////////////////////////////////////////////////////////////////////


/*	for (i=0;i<N;i++)
	{
		for (j=N;j<2*N;j++)
		{
			if (j-i==N)
			{
               ZT_S[i*2*N+j] = dx_Suf*sqrt(1+pow(df_x[i],2.0))*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm            = sqrt(pow(x[i]-x[j-N],2.0)+pow(z[i]-z[j-N],2.0));
//			   dl           = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_S[i*2*N+j] = dx_Suf*J/4.*Hankel(0,1,K*Rm);

			}			
		}
		
	}*/

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (i==j)
			{
               ZT_B[i*N+j] = dx_Suf*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm           = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl           = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_B[i*N+j]  = dx_Suf*J/4.*Hankel(0,1,K*Rm);
//cout<<ZT_B[i*N+j] <<endl;
			}			

		}
	}

//////////////////////////////////////////////////////////////////////////


 /*  for (i=N;i<2*N;i++)
   {
	   for(j=0;j<N;j++)
	   {
		   if (i-j==N)
		   {
			   ZT_S[i*2*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i-N]-x[j],2.0)+pow(z[i-N]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = (x[i-N]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i-N]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = dx_Suf*Rm*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }
	   }
   }*/

   for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   if (i==j)
		   {
			   ZT_C[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_C[i*N+j]   = dx_Suf*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }
	   }
   }

//////////////////////////////////////////////////////////////////////////



  for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   if (i==j)
		   {
              ZT_D[i*N+j] = -1.*dx_Suf*rou_0*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_D[i*N+j]   = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);

		   }
	   }
   }





//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
						       
			   Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
			   NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(z[i]-ztt[j]))/(Rm*sqrt(1+pow(df_z[j],2.0)));
		       ZT_E[i*M+j] = ds*NR*J*K1/4.*Hankel(1,1,K1*Rm); 
			   ZT_E_TE[i*M+j] = J/4.*Hankel(0,1,K1*Rm);
			   
		}
	}

//////////////////////////////////////////////////////////////////////////
 
   for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {

  			   Rm            = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = ((xtt[i]-x[j])*(-1.*df_x[j])+(ztt[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_F[i*N+j]   = dx_Suf*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   
	   }
   }


//////////////////////////////////////////////////////////////////////////

  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
		       Rm            = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_G[i*N+j]   = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);	   

	   }
   }


//////////////////////////////////////////////////////////////////////////

 


    for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
		    if (i!=j)
			{
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*sqrt(1+pow(df_z[j],2.0)));
		            ZT_H[i*M+j] = ds*NR*J*K1/4.*Hankel(1,1,K1*Rm);
					ZT_H_TE[i*M+j] = J/4.*Hankel(0,1,K1*Rm);
			}
		    else
			{
			        ZT_H[i*M+j] = 0.5+ddf_z[i]*dx_t/(4.*PI*(1.+pow(df_z[i],2.0)));
					ZT_H_TE[i*M+j] = 0.5 - ddf_z[i]*dx_t/(4.*PI*(1.+pow(df_z[i],2.0)));
			}
		   		   
		}
	}


///////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZT_S[i*2*N+j] = ZT_A[i*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZT_S[i*2*N+j] = ZT_B[i*N+j-N];
		}
	}

	for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
          ZT_S[i*2*N+j] = ZT_C[(i-N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
		  ZT_S[i*2*N+j] = ZT_D[(i-N)*N+j-N];
		}		
	}

	for (i=0;i<N;i++)
	{
		HS0[i] = H_in[i];
	}
	
	
 //////////////////////////////////////////////////////////////////////////
	int N_M = 2*N+M;
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_A[i*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_B[i*N+j-N];
		}
	}

	for (i=N;i<2*N;i++)
	{
		for (j=0;j<N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_C[(i-N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_D[(i-N)*N+j-N];
		}
		for (j=2*N;j<N_M;j++)
		{
		ZZ_MOM[i*N_M+j] = ZT_E[(i-N)*M+j-2*N];
//				ZZ_MOM[i*N_M+j] = ZT_E_TE[(i-N)*M+j-2*N];

		}
	}

	for (i=2*N;i<N_M;i++)
	{
		for (j=0;j<N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_F[(i-2*N)*N+j];
		}
		for (j=N;j<2*N;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_G[(i-2*N)*N+j-N];
		}
		for (j=2*N;j<N_M;j++)
		{
			ZZ_MOM[i*N_M+j] = ZT_H[(i-2*N)*M+j-2*N];
//            ZZ_MOM[i*N_M+j] = ZT_H_TE[(i-2*N)*M+j-2*N];
		}
	}

	for (i=0;i<N;i++)
	{
		EE_MOM[i] = H_in[i];
	}
	 
	cout<<"阻抗矩阵填充完毕。"<<endl;


//  Bi_CGSTAB(ZZ_MOM,EE_MOM,JJ_MOM);

//	Bi_CG(ZZ_MOM,EE_MOM,JJ_MOM);



//	Bi_CG(ZT_S,HS0,V12);



for (l=0;l<15;l++)
{

cout<<"进行第"<<l+1<<"次目标与环境的耦合迭代......------>>"<<endl;


	

for (m=0;m<N;m++)
{


	for (i=0;i<N;i++)
	{
		for (j=0;j<i;j++)
		{
						
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_C[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_D[i*N+j];
			
		}
	
		Jf_1[i] = ((H_in[i]-Sum0)*ZT_D[i*N+i]-(-1.*Sum1-HS_T[i])*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_D[i*N+i]-ZT_B[i*N+i]*ZT_C[i*N+i]);
		Jf_2[i] = ((-1.*Sum1-HS_T[i])*ZT_A[i*N+i]-(H_in[i]-Sum0)*ZT_C[i*N+i])/(ZT_A[i*N+i]*ZT_D[i*N+i]-ZT_B[i*N+i]*ZT_C[i*N+i]);
        

		Sum0 = zero;
		Sum1 = zero;
	}

		 
	 for (i=N-1;i>=0;i--)
	 {
		for (j=i+1;j<N;j++)
		{
		    Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_C[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_D[i*N+j];	
		}


     	Jb_1[i] = ((-1.*Sum0)*ZT_D[i*N+i]-(-1.*Sum1)*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_D[i*N+i]-ZT_B[i*N+i]*ZT_C[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*ZT_A[i*N+i]-(-1.*Sum0)*ZT_C[i*N+i])/(ZT_A[i*N+i]*ZT_D[i*N+i]-ZT_B[i*N+i]*ZT_C[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;

	}

    		 CLinearEquations num;
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta0 = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta ="<<Delta0<<endl;
			  
			 if(Delta0<=0.01)
			 {
			    cout<<"The total number of the FBM step is :"<<m+1<<"."<<"and the final error is :"<<Delta0<<endl;
		        break;
			 }
				
	          else  
		           continue;



}



	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
			Sum0 = Sum0 + ZT_F[i*N+j]*(Jf_1[j]+Jb_1[j]) + ZT_G[i*N+j]*(Jf_2[j]+Jb_2[j]);
		}

		HT_S[i] = -Sum0;
		Sum0 = zero;				
	}

   Bi_CG(ZT_H_TE,HT_S,V3);

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			Sum0 = Sum0 + ZT_E_TE[i*M+j]*V3[j]; 
		}
		HS_T[i] = Sum0;
		Sum0 = zero;
	}

         CLinearEquations num;
		  

 		     for (i=0;i<M;i++)
			 {
				 J_Target_2[i] = J_Target_1[i] - V3[i];
			 }
			 for (i=0;i<M;i++)
			 {
				 for (j=0;j<M;j++)
				 {
					  Sum_V3 = Sum_V3+ZT_H_TE[i*M+j]*J_Target_2[j];
				 }
                  Delta_V3[i] = Sum_V3;
				  Sum_V3 = zero;
			 }
			 


             Assign_Vecor(V3,J_Target_1);
             Delta1 = num.VectorNorm2(Delta_V3,M)/num.VectorNorm2(HT_S,M);

			 Delda_Matrix[l] = Delta1;

             cout<<"The error is: Delta ="<<Delta1<<endl;
			  
			 if(Delta1<=0.01)
			 {
			    cout<<"The total number of the iteration step is :"<<l+1<<"."<<"and the final error is :"<<Delta1<<endl;
		        break;
			 }
				
	          else  
		           continue;
          


}
  
/*  

   for (k=0;k<20;k++)
   {
	   
       cout<<"开始粗糙面与目标的第"<<k+1<<"次耦合迭代......------>>"<<endl;



     	for (i=0;i<N;i++)
		{
	     	for (j=0;j<M;j++)
			{
			   Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
			   NR = (x[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(z[i]-ztt[j])/Rm/sqrt(1+pow(df_z[j],2.0));
			   Sum = Sum - ds*J*K1/4.*NR*Hankel(1,1,K1*Rm)*V3[j];
			}	
			   HT_S[i] = Sum;
		   	   Sum = zero;					
		}

	   for (i=0;i<N;i++)
	   {
		   HS0[i] = H_in[i];
       	   HS0[i+N] = HT_S[i];
	   }

       Bi_CG(ZT_S,HS0,V12);

       for (i=0;i<M;i++)
		 {
	         for (j=0;j<N;j++)
			 {
                Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		        dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        NR = (xtt[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(ztt[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum + dl*NR*J*K1/4.*Hankel(1,1,K1*Rm)*V12[j];
			 }
	            HV1_T[i] = Sum;
	            Sum = zero;

		 }

        for (i=0;i<M;i++)
		{
	         for (j=0;j<N;j++)
			 {
		         Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		         dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum + rou_0*dl*J/4.*Hankel(0,1,K1*Rm)*V12[j+N];
			 }
	            HV2_T[i] = Sum;
	            Sum = zero;
		}

       for (i=0;i<M;i++)
		{
	        HT0[i] =  HV2_T[i] - HV1_T[i];
		}
 
	   Bi_CG(ZT_T,HT0,V3);

//////////////////////////////////////////////////////////////////////////

   	   
	   CLinearEquations num0;

		    for (i=0;i<2*N;i++)
		    {
				V12_2[i] = V12[i]-V12_1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',ZT_S,V12_2,Vector_temp);
              Assign_Vecor(V12,V12_1);
              Delta = num0.VectorNorm2(Vector_temp,2*N)/num0.VectorNorm2(HS0,2*N);

              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.01)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;   
   
   }*/

/////////////////////////////////////////////////////////////////////

          double R_Dsc = 100000*Lamda;
		 
		  Vec_dou Dsc;
		    Dsc.assign(N0);
		  Vec_dou Dsc_1;
		    Dsc_1.assign(N0);

		  Vec_dou Dsc_DB_1;
		    Dsc_DB_1.assign(N0);
		  Vec_dou Dsc_DB;
		    Dsc_DB.assign(N0);
		  Vec_dou x_Dsc;
		    x_Dsc.assign(N0);
		  Vec_dou z_Dsc;
		    z_Dsc.assign(N0);
		    


			for (i=0;i<N0;i++)
			{
			  Sita_Dsc[i] = (i-90)*PI/N0;
		
		      x_Dsc[i] = R_Dsc*sin(Sita_Dsc[i]);
		      z_Dsc[i] = R_Dsc*cos(Sita_Dsc[i]);
				
			}

/*			for (i=0;i<N0;i++)
			{
    			 KSR = K*(x_Dsc[i]*sin(Sita_Dsc[i])+z_Dsc[i]*cos(Sita_Dsc[i]));
				for (j=0;j<N;j++)
				{
					NKS = K*((-1.*df_x[j]/sqrt(1+pow(df_x[j],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_x[j],2.0))));                   
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					Hs[i] = Hs[i] + J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*((-1.*J)*NKS*V12[j]-V12[j+N])*exp(-1.*J*KSR)*dl;
				}
				    Hs[i] =   Hs[i]*exp(J*K*R_Dsc)/sqrt(R_Dsc);

			}*/

			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs_1[i] = Hs_1[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*JJ_MOM[j]-J/4.*Hankel(0,1,K*Rm)*JJ_MOM[j+N]);
					Hs[i] = Hs[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*Rm)*(Jf_2[j]+Jb_2[j]));


				}
		
			}
			
/*   	      for (i=0;i<N;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x[i]-x[j],2.0)+pow(H-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x[i]-x[j])+(H-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs_0[i] = Hs_0[i]-dl*(J*K1/4.*NR*Hankel(1,1,K1*Rm)*V12[j]-J/4.*rou_0*Hankel(0,1,K1*Rm)*V12[j+N]);
				}
		
			}*/

           double KG = 1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));


			for (i=0;i<N0;i++)
			{
				Dsc[i] = R_Dsc*pow(abs(Hs[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB[i] = 10.*log10(Dsc[i]);
				Dsc_1[i] = R_Dsc*pow(abs(Hs_1[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB_1[i] = 10.*log10(Dsc_1[i]);
			}

			ofstream outHs;
			         outHs.open("Buried_Target_HS.txt");
					   for (i=0;i<N0;i++)
					  {
						  outHs<<i-90<<"   "<<real(Hs[i])<<"   "<<imag(Hs[i])<<"   "<<abs(Hs[i])<<"   "<<2*PI*pow(abs(Hs[i]),2.0)<<endl;
					  }


			ofstream outDsc;
			          outDsc.open("Buried_Target_DSC.txt");
					  for (i=0;i<N0;i++)
					  {
						  outDsc<<i-90<<"   "<<Dsc[i]<<"   "<<Dsc_DB[i]<<"    "<<Dsc_DB_1[i]<<endl;
					  }
             
		    ofstream outHs_0;
			         outHs_0.open("Buried_Target_HS_0.txt");
					   for (i=0;i<N;i++)
					  {
						  outHs_0<<x[i]/Lamda<<"   "<<real(Hs_0[i])<<"   "<<imag(Hs_0[i])<<"   "<<abs(Hs_0[i])<<endl;
					  }


			 

           ofstream outV1;
	               outV1.open("Buried_Target_V1.txt");
			       for (i=0;i<N;i++)
				   {
				        outV1<<x[i]/Lamda<<"   "<<real(V12[i])<<"   "<<imag(V12[i])<<"   "<<abs(V12[i])<<"  "<<abs(JJ_MOM[i])<<"    "<<abs(Jf_1[i]+Jb_1[i])<<endl;
				   }


            ofstream outV2;
	                 outV2.open("Buried_Target_V2.txt");
			         for (i=0;i<N;i++)
					 {
				          outV2<<x[i]/Lamda<<"   "<<real(V12[i+N])<<"   "<<imag(V12[i+N])<<"   "<<abs(V12[i+N])<<"   "<<abs(JJ_MOM[i+N])<<"   "<<abs(Jf_2[i]+Jb_2[i])<<endl;
					 }

	  
             ofstream outV3;
	                  outV3.open("Buried_Target_V3.txt");
			          for (i=0;i<M;i++)
					  {
				           outV3<<i<<"   "<<real(V3[i])<<"   "<<imag(V3[i])<<"   "<<abs(V3[i])<<"   "<<abs(JJ_MOM[i+2*N])<<endl;
					  }
			ofstream out_JJ_MOM;
			        out_JJ_MOM.open("Bried_Target_JJ_MOM.txt");
			         for (i=0;i<M;i++)
			         {
						 out_JJ_MOM<<i<<"   "<<real(JJ_MOM[i+2*N])<<"   "<<imag(JJ_MOM[i+2*N])<<"   "<<abs(JJ_MOM[i+2*N])<<endl;
			         }
			ofstream out_Delta;
			        out_Delta.open("Bried_Target_Delta.txt");
					for (i=0;i<15;i++)
					{
						out_Delta<<i<<"   "<<Delda_Matrix[i]<<endl;
					}


}

void Surf_Targ_Question::Gauss_Surface_Below_Muti_Target()
{

	double wave_len=3.0e8/3e9;      //wave length
	double h=0.5*wave_len;            //rms height 
	double lc=1.*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 

  	double Epsle_0 = 8.85e-7;
//	Complex Epsle(3.7,0.13);
    Complex Epsle(2,0.);

	double xp=20*wave_len;
//	Complex rou_0= Epsle;
    Complex rou_0= 1;
		    H=-10.*wave_len;
	double a=2*wave_len;
	double b=1*wave_len;

	double g_p = 10*wave_len; 

	
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
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double e=2.71828138;
    



	Complex Sum0(0,0);
	Complex Sum1(0,0);

	double abx;
	double NR;
	double Rm;
	double dl;
	double Delta0;
	double Delta1;
	double Sita_Target_1;
	double Sita_Target_2;

	

 	int N0 = 180;

	int i,j,k,m,l;

	Vec_Com U1;U1.assign(N);	   
	Vec_Com U2;U2.assign(N);	   
	Vec_Com U3;U3.assign(M);
	Vec_Com U4;U4.assign(M);
	   
	    
    Vec_Com E_in; E_in.assign(N);
	Vec_Com E_in_1; E_in_1.assign(N);
	Vec_Com E_in_2; E_in_2.assign(N);
	  
	Vec_dou xtt_1;xtt_1.assign(M);
	Vec_dou xtt_2;xtt_2.assign(M);
	Vec_dou ztt_1;ztt_1.assign(M);
	Vec_dou ztt_2;ztt_2.assign(M);
	Vec_dou ds_Muti_1; ds_Muti_1.assign(M);
	Vec_dou ds_Muti_2; ds_Muti_2.assign(M);

	Vec_Com ZT_A_1; ZT_A_1.assign(N*N);
	Vec_Com ZT_B_1; ZT_B_1.assign(N*N);
	Vec_Com ZT_A_2; ZT_A_2.assign(N*N);
	Vec_Com ZT_B_2; ZT_B_2.assign(N*N); 
	Vec_Com ZT_A_3; ZT_A_3.assign(N*M);
	Vec_Com ZT_B_3; ZT_B_3.assign(N*M);
	Vec_Com ZT_A_4; ZT_A_4.assign(N*M);
	Vec_Com ZT_B_4; ZT_B_4.assign(N*M);
	
	Vec_Com ZT_C_1; ZT_C_1.assign(N*M);
	Vec_Com ZT_D_1; ZT_D_1.assign(N*M);
	Vec_Com ZT_C_2; ZT_C_2.assign(M*M);
	Vec_Com ZT_D_2; ZT_D_2.assign(M*M);
	Vec_Com ZT_C_3; ZT_C_3.assign(M*M);
	Vec_Com ZT_D_3; ZT_D_3.assign(M*M);

	Vec_Com Jf_1; Jf_1.assign(N);
	Vec_Com Jb_1; Jb_1.assign(N);
	Vec_Com Jf_2; Jf_2.assign(N);
	Vec_Com Jb_2; Jb_2.assign(N);
	Vec_Com Jf1,Jf2,Jf3; Jf1.assign(2*N); Jf2.assign(2*N); Jf3.assign(2*N);

	Vec_Com ES_T; ES_T.assign(N);
	Vec_Com ET1_S; ET1_S.assign(M);
	Vec_Com ET2_S; ET2_S.assign(M);

	Vec_Com E_T_S_1; E_T_S_1.assign(M);
	Vec_Com E_T_S_2; E_T_S_2.assign(M);

	Vec_Com E_T1_T2; E_T1_T2.assign(M);
	Vec_Com E_T2_T1; E_T2_T1.assign(M);

	Vec_Com ET1; ET1.assign(M);
	Vec_Com ET2; ET2.assign(M);

	Vec_Com E_S_T; E_S_T.assign(N);

	Vec_Com U4_1; U4_1.assign(M);
	Vec_Com U4_2; U4_2.assign(M);
	Vec_Com U4_3; U4_3.assign(M);

	int N_M = 2*M+2*N;

	Vec_Com ZT_MOM; ZT_MOM.assign(N_M*N_M);
	Vec_Com JJ_MOM; JJ_MOM.assign(N_M);
	Vec_Com EE_MOM; EE_MOM.assign(N_M);
    
	Vec_dou Delta_Matrix; Delta_Matrix.assign(15);


	  
	cout<<"进行目标的剖分......------>>"<<endl;
   
  	Lamda=wave_len; 
		

	ds=2.*PI*R/M;
	Sita0=30.*PI/180.;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;

		   
	for(i=0;i<M;i++)
	{
		 Sita_Target_1=i*2.*PI/M;
         Sita_Target_2=i*2.*PI/M+PI/2.;
		 xtt_1[i]=a*cos(Sita_Target_1)-xp;
		 ztt_1[i]=b*sin(Sita_Target_1)+H;
		 xtt_2[i]=b*cos(Sita_Target_2)+xp;
		 ztt_2[i]=a*sin(Sita_Target_2)+H;		 
		 ds_Muti_1[i] = sqrt(pow(a*cos(Sita_Target_1),2.0)+pow(b*sin(Sita_Target_1),2.0))*2.*PI/M;
		 ds_Muti_2[i] = sqrt(pow(b*cos(Sita_Target_2),2.0)+pow(a*sin(Sita_Target_2),2.0))*2.*PI/M;		 
	}
	
	Complex K1 = K*sqrt(Epsle);
	
    cout<<"生成粗糙面......------>>"<<endl;	
  	
    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    dx_Suf=RL/N;

	double xsz,ww1,ww2,kg,fac1,fac2;
	i_an = i_an*PI/180;


    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac1=pow((x[i]+xp+z[i]*tan(i_an))/g_p,2.0);
		fac2=pow((x[i]-xp+z[i]*tan(i_an))/g_p,2.0);
		ww1=(2.0*fac1-1.0)/kg;
		ww2=(2.0*fac2-1.0)/kg;
		E_in_1[i]=exp(i_*k_w*xsz*(1.0+ww1))*exp(-fac1);
		E_in_2[i]=exp(i_*k_w*xsz*(1.0+ww2))*exp(-fac2);		
		E_in[i] = E_in_1[i] + E_in_2[i];
	}


	cout<<" 开始填充阻抗矩阵......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{						
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   NR          = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_A_1[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);			   
			}
			else
			{
			   ZT_A_1[i*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
		}
	}



	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{						
			if (i==j)
			{
               ZT_B_1[i*N+j] = dx_Suf*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));			    
			}
			else
			{
		       Rm           = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   ZT_B_1[i*N+j]  = dx_Suf*J/4.*Hankel(0,1,K*Rm);
			}			
		}
	}

//////////////////////////////////////////////////////////////////////////

  for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   if (i==j)
		   {
			   ZT_A_2[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));
		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   NR            = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_A_2[i*N+j] = dx_Suf*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);
		   }
	   }
   }



  for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   if (i==j)
		   {
               ZT_B_2[i*N+j] = -1.*dx_Suf*rou_0*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   ZT_B_2[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);
		   }
	   }
   }
  	
///////////////////////////////////////////////////////////////////////////////////
	
  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
  		   Rm            = sqrt(pow(xtt_1[i]-x[j],2.0)+pow(ztt_1[i]-z[j],2.0));
		   NR            = ((xtt_1[i]-x[j])*(-1.*df_x[j])+(ztt_1[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
           ZT_A_3[i*N+j] = dx_Suf*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);
	   }
   }



  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
		  Rm            = sqrt(pow(xtt_1[i]-x[j],2.0)+pow(ztt_1[i]-z[j],2.0));
		  ZT_B_3[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);
	   }
   } 
  
//////////////////////////////////////////////////////////////////////////

  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
  		   Rm            = sqrt(pow(xtt_2[i]-x[j],2.0)+pow(ztt_2[i]-z[j],2.0));
		   NR            = ((xtt_2[i]-x[j])*(-1.*df_x[j])+(ztt_2[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
           ZT_A_4[i*N+j] = dx_Suf*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);
	   }
   }



  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   Rm            = sqrt(pow(xtt_2[i]-x[j],2.0)+pow(ztt_2[i]-z[j],2.0));
		   ZT_B_4[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);
	   }
   } 

//////////////////////////////////////////////////////////////////////////

    for (i=0;i<N;i++)
	{
	   for(j=0;j<M;j++)
		{
		   Rm = sqrt(pow(x[i]-xtt_1[j],2.0)+pow(z[i]-ztt_1[j],2.0));
		   ZT_C_1[i*M+j] = J/4.*Hankel(0,1,K1*Rm);		   		   
		}
	}

    for (i=0;i<N;i++)
	{
	   for(j=0;j<M;j++)
		{
		   Rm = sqrt(pow(x[i]-xtt_2[j],2.0)+pow(z[i]-ztt_2[j],2.0));
		   ZT_D_1[i*M+j] = J/4.*Hankel(0,1,K1*Rm);
		}

	}
//////////////////////////////////////////////////////////////////////////


   for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
		    if (i!=j)
			{
			   Rm = sqrt(pow(xtt_1[i]-xtt_1[j],2.0)+pow(ztt_1[i]-ztt_1[j],2.0));
			   ZT_C_2[i*M+j] = J/4.*Hankel(0,1,K1*Rm);
			}
		    else
			{
			   ZT_C_2[i*M+j] = 0.5 ;
			}		   		   
		}
	}

    for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
			 Rm = sqrt(pow(xtt_1[i]-xtt_2[j],2.0)+pow(ztt_1[i]-ztt_2[j],2.0));
	     	 ZT_D_2[i*M+j] = J/4.*Hankel(0,1,K1*Rm);		   		   
		}
	}

//////////////////////////////////////////////////////////////////////////

   for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
			Rm = sqrt(pow(xtt_2[i]-xtt_1[j],2.0)+pow(ztt_2[i]-ztt_1[j],2.0));
			ZT_C_3[i*M+j] = J/4.*Hankel(0,1,K1*Rm);		   		   
		}
	}

    for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
		    if (i!=j)
			{
			   Rm = sqrt(pow(xtt_2[i]-xtt_2[j],2.0)+pow(ztt_2[i]-ztt_2[j],2.0));
			   ZT_D_3[i*M+j] = J/4.*Hankel(0,1,K1*Rm);
			}
		    else
			{
			   ZT_D_3[i*M+j] = 0.5 ;
			}		   		   
		}
	}
////////////////////////MOM/////////////////////////////////////////////

  for (i=0;i<N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_A_1[i*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_B_1[i*N+j-N];
	  }
  }
 

 for (i=N;i<2*N;i++)
 {
	 for (j=0;j<N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_A_2[(i-N)*N+j];
	 }
	 for (j=N;j<2*N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_B_2[(i-N)*N+j-N];
	 }

	 for (j=2*N;j<2*N+M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_C_1[(i-N)*M+j-2*N];
	 }
	 for (j=2*N+M;j<N_M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_D_1[(i-N)*M+j-2*N-M];
	 }

 }

 for (i=2*N;i<2*N+M;i++)
 {
	 for (j=0;j<N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_A_3[(i-2*N)*N+j];
	 }
	 for (j=N;j<2*N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_B_3[(i-2*N)*N+j-N];
	 }

	 for (j=2*N;j<2*N+M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_C_2[(i-2*N)*M+j-2*N];
	 }
	 for (j=2*N+M;j<N_M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_D_2[(i-2*N)*M+j-2*N-M];
	 }
 }

 for (i=2*N+M;i<N_M;i++)
 {
	 for (j=0;j<N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_A_4[(i-2*N-M)*N+j];
	 }
	 for (j=N;j<2*N;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_B_4[(i-2*N-M)*N+j-N];
	 }

	 for (j=2*N;j<2*N+M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_C_3[(i-2*N-M)*M+j-2*N];
	 }
	 for (j=2*N+M;j<N_M;j++)
	 {
		 ZT_MOM[i*N_M+j] = ZT_D_3[(i-2*N-M)*M+j-2*N-M];
	 }
 }

 for (i=0;i<N;i++)
 {
	 EE_MOM[i] = E_in[i];
 }

cout<<"The matrix have been filled."<<endl;
//    GAUSS(ZT_MOM,JJ_MOM,EE_MOM,N_M);
//  Bi_CG(ZT_MOM,EE_MOM,JJ_MOM);
	  

///////////////////////////////////////////////////////////////////////
	

for (l=0;l<15;l++)
{

cout<<"Begin the "<<l+1<<"th step of the iteration."<<endl;
    for (m=0;m<N;m++)
	{


	  for (i=0;i<N;i++)
	  {
		for (j=0;j<i;j++)
		{						
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A_1[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B_1[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_A_2[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B_2[i*N+j];			
		}
	
		Jf_1[i] = ((E_in[i]-Sum0)*ZT_B_2[i*N+i]-(-1.*Sum1-E_S_T[i])*ZT_B_1[i*N+i])/(ZT_A_1[i*N+i]*ZT_B_2[i*N+i]-ZT_B_1[i*N+i]*ZT_A_2[i*N+i]);
		Jf_2[i] = ((-1.*Sum1-E_S_T[i])*ZT_A_1[i*N+i]-(E_in[i]-Sum0)*ZT_A_2[i*N+i])/(ZT_A_1[i*N+i]*ZT_B_2[i*N+i]-ZT_B_1[i*N+i]*ZT_A_2[i*N+i]);
        

		Sum0 = zero;
		Sum1 = zero;
	}

		 
	 for (i=N-1;i>=0;i--)
	 {
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A_1[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B_1[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_A_2[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B_2[i*N+j];
		}


		Jb_1[i] = ((-Sum0)*ZT_B_2[i*N+i]-(-1.*Sum1-ES_T[i])*ZT_B_1[i*N+i])/(ZT_A_1[i*N+i]*ZT_B_2[i*N+i]-ZT_B_1[i*N+i]*ZT_A_2[i*N+i]);
		Jb_2[i] = ((-1.*Sum1-ES_T[i])*ZT_A_1[i*N+i]-(-Sum0)*ZT_A_2[i*N+i])/(ZT_A_1[i*N+i]*ZT_B_2[i*N+i]-ZT_B_1[i*N+i]*ZT_A_2[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;

	}

    		 CLinearEquations num;
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta0 = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta ="<<Delta0<<endl;
			  
			 if(Delta0<=0.01)
			 {
			    cout<<"The total number of the FBM step is :"<<m+1<<"."<<"and the final error is :"<<Delta0<<endl;
		        break;
			 }
				
	          else  
		           continue;
	}    	
	
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
	      Sum0 = Sum0 + ZT_A_3[i*N+j]*(Jf_1[j]+Jb_1[j]) + ZT_B_3[i*N+j]*(Jf_2[j]+Jb_2[j]); 
		}
		E_T_S_1[i] = -Sum0;
		Sum0 = zero;
	}
	


	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
		   Sum0 = Sum0 + ZT_D_2[i*M+j]*U4[j];
		}
		E_T1_T2[i] = -Sum0;
		Sum0 = zero;
	}

	for (i=0;i<M;i++)
	{
		ET1[i] = E_T_S_1[i] + E_T1_T2[i];
	}

	    
	GAUSS(ZT_C_2,U3,ET1,M);	
	
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
	      Sum0 = Sum0 + ZT_A_4[i*N+j]*(Jf_1[j]+Jb_1[j]) + ZT_B_4[i*N+j]*(Jf_2[j]+Jb_2[j]); 
		}
		E_T_S_2[i] = -Sum0;
		Sum0 = zero;	
	}

	for (i=0;i<M;i++)
	{
		for (j=0;j<M;j++)
		{
		   Sum0 = Sum0 + ZT_C_3[i*M+j]*U3[j];
		}
		E_T2_T1[i] = -Sum0;
		Sum0 = zero;
	}	
	
	for (i=0;i<M;i++)
	{
		ET2[i] = E_T_S_2[i] + E_T2_T1[i];
	}
	
	
	GAUSS(ZT_D_3,U4,ET2,M);	

	for (i=0;i<N;i++)
	{
		for (j=0;j<M;j++)
		{
			Sum0 = Sum0 + ZT_C_1[i*M+j]*U3[j] +  ZT_D_1[i*M+j]*U4[j]; 
		}
		E_S_T[i] = Sum0;
		Sum0 = zero;
	}

	for (i=0;i<M;i++)
	{
		U4_2[i] = U4[i] - U4_1[i];
	}

              
	Assign_Vecor(U4,U4_1);

	  
	CLinearEquations num0;
 
    Delta1 = num0.VectorNorm2(U4_2,M)/num0.VectorNorm2(U4,M);

    cout<<"目标与粗糙面耦合迭代误差Delta1 ="<<Delta1<<endl;

	Delta_Matrix[l] = Delta1;

	 if (Delta1<0.01)
	{
	   cout<<"目标与粗糙面耦合迭代总次数为:"<<l+1<<"."<<endl;
	   break;
	}
				
	   else  
		  continue;   	


}


          double R_Dsc = 100000*wave_len;
		  Vec_dou Dsc; Dsc.assign(N0);		    
		  Vec_dou Dsc_DB; Dsc_DB.assign(N0);		    
		  Vec_dou x_Dsc; x_Dsc.assign(N0);		   
		  Vec_dou z_Dsc; z_Dsc.assign(N0);

		  Vec_dou Sita_Dsc; Sita_Dsc.assign(N0);

		  Vec_Com Es; Es.assign(N0);
		  Vec_Com Es1; Es1.assign(N0);		    
		  Vec_dou Dsc1; Dsc1.assign(N0);		    
		  Vec_dou Dsc_DB1; Dsc_DB1.assign(N0);		    		    
		  double	KG = 1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));

			for (i=0;i<N0;i++)
			{
			  Sita_Dsc[i] = (i-90)*PI/N0;
		
		      x_Dsc[i] = R_Dsc*sin(Sita_Dsc[i]);
		      z_Dsc[i] = R_Dsc*cos(Sita_Dsc[i]);
				
			}

			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Es[i] = Es[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*Rm)*(Jf_2[j]+Jb_2[j]));
				}
		
			}

			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Es1[i] = Es1[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(JJ_MOM[j])-J/4.*Hankel(0,1,K*Rm)*(JJ_MOM[j+N]));
				}
		
			}
			
		  	for (i=0;i<N0;i++)
			{
				Dsc[i] = R_Dsc*pow(abs(Es[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB[i] = 10.*log10(Dsc[i]);
			}

		  	for (i=0;i<N0;i++)
			{
				Dsc1[i] = R_Dsc*pow(abs(Es1[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB1[i] = 10.*log10(Dsc1[i]);
			}

			ofstream out_U1;
			out_U1.open("Buried_Muti_Target_U1.txt");
			for (i=0;i<N;i++)
			{
				out_U1<<x[i]/wave_len<<"   "<<abs(Jf_1[i]+Jb_1[i])<<"  "<<abs(JJ_MOM[i])<<endl; 
			}

			ofstream out_U2;
			out_U2.open("Buried_Muti_Target_U2.txt");
			for (i=0;i<N;i++)
			{
				out_U2<<x[i]/wave_len<<"   "<<abs(Jf_2[i]+Jb_2[i])<<"   "<<abs(JJ_MOM[i+N])<<endl; 
			}

			ofstream out_U3;
			out_U3.open("Buried_Muti_Target_U3.txt");
			for (i=0;i<M;i++)
			{
				out_U3<<i<<"   "<<abs(U3[i])<<"   "<<abs(JJ_MOM[i+2*N])<<endl; 
			}

			ofstream out_U4;
			out_U4.open("Buried_Muti_Target_U4.txt");
			for (i=0;i<M;i++)
			{
				out_U4<<i<<"   "<<abs(U4[i])<<"   "<<abs(JJ_MOM[i+2*N+M])<<endl; 
			}

			ofstream out_Es;
			out_Es.open("Buried_Muti_Target_Es.txt");
			for (i=0;i<N0;i++)
			{
				out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<Dsc[i]<<"   "<<Dsc_DB[i]<<endl; 
			}	
	
			ofstream out_Es1;
			out_Es1.open("Buried_Muti_Target_Es1.txt");
			for (i=0;i<N0;i++)
			{
				out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<Dsc1[i]<<"   "<<Dsc_DB1[i]<<endl; 
			}
			
			ofstream out_Delta;
			out_Delta.open("Buried_Muti_Target_Delta.txt");
			for (i=0;i<15;i++)
			{
				out_Delta<<i<<"   "<<Delta_Matrix[i]<<endl; 
			}
}

void Surf_Targ_Question::Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Esd()
{

	
	//一维导体粗糙面
	double wave_len=3.0e8/1.5e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.5*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 

  	double Epsle_0 = 8.85e-7;
//	Complex Epsle(3.7,0.13);
    Complex Epsle(5,0.);
	Complex rou_0= Epsle*Epsle_0/Epsle_0;
//	cout<<"rou="<<rou_0<<endl;

	double T_an = asin(abs(sin(i_an)/sqrt(Epsle)));


	
	if (T_an == atan(sqrt(abs(Epsle))))
	{
		cout<<"发生全反射。"<<endl;
	}

	double R=2.*wave_len;

	double g_p = 4.*R/cos(T_an);

//	cout<<"g_p = "<<g_p<<endl;


	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=5.*wave_len;
//	double R=5.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double e=2.71828138;
    



	Complex Sum;

	double abx;
	double NR;
	double Rm;
	double dl;
	double Delta;

	

 	int N0 = 180;

	int i,j,k;

	Vec_Com V1;
	   V1.assign(N);
	Vec_Com V2;
	   V2.assign(N);
	Vec_Com V3;
	   V3.assign(M);
	Vec_Com HT_S;
		HT_S.assign(N);

	Vec_Com HV1_T;
	    HV1_T.assign(M);
	Vec_Com HV2_T;
	    HV2_T.assign(M);
	Vec_Com HT0;
	    HT0.assign(M);
	Vec_Com H0;
	    H0.assign(M);

	Vec_Com H_in;
	   H_in.assign(N);
    Vec_Com HS0;
	   HS0.assign(2*N);
   Vec_Com HSd;
	   HSd.assign(2*N);

    Vec_Com HSd_S;
	   HSd_S.assign(N);
    Vec_Com HSd_T;
	   HSd_T.assign(M);

	Vec_Com V12;
	   V12.assign(2*N);
	Vec_Com V12_1;
	   V12_1.assign(2*N);
	Vec_Com V12_2;
	   V12_2.assign(2*N);
 	Vec_Com V12_3;
	   V12_3.assign(2*N);

	Vec_Com Vector_temp;
	   Vector_temp.assign(2*N);

	Vec_Com  Hs;
	   Hs.assign(N0);
    Vec_Com  Hs_0;
	   Hs_0.assign(N);
	Vec_Com  Hs_1;
	   Hs_1.assign(N);
	Vec_Com Hs_T;
	   Hs_T.assign(M);
	Vec_Com Hs_S;
	   Hs_S.assign(N);
	Vec_dou Sita_Dsc;
	    Sita_Dsc.assign(N0);

	Vec_Com ZT_S;
	   ZT_S.assign(2*N*2*N);
	Vec_Com ZT_T;
	   ZT_T.assign(M*M);

	cout<<"进行目标的剖分......------>>"<<endl;
   
    Cy_PoufenMatix_buried_target();

	Complex K1 = K*sqrt(Epsle);
	
    cout<<"生成粗糙面......------>>"<<endl;	
  	
    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    dx_Suf=RL/N;

	double xsz,ww,kg,fac;
	i_an = i_an*PI/180;


    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(Sita0)-z[i]*cos(Sita0);
		kg=pow(k_w*g_p*cos(Sita0),2.0);
		fac=pow((x[i]+abs(H)*tan(T_an)+z[i]*tan(Sita0))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		H_in[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}


	cout<<" 开始填充阻抗矩阵......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = (x[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = -1.*dl*(J*K/4.)*NR*Hankel(1,1,K*Rm);
			}
			else
			{
			   ZT_S[i*2*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
		}
	}

//////////////////////////////////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=N;j<2*N;j++)
		{
			if (j-i==N)
			{
               ZT_S[i*2*N+j] = dx_Suf*sqrt(1+pow(df_x[i],2.0))*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm            = sqrt(pow(x[i]-x[j-N],2.0)+pow(z[i]-z[j-N],2.0));
			   dl           = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_S[i*2*N+j] = dl*J/4.*Hankel(0,1,K*Rm);

			}			
		}
		
	}

//////////////////////////////////////////////////////////////////////////


   for (i=N;i<2*N;i++)
   {
	   for(j=0;j<N;j++)
	   {
		   if (i-j==N)
		   {
			   ZT_S[i*2*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i-N]-x[j],2.0)+pow(z[i-N]-z[j],2.0));
			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = (x[i-N]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i-N]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = dl*Rm*NR*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }
	   }
   }

//////////////////////////////////////////////////////////////////////////



  for (i=N;i<2*N;i++)
   {
	   for (j=N;j<2*N;j++)
	   {
		   if (i==j)
		   {
              ZT_S[i*2*N+j] = -1.*dx_Suf*rou_0*sqrt(1+pow(df_x[j-N],2.0))*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j-N],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i-N]-x[j-N],2.0)+pow(z[i-N]-z[j-N],2.0));
			   dl            = dx_Suf*sqrt(1+pow(df_x[j-N],2.0));
			   ZT_S[i*2*N+j] = -1.*rou_0*dl*J/4.*Hankel(0,1,K1*Rm);

	   
		   }
	   }
   }

//////////////////////////////////////////////////////////////////////////


    for (i=0;i<M;i++)
	{
	   for(j=0;j<M;j++)
		{
		    if (i!=j)
			{
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = (xtt[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(ztt[i]-ztt[j])/Rm/(1+pow(df_z[j],2.0));
		            ZT_T[i*M+j] = ds*NR*J*K1/4.*Hankel(1,1,K1*Rm); 
			}
		    else
			{
			        ZT_T[i*M+j] = 0.5+ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			}
		   		   
		}
	}



    cout<<"阻抗矩阵填充完毕。"<<endl;

    for (i=0;i<N;i++)
	{
       HS0[i] = H_in[i];
	}

    Bi_CG(ZT_S,HS0,V12_3);

	for (i=0;i<N;i++)
	{
	    for (j=0;j<N;j++)
		{
			if (i!=j)
			{
			Rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		    NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
			HSd_S[i] = HSd_S[i] - dl*(J*K1/4.*NR*Hankel(1,1,K1*Rm)*V12_3[j]-J/4.*rou_0*Hankel(0,1,K1*Rm)*V12_3[j+N]);
			}
			else
			{
			 HSd_S[i] = HSd_S[i] - dl*((0.5+ddf_x[j]*dx_Suf/(4.0*PI*(1+pow(df_x[j],2.0))) )*V12_3[j]-rou_0*dl*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j],2.0)/(2.*e)))*V12_3[j+N]);
			}
		}
		
	}


	   for (i=0;i<M;i++)
		{
			for (j=0;j<N;j++)
			{
			    Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
			    dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
				NR = (-1.*df_x[j]*(xtt[i]-x[j])+(ztt[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
				HSd_T[i] = HSd_T[i] - dl*(J*K1/4.*NR*Hankel(1,1,K1*Rm)*V12_3[j]-J/4.*rou_0*Hankel(0,1,K1*Rm)*V12_3[j+N]);
			}

			  	
		}


//////////////////////////////////////////////////////////////////////////



	
//////////////////////////////////////////////////////////////////////////
	
   for (k=0;k<50;k++)
   {

       cout<<"开始粗糙面与目标的第"<<k+1<<"次耦合迭代......------>>"<<endl;



     	for (i=0;i<N;i++)
		{
	     	for (j=0;j<M;j++)
			{
			   Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
			   NR = (x[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(z[i]-ztt[j])/Rm/sqrt(1+pow(df_z[j],2.0));
			   Sum = Sum - ds*J*K1/4.*NR*Hankel(1,1,K1*Rm)*V3[j];
			}	
			   HT_S[i] = Sum;			  
		   	   Sum = zero;					
		}



	   for (i=0;i<N;i++)
	   {
       	   HSd[i+N] = HT_S[i] + HSd_S[i];
	   }

       Bi_CG(ZT_S,HSd,V12);


      
	   
       for (i=0;i<M;i++)
		 {
	         for (j=0;j<N;j++)
			 {
                Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		        dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        NR = (xtt[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(ztt[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum + dl*NR*J*K1/4.*Hankel(1,1,K1*Rm)*V12[j];
			 }
	            HV1_T[i] = Sum;
	            Sum = zero;
		 }

        for (i=0;i<M;i++)
		{
	         for (j=0;j<N;j++)
			 {
		         Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		         dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		         Sum = Sum - rou_0*dl*J/4.*Hankel(0,1,K1*Rm)*V12[j+N];
			 }
	            HV2_T[i] = Sum;
	            Sum = zero;
		}

		for (i=0;i<M;i++)
		{
	        HT0[i] = HSd_T[i] + HV2_T[i] - HV1_T[i];
		}
 
	   Bi_CG(ZT_T,HT0,V3);

//////////////////////////////////////////////////////////////////////////

   	   
	   CLinearEquations num0;

		    for (i=0;i<2*N;i++)
		    {
				V12_2[i] = V12[i]-V12_1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',ZT_S,V12_2,Vector_temp);
              Assign_Vecor(V12,V12_1);
              Delta = num0.VectorNorm2(Vector_temp,2*N)/num0.VectorNorm2(HSd,2*N);

              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.01)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;   
   
   }

/////////////////////////////////////////////////////////////////////

          double R_Dsc = 1000*Lamda;
		  Vec_dou Dsc;
		    Dsc.assign(N0);
		  Vec_dou Dsc_DB;
		    Dsc_DB.assign(N0);
		  Vec_dou x_Dsc;
		    x_Dsc.assign(N0);
		  Vec_dou z_Dsc;
		    z_Dsc.assign(N0);
		    


			for (i=0;i<N0;i++)
			{
			  Sita_Dsc[i] = (i-90)*PI/N0;
		
		      x_Dsc[i] = R_Dsc*sin(Sita_Dsc[i]);
		      z_Dsc[i] = R_Dsc*cos(Sita_Dsc[i]);
				
			}

/*			for (i=0;i<N0;i++)
			{
    			 KSR = K*(x_Dsc[i]*sin(Sita_Dsc[i])+z_Dsc[i]*cos(Sita_Dsc[i]));
				for (j=0;j<N;j++)
				{
					NKS = K*((-1.*df_x[j]/sqrt(1+pow(df_x[j],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_x[j],2.0))));                   
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					Hs[i] = Hs[i] + J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*((-1.*J)*NKS*V12[j]-V12[j+N])*exp(-1.*J*KSR)*dl;
				}
				    Hs[i] =   Hs[i]*exp(J*K*R_Dsc)/sqrt(R_Dsc);

			}*/

/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs[i] = Hs[i]+dl*(J*K/4.*NR*Hankel(1,1,K*Rm)*V12[j]-J/4.*Hankel(0,1,K*Rm)*V12[j+N]);
				}
		
			}

 		 for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs_1[i] = Hs_1[i]+dl*(J*K/4.*NR*Hankel(1,1,K*Rm)*(V12[j]+V12_3[j])-J/4.*Hankel(0,1,K*Rm)*(V12[j+N]+V12_3[j+N]));
				}

				abx   =  1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
				Hs_1[i] = pow(abs(Hs_1[i]),2.0)*2*PI*R_Dsc;
				Hs_1[i] = 10.*log10(Hs_1[i]);

		
			}
			
   	      for (i=0;i<N;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x[i]-x[j],2.0)+pow(H-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x[i]-x[j])+(H-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs_0[i] = Hs_0[i]-dl*(J*K1/4.*NR*Hankel(1,1,K1*Rm)*V12_3[j]-J/4.*rou_0*Hankel(0,1,K1*Rm)*V12_3[j+N]);				
				}
		
			}*/




/*			for (i=0;i<N0;i++)
			{
				abx   =  1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
				Dsc[i] = pow(abs(Hs[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*abx);
				Dsc_DB[i] = 10.*log10(Dsc[i]);
			}*/

		  	for (i=0;i<N0;i++)
			{
				abx   =  1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
				Dsc[i] = pow(abs(Hs[i]),2.0)*2*PI*R_Dsc;
				Dsc_DB[i] = 10.*log10(Dsc[i]);
			}



			ofstream outHs;
			         outHs.open("Buried_Target_Esd_HS.txt");
					   for (i=0;i<N0;i++)
					  {
						  outHs<<i-90<<"   "<<real(Hs[i])<<"   "<<imag(Hs[i])<<"   "<<abs(Hs[i])<<endl;
					  }


			ofstream outDsc;
			          outDsc.open("Buried_Target_Esd_DSC.txt");
					  for (i=0;i<N0;i++)
					  {
						  outDsc<<i-90<<"   "<<Dsc[i]<<"   "<<Dsc_DB[i]<<endl;
					  }
             
		    ofstream outHs_0;
			         outHs_0.open("Buried_Target_Esd_HS_0.txt");
					   for (i=0;i<N;i++)
					  {
						  outHs_0<<x[i]/Lamda<<"   "<<real(Hs_0[i])<<"   "<<imag(Hs_0[i])<<"   "<<abs(Hs_0[i])<<endl;
					  }
	       ofstream outHs_1;
			         outHs_1.open("Buried_Target_Esd_HS_1.txt");
					   for (i=0;i<N0;i++)
					  {
						  outHs_1<<i-90<<"   "<<Hs_1[i]<<endl;
					  }




			 

           ofstream outV1;
	               outV1.open("Buried_Target_Esd_V1.txt");
			       for (i=0;i<N;i++)
				   {
				        outV1<<x[i]/Lamda<<"   "<<real(V12[i])<<"   "<<imag(V12[i])<<"   "<<abs(V12[i])<<endl;
				   }


            ofstream outV2;
	                 outV2.open("Buried_Target_Esd_V2.txt");
			         for (i=0;i<N;i++)
					 {
				          outV2<<x[i]/Lamda<<"   "<<real(V12[i+N])<<"   "<<imag(V12[i+N])<<"   "<<abs(V12[i+N])<<endl;
					 }

	  
             ofstream outV3;
	                  outV3.open("Buried_Target_Esd_V3.txt");
			          for (i=0;i<M;i++)
					  {
				           outV3<<i<<"   "<<real(V3[i])<<"   "<<imag(V3[i])<<"   "<<abs(V3[i])<<endl;
					  }

}

void Surf_Targ_Question::Gauss_Del_Surface_Con_Cylinder__Above_Target()
{

	
	//一维导体粗糙面
	double wave_len=3.0e8/3.0e8;      //wave length
	double h=0.5*wave_len;            //rms height 
	double lc=0.6*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double g_p=(6.0/pow(cos(i_an*PI/180.0),1.5))*wave_len;
//	double g_p=5.*wave_len;
	double R=1.*wave_len;
		   H=3.*wave_len;

		   
	double g_p = 4.*R/cos(i_an*PI/180.);
//		   g_p = 20.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	//rl = 25.6 * wave_len;
	//g_p = rl / 4.0;
	int N=(rl/wave_len)*10;
	//N = 256;
	cout<<"未知量 N= "<<N<<" ;"<<endl;

	int n_r=1;                 //total number of surface realizations
	int n_sc=180;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长
	//rd_bw=300;
	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;

	double e=2.71828138;
	double GAMA = 0.577216;
    
	double Epsle_0 = 8.85e-7;
	Complex Epsle(25,1.0);
	Complex Epsle_Target(3,0.03);
//	Complex rou_0= Epsle;
//	Complex rou_Target = Epsle_Target;
	Complex rou_0= 1;
    Complex rou_Target = 1;

//	cout<<"rou="<<rou_0<<endl;
	

	Complex Sum(0.,0);
	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum_V3(0.,0.);
	Complex Sum_V4(0.,0.);

	double abx;
	double NR;
	double Rm;
	double dl;
	double Delta;
	double Delta0;
	double Delta1;

	

 	int N0 = 181;
	int i,j,k,m,l;
	Vec_Com V1; V1.assign(N);	  
	Vec_Com V2; V2.assign(N);	  
	Vec_Com V3; V3.assign(2*M);	  
	Vec_Com HT_S;	HT_S.assign(M);	
	Vec_Com V4; V4.assign(M);
	Vec_Com HV1_T; HV1_T.assign(M);	   
	Vec_Com HV2_T;  HV2_T.assign(M);	  
	Vec_Com HT0;  HT0.assign(2*M);	  
	Vec_Com H0;  H0.assign(M);	  
	Vec_Com H_in; H_in.assign(N);	  
    Vec_Com HS0; HS0.assign(2*N);	  
	Vec_Com V12; V12.assign(2*N);	  
	Vec_Com V12_1;V12_1.assign(2*N);	   
	Vec_Com V12_2;V12_2.assign(2*N);	   
	Vec_Com Vector_temp;Vector_temp.assign(2*N);	   
	Vec_Com  Hs; Hs.assign(N0);	  
	Vec_Com Hs_T; Hs_T.assign(N);	  
	Vec_Com Hs_S; Hs_S.assign(N);
	Vec_dou Sita_Dsc; Sita_Dsc.assign(N0);	   
	Vec_Com ZT_S; ZT_S.assign(2*N*2*N);
	Vec_Com ZT_T; ZT_T.assign(2*M*2*M);
	   

   Vec_Com ZT_MOM;  ZT_MOM.assign((2*N+2*M)*(2*N+2*M));
   Vec_Com JJ_MOM;  JJ_MOM.assign(2*N+2*M);
   Vec_Com EE_MOM;  EE_MOM.assign(2*N+2*M);

   Vec_Com ZT_MOM_CON;  ZT_MOM_CON.assign((2*N+M)*(2*N+M));
   Vec_Com JJ_MOM_CON;  JJ_MOM_CON.assign(2*N+M);
   Vec_Com EE_MOM_CON;  EE_MOM_CON.assign(2*N+M);


   Vec_Com ZT_A; ZT_A.assign(N*N);           
   Vec_Com ZT_B; ZT_B.assign(N*N);          
   Vec_Com ZT_C; ZT_C.assign(N*M);           
   Vec_Com ZT_D; ZT_D.assign(N*M);           
   Vec_Com ZT_E; ZT_E.assign(N*N);          
   Vec_Com ZT_F; ZT_F.assign(N*N);           
   Vec_Com ZT_G; ZT_G.assign(N*M);           
   Vec_Com ZT_H; ZT_H.assign(N*M);           
   Vec_Com ZT_I; ZT_I.assign(M*M);         
   Vec_Com ZT_J; ZT_J.assign(M*M);           		   
   Vec_Com ZT_K; ZT_K.assign(M*M);
   Vec_Com ZT_L; ZT_L.assign(M*M);

    Vec_Com Jf_1; Jf_1.assign(N);
	Vec_Com Jf_2; Jf_2.assign(N);
	Vec_Com Jb_1; Jb_1.assign(N);
    Vec_Com Jb_2; Jb_2.assign(N);

	Vec_Com Jf1;  Jf1.assign(2*N);
	Vec_Com Jf2;  Jf2.assign(2*N);
	Vec_Com Jf3;  Jf3.assign(2*N);

	Vec_Com J_Target_1; J_Target_1.assign(2*M);
	Vec_Com J_Target_2; J_Target_2.assign(2*M);
	Vec_Com Delta_V3;   Delta_V3.assign(2*M);

	Vec_dou Delta_Matrix; Delta_Matrix.assign(15);


	cout<<"进行目标的剖分......------>>"<<endl;

	Lamda=wave_len; 
 
		

	ds=2.*PI*R/M;
	Gama=1.781072418;
	K=2.*PI/Lamda;
	dx_t = 4*R/M;
    Yita=120.*PI;

	for( i=0;i<M;i++)
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
   
//   Cy_PoufenMatix_above_target();

	Complex K1 = K*sqrt(Epsle);
	Complex K_Target = K*sqrt(Epsle_Target);
	
    cout<<"生成粗糙面......------>>"<<endl;	
  	
    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    dx_Suf=RL/N;
	
/*	for (i=0;i<N;i++)
	{
		z[i] = 0.;
		df_x[i] = 0.;
		ddf_x[i] = 0.;
	}*/



	double xsz,ww,kg,fac;
//	double d=4.*Lamda;

	i_an = i_an*PI/180;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+H*tan(i_an)+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		H_in[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}


	
	for (i=0;i<M;i++)
	{
		xsz=xtt[i]*sin(i_an)-ztt[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((xtt[i]+H*tan(i_an)+ztt[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
 		H0[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
 //       H0[i] = exp(i_*k_w*xsz);
		 
	}






	cout<<" 开始填充阻抗矩阵......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

/*	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
		
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = (x[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(z[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
               ZT_S[i*2*N+j] = -1.*dl*(J*K/4.)*NR*Hankel(1,1,K*Rm);
//			   cout<<ZT_S[i*2*N+j]<<endl;
			}
			else
			{
			   ZT_S[i*2*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
		}
	}*/


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (i!=j)
			{			
			   Rm          = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_A[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);
//			   cout<<ZT_S[i*2*N+j]<<endl;
			}
			else
			{
			   ZT_A[i*N+j] = 0.5-ddf_x[i]*dx_Suf/(4.*PI*(1+pow(df_x[i],2.0)));	
			}
			
		}
	}

//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
						
			if (j==i)
			{
               ZT_B[i*N+j] = dx_Suf*(J/4.)*Hankel(0,1,K*dx_Suf*sqrt(1+pow(df_x[i],2.0))/(2*e));
			    
			}
			else
			{
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_B[i*N+j]   = dx_Suf*J/4.*Hankel(0,1,K*Rm);

			}			

		}
	}
//////////////////////////////////////////////////////////////////////////

for (i=0;i<N;i++)
{
	for (j=0;j<M;j++)
	{
					       
		Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
//		NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(z[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
//		ZT_I[i*M+j] = -1.*ds*NR*J*K/4.*Hankel(1,1,K*Rm);
		ZT_D[i*M+j] = ds*J/4.*Hankel(0,1,K*Rm);
			
	}
}

for (i=0;i<N;i++)
{
	for (j=0;j<M;j++)
	{
		Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
		NR = ((x[i]-xtt[j])*(-1.*df_z[j])+(z[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		ZT_C[i*M+j] = -1.*NR*J*K/4.*ds*Hankel(1,1,K*Rm);
	}
}

//////////////////////////////////////////////////////////////////////////


   for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   		   
		   if (i==j)
		   {
			   ZT_E[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4*PI*(1+pow(df_x[j],2.0)));

		   }
		   else
		   {
  			   Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR            = ((x[i]-x[j])*(-1.*df_x[j])+(z[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_E[i*N+j]   = dx_Suf*(J*K1/4.)*Hankel(1,1,K1*Rm);

		   }

	   }
   }
   

//////////////////////////////////////////////////////////////////////////


   for (i=0;i<N;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   		   
		   if (i==j)
		   {
              ZT_F[i*N+j] = -1.*dx_Suf*rou_0*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1+pow(df_x[j],2.0))/(2*e));
		   }
		   else
		   {
		       Rm            = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_F[i*N+j]   = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*Rm);

	   
		   }

	   }
   }
//////////////////////////////////////////////////////////////////////////
for (i=0;i<M;i++)
{
	for (j=0;j<N;j++)
	{
		 			  
		       Rm          = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl          = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   NR          = ((xtt[i]-x[j])*(-1.*df_x[j])+(ztt[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
               ZT_G[i*N+j] = -1.*dx_Suf*(J*K/4.)*NR*Hankel(1,1,K*Rm);			   
		
	}
}

for (i=0;i<M;i++)
{
	for (j=0;j<N;j++)
	{			       
		       Rm            = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
//			   dl            = dx_Suf*sqrt(1+pow(df_x[j],2.0));
			   ZT_H[i*N+j]   = dx_Suf*J/4.*Hankel(0,1,K*Rm);	 			   				       		
	}
}



//////////////////////////////////////////////////////////////////////////


	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_I[i*M+j] = -1.*NR*ds*J*K/4.*Hankel(1,1,K*Rm);
					ZT_J[i*M+j] = J/4.*ds*Hankel(0,1,K*Rm);
				  }
		          else
				  {

					ZT_I[i*M+j] = 0.5-ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_J[i*M+j] = J/4.*ds*Hankel(0,1,K*dx_t*sqrt(1.+df_z[i]*df_z[i])/(2.*e));
						
//						(1.+J*2./PI*log(exp(GAMA)*ds*K/(4.*e)));
//					cout<<ZT_J[i*M+j]<<endl;
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   		          
			   if (i!=j)
				  {
			        Rm = sqrt(pow(xtt[i]-xtt[j],2.0)+pow(ztt[i]-ztt[j],2.0));
			        NR = ((xtt[i]-xtt[j])*(-1.*df_z[j])+(ztt[i]-ztt[j]))/(Rm*(1+pow(df_z[j],2.0)));
		            ZT_K[i*M+j] = ds*NR*J*K_Target/4.*Hankel(1,1,K_Target*Rm);
					ZT_L[i*M+j] = -rou_Target*ds*J/4.*Hankel(0,1,K_Target*Rm);
				  }
		          else
				  {

					ZT_K[i*M+j] = 0.5+ddf_z[i]*dx_t/(4*PI*(1+pow(df_z[i],2.0)));
			        ZT_L[i*M+j] = -rou_Target*ds*J/4.*(1.+J*2./PI*log(exp(GAMA)*ds*K_Target/(4.*e)));
                    
				  }
		   }
	   }

	   for (i=0;i<M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_I[i*M+j];
		   }
		   for (j=M;j<2*M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_J[i*M+j-M];
		   }
	   }

	   for (i=M;i<2*M;i++)
	   {
		   for (j=0;j<M;j++)
		   {
               ZT_T[i*2*M+j] = ZT_K[(i-M)*M+j];
		   }
		   for (j=M;j<2*M;j++)
		   {
			   ZT_T[i*2*M+j] = ZT_L[(i-M)*M+j-M];
		   }
	   }


//	 GAUSS(ZT_J,V3,H0,M);  

	   
/*	   for (i=0;i<M;i++)
	   {
		   HT0[i] = H0[i];
	   }*/

//	   Bi_CG(ZT_T,H0,V3);
//	   Bi_CG(ZT_L,HT0,V4);



//////////////////////////////////////////////////////////////////////////
  int N_M = 2*N+2*M;
	   
  for (i=0;i<N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_A[i*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_B[i*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_C[i*M+j-2*N];
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_D[i*M+j-2*N-M];
	  }

  }

  for (i=N;i<2*N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_E[(i-N)*N+j]; 
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_F[(i-N)*N+j-N];
	  }

  }

  for (i=2*N;i<2*N+M;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_G[(i-2*N)*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_H[(i-2*N)*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
           ZT_MOM[i*N_M+j] = ZT_I[(i-2*N)*M+j-2*N];
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_J[(i-2*N)*M+j-2*N-M];
	  }
  }
  for (i=2*N+M;i<N_M;i++)
  {
	  for (j=2*N;j<2*N+M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_K[(i-2*N-M)*M+j-2*N];		  
	  }
	  for (j=2*N+M;j<N_M;j++)
	  {
		  ZT_MOM[i*N_M+j] = ZT_L[(i-2*N-M)*M+j-2*N-M];
	  }
  }

  for (i=0;i<N;i++)
  {
	  EE_MOM[i] = H_in[i];
  }
  for (i=0;i<M;i++)
  {
	  EE_MOM[i+2*N] = H0[i];
  }

 // GAUSS(ZT_MOM,JJ_MOM,EE_MOM,N_M);
 
 //////////////////////////////////////////////////////////////// 
   int N_M_CON = 2*N+M;
	   
  for (i=0;i<N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_A[i*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_B[i*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_D[i*M+j-2*N];
	  }

  }

  for (i=N;i<2*N;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_E[(i-N)*N+j]; 
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_F[(i-N)*N+j-N];
	  }

  }

  for (i=2*N;i<2*N+M;i++)
  {
	  for (j=0;j<N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_G[(i-2*N)*N+j];
	  }
	  for (j=N;j<2*N;j++)
	  {
		  ZT_MOM_CON[i*N_M_CON+j] = ZT_H[(i-2*N)*N+j-N];
	  }
	  for (j=2*N;j<2*N+M;j++)
	  {
           ZT_MOM_CON[i*N_M_CON+j] = ZT_J[(i-2*N)*M+j-2*N];
	  }
  }
 

  for (i=0;i<N;i++)
  {
	  EE_MOM_CON[i] = H_in[i];
  }

  for (i=0;i<M;i++)
  {
	  EE_MOM_CON[i+2*N] = H0[i];
  }
   cout<<"阻抗矩阵填充完毕。"<<endl; 

//     Bi_CG(ZT_MOM,EE_MOM,JJ_MOM);
//   Bi_CGSTAB(ZT_MOM,EE_MOM,JJ_MOM);
	 GAUSS(ZT_MOM_CON,JJ_MOM_CON,EE_MOM_CON,N_M_CON);

//////////////////////////////////////////////////////////////////////////

 /*        


for (l=0;l<15;l++)
{

cout<<"Begin the "<<l+1<<"step of the iteration......------>>"<<endl;
  
   for (m=0;m<N;m++)
   {


	for (i=0;i<N;i++)
	{
		for (j=0;j<i;j++)
		{
						
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_E[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_F[i*N+j];
			
		}
	
		Jf_1[i] = ((H_in[i]-Sum0-Hs_T[i])*ZT_F[i*N+i]-(-1.*Sum1)*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Jf_2[i] = ((-1.*Sum1)*ZT_A[i*N+i]-(H_in[i]-Sum0-Hs_T[i])*ZT_E[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
        

		Sum0 = zero;
		Sum1 = zero;
	}

		 
	 for (i=N-1;i>=0;i--)
	 {
		for (j=i+1;j<N;j++)
		{
		    Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_A[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_B[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_E[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_F[i*N+j];	
		}


     	Jb_1[i] = ((-1.*Sum0)*ZT_F[i*N+i]-(-1.*Sum1)*ZT_B[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*ZT_A[i*N+i]-(-1.*Sum0)*ZT_E[i*N+i])/(ZT_A[i*N+i]*ZT_F[i*N+i]-ZT_B[i*N+i]*ZT_E[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;

	}

    		 CLinearEquations num;
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta0 = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta0 ="<<Delta0<<endl;
			  
			 if(Delta0<=0.01)
			 {
			    cout<<"The total number of the FBM step is :"<<m+1<<"."<<"and the final error is :"<<Delta0<<endl;
		        break;
			 }
				
	          else  
		           continue;



   }
   

  
  for (i=0;i<M;i++)
   {
	   for (j=0;j<N;j++)
	   {
		   Sum0 = Sum0 + ZT_G[i*N+j]*(Jf_1[j]+Jb_1[j])+ZT_H[i*N+j]*(Jf_2[j]+Jb_2[j]);
	   }
	   HT_S[i] = Sum0;
	   Sum0 = zero;
   }
   for (i=0;i<M;i++)
   {
	   HT0[i] = H0[i]-HT_S[i];
   }
 //    Bi_CG(ZT_T,HT0,V3);

    GAUSS(ZT_T,V3,HT0,2*M);
//   Bi_CG(ZT_J,HT0,V3);
	 
//    GAUSS(ZT_J,V3,H0,M);
    
   for (i=0;i<N;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   Sum0 = Sum0 + ZT_C[i*M+j]*V3[j]+ZT_D[i*M+j]*V3[j+M];
	   }
	   Hs_T[i] = Sum0;
	   Sum0 = zero;
   }
   
		  

 		     for (i=0;i<2*M;i++)
			 {
				 J_Target_2[i] = J_Target_1[i] - V3[i];
			 }
			 for (i=0;i<M;i++)
			 {
				 for (j=0;j<M;j++)
				 {
					  Sum_V3 = Sum_V3+ZT_I[i*M+j]*J_Target_2[j]+ZT_J[i*M+j]*J_Target_2[j+M];
					  Sum_V4 = Sum_V4+ZT_K[i*M+j]*J_Target_2[j]+ZT_L[i*M+j]*J_Target_2[j+M];
				 }
                  Delta_V3[i] = Sum_V3;
				  Delta_V3[i+M] = Sum_V4;
				  Sum_V3 = zero;
                  Sum_V4 = zero;
			 }
			  CLinearEquations num;

             Assign_Vecor(V3,J_Target_1);
             Delta1 = num.VectorNorm2(Delta_V3,2*M)/num.VectorNorm2(HT_S,M);

			 Delta_Matrix[l] = Delta1;

             cout<<"The error of the iteration is: Delta ="<<Delta1<<endl;
			  
			 if(Delta1<=0.01)
			 {
			    cout<<"The total number of the iteration step is :"<<l+1<<"."<<"and the final error is :"<<Delta1<<endl;
		        break;
			 }
				
	          else  
		           continue;

}
/*
   for (k=0;k<1;k++)
   {

        cout<<"开始粗糙面与目标的第"<<k+1<<"次耦合迭代......------>>"<<endl;

        for (i=0;i<N;i++)
		  {
	         for (j=0;j<M;j++)
			 {
		        Rm = sqrt(pow(x[i]-xtt[j],2.0)+pow(z[i]-ztt[j],2.0));
		        NR = (x[i]-xtt[j])/Rm*(-1.*df_z[j])/sqrt(1+pow(df_z[j],2.0))+(z[i]-ztt[j])/Rm/sqrt(1+pow(df_z[j],2.0));
		        Sum = Sum - ds*(J*K/4.)*NR*V3[j]*Hankel(1,1,K*Rm);

			 }
   
                HT_S[i] = Sum;
	            Sum = zero;
		  }     
 
          for (i=0;i<N;i++)
		  {
	          HS0[i] = H_in[i] - HT_S[i];
				 
		  }

 	Bi_CG(ZT_S,HS0,V12);
 //       resolve_FBM(ZT_S,HS0,V12,2*N);

  
         for (i=0;i<M;i++)
		 {
	         for (j=0;j<N;j++)
			 {
                Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		        dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        NR = (xtt[i]-x[j])/Rm*(-1.*df_x[j])/sqrt(1+pow(df_x[j],2.0))+(ztt[i]-z[j])/Rm/sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum - dl*NR*J*K/4.*Hankel(1,1,K*Rm)*V12[j];
			 }
	            HV1_T[i] = Sum;
	            Sum = zero;

		 }

        for (i=0;i<M;i++)
		{
	         for (j=0;j<N;j++)
			 {
		         Rm = sqrt(pow(xtt[i]-x[j],2.0)+pow(ztt[i]-z[j],2.0));
		         dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
		        Sum = Sum+dl*J/4.*Hankel(0,1,K*Rm)*V12[j+N];
			 }
	            HV2_T[i] = Sum;
	            Sum = zero;
		}

        for (i=0;i<M;i++)
		{
	        HT0[i] = H0[i] - HV1_T[i] -HV2_T[i];
		}
 
		Bi_CG(ZT_T,HT0,V3);
   	   
	   CLinearEquations num0;

		    for (i=0;i<2*N;i++)
		    {
				V12_2[i] = V12[i]-V12_1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',ZT_S,V12_2,Vector_temp);
              Assign_Vecor(V12,V12_1);
              Delta = num0.VectorNorm2(Vector_temp,2*N)/num0.VectorNorm2(HS0,2*N);

              cout<<"目标与粗糙面耦合迭代误差Delta ="<<Delta<<endl;

	          if (Delta<0.01)
			  {
				       cout<<"目标与粗糙面耦合迭代总次数为:"<<k+1<<"."<<endl;
		               break;
			  }
				
	          else  
		              continue;   
   
 }
*/
/////////////////////////////////////////////////////////////////////////////////
         
          

          double R_Dsc = 1000000*Lamda;
		 
		  Vec_dou Dsc;
		    Dsc.assign(N0);
		  Vec_dou Dsc_DB;
		    Dsc_DB.assign(N0);
		  Vec_dou x_Dsc;
		    x_Dsc.assign(N0);
		  Vec_dou z_Dsc;
		    z_Dsc.assign(N0);
		  Vec_Com Hs_MOM;  Hs_MOM.assign(N0);
		  Vec_dou Dsc_DB_MOM;  Dsc_DB_MOM.assign(N0);
		  Vec_dou Dsc_MOM;     Dsc_MOM.assign(N0);

		    


/*			for (i=0;i<N0;i++)
			{
			  Sita_Dsc[i] = (i-90)*PI/N0;
		
		      x_Dsc[i] = R_Dsc*sin(Sita_Dsc[i]);
		      z_Dsc[i] = R_Dsc*cos(Sita_Dsc[i]);
				
			}

/*			for (i=0;i<N0;i++)
			{
				 KSR = K*(x_Dsc[i]*sin(Sita_Dsc[i])+z_Dsc[i]*cos(Sita_Dsc[i]));
				for (j=0;j<N;j++)
				{
					NKS = K*((-1.*df_x[j]/sqrt(1+pow(df_x[j],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_x[j],2.0))));                   
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					Hs[i] = Hs[i] + J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*((-1.*J)*NKS*V12[j]-V12[j+N])*exp(-1.*J*KSR)*dl;
				}

				for (k=0;k<M;k++)
				{
				    NKS = K*((-1.*df_z[k]/sqrt(1+pow(df_z[k],2.0))*sin(Sita_Dsc[i])+cos(Sita_Dsc[i])/sqrt(1+pow(df_z[k],2.0))));
                    Hs[i] = Hs[i] - J/4.*sqrt(2./(PI*K))*exp(-1.*J*PI/4.)*J*NKS*V3[k]*exp(-1.*J*KSR)*ds;
                     
				}
				    Hs[i] =   Hs[i]*exp(J*K*R_Dsc)/sqrt(R_Dsc);
			}*/

/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					Rm = sqrt(pow(x_Dsc[i]-x[j],2.0)+pow(z_Dsc[i]-z[j],2.0));
					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_Dsc[i]-x[j])+(z_Dsc[i]-z[j]))/(Rm*sqrt(1+pow(df_x[j],2.0)));
					Hs[i] = Hs[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*Rm)*(Jf_2[j]+Jb_2[j]));
					Hs_MOM[i] = Hs_MOM[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*Rm)*(JJ_MOM[j])-J/4.*Hankel(0,1,K*Rm)*(JJ_MOM[j+N]));
				}
				for (k=0;k<M;k++)
				{
					Rm = sqrt(pow(x_Dsc[i]-xtt[k],2.0)+pow(z_Dsc[i]-ztt[k],2.0));
					NR = (-1.*df_z[k]*(x_Dsc[i]-xtt[k])+(z_Dsc[i]-ztt[k]))/(Rm*sqrt(1+pow(df_z[k],2.0)));
 //                   Hs[i] = Hs[i] + J*K/4.*NR*Hankel(1,1,K*Rm)*V3[k]*ds; 
                    Hs[i] = Hs[i] - ds*J/4.*Hankel(0,1,K*Rm)*V3[k]+NR*ds*J*K/4.*Hankel(1,1,K*Rm)*V3[k+M]; 
                    Hs_MOM[i] = Hs_MOM[i]-ds*J/4.*Hankel(0,1,K*Rm)*JJ_MOM[k+2*N]+NR*ds*J*K/4.*Hankel(1,1,K*Rm)*JJ_MOM[k+2*N+M];
				}		
			}

            double KG = 1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0)); 
			for (i=0;i<N0;i++)
			{
				
				Dsc[i] = R_Dsc*pow(abs(Hs[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB[i] = 10.*log10(Dsc[i]);

				Dsc_MOM[i] = R_Dsc*pow(abs(Hs_MOM[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*KG);
				Dsc_DB_MOM[i] = 10.*log10(Dsc_MOM[i]);
			}
//*2*PI*R_Dsc*/
			  ofstream outHs;
			         outHs.open("CON_HS.txt");
					   for (i=0;i<N0;i++)
					  {
						  outHs<<i-90<<"   "<<real(Hs[i])<<"   "<<imag(Hs[i])<<"   "<<abs(Hs[i])<<endl;
					  }


			  ofstream outDsc;
			          outDsc.open("Del_S_Above_CON_Target_DSC.txt");
					  for (i=0;i<N0;i++)
					  {
						  outDsc<<i-90<<"   "<<Dsc[i]<<"   "<<Dsc_DB[i]<<"   "<<Dsc_DB_MOM[i]<<endl;
					  }
			 

              ofstream outV1;
	                   outV1.open("Del_S_Above_CON_Target_V1.txt");
			        for (i=0;i<N;i++)
					{
				       outV1<<x[i]/Lamda<<"   "<<real(V12[i])<<"   "<<imag(V12[i])<<"   "<<abs(V12[i])<<"  "<<abs(Jf_1[i]+Jb_1[i])<<"   "<<abs(JJ_MOM[i])<<"   "<<abs(JJ_MOM_CON[i])<<endl;
					}


               ofstream outV2;
	                    outV2.open("Del_S_Above_CON_Target_V2.txt");
			            for (i=0;i<N;i++)
						{
				            outV2<<x[i]/Lamda<<"   "<<real(V12[i+N])<<"   "<<imag(V12[i+N])<<"   "<<abs(V12[i+N])<<"  "<<abs(Jf_2[i]+Jb_2[i])<<"   "<<abs(JJ_MOM[i+N])<<"   "<<abs(JJ_MOM_CON[i+N])<<endl;
						}

	  
                ofstream outV3;
	                     outV3.open("Del_S_Above_CON_Target_V3.txt");
			             for (i=0;i<M;i++)
						 {
				            outV3<<i<<"   "<<real(V3[i])<<"   "<<imag(V3[i])<<"   "<<abs(V3[i])<<"   "<<abs(JJ_MOM[i+2*N])<<"   "<<abs(V3[i+M])<<"   "<<abs(JJ_MOM[i+2*N+M])<<"   "<<abs(JJ_MOM_CON[i+2*N])<<endl;
						 }

               ofstream outDeDelta;
	                     outDeDelta.open("Del_S_Above_CON_Target_Delta.txt");
			             for (i=0;i<15;i++)
						 {
				            outDeDelta<<i<<"   "<<Delta_Matrix[i]<<endl;
						 }

	
}