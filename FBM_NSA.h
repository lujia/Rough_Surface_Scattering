#include "surface.h "
#include "Target.h"
//#include "Surf_Targ.h"

class FBM_SAA:public post_cp, public CylinderQuestion
{
public:
	     FBM_SAA();
		 ~FBM_SAA();
		 void Gauss_PEC_Surface_FBM();
		 void Gauss_DIL_Surface_FBM();        // The rms of the surface can not smaller than the 0.08Wavelength.

		 void Gauss_PEC_Surface_FBM_SAA();
		 void Gauss_DIL_Surface_FBM_SAA();        // The rms of the surface can not smaller than the 0.08Wavelength.
    	 void Gauss_DIL_Surface_PBTC_FBM_SAA();
		 void Gauss_Layered_Rough_surface_FBM();
};

FBM_SAA::FBM_SAA()
{

}

FBM_SAA::~FBM_SAA()
{

}
void FBM_SAA::Gauss_PEC_Surface_FBM()
{

//////////////////2008.11.01///////////////////////

///////////////////姬伟杰/////////////////////////
     int NUMBER_REALIZATION =1;
	
	//一维导体粗糙面
	 double wave_len=3.0e8/(3.0e8);      //wave length
	 h=0.5*wave_len;                    //rms height 
	 lc=1.0*wave_len;                   //correlation length 
	 double i_an=30.0;                   //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

	double R=1.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);

     g_p = 6.4*wave_len;
	
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
	int n_sc=360;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	 us=0.437198;
	 kl=100;
	 ku=4000;
	double ee=2.71828138;

	double GAMA = 0.577216;

	double euler = 1.781072418;
	double rm;

	double Bita_s ;
	double L_s;
	double R_s;
	double Q;
	double Fai_S_max;
	double Delta_0;
	double Delta_Fai;

	int N_s;

	double max_z=0.;
	double min_z=0.;


	Complex Fai_p;



	int i,j,k,m;
	//double w_Sita;

	Complex Sum;

	

	Vec_Com Ei;
	    Ei.assign(N);
	Vec_Com Js0;
	    Js0.assign(N);
    Vec_Com Jf;
	    Jf.assign(N);

	 Vec_Com Jf2;
	   Jf2.assign(N);

	 Vec_Com Jf1;
	   Jf1.assign(N);

	double Delta;

	Vec_Com Jb;
	    Jb.assign(N);
	Vec_Com Ef;
	    Ef.assign(N);
	Vec_Com Eb;
	    Eb.assign(N);

	Vec_Com Ef_w;
	    Ef_w.assign(N);
	Vec_Com Eb_w;
	    Eb_w.assign(N);

	Vec_Com Ef_s;
	    Ef_s.assign(N);
	Vec_Com Eb_s;
	    Eb_s.assign(N);

	Vec_Com ZT;
        ZT.assign(N*N);
	Vec_dou zz;
	   zz.assign(N);

    Vec_Com JJ;JJ.assign(N);
	double R_D = 10000000.*wave_len;    


	     
	int N0=180;

	
   Vec_Com Es;Es.assign(N0);		    
   Vec_dou D_Es;D_Es.assign(N0);		    
   Vec_dou DB_Es;DB_Es.assign(N0);		    	
   Vec_Com Es1;Es1.assign(N0);			    
   Vec_dou Sita_D;Sita_D.assign(N0);		    
   Vec_dou x_D;x_D.assign(N0);		     
   Vec_dou z_D;z_D.assign(N0);
		     
   K = 2.0*PI/wave_len;
   Vec_Com Vector_temp;Vector_temp.assign(N);     
   Complex zero(0.,0.);
   Complex Sum0 = zero;
   Complex Sum1 = zero;
   Complex SUM=zero;

   mean = 0.;
   vari = 1.;
   seed=123456.;
   RL=rl;
   k_w=K;

   double U195=5;

	for (int nn=0;nn<NUMBER_REALIZATION;nn++)
	{

    cout<<"Begin the "<<nn+1<<" th realization."<<endl;


	gauss_surface();
//	ocean_surface_PM(U195);
//	ocean_surface();
//    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    double dx_Suf=RL/N;

 //   i_an = i_an*PI/180;


	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an*PI/180)-z[i]*cos(i_an*PI/180);
		kg=pow(k_w*g_p*cos(i_an*PI/180),2.0);
		fac=pow((x[i]+z[i]*tan(i_an*PI/180))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}



	//The number of the strong interaction point. 
		
   for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{			
			if (j==i)
			{
                 ZT[i*N+j] = dx_Suf*J/4.*(1.+J*2./PI*(log(K*dx_Suf*euler/(4.))-1));
			}
			else
			{
				rm=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
                ZT[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
			}
		 
		}

	}
//	resolve_FBM(ZT,Ei,JJ,N);


	for(k=0;k<20;k++)
	  {
	  
           cout<<"Runing the "<<k+1<<" th step of the FBM......------->>"<<endl;

///////////////计算对电磁波前向传播有贡献的部分感应电流///////////////////////

	        for (i=0;i<N;i++)
			{
		        for (j=0;j<i;j++)
				{
			       Sum0 = Sum0+ZT[i*N+j]*(Jf[j]+Jb[j]);
			    }
			
	            Jf[i] = (Ei[i]-Sum0)/ZT[i*N+i];       
		        Sum0 = zero;	       		
			}

///////////////计算对电磁波后向传播有贡献的部分感应电流(其计算顺序为:从后向前计算，这是为了完全用更新的电流来计算)///////////////////////

	        for (i=N-1;i>=0;i--)
			{
		       for (j=i+1;j<N;j++)
			   {
			     Sum1 = Sum1+ZT[i*N+j]*(Jb[j]+Jf[j]);                 
			   }
	            Jb[i] =  -1.*(Sum1)/ZT[i*N+i];
				Sum1 = zero;	       
			}
			
///////////////////////收敛条件的判断/////////////////////////////////////////////////

		  
		    for (i=0;i<N;i++)
		    {
				Jf2[i] = Jf[i]-Jf1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',ZT,Jf2,Vector_temp);
              Assign_Vecor(Jf,Jf1);
              double DELTA = sqrt(L2_Norm(Vector_temp))/sqrt(L2_Norm(Ei));
			

           cout<<" THE error of is "<<DELTA<<endl;
	       if(DELTA<0.001)
		   {
			   cout<<"The total steps of the FBM is"<<k+1<<". and the final error is "<<DELTA<<"."<<endl;
		    	  break; 
		   }
	        else  
		          continue;
			  					
	  }


	 for(i=0;i<N;i++)
	  {
		  JJ[i] = Jf[i] + Jb[i] ;
	  
	  }
 

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
			x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);

		 }


         double EE=0.0;
 		 double abx   =  8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an*PI/180)*(1.0-((1.0+2.0*pow(tan(i_an*PI/180),2))/(2.0*pow(k_w*g_p*cos(i_an*PI/180),2))));

	     for (i=0;i<N0;i++)
		 {
			 			
			 Complex integ=0.0;

			 for (j=0;j<N;j++)
			 {
				double ksc=-k_w*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i] ));
				integ=integ+dx_Suf*JJ[j]*exp(i_*ksc);
			 }

			Es[i] = integ;

			double temp = norm(integ)/abx;

			D_Es[i] = (temp+nn*D_Es[i])/(nn+1);

        	EE=EE+D_Es[i]*PI/180;

		
			
			temp = 10.*log10(D_Es[i]);
			DB_Es[i] = (temp + nn * DB_Es[i]) / (nn + 1);
		 }

           cout<<"反射系数为："<<EE<<endl;
		}
		


	 ofstream outF_SEs;
		     outF_SEs.open("PEC_FBM_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 outF_SEs<<i-90<<"   "<<real(Es[i])<<"   "<<imag(Es[i])<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<"   "<<endl;
			 }
	 ofstream outF_SEs1;
		     outF_SEs1.open("PEC_FBM_Es1.txt");
			 for (i=0;i<N0;i++)
			 {
				 outF_SEs1<<i-90<<"   "<<real(Es1[i])<<"   "<<imag(Es1[i])<<"   "<<abs(Es1[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<"   "<<endl; 
			 } 
     

	 ofstream outJf;
		          outJf.open("PEC_FBM_Ef_Jf.txt");
                  for (i=0;i<N;i++)
                  {
					  outJf<<x[i]/wave_len<<"   "<<real(Jf[i])<<"   "<<imag(Jf[i])<<"   "<<abs(Jf[i])<<endl;
					
                  }
	             
	 ofstream outJb;
		          outJb.open("PEC_FBM_Eb_Jb.txt");
                  for (i=0;i<N;i++)
                  {
					  outJb<<x[i]/wave_len<<"   "<<real(Jb[i])<<"   "<<imag(Jb[i])<<"   "<<abs(Jb[i])<<endl;
					
                  } 
				  
     ofstream outJ;
		          outJ.open("PEC_FBM_Eb_J.txt");
                  for (i=0;i<N;i++)
                  {
					  outJ<<x[i]/wave_len<<"   "<<real(Jb[i]+Jf[i])<<"   "<<imag(Jb[i]+Jf[i])<<"   "<<abs(Jb[i]+Jf[i])<<endl;
					
                  } 
	 ofstream outJJ;
		          outJJ.open("PEC_FBM_Eb_JJ.txt");
                  for (i=0;i<N;i++)
                  {
					  outJJ<<x[i]/wave_len<<"   "<<real(JJ[i])<<"   "<<imag(JJ[i])<<"   "<<abs(JJ[i])<<endl;
					
                  } 
	 ofstream outZZ;
		          outZZ.open("PEC_FBM_ZZ.txt");
                  for (i=0;i<N;i++)
                  {
					  outZZ<<x[i]<<"   "<<z[i]<<"  "<<df_x[i]<<endl;
					
                  } 
	   			 	
}

void FBM_SAA::Gauss_PEC_Surface_FBM_SAA()
{

//////////////////2008.11.01///////////////////////

///////////////////姬伟杰/////////////////////////
     int NUMBER_REALIZATION =1;
	
	//一维导体粗糙面
	 double wave_len=3.0e8/(3.0e8);      //wave length
	 h=0.5*wave_len;                    //rms height 
	 lc=1.0*wave_len;                   //correlation length 
	 double i_an=30.0;                   //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

	double R=1.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);

     g_p = 6.4*wave_len;
	
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
	int n_sc=360;              //number of scattered angles from -90 deg to 90 deg 
	int rd_bw=int((10.0*wave_len)/(rl/N));//rd_bw 取10个波长

	cout<<"强弱相关距离 rd_bw = "<<rd_bw<<" ;"<<endl;
	
	double eps=1.0e-2;//CGM的收敛准则
	double criterion=1.0e-2;//BMIA的收敛准则
    
	double us=0.437198;
	double kl=100;
	double ku=4000;
	double ee=2.71828138;

	double GAMA = 0.577216;

	double euler = 1.781072418;
	double rm;

	double Bita_s ;
	double L_s;
	double R_s;
	double Q;
	double Fai_S_max;
	double Delta_0;
	double Delta_Fai;

	int N_s;

	double max_z=0.;
	double min_z=0.;


	Complex Fai_p;



	int i,j,k,m;
	//double w_Sita;

	Complex Sum;

	

	Vec_Com Ei;
	    Ei.assign(N);
	Vec_Com Js0;
	    Js0.assign(N);
    Vec_Com Jf;
	    Jf.assign(N);

	 Vec_Com Jf2;
	   Jf2.assign(N);

	 Vec_Com Jf1;
	   Jf1.assign(N);

	double Delta;

	Vec_Com Jb;
	    Jb.assign(N);
	Vec_Com Ef;
	    Ef.assign(N);
	Vec_Com Eb;
	    Eb.assign(N);

	Vec_Com Ef_w;
	    Ef_w.assign(N);
	Vec_Com Eb_w;
	    Eb_w.assign(N);

	Vec_Com Ef_s;
	    Ef_s.assign(N);
	Vec_Com Eb_s;
	    Eb_s.assign(N);

	Vec_Com ZT;
        ZT.assign(N*N);
	Vec_dou zz;
	   zz.assign(N);

    Vec_Com JJ;JJ.assign(N);
	double R_D = 10000000.*wave_len;    


	     
	int N0=180;

	
   Vec_Com Es;Es.assign(N0);		    
   Vec_dou D_Es;D_Es.assign(N0);		    
   Vec_dou DB_Es;DB_Es.assign(N0);		    	
   Vec_Com Es1;Es1.assign(N0);			    
   Vec_dou Sita_D;Sita_D.assign(N0);		    
   Vec_dou x_D;x_D.assign(N0);		     
   Vec_dou z_D;z_D.assign(N0);
		     
   K = 2.*PI/wave_len;
   Vec_Com Vector_temp;Vector_temp.assign(N);     
   Complex zero(0.,0.);
   Complex Sum0 = zero;
   Complex Sum1 = zero;
   Complex SUM=zero;

   mean = 0.;
   vari = 1.;
   seed=123456.;
   RL=rl;
   k_w=K;

	for (int nn=0;nn<NUMBER_REALIZATION;nn++)
	{

    cout<<"Begin the "<<nn+1<<" th realization."<<endl;


	gauss_surface();
//    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    double dx_Suf=RL/N;

 //   i_an = i_an*PI/180;


	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an*PI/180)-z[i]*cos(i_an*PI/180);
		kg=pow(k_w*g_p*cos(i_an*PI/180),2.0);
		fac=pow((x[i]+z[i]*tan(i_an*PI/180))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(i_*k_w*xsz*(1.0+ww))*exp(-fac);
	}







	

	for (i=0;i<N;i++)
	{
//		zz[i] = 0.;
		if (z[i]>max_z)
		{
			max_z=z[i];
		}
	    if (z[i]<=min_z)
	    {
			min_z=z[i];
	    }
	}
//	cout<<"max_z="<<max_z<<endl;







	L_s =10*wave_len; //10*(max_z - min_z)/4. ;                                 //The strong interaction zone.
	Bita_s = sqrt(10./(K*L_s));                                //The geometric optics lit region.
	R_s = sqrt(pow(L_s,2.0)+pow(max_z-min_z,2.0));             
	Fai_S_max = atan((max_z-min_z)/L_s);                       //The outermost possible saddle point.
	
	if (sqrt(K*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}
	
	Delta_Fai = sqrt(5./(K*R_s))/22.;                          
	Q = int(Bita_s/Delta_Fai)+5;

	N_s = int(L_s/dx_Suf);

	cout<<"N_s="<<N_s<<endl;
	//The number of the strong interaction point. 
	

	for (i=0;i<=N_s;i++)
	{
		for (j=0;j<=i+N_s;j++)
		{			

			if (j==i)
			{
			     ZT[i*N+j] = (1.+J*2./PI*(log(K*dx_Suf/(4.*ee)*sqrt(1.+pow(df_x[i],2.0)))+GAMA));
			}

			else
			{
				rm=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));

                ZT[i*N+j] = Hankel(0,1,K*rm);
			}
		 
		}

	}

	for (i=N_s+1;i<N-N_s;i++)
	{
		for (j=i-N_s;j<=i+N_s;j++)
		{			

			if (j==i)
			{
			     ZT[i*N+j] = (1.+J*2./PI*(log(K*dx_Suf/(4.*ee)*sqrt(1.+pow(df_x[i],2.0)))+GAMA));
			}

			else
			{
				rm=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));

                ZT[i*N+j] = Hankel(0,1,K*rm);
			}
		 
		}

	}

	for (i=N-N_s;i<N;i++)
	{
		for (j=N-N_s;j<N;j++)
		{			

			if (j==i)
			{
			     ZT[i*N+j] = (1.+J*2./PI*(log(K*dx_Suf/(4.*ee)*sqrt(1.+pow(df_x[i],2.0)))+GAMA));
			}

			else
			{
				rm=sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));

                ZT[i*N+j] = Hankel(0,1,K*rm);
			}
		 
		}

	}
	
	 	
//	resolve_PEC_FBM_SAA(ZT,Ei,JJ,max_z,min_z,dx_Suf,K,z);

	Complex Sum1;
 

    Vec_Com FN_f;
           FN_f.assign(2*Q+1);
	Vec_Com FN_b;
		   FN_b.assign(2*Q+1);
    Vec_Com Fai;
	       Fai.assign(2*Q+1);
    Vec_dou w_Sita;
	       w_Sita.assign(2*Q+1);
  

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);                //The window function.
         if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;		 
	}



    for (m=0;m<N;m++)
	{

	    cout<<"Runing the "<<m+1<<"th Step of FBM/SAA ......------>>"<<endl;

//////////////////////////////The strong point for the foward.////////////////////////////////////////////////////////

           
          for (i=0;i<N_s+1;i++)
		  {
	         for (j=0;j<i;j++)
			 {
                 Sum =Sum + (Jf[j]+Jb[j])*ZT[i*N+j];
			 }
               Ef_s[i] = Sum;
			   Sum = zero;
               Jf[i] = (Ei[i]-Ef_s[i])/ZT[i*N+i];	
		  } 

////////////////////////////The weak point for the foward.//////////////////////////////////////////////
		 

///////////////////////////////begin with the N_s+1 point///////////////////////////////////////////////////////
		 
       for (i=N_s+1;i<N;i++)
	   {
		     for (j=i-N_s;j<i;j++)
			 {
                Sum = Sum + (Jf[j]+Jb[j])*ZT[i*N+j];       
			 }
		        Ef_s[i] = Sum;
			    Sum = zero;
			 for (k=-Q;k<=Q;k++)
			 {   
			    Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
					  
                FN_f[k+Q] = FN_f[k+Q]*exp(J*K*dx_Suf*cos(Fai_p)) + (Jf[i-N_s-1]+Jb[i-N_s-1])*exp(J*K*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i-N_s-1]*sin(Fai_p));
  				      				      
			    Sum = Sum + exp(-1.*J*Delta_0)/PI*w_Sita[Q+k]*FN_f[k+Q]*exp(J*K*z[i]*sin(Fai_p))*Delta_Fai; 									
				  
			 }            
	            Ef_w[i] = Sum;
				Sum = zero;
	     Jf[i] = (Ei[i]-Ef_s[i]-Ef_w[i])/ZT[i*N+i];					   
	   }
	
////////////////////////////////The Strong point for the behihd.////////////////////////////////////////////////////////////
								
	  for (i=N-1;i>=N-N_s-1;i--)
	   { 
	       for (j=i+1;j<N;j++)
		   {			          
            Sum = Sum + (Jf[j]+Jb[j])*ZT[i*N+j];  
		   }
	     Eb_s[i] =Sum;	
		 Sum = zero;
		 Jb[i] = -1.*Eb_s[i]/ZT[i*N+i];
	   } 

		  

////////////////////////////The weak for the behind./////////////////////////////////////////////////////////////////



/////////////////////////////Satrt from the N-N_s-2 point./////////////////////////////////////////////


         for (i=N-N_s-2;i>=0;i--)
		  {
		     for (j=i+1;j<i+N_s+1;j++)
			 {	
                 Sum = Sum + (Jf[j]+Jb[j])*ZT[i*N+j];  		   
			 }
                  Eb_s[i] = Sum;		  
		          Sum = zero;
			for (k=-Q;k<=Q;k++)
			{   
			     Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
                 FN_b[k+Q] = FN_b[k+Q]*exp(J*K*dx_Suf*cos(Fai_p))+(Jf[i+N_s+1]+Jb[i+N_s+1])*exp(J*K*(1.+N_s)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i+N_s+1.]*sin(Fai_p));					            				 
			     Sum = Sum + exp(-1.*J*Delta_0)/PI*w_Sita[Q+k]*FN_b[k+Q]*exp(J*K*z[i]*sin(Fai_p))*Delta_Fai; 												
			}
	            Eb_w[i] = Sum;
				Sum = zero;
			    Jb[i] = -1.*(Eb_w[i]+Eb_s[i])/ZT[i*N+i];
		 }

/////////////////////////////////////////////////////////////////////////////////
	  
	        CLinearEquations num;
		  

 		    for (i=0;i<N;i++)
		    {
			  Jf2[i] = Jf[i]-Jf1[i];
		    }

              Assign_Vecor(Jf,Jf1);
              Delta = num.VectorNorm2(Jf2,N)/num.VectorNorm2(Jf,N);

              cout<<"The Delta of FBM_SAA ="<<Delta<<endl;
			  
			  if(Delta<=0.01)
			  {
				cout<<"The total steps of FBM_SAAis "<<m+1<<"."<<endl;
		        break;
			  }
				
	          else  
		        continue;

  
	}

 
    
	    

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
			x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);

		 }



 		 double abx   =  8.0*PI*k_w*g_p*sqrt(PI/2.0)*cos(i_an*PI/180)*(1.0-((1.0+2.0*pow(tan(i_an*PI/180),2))/(2.0*pow(k_w*g_p*cos(i_an*PI/180),2))));

	     for (i=0;i<N0;i++)
		 {
			 			
			 Complex integ=0.0;

			 for (j=0;j<N;j++)
			 {
				double ksc=-k_w*(x[j]*sin(Sita_D[i] )+z[j]*cos(Sita_D[i] ));
				integ=integ+4./J*(Jf[j]+Jb[j])*exp(i_*ksc);
			 }

			Es[i] = integ;

			double temp = norm(integ)/abx;
			D_Es[i] = (temp+nn*D_Es[i])/(nn+1);
			
			temp = 10.*log10(D_Es[i]);
			DB_Es[i] = (temp + nn * DB_Es[i]) / (nn + 1);
		 }

		}
				

	 ofstream outF_SEs;
		     outF_SEs.open("PEC_FBM_SAA_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 outF_SEs<<i-90<<"   "<<real(Es[i])<<"   "<<imag(Es[i])<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<"   "<<endl;
			 }
	 ofstream outF_SEs1;
		     outF_SEs1.open("PEC_FBM_SAA_Es1.txt");
			 for (i=0;i<N0;i++)
			 {
				 outF_SEs1<<i-90<<"   "<<real(Es1[i])<<"   "<<imag(Es1[i])<<"   "<<abs(Es1[i])<<"   "<<endl; 
			 } 
     

	 ofstream outJf;
		          outJf.open("PEC_FBM_SAA_Ef_Jf.txt");
                  for (i=0;i<N;i++)
                  {
					  outJf<<x[i]/wave_len<<"   "<<real(Jf[i])<<"   "<<imag(Jf[i])<<"   "<<abs(Jf[i])<<endl;
					
                  }
	             
	 ofstream outJb;
		          outJb.open("PEC_FBM_SAA_Eb_Jb.txt");
                  for (i=0;i<N;i++)
                  {
					  outJb<<x[i]/wave_len<<"   "<<real(Jb[i])<<"   "<<imag(Jb[i])<<"   "<<abs(Jb[i])<<endl;
					
                  } 
				  
     ofstream outJ;
		          outJ.open("PEC_FBM_SAA_Eb_J.txt");
                  for (i=0;i<N;i++)
                  {
					  outJ<<x[i]/wave_len<<"   "<<real(Jb[i]+Jf[i])<<"   "<<imag(Jb[i]+Jf[i])<<"   "<<abs(Jb[i]+Jf[i])<<endl;
					
                  } 
	 ofstream outJJ;
		          outJJ.open("PEC_FBM_SAA_Eb_JJ.txt");
                  for (i=0;i<N;i++)
                  {
					  outJJ<<x[i]/wave_len<<"   "<<real(JJ[i])<<"   "<<imag(JJ[i])<<"   "<<abs(JJ[i])<<endl;
					
                  } 
	 ofstream outZZ;
		          outZZ.open("PEC_FBM_SAA_ZZ.txt");
                  for (i=0;i<N;i++)
                  {
					  outZZ<<x[i]<<"   "<<z[i]<<endl;
					
                  } 
	   			 	
}

void FBM_SAA::Gauss_DIL_Surface_FBM()
{

//////////////////2008.11.03////////////////////////

///////////////////姬伟杰///////////////////////////

	
	int NUMBER_REALIZATION=1;
	int Idex_TE_TM=1;
	//一维导体粗糙面
	 wave_len=3.0e8/3e8;      //wave length
	 h=0.1*wave_len;            //rms height 
	 lc=1.0*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

	double R=1.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);
   
	double g_p = 6.4*wave_len;
	
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

	double GAMA = 0.577216;
	double rm;

	double Bita_s ;
	double L_s;
	double R_s;
	double Q;
	double Fai_S_max;
	double Delta_0;
	double Delta_Fai;
	Complex K1;
     
	Complex zero(0.,0.);
	Complex Eps_d(2.0,0.05);
	Complex rou_0 = 1;
	if (Idex_TE_TM==1)
	{
	  rou_0 = 1;
	}
	else
	{
      rou_0 = Eps_d;
	}


	int N_s;

	double max_z=0.;
	double min_z=0.;
	double dl;
	double NR;
	double Delta;


	Complex Fai_p;



	int i,j,k,m;
	//double w_Sita;

	Complex Sum0;
	Complex Sum1;
	Complex Sum2;
	Complex Sum3;

	

	Vec_Com Ei;
	    Ei.assign(N);
	Vec_Com Js0;
	    Js0.assign(N);
 
	Vec_Com Jf3;
	   Jf3.assign(2*N);  

	Vec_Com Jf2;
	   Jf2.assign(2*N);

	Vec_Com Jf1;
	   Jf1.assign(2*N);

	
    Vec_Com Jf_1;
	    Jf_1.assign(N);
	Vec_Com Jf_2;
	    Jf_2.assign(N);
	
	Vec_Com Jb_1;
	    Jb_1.assign(N);
		
	
	Vec_Com Jb_2;
	    Jb_2.assign(N);
	
		
	Vec_Com Ef;
	    Ef.assign(N);
	Vec_Com Eb;
	    Eb.assign(N);

	Vec_Com Ef_w;
	    Ef_w.assign(N);
	Vec_Com Eb_w;
	    Eb_w.assign(N);

	Vec_Com Ef_s;
	    Ef_s.assign(N);
	Vec_Com Eb_s;
	    Eb_s.assign(N);

	Vec_Com ZT_a;
        ZT_a.assign(N*N);
	Vec_Com ZT_b;
        ZT_b.assign(N*N);
	Vec_Com ZT_c;
        ZT_c.assign(N*N);
	Vec_Com ZT_d;
        ZT_d.assign(N*N);

     Vec_Com W_s;       W_s.assign(N);
		Vec_Com U_in;   U_in.assign(N);
	    Vec_Com U_s;    U_s.assign(N);
		 double R_D = 10000.*wave_len;
		 int N0=180;

		 Vec_Com Es;
		    Es.assign(N0);
		 Vec_dou D_Es;
		    D_Es.assign(N0);
    	 Vec_dou DB_Es;
		    DB_Es.assign(N0);
		
		 Vec_Com Es_d;
		    Es_d.assign(N0);	
	     Vec_dou D_Es_d;
		    D_Es_d.assign(N0);
    	 Vec_dou DB_Es_d;
		    DB_Es_d.assign(N0);

		double r_p_q=0.;
		 
		Vec_dou Rp_q;
		       Rp_q.assign(2*N0);




		 
		 Vec_Com Es1;
		    Es1.assign(N0);	
		 Vec_dou Sita_D;
		    Sita_D.assign(N0);
		 Vec_dou x_D;
		     x_D.assign(N0);
		 Vec_dou z_D;
		     z_D.assign(N0);

	    double abx;

    i_an = i_an*PI/180;
    K = 2.*PI/wave_len;
	K1 = K*sqrt(Eps_d);

//	cout<<sqrt(Eps_d)<<"  "<<K1<<endl;


		double KG = 1.0-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
        double KG1= 1.0-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
	    

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
			x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);

		 }

// g_p = rl/4.;
   mean = 0.;
   vari = 1.;
   seed=123456.;
   RL=rl;
   k_w=K;

         Complex Tou(0.0,0.0);
		 Complex Fan(0.0,0.0);
		 double  Dita_tou;
		 double  Dita_fan;


	for (int nn=0;nn<NUMBER_REALIZATION;nn++)
	{

    cout<<"Begin the "<<nn+1<<" th realization."<<endl;


	gauss_surface();

//   One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);

    double dx_Suf=RL/N;
	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*k_w*xsz*(1.0+ww))*exp(-fac);
	    U_in[i] = Ei[i]*(J*K*(1.0+ww)*(-df_x[i]*sin(i_an)-cos(i_an))+4.0*J*K*xsz*(x[i]+z[i]*tan(i_an))/pow(g_p*K*g_p*cos(i_an),2.0)*(tan(i_an)-df_x[i])+2.0*(x[i]+z[i]*tan(i_an))/pow(g_p,2.0)*(df_x[i]-tan(i_an)));	
	}


	


//////////////////////////Fill the impedance matrix////////////////////////////////////////////////

//	resolve_FBM(ZZ,EE,JJ,2*N);

/////////////////////////////////////////////////////////////////////////
	
	//The number of the strong interaction point.                          
 
	

   cout<<"Begin to fill the impedance matrix......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{			
			if (j==i)
			{
			     ZT_a[i*N+j] = 0.5-ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm);
			}		 
		}
	}


//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{			
			if(i==j)
			{
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
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
			     ZT_c[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(rm);
                 ZT_c[i*N+j] = J*K1/4.*dx_Suf*NR*Hankel(1,1,K1*rm);
			}		 
		}
	}


//////////////////////////////////////////////////////////////////////////

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{			
			if(i==j)
			{
				ZT_d[i*N+j]  =  -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_d[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*rm);
			}
		}
	}



  cout<<"The matrix has been filled."<<endl;

    Vec_Com JJ;
	    JJ.assign(2*N);
// 	resolve_DIL_FBM_SAA(ZT_a,ZT_b,ZT_c,ZT_d,Ei,JJ,max_z,min_z,dx_Suf,K,K1,rou_0,z,df_x);

/*
  Vec_Com ZZ;
       ZZ.assign(2*N*2*N);

	   for (i=0;i<N;i++)
	   {
		   for (j=0;j<N;j++)
		   {
               ZZ[i*N+j] = ZT_a[i*N+j];
			   ZZ[i*N+j+N] = ZT_b[i*N+j];
			   ZZ[(i+N)*N+j]= ZT_c[i*N+j];
			   ZZ[(i+N)*N+j+N] = ZT_d[i*N+j];
		   }
	   }
  Vec_Com EE;
       EE.assign(2*N);
	   for (i=0;i<N;i++)
	   {
		   EE[i] = Ei[i];
	   }*/

 
 
//////////////////////////////////////The window function.//////////////////////////////////////////////		   


/////////////////////////////////////Begin the FFMmethod/////////////////////////////////////////////////////////////////////

	
for (m=0;m<N;m++)
  {

	cout<<"Begin the "<<m+1<<"th step of the FBM ......------>>"<<endl;


////////////////////////////////The point which can not use only have the strong relative ,from 0 to N_s///////////////////////////////////////////////////////
	
	
	for (i=0;i<N;i++)
	{
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];
		}
	

		Jf_1[i] = ((Ei[i]-Sum0)*ZT_d[i*N+i]-(-1.*Sum1)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1)*ZT_a[i*N+i]-(Ei[i]-Sum0)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
	
	} 


/////////////////////////////////The point which can not use only have the strong relative ,from N-1 to N-N_s-1//////////////////////////////////////////////////////////


	for (i=N-1;i>=0;i--)
	{
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];

		}


    	Jb_1[i] = ((-1.*Sum0)*ZT_d[i*N+i]-(-1.*Sum1)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*ZT_a[i*N+i]-(-1.*Sum0)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
	
	}


////////////////////////////////the convergence condition ////////////////////////////////////////////////////////////

		     CLinearEquations num;
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta ="<<Delta<<endl;
			  
			 if(Delta<=0.01)
			 {
			    cout<<"The total number of the step is :"<<m+1<<"."<<"and the final error is :"<<Delta<<endl;
		        break;
			 }
				
	          else  
		           continue;

}
         

//////////////////////////////////The Es/////////////////////////////////////////////////////////////     


 


		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
			    Sum1 = Sum1 + dx_Suf*(J*K*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))*(Jf_1[j]+Jb_1[j])-(Jf_2[j]+Jb_2[j]))*exp(-J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }

			 Es[i] = Sum1;
	 
             D_Es[i]=(nn*D_Es[i]+pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG))/(nn+1);
             DB_Es[i]=(nn*DB_Es[i]+10.*log10(D_Es[i]))/(nn+1);
			 
			 Sum1 = zero;
             r_p_q = r_p_q+D_Es[i]*PI/180;
             
		 }

		 
 //       for (i=0;i<N0;i++)
//		 {
//		     Sum1=zero;
//			 double Jiaodu = (i)*PI/180.0 - PI/2.0;
			 for (j=0;j<N;j++)
			 {				     
//			    Sum1 = Sum1 + dx_Suf*(-J*K1*(df_x[j]*sin(Jiaodu)-cos(Jiaodu))/(sqrt(1.+pow(df_x[j],2.0)))*(Jf_1[j]+Jb_1[j])+(Jf_2[j]+Jb_2[j]))*exp(-J*K1*(x[j]*sin(Jiaodu)+z[j]*cos(Jiaodu)));							   				
			 
		        Tou = Tou + dx_Suf*(Jf_1[j]+Jb_1[j])*conj((Jf_2[j]+Jb_2[j]));
				Fan = Fan + dx_Suf*(Jf_1[j]+Jb_1[j]-Ei[j])*conj(Jf_2[j]+Jb_2[j]-U_in[j]);
			 }

			 cout<<Fan<<endl;

			 Dita_tou = Tou.imag()/(K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
			 Dita_fan = Fan.imag()/(K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
 
			 

  //           Es_d[i] = Sum1;
			 
   //         D_Es_d[i]=pow(abs(Es_d[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
  //          DB_Es_d[i]=10.*log10(D_Es_d[i]);			

 //          r_p_q = r_p_q+D_Es_d[i]*PI/180; 
//		 }

		 for (i=0;i<N0;i++)
		 {
			 Rp_q[i] = DB_Es[N0-i-1];
			 Rp_q[i+N0] = DB_Es_d[i];
			
		 }
	 }

		 cout<<"r_p_q="<<r_p_q<<endl;
		 cout<<"Dita_tou ="<<Dita_tou<<endl;
		 cout<<"Dita_fan ="<<Dita_fan<<endl;
		 cout<<"Dita_tou+Dita_fan="<<Dita_tou+fabs(Dita_fan)<<endl;

	



		 


  /*  	 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
				 rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);			
			 
			     Es[i] = Es[i]+J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm)*(Jf_1[j]+Jb_1[j])-(Jf_2[j]+Jb_2[j])*dx_Suf*J/4.*Hankel(0,1,K*rm);
			 }
			 
            D_Es[i]=R_D*pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

//			cout<<"Es[i] = "<<Es[i]<<endl;
             
		 }
*/


/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm*sqrt(1+pow(df_x[j],2.0)));
					Es[i] = Es[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*rm)*(Jf_2[j]+Jb_2[j]));
				}
		
			}
			


			for (i=0;i<N0;i++)
			{
				abx   =  1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
				D_Es[i] = pow(abs(Es[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*abx);
				DB_Es[i] = 10.*log10(D_Es[i]);
			}*/



	ofstream out_Jf1_Jf2;
	         out_Jf1_Jf2.open("DIL_FBM_Jf1_Jf2.txt");
			 for (i=0;i<N;i++)
			 {
               out_Jf1_Jf2<<x[i]/wave_len<<"   "<<abs(Jf_1[i])<<"   "<<abs(Jf_2[i])<<"   "<<abs(Jf_1[i]+Jb_1[i])<<endl;
			 }
	ofstream out_Jb1_Jb2;
	         out_Jb1_Jb2.open("DIL_FBM_Jb1_Jb2.txt");
			 for (i=0;i<N;i++)
			 {
               out_Jb1_Jb2<<x[i]/wave_len<<"   "<<abs(Jb_1[i])<<"   "<<abs(Jb_2[i])<<"   "<<abs(Jf_2[i]+Jb_2[i])<<endl;
			 }
/*	ofstream out_FN_b_1_2;
	         out_FN_b_1_2.open("FN_b_1_2.txt");
			 for (i=0;i<2*Q+1;i++)
			 {
               out_FN_b_1_2<<i<<"   "<<abs(FN_f_1[i])<<"   "<<abs(FN_f_2[i])<<endl;
			 }*/
	ofstream out_Es;
	         out_Es.open("DIL_FBM_Es.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es<<(i-90)<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
/*   ofstream out_Fai;
	         out_Fai.open("DIL_Fai.txt");
			 for (i=0;i<2*Q+1;i++)
			 {
               out_Fai<<i-Q<<"   "<<real(Fai[i])<<"   "<<imag(Fai[i])<<"   "<<abs(Fai[i])<<endl;
			 }*/

	ofstream out_Es_JJ;
	         out_Es_JJ.open("DIL_FBM_Es_JJ.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es_JJ<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	
    ofstream out_Es_d_JJ;
	         out_Es_d_JJ.open("DIL_FBM_Es_d.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es_d_JJ<<i-90<<"   "<<abs(Es_d[i])<<"   "<<D_Es_d[i]<<"   "<<DB_Es_d[i]<<endl;
			 }

	ofstream out_Rp_q;
	         out_Rp_q.open("DIL_FBM_Rp_q.txt");
			 for (i=0;i<2*N0;i++)
			 {
				 out_Rp_q<<i<<"      "<<Rp_q[i]<<endl;
			 }
	
}



void FBM_SAA::Gauss_DIL_Surface_FBM_SAA()
{

//////////////////2008.11.03////////////////////////

///////////////////姬伟杰///////////////////////////
	
	ofstream outtime("total_time.txt",ios::trunc);
	time_t ltime1,ltime2,ltime;
	time( &ltime1 );
	cout<<"The starting time is "<<ctime( &ltime1 )<<endl;
	//一维导体粗糙面

	int Idex_TE_TM=1;
	double wave_len=3.0e8/(3e9);      //wave length
	double h=0.5*wave_len;            //rms height 
	double lc=1.0*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

	double R=1.*wave_len;

	double g_p = 6.4*wave_len;
	
//	double g_p = 4.*R/cos(i_an*PI/180.);
   
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;

	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	int N=(rl/wave_len)*10;
	cout<<"未知量 N= "<<N<<" ;"<<endl;
	g_p = rl/4.0;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;

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

	double GAMA = 0.577216;
	double rm;

	double Bita_s ;
	double L_s;
	double R_s;
	double Q;
	double Fai_S_max;
	double Delta_0;
	double Delta_Fai;
	Complex K1;
     
	Complex zero(0.,0.);
	Complex Eps_d(6.0,0.6);
	Complex rou_0;

    if (Idex_TE_TM==1)
    {
       rou_0 = 1;
    }
	else
	{
		rou_0 = Eps_d;
	}

	int N_s;

	double max_z=0.;
	double min_z=0.;
	double dl;
	double NR;
	double Delta;


	Complex Fai_p;



	int i,j,k,m;
	//double w_Sita;

	Complex Sum0;
	Complex Sum1;
	Complex Sum2;
	Complex Sum3;

	

	Vec_Com Ei;
	    Ei.assign(N);
	Vec_Com Js0;
	    Js0.assign(N);
 
	Vec_Com Jf3;
	   Jf3.assign(2*N);  

	Vec_Com Jf2;
	   Jf2.assign(2*N);

	Vec_Com Jf1;
	   Jf1.assign(2*N);

	
    Vec_Com Jf_1;
	    Jf_1.assign(N);
	Vec_Com Jf_2;
	    Jf_2.assign(N);
	
	Vec_Com Jb_1;
	    Jb_1.assign(N);
		
	
	Vec_Com Jb_2;
	    Jb_2.assign(N);
	
		
	Vec_Com Ef;
	    Ef.assign(N);
	Vec_Com Eb;
	    Eb.assign(N);

	Vec_Com Ef_w;
	    Ef_w.assign(N);
	Vec_Com Eb_w;
	    Eb_w.assign(N);

	Vec_Com Ef_s;
	    Ef_s.assign(N);
	Vec_Com Eb_s;
	    Eb_s.assign(N);

	Vec_Com ZT_a;
        ZT_a.assign(N*N);
	Vec_Com ZT_b;
        ZT_b.assign(N*N);
	Vec_Com ZT_c;
        ZT_c.assign(N*N);
	Vec_Com ZT_d;
        ZT_d.assign(N*N);

    One_D_Conductor_Gauss(N,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);


	double dx_Suf=RL/N;

//	g_p = rl/4.;
    i_an = i_an*PI/180;
    K = 2.*PI/wave_len;
	K1 = K*sqrt(Eps_d);  cout<<"K1 = "<<K1<<endl;




	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*k_w*xsz*(1.0+ww))*exp(-fac);		
	}

/////////////////////////////////////////////////////////////////////////
		

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

	L_s = 10.0*wave_len;            //10*(max_z - min_z)/4. ;                                 //The strong interaction zone.
	
//	cout<<"L_s="<<L_s<<endl;

	Bita_s = sqrt(10./(K*L_s));                                //The geometric optics lit region.
	R_s = sqrt(pow(L_s,2.0)+pow(max_z-min_z,2.0));             
	Fai_S_max = atan((max_z-min_z)/L_s);                       //The outermost possible saddle point.

	cout<<"Fai_S_max="<<Fai_S_max<<endl;
	
	if (sqrt(K*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}

	
	Delta_Fai = sqrt(5./(K*R_s))/22.;                          
	Q = int(Bita_s/Delta_Fai)+5;

	N_s =int(L_s/dx_Suf);  
	
	cout<<"Q ="<<Q<<endl;
	cout<<"R_S = "<<R_s<<endl;
	cout<<"Fai_s_max ="<<Fai_S_max<<endl;
	cout<<"Delta_Fai="<<Delta_Fai<<endl;
	cout<<"Delta_0="<<Delta_0<<endl;
	
	cout<<"N_s="<<N_s<<endl;
	
// The number of the strong interaction point.                          
 
	

   cout<<"Begin to fill the impedance matrix......------>>"<<endl;

//////////////////////////////////////////////////////////////////////////

	for (i=0;i<=N_s;i++)
	{
		for (j=0;j<=i+N_s;j++)
		{			
			if (j==i)
			{
			     ZT_a[i*N+j] = 0.5-ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm);
			}		 
		}
	}


	for (i=N_s+1;i<N-N_s;i++)
	{
		for (j=i-N_s;j<=i+N_s;j++)
		{			
			if (j==i)
			{
			     ZT_a[i*N+j] = 0.5-ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm);
			}		 
		}
	}


	for (i=N-N_s;i<N;i++)
	{
		for (j=N-N_s;j<N;j++)
		{			
			if (j==i)
			{
			     ZT_a[i*N+j] = 0.5-ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm);
			}		 
		}
	}



//////////////////////////////////////////////////////////////////////////

	for (i=0;i<=N_s;i++)
	{
		for (j=0;j<=N_s;j++)
		{			
			if(i==j)
			{
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
			}
		}
	}


	for (i=N_s+1;i<N-N_s;i++)
	{
		for (j=i-N_s;j<=i+N_s;j++)
		{			
			if(i==j)
			{
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
			}
		}
	}

		
	for (i=N-N_s;i<N;i++)
	{
		for (j=N-N_s;j<N;j++)
		{			
			if(i==j)
			{
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
			}
		}
	}


//////////////////////////////////////////////////////////////////////////
	for (i=0;i<=N_s;i++)
	{
		for (j=0;j<=N_s;j++)
		{			
			if (j==i)
			{
			     ZT_c[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_c[i*N+j] = J*K1/4.*dx_Suf*NR*Hankel(1,1,K1*rm);
			}		 
		}
	}

	for (i=N_s+1;i<N-N_s;i++)
	{
		for (j=i-N_s;j<=i+N_s;j++)
		{			
			if (j==i)
			{
			     ZT_c[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_c[i*N+j] = J*K1/4.*dx_Suf*NR*Hankel(1,1,K1*rm);
			}		 
		}
	}

	for (i=N-N_s;i<N;i++)
	{
		for (j=N-N_s;j<N;j++)
		{			
			if (j==i)
			{
			     ZT_c[i*N+j] = 0.5+ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_c[i*N+j] = J*K1/4.*dx_Suf*NR*Hankel(1,1,K1*rm);
			}		 
		}
	}

//////////////////////////////////////////////////////////////////////////

	for (i=0;i<=N_s;i++)
	{
		for (j=0;j<=N_s;j++)
		{			
			if(i==j)
			{
				ZT_d[i*N+j]  =  -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_d[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*rm);
			}
		}
	}

	for (i=N_s+1;i<N-N_s;i++)
	{
		for (j=i-N_s;j<=i+N_s;j++)
		{			
			if(i==j)
			{
				ZT_d[i*N+j]  =  -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_d[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*rm);
			}
		}
	}


	for (i=N-N_s;i<N;i++)
	{
		for (j=N-N_s;j<N;j++)
		{			
			if(i==j)
			{
				ZT_d[i*N+j]  =  -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_d[i*N+j] = -1.*rou_0*dx_Suf*J/4.*Hankel(0,1,K1*rm);
			}
		}
	}

  cout<<"The matrix has been filled."<<endl;

    Vec_Com JJ;
	    JJ.assign(2*N);
// 	resolve_DIL_FBM_SAA(ZT_a,ZT_b,ZT_c,ZT_d,Ei,JJ,max_z,min_z,dx_Suf,K,K1,rou_0,z,df_x);

/*
  Vec_Com ZZ;
       ZZ.assign(2*N*2*N);

	   for (i=0;i<N;i++)
	   {
		   for (j=0;j<N;j++)
		   {
               ZZ[i*N+j] = ZT_a[i*N+j];
			   ZZ[i*N+j+N] = ZT_b[i*N+j];
			   ZZ[(i+N)*N+j]= ZT_c[i*N+j];
			   ZZ[(i+N)*N+j+N] = ZT_d[i*N+j];
		   }
	   }
  Vec_Com EE;
       EE.assign(2*N);
	   for (i=0;i<N;i++)
	   {
		   EE[i] = Ei[i];
	   }*/

    Vec_Com FN_f_1;FN_f_1.assign(2*Q+1);           
	Vec_Com FN_b_1;FN_b_1.assign(2*Q+1);		   
    Vec_Com FN_f_2;FN_f_2.assign(2*Q+1);           
	Vec_Com FN_b_2;FN_b_2.assign(2*Q+1);
    Vec_Com Fai;Fai.assign(2*Q+1);
	Vec_dou w_Sita;w_Sita.assign(2*Q+1);
	       
 
//////////////////////////////////////The window function.//////////////////////////////////////////////		   

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);               

		 if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;
		    Fai[k+Q]=Fai_p;
	}

/////////////////////////////////////Begin the FFM_SAA method/////////////////////////////////////////////////////////////////////

	
for (m=0;m<N;m++)
  {

	cout<<"Begin the "<<m+1<<"th step of the FBM_SAA ......------>>"<<endl;


////////////////////////////////The point which can not use only have the strong relative ,from 0 to N_s///////////////////////////////////////////////////////
	
	
	for (i=0;i<N_s+1;i++)
	{
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];
		}
		Jf_1[i] = ((Ei[i]-Sum0)*ZT_d[i*N+i]-(-1.*Sum1)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1)*ZT_a[i*N+i]-(Ei[i]-Sum0)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
	}

/////////////////////////////////The special point N_s+1//////////////////////////////////////////////////////////
 
	    i=N_s+1;
    	for (j=1;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];
		}		
		for (k=-Q;k<=Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
								  
            FN_f_1[k+Q] = ((-1.*J*K*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])+(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i-N_s-1]*sin(Fai_p));
            FN_f_2[k+Q] = ((J*K1*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])-rou_0*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K1*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K1*z[i-N_s-1]*sin(Fai_p));
			

			Sum2 = Sum2 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(J*K*z[N_s+1]*sin(Fai_p))*Delta_Fai;
        	Sum3 = Sum3 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(J*K1*z[N_s+1]*sin(Fai_p))*Delta_Fai;
		}	
		
	   Jf_1[i] = ((Ei[i]-Sum0-Sum2)*ZT_d[i*N+i]-(-1.*Sum1-Sum3)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
	   Jf_2[i] = ((-1.*Sum1-Sum3)*ZT_a[i*N+i]-(Ei[i]-Sum0-Sum2)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
       
	   
	   Sum0 = zero;
	   Sum1 = zero;
	   Sum2 = zero;
	   Sum3 = zero;	
		
//////////////////////////////////The weak  point for the forward/////////////////////////////////////////////////////////

	  for (i=N_s+2;i<N;i++)
	  {
		for (j=i-N_s;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];
		}
		for (k=-Q;k<=Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
					  
            FN_f_1[k+Q] = FN_f_1[k+Q]*exp(J*K*dx_Suf*cos(Fai_p)) + ((-1.*J*K*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])+(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i-N_s-1]*sin(Fai_p));
            FN_f_2[k+Q] = FN_f_2[k+Q]*exp(J*K1*dx_Suf*cos(Fai_p)) + ((J*K1*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])-rou_0*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K1*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K1*z[i-N_s-1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(J*K*z[i]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 									

		}
		Jf_1[i] = ((Ei[i]-Sum0-Sum2)*ZT_d[i*N+i]-(-1.*Sum1-Sum3)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1-Sum3)*ZT_a[i*N+i]-(Ei[i]-Sum0-Sum2)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);

	    Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
	}


/////////////////////////////////The point which can not use only have the strong relative ,from N-1 to N-N_s-1//////////////////////////////////////////////////////////


	for (i=N-1;i>=N-N_s-1;i--)
	{
		for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];

		}
    	Jb_1[i] = ((-1.*Sum0)*ZT_d[i*N+i]-(-1.*Sum1)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*ZT_a[i*N+i]-(-1.*Sum0)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;	
	}


////////////////////////////////The special point N-N_s-2///////////////////////////////////////////////////////////
		
	   i=N-N_s-2;
	   for (j=i+1;j<N;j++)
	   {
		   Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
		   Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];
	   }   
	   for (k=-Q;k<=Q;k++)
	   {
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
					  
            FN_b_1[k+Q] = ((-1.*J*K*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])+(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i+N_s+1]*sin(Fai_p));
            FN_b_2[k+Q] = ((J*K1*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])-rou_0*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K1*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K1*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(J*K*z[N-N_s-2]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(J*K1*z[N-N_s-2]*sin(Fai_p))*Delta_Fai; 												
		}
	   
	   Jb_1[i] = ((-1.*Sum0-Sum2)*ZT_d[i*N+i]-(-1.*Sum1-Sum3)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
	   Jb_2[i] = ((-1.*Sum1-Sum3)*ZT_a[i*N+i]-(-1.*Sum0-Sum2)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
	   Sum0 = zero;
	   Sum1 = zero;
	   Sum2 = zero;
	   Sum3 = zero;

///////////////////////////////The weak point for the behind ////////////////////////////////////////////////////////////


	 for (i=N-N_s-3;i>=0;i--)
	 {
		for (j=i+1;j<i+N_s+1;j++)
		{
		    Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*ZT_a[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*ZT_c[i*N+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N+j];	
		}

		for (k=-Q;k<=Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
					  
            FN_b_1[k+Q] = FN_b_1[k+Q]*exp(J*K*dx_Suf*cos(Fai_p)) + ((-1.*J*K*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])+(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K*z[i+N_s+1]*sin(Fai_p));
            FN_b_2[k+Q] = FN_b_2[k+Q]*exp(J*K1*dx_Suf*cos(Fai_p)) + ((J*K1*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])-rou_0*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K1*(N_s+1.)*dx_Suf*cos(Fai_p))*exp(-1.*J*K1*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(J*K*z[i]*sin(Fai_p))*Delta_Fai; 									
        	Sum3 = Sum3 + J*dx_Suf/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 												
		}

     	Jb_1[i] = ((-1.*Sum0-Sum2)*ZT_d[i*N+i]-(-1.*Sum1-Sum3)*ZT_b[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1-Sum3)*ZT_a[i*N+i]-(-1.*Sum0-Sum2)*ZT_c[i*N+i])/(ZT_a[i*N+i]*ZT_d[i*N+i]-ZT_b[i*N+i]*ZT_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
	}

////////////////////////////////the convergence condition ////////////////////////////////////////////////////////////

		     CLinearEquations num;
 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }
             Assign_Vecor(Jf3,Jf1);
             Delta = num.VectorNorm2(Jf2,2*N)/num.VectorNorm2(Jf3,2*N);

             cout<<"The error is: Delta ="<<Delta<<endl;
			  
			 if(Delta<=0.01)
			 {
			    cout<<"The total number of the step is :"<<m+1<<"."<<"and the final error is :"<<Delta<<endl;
		        break;
			 }
	          else  
		           continue;
}
         

//////////////////////////////////The Es/////////////////////////////////////////////////////////////     

         double R_D = 10000.*wave_len;
		 int N0=180;

		 Vec_Com Es;
		    Es.assign(N0);
		 Vec_dou D_Es;
		    D_Es.assign(N0);
    	 Vec_dou DB_Es;
		    DB_Es.assign(N0);
		
		 Vec_Com Es_d;
		    Es_d.assign(N0);	
	     Vec_dou D_Es_d;
		    D_Es_d.assign(N0);
    	 Vec_dou DB_Es_d;
		    DB_Es_d.assign(N0);

		double r_p_q=0.;
		 
		Vec_dou Rp_q;
		       Rp_q.assign(2*N0);




		 
		 Vec_Com Es1;
		    Es1.assign(N0);	
		 Vec_dou Sita_D;
		    Sita_D.assign(N0);
		 Vec_dou x_D;
		     x_D.assign(N0);
		 Vec_dou z_D;
		     z_D.assign(N0);

	    double abx;


		double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
 //       Complex KG1= 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));
	    

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
			x_D[i] = R_D*sin(Sita_D[i]);
			z_D[i] = R_D*cos(Sita_D[i]);
		 }

		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
			    Es[i] = Es[i] + dx_Suf*(J*K*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*(Jf_1[j]+Jb_1[j])-(Jf_2[j]+Jb_2[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);             
		 }
		 
/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
			    Es[i] = Es[i] + dx_Suf*(J*K*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*JJ[j]-JJ[j+N])*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

            r_p_q = r_p_q+D_Es[i];
		 }*/

 /*	   for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
			    Es_d[i] = Es_d[i] + dx_Suf*(J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*JJ[j]-rou_0*JJ[j+N])*exp(-1.*J*K1*(x[j]*sin(Sita_D[i])-z[j]*cos(Sita_D[i])));							   				
			 }
			 
            D_Es_d[i]=pow(abs(Es_d[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG1)/abs(rou_0);
            DB_Es_d[i]=10.*log10(D_Es_d[i]);			

            r_p_q = r_p_q+D_Es_d[i]; 
		 }


		 for (i=0;i<N0;i++)
		 {
			 Rp_q[i] = DB_Es[N0-i-1];
			 Rp_q[i+N0] = DB_Es_d[i];
			
		 }*/

		 cout<<"r_p_q="<<r_p_q<<endl;

	



		 


  /*  	 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
				 rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);			
			 
			     Es[i] = Es[i]+J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm)*(Jf_1[j]+Jb_1[j])-(Jf_2[j]+Jb_2[j])*dx_Suf*J/4.*Hankel(0,1,K*rm);
			 }
			 
            D_Es[i]=R_D*pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

//			cout<<"Es[i] = "<<Es[i]<<endl;
             
		 }
*/


/*			for (i=0;i<N0;i++)
			{
				for (j=0;j<N;j++)
				{
					rm = sqrt(pow(x_D[i]-x[j],2.0)+pow(z_D[i]-z[j],2.0));
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm*sqrt(1+pow(df_x[j],2.0)));
					Es[i] = Es[i]+dx_Suf*(J*K/4.*NR*Hankel(1,1,K*rm)*(Jf_1[j]+Jb_1[j])-J/4.*Hankel(0,1,K*rm)*(Jf_2[j]+Jb_2[j]));
				}
		
			}
			


			for (i=0;i<N0;i++)
			{
				abx   =  1.-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
				D_Es[i] = pow(abs(Es[i]),2.0)/(g_p*sqrt(PI/2.)*cos(i_an)*abx);
				DB_Es[i] = 10.*log10(D_Es[i]);
			}*/



	ofstream out_Jf1_Jf2;
	         out_Jf1_Jf2.open("DIL_FBM_SAA_Jf1_Jf2.txt");
			 for (i=0;i<N;i++)
			 {
               out_Jf1_Jf2<<x[i]/wave_len<<"   "<<abs(Jf_1[i])<<"   "<<abs(Jf_2[i])<<"   "<<abs(Jf_1[i]+Jb_1[i])<<endl;
			 }
	ofstream out_Jb1_Jb2;
	         out_Jb1_Jb2.open("DIL_FBM_SAA_Jb1_Jb2.txt");
			 for (i=0;i<N;i++)
			 {
               out_Jb1_Jb2<<x[i]/wave_len<<"   "<<abs(Jb_1[i])<<"   "<<abs(Jb_2[i])<<"   "<<abs(Jf_2[i]+Jb_2[i])<<endl;
			 }
/*	ofstream out_FN_b_1_2;
	         out_FN_b_1_2.open("FN_b_1_2.txt");
			 for (i=0;i<2*Q+1;i++)
			 {
               out_FN_b_1_2<<i<<"   "<<abs(FN_f_1[i])<<"   "<<abs(FN_f_2[i])<<endl;
			 }*/
	ofstream out_Es;
	         out_Es.open("DIL_FBM_SAA_Es.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es<<(i-90)<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
/*   ofstream out_Fai;
	         out_Fai.open("DIL_Fai.txt");
			 for (i=0;i<2*Q+1;i++)
			 {
               out_Fai<<i-Q<<"   "<<real(Fai[i])<<"   "<<imag(Fai[i])<<"   "<<abs(Fai[i])<<endl;
			 }*/

	ofstream out_Es_JJ;
	         out_Es_JJ.open("DIL_FBM_SAA_Es_JJ.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es_JJ<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	
    ofstream out_Es_d_JJ;
	         out_Es_d_JJ.open("DIL_FBM_SAA_Es_d_JJ.txt");
			 for (i=0;i<N0;i++)
			 {
               out_Es_d_JJ<<i-90<<"   "<<abs(Es_d[i])<<"   "<<D_Es_d[i]<<"   "<<DB_Es_d[i]<<endl;
			 }

	ofstream out_Rp_q;
	         out_Rp_q.open("DIL_FBM_SAA_Rp_q.txt");
			 for (i=0;i<2*N0;i++)
			 {
				 out_Rp_q<<i<<"      "<<Rp_q[i]<<endl;
			 }

	time( &ltime2 );
	cout<<"The end time is "<<ctime(&ltime2 )<<endl;
	outtime<<"The end time is "<<ctime(&ltime2 )<<endl;
	ltime = ltime2 - ltime1;
	cout<<"The total time is "<< ltime <<"seconds"<<endl;
	outtime<<"The total time is "<<ltime<<"seconds"<<endl;
	
}

void FBM_SAA::Gauss_DIL_Surface_PBTC_FBM_SAA()
{

//////////////////2008.11.8///////////////////////

///////////////////姬伟杰/////////////////////////
	
	ofstream outtime("total_time.txt",ios::trunc);
	time_t ltime1,ltime2,ltime;
	time( &ltime1 );
	cout<<"The starting time is "<<ctime( &ltime1 )<<endl;
	//一维导体粗糙面
	double wave_len=3.0e8/3e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

//	double R=1.*wave_len;
//	double g_p = 4.*R/cos(i_an*PI/180.);
    double g_p =5.*wave_len;
	
	cout<<"入射角 incident = "<<i_an<<" 度 ;"<<endl;
	cout<<"锥形参数 g = "<<g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<"粗糙面长度 RL = "<<4.0*g_p/wave_len<<"*Wavelength ;"<<endl;
	cout<<endl;
    double rl;
	cout<<"为了可以应用任意维数的 FFT"<<endl;
	cout<<"粗糙面长度 RL 的取值应为 ";
	cin>>rl;
	rl=rl*wave_len;
	int N=(rl/wave_len)*20;
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

	double GAMA = 0.577216;
	double rm;

	double Bita_s ;
	double L_s;
	double R_s;
	double Q;
	double Fai_S_max;
	double Delta_0;
	double Delta_Fai;
	double K1;
     
	double x_p;
	double z_p;
    double df_x_p;


	Complex zero(0.,0.);
	Complex Eps_d(25.,1.);
	Complex rou_0 = Eps_d;

	int N_s;

//	int n1= int(sqrt(abs(Eps_d)));
	int n1=3;
	int N_dg = 20*n1*rl/wave_len;

	double max_z=0.;
	double min_z=0.;
	double dl;
	double NR;
	double Delta;


	Complex Fai_p;



	int i,j,k,m,l;
	//double w_Sita;

	Complex Sum0(0.,0.);
	Complex Sum1(0.,0.);
	Complex Sum2(0.,0.);
	Complex Sum3(0.,0.);

	

	Vec_Com Ei;
	    Ei.assign(N);
 
	Vec_Com Jf3;
	   Jf3.assign(2*N_dg);  

	Vec_Com Jf2;
	   Jf2.assign(2*N_dg);

	Vec_Com Jf1;
	   Jf1.assign(2*N_dg);

	
    Vec_Com Jf_1;
	    Jf_1.assign(N);
	Vec_Com Jb_1;
	    Jb_1.assign(N);	
		
	Vec_Com Jf_2;
	    Jf_2.assign(N_dg);		
	Vec_Com Jb_2;
	    Jb_2.assign(N_dg);

	Vec_Com ZT_a;
        ZT_a.assign(N*N);
	Vec_Com ZT_b;
        ZT_b.assign(N*N);
	Vec_Com ZT_c;
        ZT_c.assign(N_dg*N_dg);
	Vec_Com ZT_d;
        ZT_d.assign(N_dg*N_dg);

	
cout<<"n1="<<n1<<endl;
cout<<"N_dg="<<N_dg<<endl;
    One_D_Conductor_Gauss(N_dg,h,lc,us,kl,ku,rl,wave_len,i_an,n_r,n_sc,rd_bw,g_p,eps);
	
	
    
	double dx_Suf_Ndg = RL/N_dg;
	double dx_Suf_N = RL/N;

	
    i_an = -i_an*PI/180;
    K = 2.*PI/wave_len;
	K1 = K*sqrt(abs(Eps_d));

	Vec_dou x_N;
	     x_N.assign(N);
	Vec_dou z_N;
	     z_N.assign(N);
	Vec_dou df_x_N;
	     df_x_N.assign(N);
    
	if (n1/2==0)
    {
	   for (i=0;i<N;i++)
		{
		   x_N[i] = (x[i*n1+n1/2-1]+x[i*n1+n1/2])/2.;
		   z_N[i] = (z[i*n1+n1/2-1]+z[i*n1+n1/2])/2.;
		   df_x_N[i] = (df_x[i*n1+n1/2-1]+df_x[i*n1+n1/2])/2.;

		}
    }
	else if (n1/2!=0)
	{
		for (i=0;i<N;i++)
		{
			x_N[i] = x[i*n1+(n1-1)/2];
            z_N[i] = z[i*n1+(n1-1)/2];
			df_x_N[i] = df_x[i*n1+(n1-1)/2];
		}
	}
    


	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x_N[i]*sin(i_an)-z_N[i]*cos(i_an);
		kg=pow(k_w*g_p*cos(i_an),2.0);
		fac=pow((x_N[i]+z_N[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(-1.*J*k_w*xsz*(1.0+ww))*exp(-fac);		
	}




//////////////////////////Fill the impedance matrix////////////////////////////////////////////////


    cout<<"Begin to fill the impedance matrix......------>>"<<endl;



	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{	
		    if (n1%2==0)
		    {
			    if (j==i)
				{
			     ZT_a[i*N+j] = 0.5-0.5*(ddf_x[j*n1+n1/2-1]+ddf_x[j*n1+n1/2])*dx_Suf_N/(4.*PI*(1.+pow(0.5*(df_x[j*n1+n1/2-1]+df_x[j*n1+n1/2]),2.0)));
				}

			    else
				{
				 x_p =  0.5*(x[i*n1+n1/2-1]+x[i*n1+n1/2])-0.5*(x[j*n1+n1/2-1]+x[j*n1+n1/2]);
				 z_p =  0.5*(z[i*n1+n1/2-1]+z[i*n1+n1/2])-0.5*(z[j*n1+n1/2-1]+z[j*n1+n1/2]);
				 df_x_p = 0.5*(df_x[j*n1+n1/2-1]+df_x[j*n1+n1/2]);
				 rm = sqrt(pow(x_p,2.0)+pow(z_p,2.0));
				 NR = (-1.*df_x_p*x_p+z_p)/(sqrt(1.+pow(df_x_p,2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf_N*NR*Hankel(1,1,K*rm);
				}		 				
		    }
			
			else if (n1%2!=0)
			{
		        if (j==i)
				{
			     ZT_a[i*N+j] = 0.5-ddf_x[j*n1+(n1-1)/2]*dx_Suf_N/(4.*PI*(1.+pow(df_x[j*n1+(n1-1)/2],2.0)));
				}

			    else
				{
				 x_p =  x[i*n1+(n1-1)/2]-x[j*n1+(n1-1)/2];
				 z_p =  z[i*n1+(n1-1)/2]-z[j*n1+(n1-1)/2];

				 	
				 rm = sqrt(pow(x_p,2.0)+pow(z_p,2.0));
				 NR = -1.*df_x[j*n1+(n1-1)/2]*(x_p+z_p)/(sqrt(1.+pow(df_x[j*n1+(n1-1)/2],2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf_N*NR*Hankel(1,1,K*rm);
				}		 				
			}
		}
	}

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			if (n1%2==0)
			{
			    if(i==j)
				{
				 ZT_b[i*N+j] = dx_Suf_N*J/4.*Hankel(0,1,K*dx_Suf_N*sqrt(1.+pow(0.5*(df_x[j*n1+n1/2-1]+df_x[j*n1+n1/2]),2.0))/(2.*ee));
				}
			   else
			   {
				 x_p =  0.5*(x[i*n1+n1/2-1]+x[i*n1+n1/2])-0.5*(x[j*n1+n1/2-1]+x[j*n1+n1/2]);
				 z_p =  0.5*(z[i*n1+n1/2-1]+z[i*n1+n1/2])-0.5*(z[j*n1+n1/2-1]+z[j*n1+n1/2]);
            	 rm = sqrt(pow(x_p,2.0)+pow(z_p,2.0));
				 ZT_b[i*N+j] = dx_Suf_N*J/4.*Hankel(0,1,K*rm);
			   }

			}
			else if (n1%2!=0)
			{
			    if(i==j)
				{
				 ZT_b[i*N+j] = dx_Suf_N*J/4.*Hankel(0,1,K*dx_Suf_N*sqrt(1.+pow(0.5*df_x[j*n1+(n1-1)/2],2.0))/(2.*ee));
				}
			   else
			   {
				 x_p =  x[i*n1+(n1-1)/2]-x[j*n1+(n1-1)/2];
				 z_p =  z[i*n1+(n1-1)/2]-z[j*n1+(n1-1)/2];
            	 rm = sqrt(pow(x_p,2.0)+pow(z_p,2.0));
				 ZT_b[i*N+j] = dx_Suf_N*J/4.*Hankel(0,1,K*rm);
			   }

			}

		}
	}

/* 	for (i=0;i<N_dg;i++)
	{
		for (j=0;j<N_dg;j++)
		{			
			if (j==i)
			{
			     ZT_a[i*N+j] = 0.5-ddf_x[j]*dx_Suf/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_a[i*N+j] = -1.*J*K/4.*dx_Suf*NR*Hankel(1,1,K*rm);
			}		 
		}
	}

	for (i=0;i<N_dg;i++)
	{
		for (j=0;j<N_dg;j++)
		{
			if(i==j)
			{
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*dx_Suf*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_b[i*N+j] = dx_Suf*J/4.*Hankel(0,1,K*rm);
			}
		}
	}*/


   	for (i=0;i<N_dg;i++)
	{
		for (j=0;j<N_dg;j++)
		{			
			if (j==i)
			{
			     ZT_c[i*N_dg+j] = 0.5+ddf_x[j]*dx_Suf_Ndg/(4.*PI*(1.+pow(df_x[j],2.0)));
			}

			else
			{
				 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 NR = (-1.*df_x[j]*(x[i]-x[j])+(z[i]-z[j]))/(sqrt(1.+pow(df_x[j],2.0))*rm);
                 ZT_c[i*N_dg+j] = J*K1/4.*dx_Suf_Ndg*NR*Hankel(1,1,K1*rm);
			}		 
		}
	}

	for (i=0;i<N_dg;i++)
	{
		for (j=0;j<N_dg;j++)
		{
			if(i==j)
			{
				ZT_d[i*N_dg+j]  =  -1.*rou_0*dx_Suf_Ndg*J/4.*Hankel(0,1,K1*dx_Suf_Ndg*sqrt(1.+pow(df_x[j],2.0))/(2.*ee));
			}
			else
			{
            	 rm = sqrt(pow(x[i]-x[j],2.0)+pow(z[i]-z[j],2.0));
				 ZT_d[i*N_dg+j] = -1.*rou_0*dx_Suf_Ndg*J/4.*Hankel(0,1,K1*rm);
			}
		}
	}

  cout<<"The matrix has been filled."<<endl;

        Vec_Com JJ;
	    JJ.assign(2*N);
	 

//	resolve_FBM(ZZ,EE,JJ,2*N);

/////////////////////////////////////////////////////////////////////////
		

	for (i=0;i<N_dg;i++)
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

	Vec_Com Jf_Ndg;
	       Jf_Ndg.assign(N_dg);
	Vec_Com Jf_N;
	       Jf_N.assign(N);
    


   resolve_DIL_PBTG_FBM_SAA(ZT_a,ZT_b,ZT_c,ZT_d,Ei,Jf_Ndg,Jf_N,max_z,min_z,dx_Suf_Ndg,K,K1,rou_0,z,z_N,df_x,df_x_N);

//	resolve_DIL_FBM_SAA(ZT_a,ZT_b,ZT_c,ZT_d,Ei,JJ,max_z,min_z,dx_Suf,K,K1,rou_0,z,df_x);

/*	L_s = (max_z - min_z)/4. ;                                 //The strong interaction zone.
	
//	cout<<"L_s="<<L_s<<endl;

	                                                           //The geometric optics lit region.
	R_s = sqrt(pow(L_s,2.0)+pow(max_z-min_z,2.0));             
	Bita_s = sqrt(10./(K*R_s));                                // Here the  L_s has been replaced by R_s.        
	Fai_S_max = atan((max_z-min_z)/L_s);                       //The outermost possible saddle point.
	
	if (sqrt(K*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}

//	cout<<"Delta_0="<<Delta_0<<endl;
	
	Delta_Fai = sqrt(5./(K*R_s))/22.;                          
	Q = int(Bita_s/Delta_Fai)+5;

	N_s = int(L_s/dx_Suf_Ndg);                                      //The number of the strong interaction point.                          
 
	cout<<"N_s = "<<N_s<<endl;
	cout<<"L_s = "<<L_s<<endl;

 

    Vec_Com FN_f_1;
           FN_f_1.assign(2*Q+1);
	Vec_Com FN_b_1;
		   FN_b_1.assign(2*Q+1);
    Vec_Com FN_f_2;
           FN_f_2.assign(2*Q+1);
	Vec_Com FN_b_2;
		   FN_b_2.assign(2*Q+1);

    Vec_Com Fai;
	       Fai.assign(2*Q+1);

    Vec_Com FN_s_1;
	        FN_s_1.assign(N);
    Vec_Com FN_s_2;
	        FN_s_2.assign(N_dg);
		
	
	Vec_Com FN_w_1;
	        FN_w_1.assign(N);
	Vec_Com FN_w_2;
	        FN_w_2.assign(N_dg);
  		  		   
	Vec_dou w_Sita;
	       w_Sita.assign(2*Q+1);
 
//////////////////////////////////////The window function.//////////////////////////////////////////////		   

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0); 

		 if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;
		    Fai[k+Q]=Fai_p;
	}
	
/////////////////////////////////////Begin the FFM_SAA method/////////////////////////////////////////////////////////////////////

	int i_n1,j_n1,N_s_n1;
	Complex Sum_J(0.,0.);
	Complex Sum_J2(0.,0.);
	
 for (m=0;m<N;m++)
  {

	cout<<"Begin the "<<m+1<<"th step of the FBM_SAA ......------>>"<<endl;


////////////////////////////////The point which can not use only have the strong relative ,from 0 to N_s///////////////////////////////////////////////////////
	
 
	for (i=0;i<N_s+1;i++)
	{
		
		 i_n1 = int(i/n1);
 
		      
		for (j=0;j<i;j++)
		{
             j_n1  = int(j/n1);		   
			 for (l=0;l<n1;l++)
			 {
		         Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];
			 }

			     Sum_J = Sum_J/double(n1);

		    Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];			   

				
		    Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];

			Sum_J = zero;
		}
		  Jf_1[i_n1] = ((Ei[i_n1]-Sum0)*ZT_d[i*N_dg+i]-(-1.*Sum1)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		  Jf_2[i]    = ((-1.*Sum1)*ZT_a[i_n1*N+i_n1]-(Ei[i_n1]-Sum0)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		   	
					  
		  Sum0 = zero;
		  Sum1 = zero;
		  
		  Sum_J2=zero;
		 
	
	}


	

/////////////////////////////////The special point N_s+1//////////////////////////////////////////////////////////////////////
 
	      i=N_s+1;


		  i_n1 = int(i/n1);
		
		  for (j=0;j<i;j++)
		  {
               j_n1  = int(j/n1);
			  
			   for (l=0;l<n1;l++)			  
			   {
		         Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];
			  
			   }

			   Sum_J = Sum_J/double(n1);
			   
               Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];			   
			 		        				
			   Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];
		  	   Sum_J = zero;
		  }
		        
		        N_s_n1 = int((i-N_s-1)/n1);
	  	  
	                   
               int NN = int((N_s)/n1)+1;

			                  
			   for (l=0;l<n1;l++)
				{
				   Sum_J = Sum_J + Jf_2[N_s_n1+l]+Jb_2[N_s_n1+l];
				}

			      Sum_J = Sum_J/double(n1);
	   
		      for (k=-Q;k<Q;k++)
			  {
			    Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
									        
				FN_f_1[k+Q] =  ((-1.*J*K*(-1.*sin(Fai_p)+df_x_N[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(J*K*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*J*K*z_N[N_s_n1]*sin(Fai_p));
				
				FN_f_2[k+Q] =  ((J*K1*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Sum_J)-rou_0*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K1*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*J*K1*z[i-N_s-1]*sin(Fai_p));
				
			   
				Sum2 = Sum2 + J*dx_Suf_N/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(J*K*z_N[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	    Sum3 = Sum3 + J*dx_Suf_Ndg/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 												
		 
			  }
		      Jf_1[i_n1] = ((Ei[i_n1]-Sum0-Sum2)*ZT_d[i*N_dg+i]-(-1.*Sum1-Sum3)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		      Jf_2[i]    = ((-1.*Sum1-Sum3)*ZT_a[i_n1*N+i_n1]-(Ei[i_n1]-Sum0-Sum2)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);        

 	 
		 Sum0 = zero;
		 Sum1 = zero;
		 Sum2 = zero;
		 Sum3 = zero;
		 Sum_J = zero;
		 Sum_J2=zero;
					 

//////////////////////////////////The weak  point for the foward /////////////////////////////////////////////////////////

 
	 for (i=N_s+2;i<N_dg;i++)
	  {
		 
		  i_n1 = int(i/n1);
		
		  for (j=i-N_s-1;j<i;j++)
		  {
                j_n1  = int(j/n1);
			  
				for (l=0;l<n1;l++)			  
				{
		          Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];			  
				}

			    Sum_J = Sum_J/double(n1);
				
                Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];
		   
			 		        				
				Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];
		        
				Sum_J = zero;

		  }

		          
		  N_s_n1 = int((i-N_s-1)/n1);
		 
		  int NN = int((N_s)/n1)+1;
            
		  for (l=0;l<n1;l++)
            {
				Sum_J = Sum_J + Jf_2[N_s_n1+l]+Jb_2[N_s_n1+l];
            }

			Sum_J = Sum_J/double(n1);

	 	  
	     for (k=-Q;k<Q;k++)
		 {
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);
			 
			FN_f_1[k+Q] = FN_f_1[k+Q]*exp(J*K*dx_Suf_N*cos(Fai_p)) + ((-1.*J*K*(-1.*sin(Fai_p)+df_x_N[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(J*K*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*J*K*z_N[N_s_n1]*sin(Fai_p));
            FN_f_2[k+Q] = FN_f_2[k+Q]*exp(J*K1*dx_Suf_Ndg*cos(Fai_p)) + ((J*K1*(-1.*sin(Fai_p)+df_x[i-N_s-1]*cos(Fai_p)))*(Sum_J)-rou_0*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(J*K1*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*J*K1*z[i-N_s-1]*sin(Fai_p));

			Sum2 = Sum2 + J*dx_Suf_N/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(J*K*z_N[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + J*dx_Suf_Ndg/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 												
		
		 }

		  Jf_1[i_n1] = ((Ei[i_n1]-Sum0-Sum2)*ZT_d[i*N_dg+i]-(-1.*Sum1-Sum3)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);

		  Jf_2[i]    = ((-1.*Sum1-Sum3)*ZT_a[i_n1*N+i_n1]-(Ei[i_n1]-Sum0-Sum2)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);        

		 Sum0 = zero;
		 Sum1 = zero;
		 Sum2 = zero;
		 Sum3 = zero;
		 Sum_J = zero;
	
	 	 
	 }	

        
		 

	


/////////////////////////////////The point which can not use only have the strong relative ,from N-1 to N-N_s-1//////////////////////////////////////////////////////////


	for (i=N_dg-1;i>=N_dg-N_s-1;i--)
	{
		i_n1 = int(i/n1);

		for (j=i+1;j<N_dg;j++)
		{
			  j_n1  = int(j/n1);

			  for (l=0;l<n1;l++)
			  {
		         Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];
			  }

			    Sum_J = Sum_J/double(n1);
			
             Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];
			
			 Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];
	
			 Sum_J = zero;
		}
      
		
		Jb_1[i_n1] = ((-1.*Sum0)*ZT_d[i*N_dg+i]-(-1.*Sum1)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		Jb_2[i]    = ((-1.*Sum1)*ZT_a[i_n1*N+i_n1]-(-1.*Sum0)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		
		Sum0 = zero;
		Sum1 = zero;
	
	
	}


////////////////////////////////The speclal point N-N_s-2///////////////////////////////////////////////////////////
	       i=N-N_s-2;

		   i_n1= int(i/n1);
	

		   for (j=i+1;j<N_dg;j++)
		   {
			    j_n1  = int(j/n1);

			  for (l=0;l<n1;l++)
			  {
		         Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];
			  }

			    Sum_J = Sum_J/double(n1);
			
                Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];
			
			    Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];

				Sum_J = zero;
		   }

		   		
		   N_s_n1 = int((i+N_s+1)/n1);
	       NN = int((N_s)/n1)+1;
	
		for (l=0;l<n1;l++)
		{
		     Sum_J = Sum_J + Jf_2[N_s_n1+l]+Jb_2[N_s_n1+l];
		}

			 Sum_J = Sum_J/double(n1);
        
		
	    for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);


		                					  
            FN_b_1[k+Q] =  ((-1.*J*K*(-1.*sin(Fai_p)+df_x_N[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(J*K*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*J*K*z[N_s_n1]*sin(Fai_p));
            FN_b_2[k+Q] =  ((J*K1*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Sum_J)-rou_0*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K1*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*J*K1*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + J*dx_Suf_N/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(J*K*z[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + J*dx_Suf_Ndg/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 												
		}
		
		   
		 Jb_1[i_n1] = ((-1.*Sum0-Sum2)*ZT_d[i*N_dg+i]-(-1.*Sum1-Sum3)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);

		 Jb_2[i]    = ((-1.*Sum1-Sum3)*ZT_a[i_n1*N+i_n1]-(-1.*Sum0-Sum2)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
	    
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
		Sum_J = zero;
	
			
	
///////////////////////////////The weak point for the behind ////////////////////////////////////////////////////////////


	 for (i=N_dg-N_s-3;i>=0;i--)
	 {
		 
		 i_n1 = int(i/n1);

	    for (j=i+1;j<i+N_s+1;j++)
		{		   		  
			j_n1  = int(j/n1);
			  
			for (l=0;l<n1;l++)
			  {
		         Sum_J = Sum_J + Jf_2[j_n1*n1+l]+Jb_2[j_n1*n1+l];
			  }

			    Sum_J = Sum_J/double(n1);			
            Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*ZT_a[i_n1*N+j_n1]+Sum_J*ZT_b[i_n1*N+j_n1];
			
			Sum1 = Sum1 + Sum_J*ZT_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*ZT_d[i*N_dg+j];

			Sum_J = zero;
		 
		}
		    
		
		int NN = int((N_s)/n1)+1;

		N_s_n1 = int((i+N_s+1)/n1);
        
		
		for (l=0;l<n1;l++)
		{
		      Sum_J = Sum_J + Jf_2[N_s_n1+l]+Jb_2[N_s_n1+l];
		}

			   Sum_J = Sum_J/double(n1);

	    for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*J*Delta_0);

 
		    
					  
            FN_b_1[k+Q] = FN_b_1[k+Q]*exp(J*K*dx_Suf_N*cos(Fai_p)) + ((-1.*J*K*(-1.*sin(Fai_p)+df_x_N[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(J*K*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*J*K*z[N_s_n1]*sin(Fai_p));
            FN_b_2[k+Q] = FN_b_2[k+Q]*exp(J*K1*dx_Suf_Ndg*cos(Fai_p)) + ((J*K1*(-1.*sin(Fai_p)+df_x[i+N_s+1]*cos(Fai_p)))*(Sum_J)-rou_0*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(J*K1*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*J*K1*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + J*dx_Suf_N/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(J*K*z[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + J*dx_Suf_Ndg/(4.*PI)*exp(-1.*J*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(J*K1*z[i]*sin(Fai_p))*Delta_Fai; 												
		}
		 
		
		Jb_1[i_n1] = ((-Sum0-Sum2)*ZT_d[i*N_dg+i]-(-1.*Sum1-Sum3)*ZT_b[i_n1*N+i_n1])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);

		Jb_2[i]    = ((-1.*Sum1-Sum3)*ZT_a[i_n1*N+i_n1]-(-Sum0-Sum2)*ZT_c[i*N_dg+i])/(ZT_a[i_n1*N+i_n1]*ZT_d[i*N_dg+i]-ZT_b[i_n1*N+i_n1]*ZT_c[i*N_dg+i]);
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
		Sum_J = zero;
	
	
	}

////////////////////////////////the convergence condition ////////////////////////////////////////////////////////////

		     CLinearEquations num;
		  

 		     for (i=0;i<N_dg;i++)
			 {
			    Jf2[i]   = (Jf_2[i]+Jb_2[i])-Jf1[i];
			    Jf3[i]   =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
             Delta = num.VectorNorm2(Jf2,N_dg)/num.VectorNorm2(Jf3,N_dg);

             cout<<"The error is: Delta ="<<Delta<<endl;
			  
			 if(Delta<=0.01)
			 {
			    cout<<"The total number of the step is :"<<m+1<<"."<<"and the final error is :"<<Delta<<endl;
		        break;
			 }
				
	          else  
		           continue;

}

*/
 
		 int N0=181;

		 Vec_Com Es_u;
		    Es_u.assign(N0);
		 Vec_dou D_Es_u;
		    D_Es_u.assign(N0);
    	 Vec_dou DB_Es_u;
		    DB_Es_u.assign(N0);
		
		 Vec_Com Es_d;
		    Es_d.assign(N0);	
	     Vec_dou D_Es_d;
		    D_Es_d.assign(N0);
    	 Vec_dou DB_Es_d;
		    DB_Es_d.assign(N0);

		double r_p_q=0.;
		 
		Vec_dou Rp_q;
		       Rp_q.assign(2*N0);




		 
		 Vec_Com Es1;
		    Es1.assign(N0);	
		 Vec_dou Sita_D;
		    Sita_D.assign(N0);



		double KG = 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K*g_p*cos(i_an),2.0));
        double KG1= 1-(1.+2.*pow(tan(i_an),2.0))/(2.*pow(K1*g_p*cos(i_an),2.0));
	    

		 for (i=0;i<N0;i++)
		 {
            Sita_D[i] = i*PI/180-PI/2.;
		 }


/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {				     
			    Es[i] = Es[i] + dx_Suf*(J*K*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*(Jf_1[j]+Jb_1[j])-(Jf_2[j]+Jb_2[j]))*exp(-1.*J*K*(x[j]*sin(Sita_D[i])+z[j]*cos(Sita_D[i])));							   				
			 }
			 
            D_Es[i]=pow(abs(Es[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG);
            DB_Es[i]=10.*log10(D_Es[i]);

             
		 }*/
		 

/*		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {	
				 for (l=0;l<n1;l++)
				 {
					 Sum_J = Sum_J + Jf_2[j*n1+l]+Jb_2[j*n1+l];
				 }
				     Sum_J = Sum_J/double(n1);
				 
				 Es_u[i] = Es_u[i] + dx_Suf_N*(J*K*(df_x_N[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x_N[j],2.0)))*(Jf_1[j]+Jb_1[j])-Sum_J)*exp(-1.*J*K*(x_N[j]*sin(Sita_D[i])+z_N[j]*cos(Sita_D[i])));							   				
			     Sum_J = zero;
			 }
			     
			 
            D_Es_u[i]=pow(abs(Es_u[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG)/abs(rou_0);
            DB_Es_u[i]=10.*log10(D_Es_u[i]);

            r_p_q = r_p_q+D_Es_u[i];
		 }

 	    for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N_dg;j++)
			 {	
				j_n1 = int(j/n1);

			    Es_d[i] = Es_d[i] + dx_Suf_Ndg*(J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*(Jf_1[j_n1]+Jb_1[j_n1])-rou_0*(Jf_2[j]+Jb_2[j]))*exp(-1.*J*K1*(x[j]*sin(Sita_D[i])-z[j]*cos(Sita_D[i])));							   				
			       
			 }
			 
               D_Es_d[i]=pow(abs(Es_d[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG1)/abs(rou_0);
               DB_Es_d[i]=10.*log10(D_Es_d[i]);			

            r_p_q = r_p_q+D_Es_d[i]; 
		 }

		 for (i=0;i<N0;i++)
		 {
			 Rp_q[i] = DB_Es_u[N0-i-1];
			 Rp_q[i+N0] = DB_Es_d[i];
			
		 }

		 cout<<"r_p_q="<<r_p_q<<endl;*/

		 Complex Sum_J(0.,0.);
		 int j_n1;
		 for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N;j++)
			 {	
				 for (l=0;l<n1;l++)
				 {
					 Sum_J = Sum_J + Jf_Ndg[j*n1+l];
				 }
				     Sum_J = Sum_J/double(n1);
				 
				 Es_u[i] = Es_u[i] + dx_Suf_N*(J*K*(df_x_N[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x_N[j],2.0)))*(Jf_N[j])-Sum_J)*exp(-1.*J*K*(x_N[j]*sin(Sita_D[i])+z_N[j]*cos(Sita_D[i])));							   				
			     Sum_J = zero;
			 }
			     
			 
            D_Es_u[i]=pow(abs(Es_u[i]),2.0)/(8.*PI*K*sqrt(PI/2.)*g_p*cos(i_an)*KG)/abs(rou_0);
            DB_Es_u[i]=10.*log10(D_Es_u[i]);

            r_p_q = r_p_q+D_Es_u[i];
		 }

 	    for (i=0;i<N0;i++)
		 {
			 for (j=0;j<N_dg;j++)
			 {	
				j_n1 = int(j/n1);

			    Es_d[i] = Es_d[i] + dx_Suf_Ndg*(J*K1*(df_x[j]*sin(Sita_D[i])-cos(Sita_D[i]))/(sqrt(1.+pow(df_x[j],2.0)))*(Jf_N[j_n1])-rou_0*Jf_Ndg[j])*exp(-1.*J*K1*(x[j]*sin(Sita_D[i])-z[j]*cos(Sita_D[i])));							   				
			       
			 }
			 
               D_Es_d[i]=pow(abs(Es_d[i]),2.0)/(8.*PI*K1*sqrt(PI/2.)*g_p*cos(i_an)*KG1)/abs(rou_0);
               DB_Es_d[i]=10.*log10(D_Es_d[i]);			

            r_p_q = r_p_q+D_Es_d[i]; 
		 }

		 for (i=0;i<N0;i++)
		 {
			 Rp_q[i] = DB_Es_u[N0-i-1];
			 Rp_q[i+N0] = DB_Es_d[i];
			
		 }

		 cout<<"r_p_q="<<r_p_q<<endl;


 /*  ofstream out_ZT_a;
           out_ZT_a.open("ZT_a.txt");
		   for (i=0;i<N;i++)
		   {
			   for (j=0;j<N;j++)
			   {
				   out_ZT_a<<ZT_a[i*N+j]<<"   ";
				   if (i%N==0)
				   {
					   out_ZT_a<<endl;
				   }
			   }
		   }
   ofstream out_ZT_c;
           out_ZT_c.open("ZT_c.txt");
		   for (i=0;i<N_dg;i++)
		   {
			   for (j=0;j<N_dg;j++)
			   {
				   out_ZT_c<<ZT_c[i*N_dg+j]<<"   ";
				   if (i%N==0)
				   {
					   out_ZT_c<<endl;
				   }
			   }
		   }   
	ofstream out_ZT_d;
           out_ZT_d.open("ZT_d.txt");
		   for (i=0;i<N_dg;i++)
		   {
			   for (j=0;j<N_dg;j++)
			   {
				   out_ZT_d<<ZT_d[i*N_dg+j]<<"   ";
				   if (i%N==0)
				   {
					   out_ZT_d<<endl;
				   }
			   }
		   }  


		   ofstream out_J1;
		            out_J1.open("PBTG_J1.txt");
						for (i=0;i<N;i++)
						{
							out_J1<<x_N[i]/wave_len<<"   "<<abs(Jf_1[i])<<"   "<<abs(Jb_1[i])<<"   "<<abs(Jf_1[i]+Jb_1[i])<<endl;
						}

           ofstream out_J2;
		            out_J2.open("PBTG_J2.txt");
						for (i=0;i<N_dg;i++)
						{
							out_J2<<x[i]/wave_len<<"   "<<abs(Jf_2[i])<<"     "<<abs(Jb_2[i])<<"    "<<abs(Jf_2[i]+Jb_2[i])<<endl;
						}

		   ofstream out_FN_1;
		            out_FN_1.open("PBTG_FN_1.txt");
						for (i=0;i<2*Q+1;i++)
						{
							out_FN_1<<i-Q<<"   "<<abs(FN_f_1[i])<<endl;
						}

	      ofstream out_FN_2;
		            out_FN_2.open("PBTG_FN_2.txt");
						for (i=0;i<2*Q+1;i++)
						{
							out_FN_2<<i-Q<<"   "<<abs(FN_f_2[i])<<endl;

						}
		  ofstream out_Fai;
			         out_Fai.open("PBTG_Fai.txt");
						 for (i=0;i<2*Q+1;i++)
						 {
							 out_Fai<<i-Q<<"   "<<abs(Fai[i])<<endl;
						 }

		  ofstream out_w_Sita;
			         out_w_Sita.open("PBTG_w_Sita.txt");
						 for (i=0;i<2*Q+1;i++)
						 {
							 out_w_Sita<<i<<"   "<<w_Sita[i]<<endl;
						 }*/
		  ofstream out_Es_u;
	                out_Es_u.open("PBTG_Es_u.txt");
			        for (i=0;i<N0;i++)
					{
                        out_Es_u<<i-90<<"   "<<abs(Es_u[i])<<"   "<<D_Es_u[i]<<"   "<<DB_Es_u[i]<<endl;
					}
	
          ofstream out_Es_d;
	               out_Es_d.open("PBTGL_Es_d.txt");
			       for (i=0;i<N0;i++)
				   { 
                      out_Es_d<<i-90<<"   "<<abs(Es_d[i])<<"   "<<D_Es_d[i]<<"   "<<DB_Es_d[i]<<endl;
				   } 

	      ofstream out_Rp_q;
	               out_Rp_q.open("PBTG_Rp_q.txt");
			      for (i=0;i<2*N0;i++)
				  {
				        out_Rp_q<<i<<"      "<<Rp_q[i]<<endl;
				  }

          ofstream out_J_Ndg;
		           out_J_Ndg.open("PBTG_J_Ndg.txt");
                     for (i=0;i<N_dg;i++)
                     {
						 out_J_Ndg<<x[i]/wave_len<<"   "<<abs(Jf_Ndg[i])<<endl;
                     }
          ofstream out_J_N;
		           out_J_N.open("PBTG_J_N.txt");
                     for (i=0;i<N;i++)
                     {
						 out_J_N<<x[i]/wave_len<<"   "<<abs(Jf_N[i])<<endl;
                     } 




 

	time( &ltime2 );
	cout<<"The end time is "<<ctime(&ltime2 )<<endl;
	outtime<<"The end time is "<<ctime(&ltime2 )<<endl;
	ltime = ltime2 - ltime1;
	cout<<"The total time is "<< ltime <<"seconds"<<endl;
	outtime<<"The total time is "<<ltime<<"seconds"<<endl;

}

void FBM_SAA::Gauss_Layered_Rough_surface_FBM()
{
//////////////////2008.11.13///////////////////////
//////////////////////////////////////////////////////////////////////////
/////        如果媒质介电常数虚数部分过大则TM波收敛性更好           //////
//////////////////////////////////////////////////////////////////////////

///////////////////姬伟杰/////////////////////////
	
	ofstream outtime("total_time.txt",ios::trunc);
	time_t ltime1,ltime2,ltime;
	time( &ltime1 );
	cout<<"The starting time is "<<ctime( &ltime1 )<<endl;
	//一维导体粗糙面
	double wave_len=3.0e8/3e9;      //wave length
	double h=0.1*wave_len;            //rms height 
	double lc=1.*wave_len;           //correlation length 
	double i_an=30.0;                 //incident angle in degree 
	
	//锥形波的宽度参数g的选择采用 叶红霞 博士论文的结果

    double g_p =4.*wave_len;
	
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
	double ee=2.71828183;

	double d=5.*wave_len;

	double GAMA = 0.577216;
	double rm;
	double dl;
	double Delta;


	Complex NR(0.,0.);
	
	Complex er2(2.5,0.05);
	Complex er3(250,0);

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

//	cout<<"dx_Suf="<<dx_Suf<<endl;
	i_an = i_an*PI/180.;

	K1 = 2.0*PI/wave_len;
	K2 = K1*sqrt(er2);
	K3 = K1*sqrt(er3);
/*    for (i=0;i<N;i++)
    {
		z[i] = 0.;
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

//	double Yita = 120.*PI;

	double xsz,ww,kg,fac;
    for(i=0;i<N;i++)
	{
		xsz=x[i]*sin(i_an)-z[i]*cos(i_an);
		kg=pow(K1*g_p*cos(i_an),2.0);
		fac=pow((x[i]+z[i]*tan(i_an))/g_p,2.0);
		ww=(2.0*fac-1.0)/kg;
		Ei[i]=exp(J*K1*xsz*(1.0+ww))*exp(-fac);	
	}

	cout<<" Filling the impedance matrix......------>> "<<endl;

//////////////////////////////////////1_1_1///////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
//			dl = sqrt(1.+pow(df_1[j],2.0));
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
//				cout<<abs(A_2_1_1[i*N+j])<<endl;
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
	}




	cout<<"The matrix has been filled."<<endl;


//	Bi_CG(ZT,EE,Jf);                     //The MOM method


CLinearEquations linear;

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

		


			linear.GAUSS(A,X,V,4);

		

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

			
				
			linear.GAUSS(A,X,V,4);


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
	  
	  cout<<"The error is: Delta ="<<Delta<<endl;
			  
	  if(Delta<=0.01)
	  {
	    cout<<"The total number of the step is :"<<l+1<<"."<<"and the final error is :"<<Delta<<endl;
		  break;
	  }
				
	     else  
		     continue;

}

    	 int N0=181;
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
			
	double max_z;
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
//					dl = dx_Suf*sqrt(1+pow(df_x[j],2.0));
					NR = (-1.*df_x[j]*(x_D[i]-x[j])+(z_D[i]-z[j]))/(rm);
					Es[i] = Es[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(W1_f[j]+W1_b[j])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(U1_f[j]+U1_b[j]);
				}

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
					Es1[i] = Es1[i]+dx_Suf*J*K1/4.*NR*Hankel(1,1,K1*rm)*(Jf[j+N])-dx_Suf*J/4.*Hankel(0,1,K1*rm)*(Jf[j]);
				}

				 D_Es1[i]=R_D*pow(abs(Es1[i]),2.0)/(sqrt(PI/2.)*g_p*cos(i_an)*KG);
                 DB_Es1[i]=10.*log10(D_Es1[i]);
		
			}


		 
    

     

     


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
	         out_U1.open("Layered_U1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U1<<x[i]/wave_len<<"   "<<abs(U1_f[i])<<"   "<<abs(U1_b[i])<<"   "<<abs(U1_f[i]+U1_b[i])<<"   "<<abs(Jf[i])<<endl;
			 }
   
    ofstream out_W1;
	         out_W1.open("Layered_W1.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W1<<x[i]/wave_len<<"   "<<abs(W1_f[i])<<"   "<<abs(W1_b[i])<<"   "<<abs(W1_b[i]+W1_f[i])<<"   "<<abs(Jf[i+N])<<endl;
			 }

	ofstream out_U2;
	         out_U2.open("Layered_U2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_U2<<x[i]/wave_len<<"   "<<abs(U2_f[i])<<"   "<<abs(U2_b[i])<<"   "<<abs(U2_f[i]+U2_b[i])<<"   "<<abs(Jf[i+2*N])<<endl;
			 }
	ofstream out_W2;
	         out_W2.open("Layered_W2.txt");
			 for (i=0;i<N;i++)
			 {
				 out_W2<<x[i]/wave_len<<"   "<<abs(W2_f[i])<<"   "<<abs(W2_b[i])<<"   "<<abs(W2_b[i]+W2_f[i])<<"   "<<abs(Jf[i+3*N])<<endl;
			 }
	ofstream out_Es;
	         out_Es.open("Layered_Es.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es<<i-90<<"   "<<abs(Es[i])<<"   "<<D_Es[i]<<"   "<<DB_Es[i]<<endl;
			 }
	ofstream out_Es1;
	         out_Es1.open("Layered_Es1.txt");
			 for (i=0;i<N0;i++)
			 {
				 out_Es1<<i-90<<"   "<<abs(Es1[i])<<"   "<<D_Es1[i]<<"   "<<DB_Es1[i]<<endl;
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

	

	time( &ltime2 );
	cout<<"The end time is "<<ctime(&ltime2 )<<endl;
	outtime<<"The end time is "<<ctime(&ltime2 )<<endl;
	ltime = ltime2 - ltime1;
	cout<<"The total time is "<< ltime <<"seconds"<<endl;
	outtime<<"The total time is "<<ltime<<"seconds"<<endl;

}






