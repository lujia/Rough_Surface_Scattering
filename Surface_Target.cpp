//#include "post_computer.h"
//#include "target.h"
//#include "LinearEquations.h"
//#include "Surface.h"

//#include "Surf_Targ.h"
#include "FBM_NSA.h"
//#include "Lyaered_FBM_SAA.h"

ofstream outtime("total_time.txt",ios::trunc);

void main()
{
    ofstream outtime("total_time.txt",ios::trunc);
	time_t ltime1,ltime2,ltime;
	time( &ltime1 );
	cout<<"The starting time is "<<ctime( &ltime1 )<<endl;

/*
	CylinderQuestion oo;
	oo.Cy_ResultofJ();

*/
// 	  Surf_Targ_Question o;
      FBM_SAA ooo;
 //   Lyaered_FBM_SAA oooo;


//	o.Gauss_Surface_Cylinder_Target_J();

//	o.Gauss_Surface_Cylinder_Target_Jsd1_Jsd2();

//	o.Gauss_Surface_Cylinder_Target_Muti_Target();

//	o.Gauss_Surface_medium_Cylinder_Target_J();
	
//	o.Gauss_Del_Surface_Con_Cylinder__Above_Target();

//	o.Gauss_Del_Surface_Del_Cylinder__Above_Target();

//  o.Gauss_Surface_Below_Muti_Target();

//	o.Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Es1();

//	o.Gauss_Del_Surface_Conduct_Cylinder_Buried_Target_Esd();

	ooo.Gauss_PEC_Surface_FBM();

// 	ooo.Gauss_PEC_Surface_FBM_SAA();

//	 	ooo.Gauss_DIL_Surface_FBM();

 //  ooo.Gauss_DIL_Surface_FBM_SAA();

// ooo.Gauss_DIL_Surface_PBTC_FBM_SAA();
    
//	oooo.Gauss_Layered_Rough_surface_FBM();

// oooo.Gauss_Layered_Rough_surface_FBM_Bried_Target();

//	oooo.Gauss_Layered_Rough_surface_FBM_SAA();

//  oooo.Gauss_Layered_Rough_surface_FBM_Above_Target();
	
//  oooo.Gauss_Layered_Rough_surface_FBM_Bried_conduct_Target();	
	
	
	
	
	
	
	
	time( &ltime2 );
	cout<<"The end time is "<<ctime( &ltime2 )<<endl;
	outtime<<"The end time is "<<ctime( &ltime2 )<<endl;
	ltime = ltime2 - ltime1;
	cout<<" The total time is "<< ltime <<" seconds"<<endl;
	outtime<<" The total time is "<<ltime<<" seconds"<<endl;	   
}