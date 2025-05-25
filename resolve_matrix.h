#include "Mathe_Functions.h"

ofstream out_resual("Step_Resual.txt",ios::trunc);

//该程序解带状矩阵，且所有的带状矩阵元素都存储了，不考虑非零元素

class resolve_matrix: public Mathe_Functions
{
public:
	resolve_matrix();
	~resolve_matrix();
	
	//////////////////////////////////////////////////////
	//FBM迭代解法
	void resolve_FBM(Vec_Com &zm,Vec_Com &vm, Vec_Com &Jm, int n);
    void resolve_PEC_FBM_SAA(Vec_Com &zm,Vec_Com &vm, Vec_Com &Jm,double  max_z_i, double min_z_i,double dx_suf_i,double K_i,Vec_dou z);
    void resolve_DIL_FBM_SAA(Vec_Com &zm_a,Vec_Com &zm_b,Vec_Com &zm_c,Vec_Com &zm_d,Vec_Com &vm, Vec_Com &Jm,double  max_z_i, double min_z_i,double dx_suf_i,double K_i,double K1_i,Complex rou_i,Vec_dou z,Vec_dou df_x_i);
    ////////////////////Maybe there are also other questions for the PBTG.The backword current is not so well convergence////////////////////////////////////////////////
	void resolve_DIL_PBTG_FBM_SAA(Vec_Com &zm_a,Vec_Com &zm_b,Vec_Com &zm_c,Vec_Com &zm_d,Vec_Com &vm, Vec_Com &Jm_Ndg,Vec_Com &Jm_N, double  max_z_i, double min_z_i,double dx_suf_Ndg_i,double K_i,double K1_i,Complex rou_i,Vec_dou &z_Ndg_i,Vec_dou &z_N_i,Vec_dou &df_x_Ndg_i,Vec_dou &df_x_N_i);

	void Bi_CGSTAB(Vec_Com &zs,Vec_Com &vmn,Vec_Com &rex);
//	void Multiply_Full_Matrix_Vector(char a, Vec_Com &zm, Vec_Com &vec, Vec_Com &result );
    void Bi_CG(Vec_Com &zs,Vec_Com &vmn,Vec_Com &rex);   
	
	
	
	
    int   GAUSS(Vec_Com &A,Vec_Com &x,Vec_Com &b,int n);	
	
	
	//带状矩阵解法
	void resolve_bm(Vec_Com &zm,Vec_Com &vm,int n,int bw);
	//高斯消元解法
	void resolve_m(Vec_Com &a,Vec_Com &b,Vec_Com &re);
	void resolve_m(Vec_Com &a,Vec_Com &b,int m);
	//CGM解带状矩阵
	void resolve_cgm(Complex *am,Complex *bm,int n,int bw);
	//CGM结合PBTG解块带状矩阵
	void resolve_pbtg_cgm(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
    //CGM结合PBTG解块带状矩阵
	void CG(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//胡家赣CGM结合PBTG解块带状矩阵
	void resolve_pbtg_HJG_cgm(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *vm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//ILUCGM结合BMIA\CAG解带状矩阵
	void resolve_bmia_cag_ilucg(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *vm,int ndg,int bw,double ebs);
	//LSQR结合PBTG解块带状矩阵,效果一般
	void LSQR(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//CGLS结合PBTG解块带状矩阵，效果较好
	void CGLS(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//GMRES(k)，其中k仍然没有较好地确定
	Complex *GMRES_K(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw,double tol);
	//GMRES(k)中：解上三角矩阵，矩阵是按列存储的
	Complex *Resolve_Tri_Matrix_GMRES(Vec_Com &am,Vec_Com &bm,int m);
	//GMRES(k)中：(m,n)的矩阵与(n,1)的矩阵相乘，(m,n)的矩阵是按列存储的
	Complex *multiply_mn_GMRES(Vec_Com &am,Complex *bm,int m,int n);
	//Flexible GMRES(k)
	void F_GMRES_K(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw,double eps_m,double eps_p);
	//预处理的双共轭梯度法：不收敛
	void Bi_CG_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//CGS:P104
	void CGS(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//Bi-CGSTAB algorithm with preconditioning :P138
	void BI_CGSTAB_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//The GPBi-CG algorithm without preconditioning :P143
	void GPBI_CG(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	//The GPBi-CG algorithm with preconditioning :P143
	void GPBI_CG_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
	
	//矩阵与列向量相乘考虑了矩阵的带状性
	//具有四个子带状矩阵的矩阵与列向量相乘
	Complex *Multi_A_R(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *VR,int ndg,int bw);
	//具有四个子带状矩阵的矩阵的转置与列向量相乘，其中还有插值
	Complex *Multi_Tran_M_V(Complex *asmnt,Complex *bsmnt,Complex *a1mnt,Complex *b1mnt,Complex *aabmnt,Complex *babmnt,Complex *VR,Vec_Dou &CX,Vec_Dou &DX,int fbw,int lbw,int dbw,int ndg,int ncg);
	//具有四个子带状矩阵的矩阵与列向量相乘，其中还有插值
	Complex *Multi_M_V(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *aabmn,Complex *babmn,Complex *VR,Vec_Dou &CX,Vec_Dou &DX,int fbw,int lbw,int dbw,int ndg,int ncg);
	//压缩存储的带状矩阵与列向量相乘
	Complex *BandMatrix_multiply(Complex *a,Complex *b,int n,int bw);
	//PBTG中，粗格中的矩阵向量积
	Complex *Matrix_multiply(Complex *a,Complex *b,int fbw,int dbw,int n,int nn1);
	//压缩存储的带状矩阵的共轭转置
	Complex *Transpose_BandMatrix(Complex *am,int n,int bm);
	////////////////////不完全 Cholesky 分解//////////////////
	Complex *Creat_Diag(Complex *asmn,Complex *a1mn,Complex *b1mn,int ndg,int lbw,int fbw);//对角阵
	Complex *Creat_triang(Complex *M,Complex *Dm,int ndg,int bw);//生成三角阵U(Dm只为对角阵的一半)
	Complex *Creat_Inv_Diag(Complex *Dm,int ndg);//对角阵的求逆且开平方
	Complex *Multi_HalfBM_Diag(Complex *Half_BM,Complex *Dm_U_D,int ndg,int bw);//带状矩阵的下三角阵与对角阵相乘
	Complex *Multi_BM_Diag(Complex *BM,Complex *Dm_U_D,int ndg,int bw);//带状矩阵与对角阵相乘
	Complex *Tria_Decom(Complex *TM,int ndg,int bw);//带状三角阵求逆，生成三角阵
	Complex *Multi_BM_TM(Complex *BM,Complex *TM,int ndg,int bw);//带状矩阵与下三角阵相乘
	Complex *Multi_TM_LM(Complex *TM,Complex *LM,int ndg,int bw);//下三角阵与上三角阵为带状结构的矩阵相乘
	//乔列斯基分解的下三角矩阵的逆与列向量相乘，得到列向量
	Complex *Multi_CN_Vector(Complex *CN_Left_up,Complex *CN_Light_Down,Complex *CN_Left_Down,Complex *Vector,int ndg,int bw);
	//乔列斯基分解的下三角矩阵的转置与列向量相乘，得到列向量
	Complex *Multi_C_Trans_Vector(Complex *C_Left_up,Complex *C_Light_Down,Complex *C_Left_Down,Complex *Vector,int ndg,int bw);
	//乔列斯基分解的下三角矩阵的逆的转置与列向量相乘，得到列向量
	Complex *Multi_CN_Trans_Vector(Complex *CN_Left_up,Complex *CN_Light_Down,Complex *CN_Left_Down,Complex *Vector,int ndg,int bw);
	Complex *Multi_UALR(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *CN_L_U,Complex *CN_Li_D,Complex *CN_L_D,Complex *VR,int ndg,int bw);

	////////////////预处理矩阵////////////////////////
	Complex *Preconditoning_Matrix_Diag_Inverse(Complex *asmn,Complex *b1mn,int ndg,int lbw,int fbw);//对角阵的逆，预处理矩阵为对角阵
	Complex *Resolve_Cholesky_LU(Complex *C_Left_Up,Complex *C_Left_Down,Complex *C_Light_Down,Complex *Bm,int ndg,int lbw,int fbw);
	//预处理矩阵为对角阵的共轭梯度法
	void Preconditioning_CGM(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw);
//	Complex *Preconditioning_ICholesky();
	
	//在 GPBi-CG 中应用预处理矩阵
	Complex *Multi_DiagM_BandM(Complex *DiagM,Complex *BandM,int m,int bw);
	Complex *Multi_DiagM_ABM(Complex *DiagM,Complex *ABM,int md,int mc,int rd_bw,int rf_bw);

	//二维粗糙面的方程解法
	//二维标量波
	void Two_D_Cond_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs,Vec_Com &zfs,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &rex);
	void Two_D_Cond_CGM(Vec_Com &zs,Vec_Com &zfs,Vec_Com zs_inver,Vec_Com zfs_inver,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &rex);
	void Two_D_Cond_CGS_Pre(Vec_Com &zs,Vec_Com &zfs,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &resx);
	//对称带状矩阵向量积
	void Multi_Zs_Vec(Vec_Com &Zs,Vec_Com &vec,int RD_num,int nxy,Vec_Com &result);
	//二维标量波入射导体粗糙面弱相关矩阵向量积
	void Two_D_Conductor_Scalar_Multi_Zw_Vec(int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re);
    Complex Para_Two_D_Conductor_Scalar(int Bool,double R,Complex k);
	Complex G_0_fs_R(double R,Complex k);

	//二维矢量波
	void Two_D_Cond_Vector_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs_lu,vector <Zmn_Cell> &zs_ru,vector <Zmn_Cell> &zs_ld,vector <Zmn_Cell> &zs_rd,Vec_Com &zfs_lu,Vec_Com &zfs_ru,Vec_Com &zfs_ld,Vec_Com &zfs_rd,Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d,int nxy,int rd_num,double eeps,Vec_Com &rex);
   	void Two_D_Cond_Vector_Bi_CGSTAB_Pre(Vec_Com &zs_lu,Vec_Com &zs_ru,Vec_Com &zs_ld,Vec_Com &zs_rd,Vec_Com &zfs_lu,Vec_Com &zfs_ru,Vec_Com &zfs_ld,Vec_Com &zfs_rd,Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d,int nxy,int rd_num,double eeps,Vec_Com &rex);
	void Two_D_Cond_Vector_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs_lu,vector <Zmn_Cell> &zs_ru,vector <Zmn_Cell> &zs_ld,vector <Zmn_Cell> &zs_rd, Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d, Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn, Vec_Dou &dfymn,Complex kw,double dxy,int nxy,int rd_num,double eeps, Vec_Com &rex);

	//二维标量波入射近场相关矩阵向量积
	void NearMultiplyVector(char a, vector <Zmn_Cell> &zs, vector<Complex> &vec, vector<Complex> &result);
	//二维矢量波入射场相关矩阵向量积
	void Two_D_Vector_NearMultiplyVector(char a, vector <Zmn_Cell> &zs, vector<Complex> &vec, vector<Complex> &result);

	//二维矢量波入射弱矩阵向量积
	void Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re);
	void Two_D_Conductor_Vector_Multi_Zw_Vec_Child(int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re);
	Complex Para_Two_D_Conductor_Vector(int Bool,double R,Complex k);
    //二维矢量波入射平面矩阵向量积的计算：相当于零阶泰勒级数展开，运用2-D FFT计算

	//满阵与列向量相乘
	void Multiply_Full_Matrix_Vector(char a, Vec_Com &zm, Vec_Com &vec, Vec_Com &result );
	
	//二维介质粗糙面电磁散射仿真计算的弱矩阵向量积
	Complex Para_Two_D_Dielectric_G_a(int Bool,double R,Complex k);
	Complex Para_Two_D_Dielectric_g_b(int Bool,double R,Complex k);
	void Two_D_Diel_Multi_Zw_Vec_g_b_Child(char a, int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex kk,double dxy,int nxy,int rd_num,Vec_Com &re);
	void Two_D_Diel_Multi_Zw_Vec_G_a_Child(char a, int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex kk,double dxy,int nxy,int rd_num,Vec_Com &re);
	void Two_D_Diel_Multi_Zw_Vec(int Para_G_a_Num,int Para_g_b_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6]);
    void Two_D_Diel_Multi_Zw_Vector(int Para_G_a_Num,int Para_g_b_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6]);

	void Multi_vec_vec(Vec_Com &va, Vec_Com &vb, Vec_Com &re );
	void Multi_vec_num(Vec_Com &va, Complex nb, Vec_Com &re );

	void Two_D_Diel_Bi_CGSTAB_Pre(Vec_Zmn_Cell zs[6][6],Vec_Com zs_diag[6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,int nxy,int rd_num,double eeps,Vec_Com rex[6]);
    //不收敛
	void Two_D_Diel_GPBI_CG(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,int nxy,int rd_num,double eeps,Vec_Com rex[6]);
	//收敛速度慢
	void Two_D_Diel_CGM(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,int nxy,int rd_num,double eeps,Vec_Com rex[6]);
    //收敛速度习快于CGM，但仍然较慢
	void Two_D_Diel_CGLS(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,int nxy,int rd_num,double eps,Vec_Com rex[6]);
	//
	void Two_D_Diel_CGM(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],int nxy,double eps,Vec_Com rex[6]);

	//一般矩阵向量积：当a为' '时为向量积，当a为其它字符时为转置向量积
	void Two_D_Diel_Matrix_Vector_Multiply(char a, Vec_Zmn_Cell &zm, Vec_Com &vm, Vec_Com &re);
	//强相关矩阵向量积和其转置向量积：当a为' '时为向量积，当a为其它字符时为转置向量积
	void Two_D_Diel_Zs_Vector_Multiply(char a, Vec_Zmn_Cell zs[6][6], Vec_Com vm[6], Vec_Com re[6]);
	//弱相关矩阵向量积和其转置向量积：当a为' '时为向量积，当a为其它字符时为转置向量积
    void Two_D_Diel_Multi_Zfs_Vec(char a, Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6]);

	//////////////////////////////////////////////////////////////////////////
	//ILU(p,t)：不完全LU分解
	void ILU_p_t( int n,Vec_Zmn_Cell zs[6][6],int p,double tau,Vec_Zmn_Cell *M_L,Vec_Zmn_Cell *M_U );
	void Test_ILU_p_t( int nn, Vec_Zmn_Cell *WW,int P,double tau );
	Complex Get_Zmn(int j,Vec_Zmn_Cell &a);

	void ILU_0( int n,Vec_Zmn_Cell *WW, Vec_Zmn_Cell *M_L, Vec_Zmn_Cell *M_U );
	void Resolve_ILU( Vec_Zmn_Cell *M_L, Vec_Zmn_Cell *M_U , Vec_Com &Bm, Vec_Com &result );

	//排序:对元素为Vec_Zmn_Cell类型的向量排序
	void sort_vec_zmn_cell(Vec_Zmn_Cell &vec_cell);
	void test_zmn_cell();

	Complex *rex;
	double eps;
	double criterion;

};

resolve_matrix::resolve_matrix()
{
	rex = NULL;
}
resolve_matrix::~resolve_matrix()
{
	delete [] rex;
}

//////////////////////////////////////////////////////////////////////////


int resolve_matrix::GAUSS(Vec_Com &A,Vec_Com &x,Vec_Com &b,int n)
 // double ar[],ai[],br[],bi[];
{
	double *ar,*ai,*br,*bi;
	ar=new double[n*n];  ai=new double[n*n];
	br=new double[n];    bi=new double[n];
	for(int m=0;m<n*n;m++)
	{
		ar[m]=A[m].real();
		ai[m]=A[m].imag();
	}
	for( m=0;m<n;m++)
	{
		br[m]=b[m].real();
		bi[m]=b[m].imag();
	}
	int *js,l,k,i,j,is,u,v;
    double p,q,s,d;
    //js=malloc(n*sizeof(int));
	js=new int[n];
    for (k=0;k<=n-2;k++)
      { d=0.0;
        for (i=k;i<=n-1;i++)
        for (j=k;j<=n-1;j++)
          { u=i*n+j;
            p=ar[u]*ar[u]+ai[u]*ai[u];
            if (p>d) {d=p;js[k]=j;is=i;}
          }
        if (d+1.0==1.0)
          { 
//			delete[js];
			cout<<"err**fail1"<<endl;
			//free(js); printf("err**fail\n");
            return(0);
          }
        if (is!=k)
          { for (j=k;j<=n-1;j++)
              { u=k*n+j; v=is*n+j;
                p=ar[u]; ar[u]=ar[v]; ar[v]=p;
                p=ai[u]; ai[u]=ai[v]; ai[v]=p;
              }
            p=br[k]; br[k]=br[is]; br[is]=p;
            p=bi[k]; bi[k]=bi[is]; bi[is]=p;
          }
        if (js[k]!=k)
          for (i=0;i<=n-1;i++)
            { u=i*n+k; v=i*n+js[k];
              p=ar[u]; ar[u]=ar[v]; ar[v]=p;
              p=ai[u]; ai[u]=ai[v]; ai[v]=p;
            }
        v=k*n+k;
        for (j=k+1;j<=n-1;j++)
          { u=k*n+j;
            p=ar[u]*ar[v]; q=-ai[u]*ai[v];
            s=(ar[v]-ai[v])*(ar[u]+ai[u]);
            ar[u]=(p-q)/d; ai[u]=(s-p-q)/d;
          }
        p=br[k]*ar[v]; q=-bi[k]*ai[v];
        s=(ar[v]-ai[v])*(br[k]+bi[k]);
        br[k]=(p-q)/d; bi[k]=(s-p-q)/d;
        for (i=k+1;i<=n-1;i++)
          { u=i*n+k;
            for (j=k+1;j<=n-1;j++)
              { v=k*n+j; l=i*n+j;
                p=ar[u]*ar[v]; q=ai[u]*ai[v];
                s=(ar[u]+ai[u])*(ar[v]+ai[v]);
                ar[l]=ar[l]-p+q;
                ai[l]=ai[l]-s+p+q;
              }
            p=ar[u]*br[k]; q=ai[u]*bi[k];
            s=(ar[u]+ai[u])*(br[k]+bi[k]);
            br[i]=br[i]-p+q; bi[i]=bi[i]-s+p+q;
          }
      }
    u=(n-1)*n+n-1;
    d=ar[u]*ar[u]+ai[u]*ai[u];
    if (d+1.0==1.0)
      { 
//		delete[js];
	    cout<<"err**fail2"<<endl;
		//free(js); printf("err**fail\n");
        return(0);
      }
    p=ar[u]*br[n-1]; q=-ai[u]*bi[n-1];
    s=(ar[u]-ai[u])*(br[n-1]+bi[n-1]);
    br[n-1]=(p-q)/d; bi[n-1]=(s-p-q)/d;
    for (i=n-2;i>=0;i--)
    for (j=i+1;j<=n-1;j++)
      { u=i*n+j;
        p=ar[u]*br[j]; q=ai[u]*bi[j];
        s=(ar[u]+ai[u])*(br[j]+bi[j]);
        br[i]=br[i]-p+q;
        bi[i]=bi[i]-s+p+q;
      }
    js[n-1]=n-1;
    for (k=n-1;k>=0;k--)
      if (js[k]!=k)
        { p=br[k]; br[k]=br[js[k]]; br[js[k]]=p;
          p=bi[k]; bi[k]=bi[js[k]]; bi[js[k]]=p;
        }
//	  delete[js];
    //free(js);
   for(m=0;m<n;m++)
   {
	   x[m]=br[m]+i_*bi[m];
   }
    return(1);
}
//////////////////////////////////////////////////////////////////////////
void resolve_matrix::Bi_CGSTAB(Vec_Com &zs,Vec_Com &vmn,Vec_Com &rex)
{
	int i;
	int Num = vmn.size();
	rex.assign(Num);
	Vec_Com vector_r,vector_p,vector_pt,vector_r0p,vector_s,vector_st,vector_t,vector_v;
	vector_r.assign(Num);
	vector_p.assign(Num);
	vector_pt.assign(Num);
	vector_r0p.assign(Num);
	vector_s.assign(Num);
	vector_st.assign(Num);
	vector_t.assign(Num);
	vector_v.assign(Num);
	
/*
	for ( i=0; i<zs.size(); i++ )
	{
		out_resual<<zs[i]<<endl;
	}*/
	
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	double norm_0=sqrt(L2_Norm(vmn));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(Num);
	for ( i=0; i<Num; i++ )
	{
		Diag_Matrix_Inverse[i] = 1. / zs[i*Num+i];
	}
	
	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r);
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<Num; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<Num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
        
		Multiply_Full_Matrix_Vector(' ',zs,vector_pt,vector_v);
		
		temp1=Multi_vt_v(vector_r0p,vector_v);
	/*	for (i=0; i<Num; i++)
		{
			out_resual<<zs[i]<<"  "<<vector_v[i]<<"  "<<vector_pt[i]<<endl;
		}*/

		alfa=rou/temp1;
		for (i=0; i<Num; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
		}

		Norm=sqrt(L2_Norm(vector_s));
	
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<Num; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<Num; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}

		Multiply_Full_Matrix_Vector(' ',zs,vector_st,vector_t);
	/*	for ( i=0; i<Num; i++ )
		{
			cout<<vector_st[i]<<endl;
		}*/
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for (i=0; i<Num; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < 1.e-4)
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	vector_r.clear();
	vector_p.clear();
	vector_pt.clear();
	vector_r0p.clear();
	vector_s.clear();
	vector_st.clear();
	vector_t.clear();
	vector_v.clear();
	Diag_Matrix_Inverse.clear();
}
/*
void resolve_matrix::Multiply_Full_Matrix_Vector(char a, Vec_Com &zm, Vec_Com &vec, Vec_Com &result )
{
	int i,j,Num;
	Num = vec.size();
	result.assign(Num);

	if(a == ' ')
	{
		for ( i=0; i<Num; i++ )
		{
			for ( j=0; j<Num; j++ )
			{
				result[i] += zm[i*Num+j]*vec[j];
			}
		}
	}
	else
	{
		for ( i=0; i<Num; i++ )
		{
			for ( j=0; j<Num; j++ )
			{
				result[i] += conj(zm[j*Num+i])*vec[j];
			}
		}
	}
}
*/
void resolve_matrix::Bi_CG(Vec_Com &zs,Vec_Com &vmn,Vec_Com &rex)  
{
//////////////////2008.10.25//////////////////////

///////////////////姬伟杰/////////////////////////
	
	int i;
	int Num = vmn.size();
	      rex.assign(Num);
	Vec_Com vector_r,vector_rt,vector_r0,vector_r0t,
		    vector_p,vector_pt,vector_temp;
	vector_r.assign(Num);
	vector_rt.assign(Num);
	vector_r0.assign(Num);
	vector_r0t.assign(Num);
	vector_p.assign(Num);
	vector_pt.assign(Num);
	vector_temp.assign(Num);
	
	Complex alfa,beta;
	Complex Norm;
	Complex norm_0 = Multi_vt_v(vmn,vmn);

	Assign_Vecor(vmn,vector_p);	
	Assign_Vecor(vmn,vector_pt);
	Assign_Vecor(vmn,vector_r0);
	Assign_Vecor(vmn,vector_r0t);
	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_rt);

	
    for (int k=0;k<Num;k++)
	{

		Multiply_Full_Matrix_Vector(' ',zs,vector_p,vector_temp);
		alfa = Multi_vt_v(vector_rt,vector_r) / Multi_vt_v(vector_pt,vector_temp);

		for (i=0; i<Num; i++)
		{
			rex[i] = rex[i] + alfa * vector_p[i];
			vector_r[i] -= alfa * vector_temp[i]; 
		}

		Multiply_Full_Matrix_Vector('*',zs,vector_pt,vector_temp);

		for (i=0; i<Num; i++)
		{
			vector_rt[i] -= conj(alfa) * vector_temp[i]; 
		}

		beta = Multi_vt_v(vector_rt,vector_r) / Multi_vt_v(vector_r0t,vector_r0);

		for (i=0; i<Num; i++)
		{
			vector_p[i] = vector_r[i] + beta * vector_p[i]; 
			vector_pt[i] = vector_rt[i] + conj(beta) * vector_pt[i]; 
//			out_resual<<vector_rt[i]<<"  "<<vector_r[i]<<endl;
		}		
		
		Norm = Multi_vt_v(vector_rt,vector_r) / norm_0;
		cout<<"Bi-CG steps "<<k+1<<",  res = "<<abs(Norm)<<endl;
		if ( abs(Norm) < 1.e-4)
		{
			break;
		}

		Assign_Vecor(vector_r,vector_r0);
		Assign_Vecor(vector_rt,vector_r0t);

	}
	cout<<"Total Bi-CG Steps = "<<k+1<<", to achieve "<<abs(Norm)<<endl;
	
}










void resolve_matrix::resolve_FBM(Vec_Com &zm,Vec_Com &vm, Vec_Com &Jm,int n)
{
//////////////////2008.10.20//////////////////////

///////////////////姬伟杰/////////////////////////

   cout<<"FBM法解粗糙面方程组......------>>"<<endl;
	
  int i,j,k;
   Complex zero(0.,0.);
   Complex Sum0 = zero;
   Complex Sum1 = zero;
   Complex SUM=zero;
   Vec_Com zf;

   Vec_Com Jf;
     Jf.assign(n);  
   Vec_Com Jf1;
     Jf1.assign(n);
  Vec_Com Jf2;
     Jf2.assign(n); 
   Vec_Com Jb1;
     Jb1.assign(n);
   Vec_Com Jb;
     Jb.assign(n);
  Vec_Com Vector_temp;
     Vector_temp.assign(n);

   
   double DELTA=0.;
   Complex max =zero;
   Complex max1=zero;
	  

	  for(k=0;k<20;k++)
	  {
	  
           cout<<"Runing the"<<k+1<<"th step of the FBM......------->>"<<endl;

///////////////计算对电磁波前向传播有贡献的部分感应电流///////////////////////

	        for (i=0;i<n;i++)
			{
		        for (j=0;j<i;j++)
				{
			       Sum0 = Sum0+zm[i*n+j]*(Jf[j]+Jb[j]);
			    }
			
	            Jf[i] = (vm[i]-Sum0)/zm[i*n+i];       
		        Sum0 = zero;	       		
			}

///////////////计算对电磁波后向传播有贡献的部分感应电流(其计算顺序为:从后向前计算，这是为了完全用更新的电流来计算)///////////////////////

	        for (i=n-1;i>=0;i--)
			{
		       for (j=i+1;j<n;j++)
			   {
			     Sum1 = Sum1+zm[i*n+j]*(Jb[j]+Jf[j]);                 
			   }
	            Jb[i] =  -1.*(Sum1)/zm[i*n+i];
				Sum1 = zero;	       
			}
			
///////////////////////收敛条件的判断/////////////////////////////////////////////////

		  
		    for (i=0;i<n;i++)
		    {
				Jf2[i] = Jf[i]-Jf1[i];
		    }


			  Multiply_Full_Matrix_Vector(' ',zm,Jf2,Vector_temp);
              Assign_Vecor(Jf,Jf1);
              DELTA = sqrt(L2_Norm(Vector_temp))/sqrt(L2_Norm(vm));
			

           cout<<" THE error of is "<<DELTA<<endl;
	       if(DELTA<0.01)
		   {
			   cout<<"The total steps of the FBM is"<<k+1<<". and the final error is "<<DELTA<<"."<<endl;
		    	  break; 
		   }
	        else  
		          continue;
			  					
	  }

//////////////////////输出总电流/////////////////////////////////////
	 
	  for(i=0;i<n;i++)
	  {
		  Jm[i] = Jf[i] + Jb[i] ;
	  
	  }
	     
}


void resolve_matrix::resolve_PEC_FBM_SAA(Vec_Com &zm,Vec_Com &vm, Vec_Com &Jm,double  max_z_i, double min_z_i,double dx_suf_i,double K_i,Vec_dou z)
{

//////////////////2008.11.3///////////////////////

///////////////////姬伟杰/////////////////////////
	
	int i,j,k,m;
	double Delta_0;
	
	double Delta;
	
	int N = vm.size();

	Complex JJ(0.,1.);
	Complex zero(0.,0.);
    Complex Fai_p;

	double L_s = (max_z_i - min_z_i)/4.;                                 //The strong interaction zone.
	double Bita_s = sqrt(10./(K_i*L_s));                                //The geometric optics lit region.
	double R_s = sqrt(pow(L_s,2.0)+pow(max_z_i-min_z_i,2.0));             
	double Fai_S_max = atan((max_z_i-min_z_i)/L_s);                       //The outermost possible saddle point.
	
	if (sqrt(K_i*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K_i*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}
	
	double Delta_Fai = sqrt(5./(K_i*R_s))/22.;                          
	int  Q = int(Bita_s/Delta_Fai)+5;

	int N_s = 1.+int(L_s/dx_suf_i);                                      //The number of the strong interaction point.                           

	cout<<"L_s = "<<L_s<<endl;
	cout<<"N_s = "<<N_s<<endl;
	Complex Sum,Sum1;
    


    Vec_Com FN_f;
           FN_f.assign(2*Q+1);
	Vec_Com FN_b;
		   FN_b.assign(2*Q+1);
    Vec_dou w_Sita;
	       w_Sita.assign(2*Q+1);
    

     Vec_Com Jf2;
	   Jf2.assign(N);

	 Vec_Com Jf1;
	   Jf1.assign(N);

	 Vec_Com Jf;
	   Jf.assign(N);



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
  

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);                //The window function.
         if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;		 
	}



    for (m=0;m<N;m++)
	{

	    cout<<"Runing the "<<m+1<<"th Step of FBM/SAA ......------>>"<<endl;

//////////////////////////////The strong point for the forward.////////////////////////////////////////////////////////

           
          for (i=0;i<N_s+1;i++)
		  {
	         for (j=0;j<i;j++)
			 {
                 Sum =Sum + (Jf[j]+Jb[j])*zm[i*N+j];
			 }
			  Jf[i] = (vm[i]-Sum)/zm[i*N+i];
 //           Ef_s[i] = Sum;
			  Sum = zero;
              	
		  } 
		  

////////////////////////////The weak point for the forward.//////////////////////////////////////////////

  			 
		   for (j=0;j<N_s+1;j++)
		   {
              Sum =Sum + (Jf[j]+Jb[j])*zm[(N_s+1)*N+(N_s+1)];
		   }

		   for (k=-Q;k<=Q;k++)
			 {   
			    Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
                FN_f[k+Q] =  (Jf[0]+Jb[0])*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[0]*sin(Fai_p));
  				      				      
			    Sum1 = Sum1 + exp(-1.*JJ*Delta_0)/PI*w_Sita[Q+k]*FN_f[k+Q]*exp(JJ*K_i*z[N_s+1]*sin(Fai_p))*Delta_Fai; 									
				  
			 }
			    Jf[N_s+1] = (vm[N_s+1]-Sum1-Sum)/zm[(N_s+1)*N+(N_s+1)];	
	            Ef_w[i] = Sum;
				Sum = zero;
				Sum1 = 0.;
				  



///////////////////////////////begin with the N_s+2 point///////////////////////////////////////////////////////
		 
       for (i=N_s+2;i<N;i++)
	   {
		     for (j=i-N_s;j<i;j++)
			 {
                Sum1 = Sum1 + (Jf[j]+Jb[j])*zm[i*N+j];       
			 }
//		        Ef_s[i] = Sum;
//			    Sum = zero;
			 for (k=-Q;k<=Q;k++)
			 {   
			    Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
                FN_f[k+Q] = FN_f[k+Q]*exp(JJ*K_i*dx_suf_i*cos(Fai_p)) + (Jf[i-N_s-1]+Jb[i-N_s-1])*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i-N_s-1]*sin(Fai_p));
  				      				      
			    Sum = Sum + exp(-1.*JJ*Delta_0)/PI*w_Sita[Q+k]*FN_f[k+Q]*exp(JJ*K_i*z[i]*sin(Fai_p))*Delta_Fai; 									
				  
			 }
			    Jf[i] = (vm[i]-Sum1-Sum)/zm[i*N+i];	
//	            Ef_w[i] = Sum;
				Sum = zero;
				Sum1 = zero;
					     				   
	   }
	
////////////////////////////////The Strong point for the behihd.////////////////////////////////////////////////////////////
								
	  for (i=N-1;i>=N-N_s-1;i--)
	   { 
	       for (j=i+1;j<N;j++)
		   {			          
            Sum = Sum + (Jf[j]+Jb[j])*zm[i*N+j];  
		   }
	  	 Jb[i] = -1.*Sum/zm[i*N+i];
//	     Eb_s[i] =Sum;	
		 Sum = zero;
	
	   } 

		  

////////////////////////////The weak for the behind./////////////////////////////////////////////////////////////////
				  
	  
	        for (j=N-N_s-2;j<N;j++)
			{			          
                Sum = Sum + (Jf[j]+Jb[j])*zm[(N-N_s-2)*N+(N-N_s-2)];  
			}

	        for (k=-Q;k<=Q;k++)
			{   
			     Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
                 
				 FN_b[k+Q] = (Jf[N-1]+Jb[N-1])*exp(JJ*K_i*(1.+N_s)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[N-1]*sin(Fai_p));					            				 
			     
				 Sum1 = Sum1 + exp(-1.*JJ*Delta_0)/PI*w_Sita[Q+k]*FN_b[k+Q]*exp(JJ*K_i*z[N-N_s-2]*sin(Fai_p))*Delta_Fai; 												
			}
			    Jb[N-N_s-2] = -1.*(Sum1+Sum)/zm[(N-N_s-2)*N+(N-N_s-2)];
//	            Eb_w[i] = Sum;
				Sum1 = zero;
				Sum=zero;
		

/////////////////////////////Start from the N-N_s-3 point./////////////////////////////////////////////


         for (i=N-N_s-3;i>=0;i--)
		  {
		     for (j=i+1;j<i+N_s+1;j++)
			 {	
                 Sum1 = Sum1 + (Jf[j]+Jb[j])*zm[i*N+j];  		   
			 }
//                Eb_s[i] = Sum;		  
//		          Sum = zero;
			for (k=-Q;k<=Q;k++)
			{   
			     Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
                 
				 FN_b[k+Q] = FN_b[k+Q]*exp(JJ*K_i*dx_suf_i*cos(Fai_p))+(Jf[i+N_s+1.]+Jb[i+N_s+1.])*exp(JJ*K_i*(1.+N_s)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i+N_s+1]*sin(Fai_p));					            				 
			     
				 Sum = Sum + exp(-1.*JJ*Delta_0)/PI*w_Sita[Q+k]*FN_b[k+Q]*exp(JJ*K_i*z[i]*sin(Fai_p))*Delta_Fai; 												
			}
			    Jb[i] = -1.*(Sum+Sum1)/zm[i*N+i];
//	            Eb_w[i] = Sum;
				Sum = zero;
				Sum1= zero;
			  
		 }

//////////////////////////判断收敛条件///////////////////////////////////////////////////////
	  
	      
 		    
		 
		      for (i=0;i<N;i++)
			  {
			     Jf2[i] = Jf[i]-Jf1[i];
			  }

              
			  Delta = sqrt(L2_Norm(Jf2))/sqrt(L2_Norm(Jf));
              Assign_Vecor(Jf,Jf1);
              cout<<"The error of the FBM_SAA is Delta ="<<Delta<<"."<<endl;
			  
			  if(Delta<=0.01)
			  {
	    	   cout<<"The total steps of the FBM_SAA is"<<m<<". and the final error is "<<Delta<<"."<<endl;
		        break;
			  }
				
	          else  
		        continue;
	}

/////////////////////输出最终结果/////////////////////////////////////////////////////


	for (i=0;i<N;i++)
	{
		Jm[i] = Jf[i] + Jb[i];
	}

}

void resolve_matrix::resolve_DIL_FBM_SAA(Vec_Com &zm_a,Vec_Com &zm_b,Vec_Com &zm_c,Vec_Com &zm_d,Vec_Com &vm, Vec_Com &Jm,double  max_z_i, double min_z_i,double dx_suf_i,double K_i,double K1_i,Complex rou_i,Vec_dou z,Vec_dou df_x_i)
{

//////////////////2008.11.8///////////////////////

///////////////////姬伟杰/////////////////////////	
	
	int i,j,k,m;
	double Delta_0;
	
	double Delta;
	
	int N = vm.size();

	Complex JJ(0.,1.);
	Complex zero(0.,0.);
    Complex Fai_p;

	double L_s = (max_z_i - min_z_i)/4.;                                 //The strong interaction zone.
	double Bita_s = sqrt(10./(K_i*L_s));                                //The geometric optics lit region.
	double R_s = sqrt(pow(L_s,2.0)+pow(max_z_i-min_z_i,2.0));             
	double Fai_S_max = atan((max_z_i-min_z_i)/L_s);                       //The outermost possible saddle point.
	
	if (sqrt(K_i*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K_i*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}
	
	double Delta_Fai = sqrt(5./(K_i*R_s))/22.;                          
	int  Q = int(Bita_s/Delta_Fai)+5;

	int N_s = 1+int(L_s/dx_suf_i);                                      //The number of the strong interaction point.                           

	Complex Sum0,Sum1,Sum2,Sum3;
    


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
		   
	
		   
	Vec_dou w_Sita;
	       w_Sita.assign(2*Q+1);
    
  

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);                //The window function.
         if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;		 
	}

  for (m=0;m<N;m++)
  {

	cout<<"Begin the "<<m+1<<"th step of the FBM_SAA ......------>>"<<endl;


////////////////////////////////The point which can not use only have the strong relative ,from 0 to N_s///////////////////////////////////////////////////////
	
	
	for (i=0;i<N_s+1;i++)
	{
		for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];
		}
	

		Jf_1[i] = ((vm[i]-Sum0)*zm_d[i*N+i]-(-1.*Sum1)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1)*zm_a[i*N+i]-(vm[i]-Sum0)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
	
	}

/////////////////////////////////The special point N_s+1//////////////////////////////////////////////////////////
 
	    i=N_s+1;

    	for (j=0;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];
		}	
	
	   for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
            FN_f_1[k+Q] = ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_i[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])+(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i-N_s-1]*sin(Fai_p));
            FN_f_2[k+Q] = ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_i[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])-rou_i*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K1_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K1_i*z[i-N_s-1]*sin(Fai_p));
			
			
			Sum2 = Sum2 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(JJ*K_i*z[N_s+1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(JJ*K1_i*z[N_s+1]*sin(Fai_p))*Delta_Fai; 												
		}
		Jf_1[i] = ((vm[i]-Sum0-Sum2)*zm_d[i*N+i]-(-1.*Sum1-Sum3)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1-Sum3)*zm_a[i*N+i]-(vm[i]-Sum0-Sum2)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
        
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;		

		
//////////////////////////////////The weak  point for the forward /////////////////////////////////////////////////////////

	for (i=N_s+2;i<N;i++)
	  {
		for (j=i-N_s-1;j<i;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];
		}
         
		for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
            FN_f_1[k+Q] = FN_f_1[k+Q]*exp(JJ*K_i*dx_suf_i*cos(Fai_p)) + ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_i[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])+(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i-N_s-1]*sin(Fai_p));
            FN_f_2[k+Q] = FN_f_2[k+Q]*exp(JJ*K1_i*dx_suf_i*cos(Fai_p)) + ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_i[i-N_s-1]*cos(Fai_p)))*(Jf_1[i-N_s-1]+Jb_1[i-N_s-1])-rou_i*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K1_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K1_i*z[i-N_s-1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(JJ*K_i*z[i]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(JJ*K1_i*z[i]*sin(Fai_p))*Delta_Fai; 									

		}

		Jf_1[i] = ((vm[i]-Sum0-Sum2)*zm_d[i*N+i]-(-1.*Sum1-Sum3)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jf_2[i] = ((-1.*Sum1-Sum3)*zm_a[i*N+i]-(vm[i]-Sum0-Sum2)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
        

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
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];

		}


    	Jb_1[i] = ((-1.*Sum0)*zm_d[i*N+i]-(-1.*Sum1)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1)*zm_a[i*N+i]-(-1.*Sum0)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
	
	}


////////////////////////////////The speczial point N-N_s-2///////////////////////////////////////////////////////////
		i=N-N_s-2;
	
	   for (j=i+1;j<N;j++)
		{
			Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];

		}

	   for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
            FN_b_1[k+Q] = ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_i[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])+(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i+N_s+1]*sin(Fai_p));
            FN_b_2[k+Q] = ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_i[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])-rou_i*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K1_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K1_i*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(JJ*K_i*z[N-N_s-2]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(JJ*K1_i*z[N-N_s-2]*sin(Fai_p))*Delta_Fai; 												
		}

     	Jb_1[i] = ((-1.*Sum0-Sum2)*zm_d[i*N+i]-(-1.*Sum1-Sum3)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1-Sum3)*zm_a[i*N+i]-(-1.*Sum0-Sum2)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;				
	
///////////////////////////////The weak point for the behind ////////////////////////////////////////////////////////////


	 for (i=N-N_s-3;i>=0;i--)
	 {
		for (j=i+1;j<i+N_s+1;j++)
		{
		    Sum0 = Sum0 + (Jf_1[j]+Jb_1[j])*zm_a[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_b[i*N+j];
			Sum1 = Sum1 + (Jf_1[j]+Jb_1[j])*zm_c[i*N+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N+j];	
		}

		for (k=-Q;k<Q;k++)
		{
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
					  
            FN_b_1[k+Q] = FN_b_1[k+Q]*exp(JJ*K_i*dx_suf_i*cos(Fai_p)) + ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_i[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])+(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K_i*z[i+N_s+1]*sin(Fai_p));
            FN_b_2[k+Q] = FN_b_2[k+Q]*exp(JJ*K1_i*dx_suf_i*cos(Fai_p)) + ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_i[i+N_s+1]*cos(Fai_p)))*(Jf_1[i+N_s+1]+Jb_1[i+N_s+1])-rou_i*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K1_i*(N_s+1.)*dx_suf_i*cos(Fai_p))*exp(-1.*JJ*K1_i*z[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(JJ*K_i*z[i]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_suf_i/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(JJ*K1_i*z[i]*sin(Fai_p))*Delta_Fai; 												
		}

     	Jb_1[i] = ((-1.*Sum0-Sum2)*zm_d[i*N+i]-(-1.*Sum1-Sum3)*zm_b[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Jb_2[i] = ((-1.*Sum1-Sum3)*zm_a[i*N+i]-(-1.*Sum0-Sum2)*zm_c[i*N+i])/(zm_a[i*N+i]*zm_d[i*N+i]-zm_b[i*N+i]*zm_c[i*N+i]);
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
	}

////////////////////////////////the convergence condition ////////////////////////////////////////////////////////////

		   
		  

 		     for (i=0;i<N;i++)
			 {
			    Jf2[i]   = (Jf_1[i]+Jb_1[i])-Jf1[i];
        	    Jf2[i+N] = (Jf_2[i]+Jb_2[i])-Jf1[i+N];
			    Jf3[i]   =  Jf_1[i]+Jb_1[i];
                Jf3[i+N] =  Jf_2[i]+Jb_2[i];
			 }

	          
			 Delta = sqrt(L2_Norm(Jf2))/sqrt(L2_Norm(Jf3));

             Assign_Vecor(Jf3,Jf1);


             cout<<"The error is: Delta ="<<Delta<<endl;
			  
			 if(Delta<=0.01)
			 {
			    cout<<"The total number of the step is :"<<m+1<<"."<<"and the final error is :"<<Delta<<endl;
		        break;
			 }
				
	          else  
		           continue;
	}

	
/////////////////////输出最终结果/////////////////////////////////////////////////////


	for (i=0;i<2*N;i++)
	{
		Jm[i] = Jf3[i] ;
	}



}


void resolve_matrix::resolve_DIL_PBTG_FBM_SAA(Vec_Com &zm_a,Vec_Com &zm_b,Vec_Com &zm_c,Vec_Com &zm_d,Vec_Com &vm, Vec_Com &Jm_Ndg,Vec_Com &Jm_N, double  max_z_i, double min_z_i,double dx_suf_Ndg_i,double K_i,double K1_i,Complex rou_i,Vec_dou &z_Ndg_i,Vec_dou &z_N_i,Vec_dou &df_x_Ndg_i,Vec_dou &df_x_N_i)
{
	int i,j,k,m,l;
	int i_n1,j_n1;
	int N_s_n1;
	int NN;
	

	
	double Delta_0;
	
	double Delta;
//	int n1 = sqrt(Epsle_d);
	int N_dg = Jm_Ndg.size();
    int n1=3;                                                         //Usually it can been const.
	int N = N_dg/n1;

//	Vec_dou z;
//	       z.assign(N_dg);
//	Vec_dou z_N;
//	        z_N.assign(N);
	
	double dx_Suf_Ndg = dx_suf_Ndg_i;
	double dx_Suf_N = n1*dx_suf_Ndg_i;


	Complex JJ(0.,1.);
	Complex zero(0.,0.);
    Complex Fai_p;

	double L_s = (max_z_i - min_z_i)/4.;                                 //The strong interaction zone.
	double R_s = sqrt(pow(L_s,2.0)+pow(max_z_i-min_z_i,2.0));  
	double Bita_s = sqrt(10./(K_i*R_s));                                //The geometric optics lit region.The L_s has been replaced by R_s.
	           
	double Fai_S_max = atan((max_z_i-min_z_i)/L_s);                       //The outermost possible saddle point.
	
	if (sqrt(K_i*R_s/20.)*Fai_S_max-1>1)
	{
		Delta_0 = atan(1./(sqrt(K_i*R_s/20.)*Fai_S_max-1.));     
	}
	else
	{
		Delta_0 = PI/4.;
	}
	
	double Delta_Fai = sqrt(5./(K_i*R_s))/22.;                          
	int  Q = int(Bita_s/Delta_Fai)+5;

	int N_s = int(L_s/dx_suf_Ndg_i);                                      //The number of the strong interaction point.                           

	Complex Sum0,Sum1,Sum2,Sum3;
    

    
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
		   
	
		   
	Vec_dou w_Sita;
	       w_Sita.assign(2*Q+1);


    Vec_Com FN_f_1;
           FN_f_1.assign(2*Q+1);
	Vec_Com FN_b_1;
		   FN_b_1.assign(2*Q+1);
    Vec_Com FN_f_2;
           FN_f_2.assign(2*Q+1);
	Vec_Com FN_b_2;
		   FN_b_2.assign(2*Q+1);
    
  

	for (k=-Q;k<=Q;k++)
	{   
	     Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);                //The window function.
         if (abs(Fai_p)<=Bita_s)				  
			w_Sita[Q+k] = 1.;				  
		 else				  
		    w_Sita[Q+k] = 0.;		 
	}


	Complex Sum_J(0.,0.);

	
  for (m=0;m<N;m++)
  {

	cout<<"Begin the "<<m+1<<"th step of the PBTG_FBM_SAA ......------>>"<<endl;


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

		         Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];			   

				
		         Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];

			     Sum_J = zero;
		}
		  Jf_1[i_n1] = ((vm[i_n1]-Sum0)*zm_d[i*N_dg+i]-(-1.*Sum1)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
		  Jf_2[i]    = ((-1.*Sum1)*zm_a[i_n1*N+i_n1]-(vm[i_n1]-Sum0)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
		   	
	

		  Sum0 = zero;
		  Sum1 = zero;
		  

	
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
			   
               Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];			   
			 		        				
			   Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];
		  	   Sum_J = zero;
		  }
		        
		        N_s_n1 = int((i-N_s-1)/n1);
	  	  
	                   
               NN = int((N_s)/n1)+1;

			                  
			   for (l=0;l<n1;l++)
				{
				   Sum_J = Sum_J + Jf_2[N_s_n1+l]+Jb_2[N_s_n1+l];
				}

			      Sum_J = Sum_J/double(n1);
	   
		      for (k=-Q;k<Q;k++)
			  {
			    Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
									        
			    
			    FN_f_1[k+Q] =  ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_N_i[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(JJ*K_i*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*JJ*K_i*z_N_i[N_s_n1]*sin(Fai_p));
                FN_f_2[k+Q] =  ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_Ndg_i[i-N_s-1]*cos(Fai_p)))*((Jf_1[j_n1]+Jb_1[j_n1]))-rou_i*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K1_i*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*JJ*K1_i*z_Ndg_i[i-N_s-1]*sin(Fai_p));
				
			 	  
				Sum2 = Sum2 + JJ*dx_Suf_N/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(JJ*K_i*z_N_i[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	    Sum3 = Sum3 + JJ*dx_Suf_Ndg/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(JJ*K1_i*z_Ndg_i[i]*sin(Fai_p))*Delta_Fai; 												
		 
			  }
		      Jf_1[i_n1] = ((vm[i_n1]-Sum0-Sum2)*zm_d[i*N_dg+i]-(-1.*Sum1-Sum3)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
		      Jf_2[i]    = ((-1.*Sum1-Sum3)*zm_a[i_n1*N+i_n1]-(vm[i_n1]-Sum0-Sum2)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);        

        	
		 Sum0 = zero;
		 Sum1 = zero;
		 Sum2 = zero;
		 Sum3 = zero;
		 Sum_J = zero;

					 

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
				
                Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];
		   
			 		        				
				Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];
		        
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
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
			 
			FN_f_1[k+Q] = FN_f_1[k+Q]*exp(JJ*K_i*dx_Suf_N*cos(Fai_p)) + ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_N_i[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(JJ*K_i*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*JJ*K_i*z_N_i[N_s_n1]*sin(Fai_p));
            FN_f_2[k+Q] = FN_f_2[k+Q]*exp(JJ*K1_i*dx_Suf_Ndg*cos(Fai_p)) + ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_Ndg_i[i-N_s-1]*cos(Fai_p)))*((Jf_1[j_n1]+Jb_1[j_n1]))-rou_i*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K1_i*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*JJ*K1_i*z_Ndg_i[i-N_s-1]*sin(Fai_p));

			Sum2 = Sum2 + JJ*dx_Suf_N/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(JJ*K_i*z_N_i[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_Suf_Ndg/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(JJ*K1_i*z_Ndg_i[i]*sin(Fai_p))*Delta_Fai; 												
		
		 }

		  Jf_1[i_n1] = ((vm[i_n1]-Sum0-Sum2)*zm_d[i*N_dg+i]-(-1.*Sum1-Sum3)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);

		  Jf_2[i]    = ((-1.*Sum1-Sum3)*zm_a[i_n1*N+i_n1]-(vm[i_n1]-Sum0-Sum2)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);        

		 Sum0 = zero;
		 Sum1 = zero;
		 Sum2 = zero;
		 Sum3 = zero;
		 Sum_J = zero;
	
	 	 
	 }	

 
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
				
                Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];
		   
			 		        				
				Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];
		        
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
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);
			 	        
			FN_f_1[k+Q] = FN_f_1[k+Q]*exp(JJ*K_i*dx_Suf_N*cos(Fai_p)) + ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_N_i[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(JJ*K_i*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*JJ*K_i*z_N_i[N_s_n1]*sin(Fai_p));
            FN_f_2[k+Q] = FN_f_2[k+Q]*exp(JJ*K1_i*dx_Suf_Ndg*cos(Fai_p)) + ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_Ndg_i[i-N_s-1]*cos(Fai_p)))*((Jf_1[j_n1]+Jb_1[j_n1]))-rou_i*(Jf_2[i-N_s-1]+Jb_2[i-N_s-1]))*exp(JJ*K1_i*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*JJ*K1_i*z_Ndg_i[i-N_s-1]*sin(Fai_p));
				
			Sum2 = Sum2 + JJ*dx_Suf_N/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_1[k+Q]*exp(JJ*K_i*z_N_i[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_Suf_Ndg/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_f_2[k+Q]*exp(JJ*K1_i*z_Ndg_i[i]*sin(Fai_p))*Delta_Fai; 												
		
		 }

		  Jf_1[i_n1] = ((vm[i_n1]-Sum0-Sum2)*zm_d[i*N_dg+i]-(-1.*Sum1-Sum3)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);

		  Jf_2[i]    = ((-1.*Sum1-Sum3)*zm_a[i_n1*N+i_n1]-(vm[i_n1]-Sum0-Sum2)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);        

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
			
             Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];
			
			 Sum1 = Sum1 + ((Jf_1[j_n1]+Jb_1[j_n1]))*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];
	
			 Sum_J = zero;
		}
      
		
		Jb_1[i_n1] = ((-1.*Sum0)*zm_d[i*N_dg+i]-(-1.*Sum1)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
		Jb_2[i]    = ((-1.*Sum1)*zm_a[i_n1*N+i_n1]-(-1.*Sum0)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
		

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
			
                Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];
			
			    Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];

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
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);


		                					  
            FN_b_1[k+Q] =  ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_N_i[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(JJ*K_i*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*JJ*K_i*z_N_i[N_s_n1]*sin(Fai_p));
            FN_b_2[k+Q] =  ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_Ndg_i[i+N_s+1]*cos(Fai_p)))*((Jf_1[j_n1]+Jb_1[j_n1]))-rou_i*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K1_i*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*JJ*K1_i*z_Ndg_i[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + JJ*dx_Suf_N/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(JJ*K_i*z_N_i[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_Suf_Ndg/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(JJ*K1_i*z_Ndg_i[i]*sin(Fai_p))*Delta_Fai; 												
		}
		
		   
		 Jb_1[i_n1] = ((-1.*Sum0-Sum2)*zm_d[i*N_dg+i]-(-1.*Sum1-Sum3)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);

		 Jb_2[i]    = ((-1.*Sum1-Sum3)*zm_a[i_n1*N+i_n1]-(-1.*Sum0-Sum2)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
			
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
            Sum0 = Sum0 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_a[i_n1*N+j_n1]+Sum_J*zm_b[i_n1*N+j_n1];
			
			Sum1 = Sum1 + (Jf_1[j_n1]+Jb_1[j_n1])*zm_c[i*N_dg+j]+(Jf_2[j]+Jb_2[j])*zm_d[i*N_dg+j];

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
			Fai_p = k*Delta_Fai*exp(-1.*JJ*Delta_0);

 
		    
					  
            FN_b_1[k+Q] = FN_b_1[k+Q]*exp(JJ*K_i*dx_Suf_N*cos(Fai_p)) + ((-1.*JJ*K_i*(-1.*sin(Fai_p)+df_x_N_i[N_s_n1]*cos(Fai_p)))*(Jf_1[N_s_n1]+Jb_1[N_s_n1])+Sum_J)*exp(JJ*K_i*double(NN)*dx_Suf_N*cos(Fai_p))*exp(-1.*JJ*K_i*z_N_i[N_s_n1]*sin(Fai_p));
            FN_b_2[k+Q] = FN_b_2[k+Q]*exp(JJ*K1_i*dx_Suf_Ndg*cos(Fai_p)) + ((JJ*K1_i*(-1.*sin(Fai_p)+df_x_Ndg_i[i+N_s+1]*cos(Fai_p)))*((Jf_1[j_n1]+Jb_1[j_n1]))-rou_i*(Jf_2[i+N_s+1]+Jb_2[i+N_s+1]))*exp(JJ*K1_i*(N_s+1.)*dx_Suf_Ndg*cos(Fai_p))*exp(-1.*JJ*K1_i*z_Ndg_i[i+N_s+1]*sin(Fai_p));
	      				      
			Sum2 = Sum2 + JJ*dx_Suf_N/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_1[k+Q]*exp(JJ*K_i*z_N_i[i_n1]*sin(Fai_p))*Delta_Fai; 									
         	Sum3 = Sum3 + JJ*dx_Suf_Ndg/(4.*PI)*exp(-1.*JJ*Delta_0)*w_Sita[Q+k]*FN_b_2[k+Q]*exp(JJ*K1_i*z_Ndg_i[i]*sin(Fai_p))*Delta_Fai; 												
		}
		 
		
		Jb_1[i_n1] = ((-Sum0-Sum2)*zm_d[i*N_dg+i]-(-1.*Sum1-Sum3)*zm_b[i_n1*N+i_n1])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);

		Jb_2[i]    = ((-1.*Sum1-Sum3)*zm_a[i_n1*N+i_n1]-(-Sum0-Sum2)*zm_c[i*N_dg+i])/(zm_a[i_n1*N+i_n1]*zm_d[i*N_dg+i]-zm_b[i_n1*N+i_n1]*zm_c[i*N_dg+i]);
			
		Sum0 = zero;
		Sum1 = zero;
		Sum2 = zero;
		Sum3 = zero;
		Sum_J = zero;
	
	
	 }

////////////////////////////////the convergence condition ////////////////////////////////////////////////////////////


		  

 		     for (i=0;i<N_dg;i++)
			 {
			    Jf2[i]   = (Jf_2[i]+Jb_2[i])-Jf1[i];
			    Jf3[i]   =  Jf_2[i]+Jb_2[i];
			 }

             Assign_Vecor(Jf3,Jf1);
			 Delta = sqrt(L2_Norm(Jf2))/sqrt(L2_Norm(Jf3));
             cout<<"The error is: Delta ="<<Delta<<endl;
			  
			 if(Delta<=0.01)
			 {
			    cout<<"The total number of the step is :"<<m+1<<"."<<"and the final error is :"<<Delta<<endl;
		        break;
			 }
				
	          else  
		           continue;

}

   

	for (i=0;i<N_dg;i++)
	{
		Jm_Ndg[i] = Jf3[i];
	}

	for (i=0;i<N;i++)
	{
		Jm_N[i] = Jf_1[i] + Jb_1[i];
	}


}










Complex resolve_matrix::Get_Zmn(int j,Vec_Zmn_Cell &a)
{//对稀疏存储的向量，只要给定行数和列数，即可找到对应的元素
	Complex result;
	result = 0.;
	for ( int ii=0; ii<a.size(); ii++ )
	{
		while ( a[ii].Zmn_j == j )
		{
			result = a[ii].Zmn_z;
			break;
		}
	}
	return result;
}
void resolve_matrix::Resolve_ILU( Vec_Zmn_Cell *M_L, Vec_Zmn_Cell *M_U , Vec_Com &Bm, Vec_Com &result )
{
	int ii,jj;
	int n = Bm.size();
	result.assign(n);
	for ( ii=0; ii<n; ii++ )
	{
		for ( jj=0; jj<M_L[ii].size(); jj++ )
		{
//			result[ii] = Bm[ii] - 
		}
	}
}
void resolve_matrix::test_zmn_cell()
{
//	b = a;
	int i,j;
	//////////////////////////////////////////////////////////////////////////
	Complex temp_matrix[8][8] =
	{
		{0.0737,         0,    0.6868,    0.7012,    0.7564,         0,         0,    0.5045},
		{0     ,    0.5007,    0.5890,    0.9103,    0.9910,    0.6316,         0,    0.5163},
		{0.6789,    0.5297,    1.3146,    0.7622,         0,    0.8847,         0,         0},
		{0.6793,    0.6711,    0.8462,    0.2771,         0,         0,    0.4865,    0.9866},
		{0.9347,         0,    0.5269,         0,    1.8964,    0.4364,    0.8977,    0.4940},
        {0     ,         0,         0,    0.7361,    0.7227,    1.2962,    0.9092,         0},
		{0.5194,         0,    0.6539,         0,    0.7534,    0.4777,    0.4644,         0},
		{0.8310,    0.4175,    0.4160,    0.6326,    0.6515,         0,    0.9047,    1.8887}
	};
	Vec_Zmn_Cell WW[8];
	Zmn_Cell temp_zmn;
	for ( i=0; i<8; i++ )
	{
		for ( j=0; j<8; j++ )
		{
			if ( abs(temp_matrix[i][j]) > 0. )
			{
				temp_zmn.Zmn_i = i;
				temp_zmn.Zmn_j = j;
				temp_zmn.Zmn_z = temp_matrix[i][j];
				WW[i].push_back(temp_zmn);
			}
		}
	}
/*	for ( i=0; i<8; i++ )
	{
		for ( j=0; j<WW[i].size(); j++ )
		{
			cout<<WW[i][j].Zmn_i<<" "<<WW[i][j].Zmn_j<<" "<<WW[i][j].Zmn_z<<endl;
		}
	}*/
	//////////////////////////////////////////////////////////////////////////
//	Vec_Zmn_Cell M_L[8], M_U[8];
//	Vec_Zmn_Cell *input;
//	Vec_Zmn_Cell *outp_l,*outp_u;
//	input = &WW[0];
//	outp_l = &M_L[0];
//	outp_u = &M_U[0];
/*	sort_vec_zmn_cell(WW[3]);
	for ( i=0; i<WW[3].size(); i++ )
	{
		cout<<WW[3][i].Zmn_i<<" "<<WW[3][i].Zmn_j<<" "<<WW[3][i].Zmn_z<<endl;
	}*/
	
	Vec_Zmn_Cell *M_L,*M_U;
	int Num = 8;
	M_L = new Vec_Zmn_Cell[Num];
	M_U = new Vec_Zmn_Cell[Num];
//	Test_ILU_p_t( 8, WW, 4, 0.05 );
	ILU_0( 8, WW, M_L, M_U );

	for ( int ii=0; ii<Num; ii++ )
	{
		for ( int jj=0; jj<M_L[ii].size(); jj++ )
		{
			cout<<M_L[ii][jj].Zmn_z<<"  ";
		}
		cout<<endl;
	}
}
//ILU(p,t):Matrix Preconditioning Techniques and Applications(P175)
void resolve_matrix::Test_ILU_p_t( int nn, Vec_Zmn_Cell *WW,int P,double tau )
{//仍然有问题
//	cout<<WW[0][1].Zmn_i<<" "<<WW[0][1].Zmn_j<<" "<<WW[0][1].Zmn_z<<endl;
	int ii,jj,k;
	double TOL = 1.e-12;
	int Num = nn;
	Vec_Zmn_Cell *B = new Vec_Zmn_Cell[Num];
	Vec_Com *Value_WW = new Vec_Com[Num];
	Vec_Zmn_Cell *M_L,*M_U;
	M_L = new Vec_Zmn_Cell[Num];
	M_U = new Vec_Zmn_Cell[Num];
	Zmn_Cell zsmn_t;

	//Step1
	for ( ii=0; ii<Num; ii++ )
	{
		for ( jj=0; jj<WW[ii].size(); jj++ )
		{
			Value_WW[ii].push_back(WW[ii][jj].Zmn_z);
//			cout<<WW[ii][jj].Zmn_i<<" "<<WW[ii][jj].Zmn_j<<" "<<WW[ii][jj].Zmn_z<<endl;
		}
//		cout<<WW[ii].size()<<endl;
	}
	//Step2
	Vec_Dou VV;
	VV.assign(Num);
	for ( ii=0; ii<Num; ii++ )
	{
		VV[ii] = tau * L2_Norm_Inf(Value_WW[ii]);
	}
	//Step3
	int JL,JU;
	Vec_Zmn_Cell temp_ww1,temp_ww2;
	vector <int> jj_ii,W_I_N;
	jj_ii.assign(Num);
	W_I_N.assign(Num);
	Complex tmp,Value_temp;
	//计算每行下三角的非零个数
	for ( ii=0; ii<Num; ii++ )
	{
		W_I_N[ii] = 0;
		for ( jj=0; jj<WW[ii].size(); jj++ )
		{
			if ( WW[ii][jj].Zmn_j < ii )
			{
				W_I_N[ii] += 1;
			}
		}
//		cout<<W_I_N[ii]<<endl;
	}

//	cout<<WW[0][1].Zmn_i<<" "<<WW[0][1].Zmn_j<<" "<<WW[0][1].Zmn_z<<endl;

	//Step
	
//		cout<<JL<<"  "<<JU<<endl;

	for ( ii=0; ii<Num; ii++ )
	{
		for ( k=0; k<W_I_N[ii]; k++ )
		{
//			WW[ii][k].Zmn_z = WW[ii][k].Zmn_z / M_U[WW[ii][k].Zmn_j][0].Zmn_z;
			WW[ii][k].Zmn_z = WW[ii][k].Zmn_z / B[WW[ii][k].Zmn_j][jj_ii[WW[ii][k].Zmn_j]].Zmn_z;
			if ( abs(WW[ii][k].Zmn_z) < VV[ii] )
			{
				WW[ii][k].Zmn_z = 0.;
			}
//			cout<<B[WW[ii][k].Zmn_j][jj_ii[WW[ii][k].Zmn_j]].Zmn_z<<endl;
//			cout<<WW[ii][k].Zmn_i<<" "<<WW[ii][k].Zmn_j<<" "<<WW[ii][k].Zmn_z<<endl;
/*			for ( jj=0; jj<W_I_N[ii]; jj++ )
			{
				cout<<WW[ii][jj].Zmn_i<<" "<<WW[ii][jj].Zmn_j<<" "<<WW[ii][jj].Zmn_z<<endl;
			}*/

			if ( abs(WW[ii][k].Zmn_z) >= VV[ii] )
			{				
				for ( jj=k+1; jj<WW[ii].size(); jj++ )
				{
//					cout<<Get_Zmn(7,WW[0])<<endl;
/*					if ( WW[ii][jj].Zmn_j < ii )
					{
						tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*Get_Zmn(WW[ii][jj].Zmn_j,M_L[WW[ii][k].Zmn_j]);
//						cout<<Get_Zmn(WW[ii][jj].Zmn_j,M_L[WW[ii][k].Zmn_j])<<endl;
					}
					else
					{
						tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*Get_Zmn(WW[ii][jj].Zmn_j,M_U[WW[ii][k].Zmn_j]);
//						cout<<Get_Zmn(WW[ii][jj].Zmn_j,M_U[WW[ii][k].Zmn_j])<<endl;
					}*/
					//下面这步是对所有的元素进行处理，程序还有问题
/*					Value_temp = Get_Zmn(WW[ii][jj].Zmn_j,B[WW[ii][k].Zmn_j]);
					if ( Value_temp != 0. )
					{
						tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*Value_temp;
					}*/				
//					cout<<Value_temp<<endl;
					tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*Get_Zmn(WW[ii][jj].Zmn_j,B[WW[ii][k].Zmn_j]);
//					cout<<Get_Zmn(WW[ii][jj].Zmn_j,B[WW[ii][k].Zmn_j])<<endl;
					if ( abs(tmp) <= VV[ii] )
					{
						tmp = 0.;
					}
					WW[ii][jj].Zmn_z = tmp;
				}
			}
		}

/*		for ( jj=0; jj<WW[ii].size(); jj++ )
		{
			cout<<WW[ii][jj].Zmn_i<<" "<<WW[ii][jj].Zmn_j<<" "<<WW[ii][jj].Zmn_z<<endl;
		}*/		
		
		JL = 0;	JU = 0;
		for ( jj=0; jj<WW[ii].size(); jj++ )
		{
			if ( WW[ii][jj].Zmn_j < ii )
			{
				JL+=1;
			}
			else if ( WW[ii][jj].Zmn_j == ii )
			{
				jj_ii[ii] = jj;
			}
			else
			{
				JU+=1;
			}
		}
//		cout<<JL<<"  "<<JU<<endl;
//		cout<<jj_ii[ii]<<endl;
		
		//Step
		//前提条件是每行元素都是按行数从小到大的顺序排列
		if ( JL > P )
		{
			//存储
			for ( jj=0; jj<JL; jj++ )
			{
				temp_ww1.push_back(WW[ii][jj]);				
			}
			//排序（对元素为Vec_Zmn_Cell类型的向量排序）
			sort_vec_zmn_cell(temp_ww1);			
		}

		if ( JU > P-1 )
		{
			for ( jj=0; jj<JU; jj++ )
			{
				temp_ww2.push_back(WW[ii][JL+1+jj]);
			}
			sort_vec_zmn_cell(temp_ww2);
		}

		//B
		if ( JL > P )
		{
			//取前P个值
			for ( jj=0; jj<P; jj++ )
			{
				B[ii].push_back(temp_ww1[jj]);
			}
		}
		else
		{			
			for ( jj=0; jj<W_I_N[ii]; jj++ )
			{
				B[ii].push_back(WW[ii][jj]);
			}
		}
		
		B[ii].push_back(WW[ii][jj_ii[ii]]);
		if ( JU > P-1 )
		{
			for ( jj=0; jj<P-1; jj++ )
			{
				B[ii].push_back(temp_ww2[jj]);
			}
		}
		else
		{
			for ( jj=JL+1; jj<WW[ii].size(); jj++ )
			{
				B[ii].push_back(WW[ii][jj]);
			}
		}
		
		
/*		if ( JL > P )
		{
			//取前P个值
			for ( jj=0; jj<P; jj++ )
			{
				M_L[ii].push_back(temp_ww1[jj]);
			}
		}
		else
		{
			for ( jj=0; jj<W_I_N[ii]; jj++ )
			{
				M_L[ii].push_back(WW[ii][jj]);
			}
		}
		
		M_U[ii].push_back(WW[ii][jj_ii[ii]]);
		if ( JU > P-1 )
		{
			for ( jj=0; jj<P-1; jj++ )
			{
				M_U[ii].push_back(temp_ww2[jj]);
			}
		}
		else
		{
			for ( jj=JL+1; jj<WW[ii].size(); jj++ )
			{
				M_U[ii].push_back(WW[ii][jj]);
			}
		}*/

		//
//		cout<<temp_ww1.size()<<"  "<<temp_ww2.size()<<endl;
		temp_ww1.clear();
		temp_ww2.clear();
	}

	for ( ii=0; ii<Num; ii++ )
	{
		zsmn_t.Zmn_i = ii;
		zsmn_t.Zmn_j = ii;
		zsmn_t.Zmn_z = 1.0;
		M_L[ii].push_back(zsmn_t);
		for ( jj=0; jj<B[ii].size(); jj++ )
		{
			if ( jj < W_I_N[ii] )
			{
				M_L[ii].push_back(B[ii][jj]);
			}
			else
			{
				M_U[ii].push_back(B[ii][jj]);
			}
		}
	}

	for ( ii=0; ii<Num; ii++ )
	{
		for ( jj=0; jj<B[ii].size(); jj++ )
		{
			cout<<B[ii][jj].Zmn_z<<"  ";
		}
		cout<<endl;
	}

	for ( ii=0; ii<Num; ii++ )
	{
		WW[ii].clear();
		Value_WW[ii].clear();
	}
	jj_ii.clear();
}
void resolve_matrix::ILU_p_t( int nn, Vec_Zmn_Cell zs[6][6],int P,double tau,Vec_Zmn_Cell *M_L,Vec_Zmn_Cell *M_U )
{
	int ii,jj,pp,qq,k;
	double TOL = 1.e-12;
	int Num = nn*6;
	Vec_Zmn_Cell *WW = new Vec_Zmn_Cell[Num];
	Vec_Com *Value_WW = new Vec_Com[Num];
	M_L = new Vec_Zmn_Cell[Num];
	M_U = new Vec_Zmn_Cell[Num];
	Zmn_Cell zsmn_t;

	//Step1
	for ( ii=0; ii<6; ii++ )
	{
		for ( jj=0; jj<6; jj++ )
		{
			for ( pp=0; pp<zs[ii][jj].size(); pp++ )
			{
				for ( qq=0; qq<Num; qq++ )
				{
					if ( zs[ii][jj][pp].Zmn_i == qq )
					{
						zsmn_t.Zmn_i = qq;
						zsmn_t.Zmn_j = zs[ii][jj][pp].Zmn_j;
						zsmn_t.Zmn_z = zs[ii][jj][pp].Zmn_z;
						WW[qq].push_back(zsmn_t);
						Value_WW[qq].push_back(zsmn_t.Zmn_z);
					}
				}
			}
		}
	}
	//Step2
	Vec_Dou VV;
	VV.assign(Num);
	for ( ii=0; ii<Num; ii++ )
	{
		VV[ii] = tau * L2_Norm_Inf(Value_WW[ii]);
	}
	//Step3
	int JL,JU;
	Vec_Zmn_Cell temp_ww1,temp_ww2;
	vector <int> jj_ii;
	Complex tmp;

	for ( ii=0; ii<Num; ii++ )
	{
		k = 0;
		if ( WW[ii][k].Zmn_j < ii )
		{
			if ( abs( WW[ii][k].Zmn_z ) > TOL )
			{
				WW[ii][k].Zmn_z = WW[ii][k].Zmn_z / M_U[k][jj_ii[k]].Zmn_z;
				if ( abs(WW[ii][k].Zmn_z) < VV[ii] )
				{
					WW[ii][k].Zmn_z = 0.;
				}
			}

			if ( abs(WW[ii][k].Zmn_z) >= VV[ii] )
			{				
				for ( jj=0; jj<WW[ii].size(); jj++ )
				{
					if ( WW[ii][jj].Zmn_j >= k+1 )
					{
						if ( WW[ii][jj].Zmn_j < ii )
						{
							tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*M_L[k][jj].Zmn_z;
						}
						else
						{
							tmp = WW[ii][jj].Zmn_z - WW[ii][k].Zmn_z*M_U[k][jj].Zmn_z;
						}						
						if ( abs(tmp) <= VV[ii] )
						{
							tmp = 0.;
						}
						WW[ii][jj].Zmn_z = tmp;
					}					
				}
			}

			k+=1;
		}			
		
		//Step
		JL = 0;	JU = 0;
		for ( jj=0; jj<WW[ii].size(); jj++ )
		{			
			if ( WW[ii][jj].Zmn_j < ii )
			{
				JL+=1;
			}
			else if ( WW[ii][jj].Zmn_j = ii )
			{
				jj_ii.push_back(jj);
			}
			else
			{
				JU+=1;
			}
		}
		
		//Step
		if ( JL > P )
		{
			//存储
			for ( jj=0; jj<JL; jj++ )
			{
//				zsmn_t.Zmn_i = WW[ii][jj].Zmn_i;
//				zsmn_t.Zmn_j = WW[ii][jj].Zmn_j;
//				zsmn_t.Zmn_z = WW[ii][jj].Zmn_z;
				temp_ww1.push_back(WW[ii][jj]);				
			}
			//排序（对元素为Vec_Zmn_Cell类型的向量排序）
			sort_vec_zmn_cell(temp_ww1);			
		}

		if ( JU > P-1 )
		{
			for ( jj=0; jj<JU; jj++ )
			{
//				zsmn_t.Zmn_i = WW[ii][JL+1+jj].Zmn_i;
//				zsmn_t.Zmn_j = WW[ii][JL+1+jj].Zmn_j;
//				zsmn_t.Zmn_z = WW[ii][JL+1+jj].Zmn_z;
				temp_ww2.push_back(WW[ii][JL+1+jj]);
			}
			sort_vec_zmn_cell(temp_ww2);
		}		

		//
		zsmn_t.Zmn_i = ii;
		zsmn_t.Zmn_j = ii;
		zsmn_t.Zmn_z = 1.0;
		M_L[ii].push_back(zsmn_t);
		if ( JL > P )
		{
			//取前P个值
			for ( jj=0; jj<P; jj++ )
			{
//				zsmn_t.Zmn_i = temp_ww1[jj].Zmn_i;
//				zsmn_t.Zmn_j = temp_ww1[jj].Zmn_j;
//				zsmn_t.Zmn_z = temp_ww1[jj].Zmn_z;
				M_L[ii].push_back(temp_ww1[jj]);
			}
		}
		else
		{
			jj=0;
			if ( WW[ii][jj].Zmn_j < ii )
			{
//				zsmn_t.Zmn_i = WW[ii][jj].Zmn_i;
//				zsmn_t.Zmn_j = WW[ii][jj].Zmn_j;
//				zsmn_t.Zmn_z = WW[ii][jj].Zmn_z;
				M_L[ii].push_back(WW[ii][jj]);
				jj+=1;
			}
		}
		
		M_U[ii].push_back(WW[ii][jj_ii[ii]]);
		if ( JU > P-1 )
		{
			for ( jj=0; jj<P-1; jj++ )
			{
//				zsmn_t.Zmn_i = temp_ww2[jj].Zmn_i;
//				zsmn_t.Zmn_j = temp_ww2[jj].Zmn_j;
//				zsmn_t.Zmn_z = temp_ww2[jj].Zmn_z;
				M_U[ii].push_back(temp_ww2[jj]);
			}
		}
		else
		{
			for ( jj=JL+1; jj<WW[ii].size(); jj++ )
			{
//				zsmn_t.Zmn_i = WW[ii][jj].Zmn_i;
//				zsmn_t.Zmn_j = WW[ii][jj].Zmn_j;
//				zsmn_t.Zmn_z = WW[ii][jj].Zmn_z;
				M_U[ii].push_back(WW[ii][jj]);
			}
		}

		//
		temp_ww1.clear();
		temp_ww2.clear();
	}

	for ( ii=0; ii<Num; ii++ )
	{
		WW[ii].clear();
		Value_WW[ii].clear();
	}
	jj_ii.clear();
}
void resolve_matrix::sort_vec_zmn_cell(Vec_Zmn_Cell &vec_cell)
{
	int n = vec_cell.size();
	int m,k,j,i;
	Zmn_Cell zmn_t;
//  Complex Z_z; int Z_I,Z_J;
    k=0;
	m=n-1;
    while ( k<m )
	{
		j=m-1;
		m=0;
        for ( i=k; i<=j; i++ )
		{
			if ( abs(vec_cell[i].Zmn_z) < abs(vec_cell[i+1].Zmn_z) )
            {
/*				Z_z = vec_cell[i].Zmn_z;
				vec_cell[i].Zmn_z = vec_cell[i+1].Zmn_z;
				vec_cell[i+1].Zmn_z = Z_z;

				Z_I = vec_cell[i].Zmn_i;
				vec_cell[i].Zmn_i = vec_cell[i+1].Zmn_i;
				vec_cell[i+1].Zmn_i = Z_I;

				Z_J = vec_cell[i].Zmn_j;
				vec_cell[i].Zmn_j = vec_cell[i+1].Zmn_j;
				vec_cell[i+1].Zmn_j = Z_J;*/

				zmn_t = vec_cell[i];
				vec_cell[i] = vec_cell[i+1];
				vec_cell[i+1] = zmn_t;
				
				m=i;
			}
		}
		j=k+1;
		k=0;
		for ( i=m; i>=j; i-- )
		{
			if ( abs(vec_cell[i-1].Zmn_z) < abs(vec_cell[i].Zmn_z) )
			{ 
/*				Z_z = vec_cell[i].Zmn_z; 
				vec_cell[i].Zmn_z = vec_cell[i-1].Zmn_z; 
				vec_cell[i-1].Zmn_z = Z_z;
				
				Z_I = vec_cell[i].Zmn_i;
				vec_cell[i].Zmn_i = vec_cell[i-1].Zmn_i;
				vec_cell[i-1].Zmn_i = Z_I;

				Z_J = vec_cell[i].Zmn_j;
				vec_cell[i].Zmn_j = vec_cell[i-1].Zmn_j;
				vec_cell[i-1].Zmn_j = Z_J;*/

				zmn_t = vec_cell[i];
				vec_cell[i] = vec_cell[i-1];
				vec_cell[i-1] = zmn_t;
				
				k=i;
			}
		}
	}
}
void resolve_matrix::ILU_0( int n,Vec_Zmn_Cell *WW, Vec_Zmn_Cell *M_L, Vec_Zmn_Cell *M_U )
{//注意：Matlab存储数据位数很长，所以数据运算有误差
	int ii,jj,k;
	int Num = n;
	
//	M_L = new Vec_Zmn_Cell[Num];
//	M_U = new Vec_Zmn_Cell[Num];
	Vec_Zmn_Cell *B = new Vec_Zmn_Cell[Num];
	Zmn_Cell zsmn_t;

	//赋给B后不会改变原矩阵
	for ( ii=0; ii<Num; ii++ )
	{
		for ( jj=0; jj<WW[ii].size(); jj++ )
		{
			B[ii].push_back(WW[ii][jj]);
//			cout<<B[ii][jj].Zmn_z<<" ";
		}
//		cout<<endl;
	}

	//Step3
	vector <int> jj_ii,W_I_N;
	jj_ii.assign(Num);
	W_I_N.assign(Num);
	Complex tmp;
	//计算每行下三角的非零个数
	for ( ii=0; ii<Num; ii++ )
	{
		W_I_N[ii] = 0;
		for ( jj=0; jj<B[ii].size(); jj++ )
		{
			if ( B[ii][jj].Zmn_j < ii )
			{
				W_I_N[ii] += 1;
			}
			if ( B[ii][jj].Zmn_j == ii )
			{
				jj_ii[ii] = jj;
			}
//			cout<<Get_Zmn(WW[ii][jj].Zmn_j,WW[ii])<<" ";
		}
//		cout<<W_I_N[ii]<<endl;
//		cout<<jj_ii[ii]<<endl;
//		cout<<endl;
	}

	for ( ii=1; ii<Num; ii++ )
	{
		for ( k=0; k<W_I_N[ii]; k++ )
		{
			B[ii][k].Zmn_z = B[ii][k].Zmn_z / B[B[ii][k].Zmn_j][jj_ii[B[ii][k].Zmn_j]].Zmn_z;
//			cout<<B[ii][k].Zmn_i<<" "<<B[ii][k].Zmn_j<<" "<<B[ii][k].Zmn_z<<endl;

			for ( jj=k+1; jj<B[ii].size(); jj++ )
			{
				tmp = Get_Zmn(B[ii][jj].Zmn_j,B[B[ii][k].Zmn_j]);
				if ( tmp != 0. )
				{
					B[ii][jj].Zmn_z -= B[ii][k].Zmn_z*tmp;
				}				
//				cout<<tmp<<endl;
			}
		}
	}
	
	//Step
	for ( ii=0; ii<Num; ii++ )
	{
		zsmn_t.Zmn_i = ii;
		zsmn_t.Zmn_j = ii;
		zsmn_t.Zmn_z = 1.0;
		M_L[ii].push_back(zsmn_t);
		for ( jj=0; jj<B[ii].size(); jj++ )
		{
			if ( jj < W_I_N[ii] )
			{
				M_L[ii].push_back(B[ii][jj]);
			}
			else
			{
				M_U[ii].push_back(B[ii][jj]);
			}
		}
	}
	
/*	for ( ii=0; ii<Num; ii++ )
	{
		for ( jj=0; jj<B[ii].size(); jj++ )
		{
			cout<<B[ii][jj].Zmn_z<<"  ";
		}
		cout<<endl;
	}*/

	for ( ii=0; ii<Num; ii++ )
	{
		B[ii].clear();
	}
	jj_ii.clear();
	W_I_N.clear();
}
//////////////////////////////////////////////////////////////////////////
//二维介质粗糙面电磁散射仿真计算的弱矩阵向量积的计算及其系数等的计算
void resolve_matrix::Two_D_Diel_Zs_Vector_Multiply(char a, Vec_Zmn_Cell zs[6][6], Vec_Com vm[6], Vec_Com re[6])
{
	int i,j,k;
	int nxy = vm[0].size();
	for ( i=0; i<6; i++ )
	{
		re[i].assign(nxy);
	}
	Vec_Com MV[6][6];

	if ( a == ' ' )
	{
		for ( i=0; i<6; i++ )
		{
			for ( j=0; j<6; j++ )
			{
				if ( zs[i][j].size() != 0 )
				{
					Two_D_Diel_Matrix_Vector_Multiply(' ',zs[i][j],vm[j],MV[i][j]);
					for ( k=0; k<nxy; k++ )
					{
						re[i][k] += MV[i][j][k];
					}
				}
			}
		}
	}
	else//转置矩阵向量积
	{
		for ( i=0; i<6; i++ )
		{
			for ( j=0; j<6; j++ )
			{
				if ( zs[j][i].size() != 0 )
				{
					Two_D_Diel_Matrix_Vector_Multiply('*',zs[j][i],vm[j],MV[i][j]);
					for ( k=0; k<nxy; k++ )
					{
						re[i][k] += MV[i][j][k];
					}
				}
			}
		}
	}
	
	for ( i=0; i<6; i++ )
	{
		for ( j=0; j<6; j++ )
		{
			MV[i][j].clear();
		}
	}
}
void resolve_matrix::Two_D_Diel_Matrix_Vector_Multiply(char a, Vec_Zmn_Cell &zm, Vec_Com &vm, Vec_Com &re)
{
	int nxy = vm.size();
	re.assign(nxy);
	int zm_num;
	
	if( a == ' ' )
	{
		for( zm_num=0; zm_num<zm.size(); zm_num++)
		{
			re[zm[zm_num].Zmn_i] += zm[zm_num].Zmn_z*vm[zm[zm_num].Zmn_j];
		}
	}
	else
	{//矩阵转置共轭的矩阵向量积
		for( zm_num=0; zm_num<zm.size(); zm_num++ )
		{
			re[zm[zm_num].Zmn_j] += conj(zm[zm_num].Zmn_z)*vm[zm[zm_num].Zmn_i];
		}
	}
}
void resolve_matrix::Two_D_Diel_Multi_Zw_Vec(int Para_G_a_Num,int Para_g_b_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6])
{
	int i;
	int Num = nxy*nxy;
	for ( i=0; i<6; i++ )
	{
		re[i].assign(Num);
	}

	Vec_Com Tre[6],Tvec[5];
	for ( i=0; i<5; i++ )
	{
		Tvec[i].assign(Num);
	}

	//////////////////////////////////////////////////////////////////////////
	//自由空间中的积分方程
	//第一个矩阵方程
	//re0
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[0][i]*ymn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i];
		Tvec[2][i] = vec[0][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[0][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
	for ( i=0; i<Num; i++ )
	{
		re[0][i] = dfymn[i]*(ymn[i]*Tre[0][i]-Tre[1][i])+(xmn[i]*Tre[2][i]-Tre[3][i])-(zmn[i]*Tre[0][i]-Tre[4][i]);
	}
	for ( i=0; i<5; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[1][i]*xmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[0][i] += -dfymn[i]*(xmn[i]*Tre[0][i]-Tre[1][i])+(xmn[i]*Tre[2][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}

	//re2
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[2][i]*zmn[i];
		Tvec[1][i] = vec[2][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	for ( i=0; i<Num; i++ )
	{
		re[0][i] += dfymn[i]*(zmn[i]*Tre[0][i]-Tre[1][i])+(ymn[i]*Tre[0][i]-Tre[2][i]);
	}
	for ( i=0; i<3; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[0]);
	for ( i=0; i<Num; i++ )
	{
		re[0][i] += -i_*(k1/wave_impe1)*dfymn[i]*Tre[0][i];
	}
	for ( i=0; i<1; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[4],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[0][i] += -i_*(k1/wave_impe1)*(dfymn[i]*Tre[1][i]+Tre[0][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//第二个积分方程
	//re0
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[0][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[1][i] = (ymn[i]*Tre[1][i]-Tre[2][i])-dfxmn[i]*(ymn[i]*Tre[0][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*zmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[1][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
	for ( i=0; i<Num; i++ )
	{
		re[1][i] += -(zmn[i]*Tre[0][i]-Tre[1][i])+(ymn[i]*Tre[2][i]-Tre[3][i])+dfxmn[i]*(xmn[i]*Tre[0][i]-Tre[4][i]);
	}
	for ( i=0; i<5; i++ )
	{
		Tre[i].clear();
	}

	//re2
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[2][i]*zmn[i];
		Tvec[1][i] = vec[2][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	for ( i=0; i<Num; i++ )
	{
		re[1][i] += -dfxmn[i]*(zmn[i]*Tre[0][i]-Tre[1][i])-(xmn[i]*Tre[0][i]-Tre[2][i]);
	}
	for ( i=0; i<3; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[3],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[1][i] += -i_*(k1/wave_impe1)*((-1.-dfxmn[i])*Tre[0][i]-Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[0]);
	for ( i=0; i<Num; i++ )
	{
		re[1][i] += i_*(k1/wave_impe1)*dfxmn[i]*Tre[0][i];
	}
	for ( i=0; i<1; i++ )
	{
		Tre[i].clear();
	}

	//第三个积分方程
	//re0
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[2][i] = i_*(k1*wave_impe1)*(dfxmn[i]*Tre[0][i]-Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[2][i] += i_*(k1*wave_impe1)*(dfymn[i]*Tre[0][i]-Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[3][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[3][i]*zmn[i];
		Tvec[4][i] = vec[3][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,Tre[5]);
	for ( i=0; i<Num; i++ )
	{
		re[2][i] += -dfxmn[i]*(ymn[i]*Tre[1][i]-Tre[2][i])+dfymn[i]*((xmn[i]*Tre[1][i]-Tre[3][i])-(zmn[i]*Tre[0][i]-Tre[4][i]))-(ymn[i]*Tre[0][i]-Tre[5][i]);
	}
	for ( i=0; i<6; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*zmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[4][i]*dfymn[i]*xmn[i];
		Tvec[4][i] = vec[4][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,Tre[5]);
	for ( i=0; i<Num; i++ )
	{
		re[2][i] += dfxmn[i]*((zmn[i]*Tre[0][i]-Tre[1][i])-(ymn[i]*Tre[2][i]-Tre[3][i]))+dfymn[i]*(xmn[i]*Tre[2][i]-Tre[4][i])+(xmn[i]*Tre[0][i]-Tre[5][i]);
	}
	for ( i=0; i<6; i++ )
	{
		Tre[i].clear();
	}

	//re5
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*xmn[i];
		Tvec[1][i] = vec[5][i]*ymn[i];
		Tvec[2][i] = vec[5][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k1,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[2][i] += dfxmn[i]*(xmn[i]*Tre[0][i]-Tre[1][i])+dfymn[i]*(ymn[i]*Tre[0][i]-Tre[2][i])-(zmn[i]*Tre[0][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}
	
	//////////////////////////////////////////////////////////////////////////
	//介质空间中的积分方程
	//第四个积分方程
	//re0
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[0]);
	for ( i=0; i<Num; i++ )
	{
		re[3][i] = i_*(k1*wave_impe1)*dfymn[i]*Tre[0][i];//原来有错误：dfymn写成了dfxmn
	}
	for ( i=0; i<1; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[1],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[3][i] += i_*(k1*wave_impe1)*(dfymn[i]*Tre[1][i]+Tre[0][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*ymn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i];
		Tvec[2][i] = vec[3][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[3][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
	for ( i=0; i<Num; i++ )
	{
		re[3][i] += dfymn[i]*(ymn[i]*Tre[0][i]-Tre[1][i])+(xmn[i]*Tre[2][i]-Tre[3][i])-(zmn[i]*Tre[0][i]-Tre[4][i]);
	}
	for ( i=0; i<5; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*xmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[3][i] += -dfymn[i]*(xmn[i]*Tre[0][i]-Tre[1][i])+(xmn[i]*Tre[2][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}

	//re5
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*zmn[i];
		Tvec[1][i] = vec[5][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	for ( i=0; i<Num; i++ )
	{
		re[3][i] += (1./permi)*(dfymn[i]*(zmn[i]*Tre[0][i]-Tre[1][i])+(ymn[i]*Tre[0][i]-Tre[2][i]));
	}
	for ( i=0; i<3; i++ )
	{
		Tre[i].clear();
	}
	
	//第五个积分方程
	//re0
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[0],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[4][i] = i_*(k1*wave_impe1)*(-Tre[0][i]-dfxmn[i]*Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[0]);
	for ( i=0; i<Num; i++ )
	{
		re[4][i] += i_*(k1*wave_impe1)*(-dfxmn[i]*Tre[0][i]);
	}
	for ( i=0; i<1; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[3][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[4][i] += (ymn[i]*Tre[1][i]-Tre[2][i])-dfxmn[i]*(ymn[i]*Tre[0][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*zmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[4][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
	for ( i=0; i<Num; i++ )
	{
		re[4][i] += -(zmn[i]*Tre[0][i]-Tre[1][i])+(ymn[i]*Tre[2][i]-Tre[3][i])+dfxmn[i]*(xmn[i]*Tre[0][i]-Tre[4][i]);
	}
	for ( i=0; i<5; i++ )
	{
		Tre[i].clear();
	}

	//re5
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*zmn[i];
		Tvec[1][i] = vec[5][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	for ( i=0; i<Num; i++ )
	{
		re[4][i] += (1./permi)*(-dfxmn[i]*(zmn[i]*Tre[0][i]-Tre[1][i])+(xmn[i]*Tre[0][i]-Tre[2][i]));
	}
	for ( i=0; i<3; i++ )
	{
		Tre[i].clear();
	}

	//第六个积分方程
	//re0
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[0][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[0][i]*zmn[i];
		Tvec[4][i] = vec[0][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,Tre[5]);
	for ( i=0; i<Num; i++ )
	{
		re[5][i] = -dfxmn[i]*(ymn[i]*Tre[1][i]-Tre[2][i])+dfymn[i]*((xmn[i]*Tre[1][i]-Tre[3][i])-(zmn[i]*Tre[0][i]-Tre[4][i]))-(ymn[i]*Tre[0][i]-Tre[5][i]);
	}
	for ( i=0; i<6; i++ )
	{
		Tre[i].clear();
	}

	//re1
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*zmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*xmn[i];
		Tvec[3][i] = vec[1][i]*dfymn[i]*ymn[i];
		Tvec[4][i] = vec[1][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,Tre[5]);
	for ( i=0; i<Num; i++ )
	{
		re[5][i] += dfxmn[i]*(zmn[i]*Tre[0][i]-Tre[1][i])+dfymn[i]*(xmn[i]*Tre[2][i]-Tre[3][i])-dfxmn[i]*(ymn[i]*Tre[2][i]-Tre[4][i])+(xmn[i]*Tre[0][i]-Tre[5][i]);
	}
	for ( i=0; i<6; i++ )
	{
		Tre[i].clear();
	}

	//re2
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[2][i]*xmn[i];
		Tvec[1][i] = vec[2][i]*ymn[i];
		Tvec[2][i] = vec[2][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
	for ( i=0; i<Num; i++ )
	{
		re[5][i] += dfxmn[i]*(xmn[i]*Tre[0][i]-Tre[1][i])+dfymn[i]*(ymn[i]*Tre[0][i]-Tre[2][i])-(zmn[i]*Tre[0][i]-Tre[3][i]);
	}
	for ( i=0; i<4; i++ )
	{
		Tre[i].clear();
	}

	//re3
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[5][i] += -i_*(k1*wave_impe1*permi)*(dfxmn[i]*Tre[0][i]-Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//re4
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
	for ( i=0; i<Num; i++ )
	{
		re[5][i] += -i_*(k1*wave_impe1*permi)*(dfymn[i]*Tre[0][i]-Tre[1][i]);
	}
	for ( i=0; i<2; i++ )
	{
		Tre[i].clear();
	}

	//////////////////////////////////////////////////////////////////////////
	
	for ( i=0; i<5; i++ )
	{
		Tvec[i].clear();
	}
}
void resolve_matrix::Two_D_Diel_Multi_Zfs_Vec(char a, Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6])
{
	int i;
	int Num = nxy*nxy;
	for ( i=0; i<6; i++ )
	{
		re[i].assign(Num);
	}

	Vec_Com Tre[7],Tvec[6];
	for ( i=0; i<6; i++ )
	{
		Tvec[i].assign(Num);
	}

	if ( a == ' ' )
	{
		Two_D_Diel_Multi_Zw_Vec(0,0,xmn,ymn,zmn,dfxmn,dfymn,vec,k1,k2,wave_impe1,permi,dxy,nxy,rd_num,re);
	}
	else
	{
		//////////////////////////////////////////////////////////////////////////
		//第一行
		//re0_0
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[0][i]*dfymn[i]*ymn[i];
			Tvec[1][i] = vec[0][i]*dfymn[i];
			Tvec[2][i] = vec[0][i]*xmn[i];
			Tvec[3][i] = vec[0][i]*zmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] = (Tre[1][i]-ymn[i]*Tre[2][i])+dfxmn[i]*(Tre[3][i]-xmn[i]*Tre[0][i])-(Tre[4][i]-zmn[i]*Tre[0][i]);
		}
		for ( i=0; i<5; i++ )
		{
			Tre[i].clear();
		}
		
		//re0_1
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[1][i]*ymn[i];
			Tvec[1][i] = vec[1][i]*dfxmn[i]*ymn[i];
			Tvec[2][i] = vec[1][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] += dfxmn[i]*(Tre[1][i]-ymn[i]*Tre[0][i])-(Tre[2][i]-ymn[i]*Tre[3][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}
		
		//re0_2
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[2][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] += conj(i_)*(k1*wave_impe1)*(Tre[1][i]-dfxmn[i]*Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//re0_3
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[3][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[0]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] += conj(i_)*(k1*wave_impe1)*dfxmn[i]*Tre[0][i];//原来有错误：dfymn写成了dfxmn
		}
		for ( i=0; i<1; i++ )
		{
			Tre[i].clear();
		}

		//re0_4
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[4][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] += conj(i_)*(k1*wave_impe1)*(-Tre[0][i]-dfxmn[i]*Tre[1][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//re0_5
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[5][i]*dfxmn[i]*ymn[i];
			Tvec[1][i] = vec[5][i]*dfxmn[i];
			Tvec[2][i] = vec[5][i]*dfymn[i]*xmn[i];
			Tvec[3][i] = vec[5][i]*dfymn[i];
			Tvec[4][i] = vec[5][i]*dfymn[i]*zmn[i];
			Tvec[5][i] = vec[5][i]*ymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,Tre[5]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[5],k2,dxy,nxy,rd_num,Tre[6]);
		for ( i=0; i<Num; i++ )
		{
			re[0][i] += -dfxmn[i]*(Tre[1][i]-ymn[i]*Tre[2][i])+(dfxmn[i]*(Tre[3][i]-xmn[i]*Tre[4][i])-(Tre[5][i]-zmn[i]*Tre[4][i]))-(Tre[6][i]-ymn[i]*Tre[0][i]);
		}
		for ( i=0; i<7; i++ )
		{
			Tre[i].clear();
		}

		//第二行
		//re1_0
		for (i=0; i<Num; i++)
		{
			Tvec[0][i] = vec[0][i]*dfymn[i]*xmn[i];
			Tvec[1][i] = vec[0][i]*dfymn[i];
			Tvec[2][i] = vec[0][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] = -(Tre[1][i]-xmn[i]*Tre[2][i])+dfymn[i]*(Tre[3][i]-xmn[i]*Tre[0][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re1_1
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[1][i]*zmn[i];
			Tvec[1][i] = vec[1][i]*ymn[i];
			Tvec[2][i] = vec[1][i]*dfxmn[i]*xmn[i];
			Tvec[3][i] = vec[1][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] += -(Tre[1][i]-zmn[i]*Tre[0][i])+dfymn[i]*(Tre[2][i]-ymn[i]*Tre[0][i])+(Tre[3][i]-xmn[i]*Tre[4][i]);
		}
		for ( i=0; i<5; i++ )
		{
			Tre[i].clear();
		}	

		//re1_2
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[2][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] += conj(i_)*(k1*wave_impe1)*(Tre[1][i]-dfymn[i]*Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}	

		//re1_3
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[3][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] += conj(i_)*(k1*wave_impe1)*(dfymn[i]*Tre[1][i]+Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}	

		//re1_4
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[4][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[0]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] += conj(i_)*(k1*wave_impe1)*(-dfymn[i]*Tre[0][i]);
		}
		for ( i=0; i<1; i++ )
		{
			Tre[i].clear();
		}

		//re1_5
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[5][i]*dfxmn[i]*zmn[i];
			Tvec[1][i] = vec[5][i]*dfxmn[i];
			Tvec[2][i] = vec[5][i]*dfymn[i]*xmn[i];
			Tvec[3][i] = vec[5][i]*dfymn[i];
			Tvec[4][i] = vec[5][i]*dfxmn[i]*ymn[i];
			Tvec[5][i] = vec[5][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,Tre[5]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[5],k2,dxy,nxy,rd_num,Tre[6]);
		for ( i=0; i<Num; i++ )
		{
			re[1][i] += (Tre[1][i]-zmn[i]*Tre[2][i])+dfymn[i]*(Tre[3][i]-xmn[i]*Tre[4][i])-dfymn[i]*(Tre[5][i]-ymn[i]*Tre[2][i])+(Tre[6][i]-xmn[i]*Tre[0][i]);
		}
		for ( i=0; i<7; i++ )
		{
			Tre[i].clear();
		}

		//第三行
		//re2_0
		for (i=0; i<Num; i++)
		{
			Tvec[0][i] = vec[0][i]*dfymn[i]*zmn[i];
			Tvec[1][i] = vec[0][i]*dfymn[i];
			Tvec[2][i] = vec[0][i]*ymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[2][i] = (Tre[1][i]-zmn[i]*Tre[2][i])+(Tre[3][i]-ymn[i]*Tre[0][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re2_1
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[1][i]*dfxmn[i]*zmn[i];
			Tvec[1][i] = vec[1][i]*dfxmn[i];
			Tvec[2][i] = vec[1][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[2][i] += -(Tre[1][i]-zmn[i]*Tre[2][i])-(Tre[3][i]-xmn[i]*Tre[0][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re2_5
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[5][i]*dfxmn[i]*xmn[i];
			Tvec[1][i] = vec[5][i]*dfxmn[i];
			Tvec[2][i] = vec[5][i]*dfymn[i]*ymn[i];
			Tvec[3][i] = vec[5][i]*dfymn[i];
			Tvec[4][i] = vec[5][i]*zmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,Tre[5]);
		for ( i=0; i<Num; i++ )
		{
			re[2][i] += (Tre[1][i]-xmn[i]*Tre[2][i])+(Tre[3][i]-ymn[i]*Tre[4][i])-(Tre[5][i]-zmn[i]*Tre[0][i]);
		}
		for ( i=0; i<6; i++ )
		{
			Tre[i].clear();
		}

		//第四行
		//re3_0
		for (i=0; i<Num; i++)
		{
			Tvec[0][i] = vec[0][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[0]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] = -conj(i_)*(k1/wave_impe1)*dfxmn[i]*Tre[0][i];
		}
		for ( i=0; i<1; i++ )
		{
			Tre[i].clear();
		}

		//re3_1
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[1][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] += -conj(i_)*(k1/wave_impe1)*((-1.-dfxmn[i])*Tre[0][i]-Tre[1][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//re3_2
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[2][i]*dfxmn[i]*ymn[i];
			Tvec[1][i] = vec[2][i]*dfxmn[i];
			Tvec[2][i] = vec[2][i]*dfymn[i]*xmn[i];
			Tvec[3][i] = vec[2][i]*dfymn[i];
			Tvec[4][i] = vec[2][i]*dfymn[i]*zmn[i];
			Tvec[5][i] = vec[2][i]*ymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,Tre[5]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[5],k1,dxy,nxy,rd_num,Tre[6]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] += -dfxmn[i]*(Tre[1][i]-ymn[i]*Tre[2][i])+(dfxmn[i]*(Tre[3][i]-xmn[i]*Tre[4][i])-(Tre[5][i]-zmn[i]*Tre[4][i]))-(Tre[6][i]-ymn[i]*Tre[0][i]);
		}
		for ( i=0; i<7; i++ )
		{
			Tre[i].clear();
		}

		//re3_3
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[3][i]*dfymn[i]*ymn[i];
			Tvec[1][i] = vec[3][i]*dfymn[i];
			Tvec[2][i] = vec[3][i]*xmn[i];
			Tvec[3][i] = vec[3][i]*zmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] += (Tre[1][i]-ymn[i]*Tre[2][i])+dfxmn[i]*(Tre[3][i]-xmn[i]*Tre[0][i])-(Tre[4][i]-zmn[i]*Tre[0][i]);
		}
		for ( i=0; i<5; i++ )
		{
			Tre[i].clear();
		}

		//re3_4
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[4][i]*ymn[i];
			Tvec[1][i] = vec[4][i]*dfxmn[i]*ymn[i];
			Tvec[2][i] = vec[4][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] += dfxmn[i]*(Tre[1][i]-ymn[i]*Tre[0][i])-(Tre[2][i]-ymn[i]*Tre[3][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re3_5
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[5][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[3][i] += -conj(i_)*(k1*wave_impe1)*conj(permi)*(Tre[1][i]-dfxmn[i]*Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//第五行
		//re4_0
		for (i=0; i<Num; i++)
		{
			Tvec[0][i] = vec[0][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[4][i] = -conj(i_)*(k1/wave_impe1)*(dfymn[i]*Tre[1][i]+Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//re4_1
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[1][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[0]);
		for ( i=0; i<Num; i++ )
		{
			re[4][i] += conj(i_)*(k1/wave_impe1)*dfymn[i]*Tre[0][i];
		}
		for ( i=0; i<1; i++ )
		{
			Tre[i].clear();
		}

		//re4_2
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[2][i]*dfxmn[i]*zmn[i];
			Tvec[1][i] = vec[2][i]*dfxmn[i];
			Tvec[2][i] = vec[2][i]*dfxmn[i]*ymn[i];
			Tvec[3][i] = vec[2][i]*dfymn[i]*xmn[i];
			Tvec[4][i] = vec[2][i]*dfymn[i];
			Tvec[5][i] = vec[2][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,Tre[5]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[5],k1,dxy,nxy,rd_num,Tre[6]);		
		for ( i=0; i<Num; i++ )
		{
			re[4][i] += ((Tre[1][i]-zmn[i]*Tre[2][i])-dfymn[i]*(Tre[3][i]-ymn[i]*Tre[2][i]))+dfymn[i]*(Tre[4][i]-xmn[i]*Tre[5][i])+(Tre[6][i]-xmn[i]*Tre[0][i]);
		}
		for ( i=0; i<7; i++ )
		{
			Tre[i].clear();
		}

		//re4_3
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[3][i]*dfymn[i]*xmn[i];
			Tvec[1][i] = vec[3][i]*dfymn[i];
			Tvec[2][i] = vec[3][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);		
		for ( i=0; i<Num; i++ )
		{
			re[4][i] += -(Tre[1][i]-xmn[i]*Tre[2][i])+dfymn[i]*(Tre[3][i]-xmn[i]*Tre[0][i]);
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re4_4
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[4][i]*zmn[i];
			Tvec[1][i] = vec[4][i]*ymn[i];
			Tvec[2][i] = vec[4][i]*dfxmn[i]*xmn[i];
			Tvec[3][i] = vec[4][i]*dfxmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,Tre[4]);
		for ( i=0; i<Num; i++ )
		{
			re[4][i] += -(Tre[1][i]-zmn[i]*Tre[0][i])+dfymn[i]*(Tre[2][i]-ymn[i]*Tre[0][i])+(Tre[3][i]-xmn[i]*Tre[4][i]);
		}
		for ( i=0; i<5; i++ )
		{
			Tre[i].clear();
		}

		//re4_5
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[5][i]*dfymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_g_b_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		for ( i=0; i<Num; i++ )
		{
			re[4][i] += -conj(i_)*(k1*wave_impe1)*conj(permi)*(Tre[1][i]-dfymn[i]*Tre[0][i]);
		}
		for ( i=0; i<2; i++ )
		{
			Tre[i].clear();
		}

		//第六行
		//re5_2
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[2][i]*dfxmn[i]*xmn[i];
			Tvec[1][i] = vec[2][i]*dfxmn[i];
			Tvec[2][i] = vec[2][i]*dfymn[i]*ymn[i];
			Tvec[3][i] = vec[2][i]*dfymn[i];
			Tvec[4][i] = vec[2][i]*zmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,Tre[3]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,Tre[4]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,Tre[5]);
		for ( i=0; i<Num; i++ )
		{
			re[5][i] = (Tre[1][i]-xmn[i]*Tre[2][i])+(Tre[3][i]-ymn[i]*Tre[4][i])-(Tre[5][i]-zmn[i]*Tre[0][i]);
		}
		for ( i=0; i<6; i++ )
		{
			Tre[i].clear();
		}

		//re5_3
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[3][i]*dfymn[i]*zmn[i];
			Tvec[1][i] = vec[3][i]*dfymn[i];
			Tvec[2][i] = vec[3][i]*ymn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[5][i] += conj(1./permi)*((Tre[1][i]-zmn[i]*Tre[2][i])+(Tre[3][i]-ymn[i]*Tre[0][i]));
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//re5_4
		for ( i=0; i<Num; i++ )
		{
			Tvec[0][i] = vec[4][i]*dfxmn[i]*zmn[i];
			Tvec[1][i] = vec[4][i]*dfxmn[i];
			Tvec[2][i] = vec[4][i]*xmn[i];
		}
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,Tre[0]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,Tre[1]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,Tre[2]);
		Two_D_Diel_Multi_Zw_Vec_G_a_Child('*',0,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,Tre[3]);
		for ( i=0; i<Num; i++ )
		{
			re[5][i] += conj(1./permi)*(-(Tre[1][i]-zmn[i]*Tre[2][i])+(Tre[3][i]-xmn[i]*Tre[0][i]));
		}
		for ( i=0; i<4; i++ )
		{
			Tre[i].clear();
		}

		//////////////////////////////////////////////////////////////////////////

		for ( i=0; i<6; i++ )
		{
			Tvec[i].clear();
		}
	}	
}
void resolve_matrix::Two_D_Diel_Multi_Zw_Vector(int Para_G_a_Num,int Para_g_b_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com vec[6],double k1,Complex k2,double wave_impe1,Complex permi,double dxy,int nxy,int rd_num,Vec_Com re[6])
{
	int i;
	int Num = nxy*nxy;
	for ( i=0; i<6; i++ )
	{
		re[i].assign(Num);
	}

	Vec_Com Tvec[5],Ti[6];
	for ( i=0; i<5; i++ )
	{
		Tvec[i].assign(Num);
	}
	for ( i=0; i<6; i++ )
	{
		Ti[i].assign(Num);
	}

	//////////////////////////////////////////////////////////////////////////
	//自由空间中的积分方程
	//第一个矩阵方程
	//re0
	Vec_Com re000,re001,re002,re003,re004;
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[0][i]*ymn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i];
		Tvec[2][i] = vec[0][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[0][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,re000);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re001);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re002);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re003);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,re004);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = dfymn[i]*(ymn[i]*re000[i]-re001[i])+(xmn[i]*re002[i]-re003[i])-(zmn[i]*re000[i]-re004[i]);
	}
	re000.clear();
	re001.clear();
	re002.clear();
	re003.clear();
	re004.clear();

	//re1
	Vec_Com re010,re011,re012,re013;
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[1][i]*xmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,re010);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re011);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re012);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re013);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = -dfymn[i]*(xmn[i]*re010[i]-re011[i])+(xmn[i]*re012[i]-re013[i]);
	}
	re010.clear();
	re011.clear();
	re012.clear();
	re013.clear();

	//re2
	Vec_Com re020,re021,re022;
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[2][i]*zmn[i];
		Tvec[1][i] = vec[2][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,re020);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re021);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re022);
	for ( i=0; i<Num; i++ )
	{
		Ti[2][i] = dfymn[i]*(zmn[i]*re020[i]-re021[i])+(ymn[i]*re020[i]-re022[i]);
	}
	re020.clear();
	re021.clear();
	re022.clear();

	//re3
	Vec_Com re030;
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re030);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = -i_*(k1/wave_impe1)*dfymn[i]*re030[i];
	}
	re030.clear();

	//re4
	Vec_Com re040,re041;
	for (i=0; i<Num; i++)
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[4],k1,dxy,nxy,rd_num,re040);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re041);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = -i_*(k1/wave_impe1)*(dfymn[i]*re041[i]+re040[i]);
	}
	re040.clear();
	re041.clear();

	for ( i=0; i<Num; i++ )
	{
		re[0][i] = Ti[0][i] + Ti[1][i] + Ti[2][i] + Ti[3][i] + Ti[4][i];
	}

	//第二个积分方程
	//re0
	Vec_Com re100,re101,re102,re103;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[0][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,re100);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re101);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re102);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re103);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = (ymn[i]*re101[i]-re102[i])-dfxmn[i]*(ymn[i]*re100[i]-re103[i]);
	}
	re100.clear();
	re101.clear();
	re102.clear();
	re103.clear();

	//re1
	Vec_Com re110,re111,re112,re113,re114;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*zmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[1][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,re110);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re111);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re112);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re113);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,re114);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = -(zmn[i]*re110[i]-re111[i])+(ymn[i]*re112[i]-re113[i])+dfxmn[i]*(xmn[i]*re110[i]-re114[i]);
	}
	re110.clear();
	re111.clear();
	re112.clear();
	re113.clear();
	re114.clear();

	//re2
	Vec_Com re120,re121,re122;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[2][i]*zmn[i];
		Tvec[1][i] = vec[2][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k1,dxy,nxy,rd_num,re120);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re121);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re122);
	for ( i=0; i<Num; i++ )
	{
		Ti[2][i] = -dfxmn[i]*(zmn[i]*re120[i]-re121[i])-(xmn[i]*re120[i]-re122[i]);
	}
	re120.clear();
	re121.clear();
	re122.clear();

	//re3
	Vec_Com re130,re131;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[3],k1,dxy,nxy,rd_num,re130);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re131);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = -i_*(k1/wave_impe1)*((-1.-dfxmn[i])*re130[i]-re131[i]);
	}
	re130.clear();
	re131.clear();

	//re4
	Vec_Com re140;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re140);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = i_*(k1/wave_impe1)*dfxmn[i]*re140[i];
	}
	re140.clear();

	for ( i=0; i<Num; i++ )
	{
		re[1][i] = Ti[0][i] + Ti[1][i] + Ti[2][i] + Ti[3][i] + Ti[4][i];
	}

	//第三个积分方程
	//re0
	Vec_Com re200,re201;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[0],k1,dxy,nxy,rd_num,re200);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re201);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = i_*(k1*wave_impe1)*(dfxmn[i]*re200[i]-re201[i]);
	}
	re200.clear();
	re201.clear();

	//re1
	Vec_Com re210,re211;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[1],k1,dxy,nxy,rd_num,re210);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re211);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = i_*(k1*wave_impe1)*(dfymn[i]*re210[i]-re211[i]);
	}
	re210.clear();
	re211.clear();
	
	//re2
	Vec_Com re220,re221,re222,re223,re224,re225;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[3][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[3][i]*zmn[i];
		Tvec[4][i] = vec[3][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k1,dxy,nxy,rd_num,re220);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re221);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re222);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re223);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,re224);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,re225);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = -dfxmn[i]*(ymn[i]*re221[i]-re222[i])+dfymn[i]*((xmn[i]*re221[i]-re223[i])-(zmn[i]*re220[i]-re224[i]))-(ymn[i]*re220[i]-re225[i]);
	}
	re220.clear();
	re221.clear();
	re222.clear();
	re223.clear();
	re224.clear();
	re225.clear();

	//re4
	Vec_Com re240,re241,re242,re243,re244,re245;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*zmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[4][i]*dfymn[i]*xmn[i];
		Tvec[4][i] = vec[4][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k1,dxy,nxy,rd_num,re240);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re241);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re242);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re243);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k1,dxy,nxy,rd_num,re244);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k1,dxy,nxy,rd_num,re245);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = dfxmn[i]*((zmn[i]*re240[i]-re241[i])-(ymn[i]*re242[i]-re243[i]))+dfymn[i]*(xmn[i]*re242[i]-re244[i])+(xmn[i]*re240[i]-re245[i]);
	}
	re240.clear();
	re241.clear();
	re242.clear();
	re243.clear();
	re244.clear();
	re245.clear();

	//re5
	Vec_Com re250,re251,re252,re253;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*xmn[i];
		Tvec[1][i] = vec[5][i]*ymn[i];
		Tvec[2][i] = vec[5][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k1,dxy,nxy,rd_num,re250);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k1,dxy,nxy,rd_num,re251);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k1,dxy,nxy,rd_num,re252);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k1,dxy,nxy,rd_num,re253);
	for ( i=0; i<Num; i++ )
	{
		Ti[5][i] = dfxmn[i]*(xmn[i]*re250[i]-re251[i])+dfymn[i]*(ymn[i]*re250[i]-re252[i])-(zmn[i]*re250[i]-re253[i]);
	}
	re250.clear();
	re251.clear();
	re252.clear();
	re253.clear();
	
	for ( i=0; i<Num; i++ )
	{
		re[2][i] = Ti[0][i] + Ti[1][i] + Ti[3][i] + Ti[4][i] + Ti[5][i];
	}
	//////////////////////////////////////////////////////////////////////////
	//介质空间中的积分方程
	//第四个积分方程
	//re0
	Vec_Com re300;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re300);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = i_*(k1*wave_impe1)*dfxmn[i]*re300[i];
	}
	re300.clear();

	//re1
	Vec_Com re310,re311;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[1],k2,dxy,nxy,rd_num,re310);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re311);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = i_*(k1*wave_impe1)*(dfymn[i]*re311[i]+re310[i]);
	}
	re310.clear();
	re311.clear();

	//re3
	Vec_Com re330,re331,re332,re333,re334;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*ymn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i];
		Tvec[2][i] = vec[3][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[3][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,re330);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re331);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re332);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re333);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,re334);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = dfymn[i]*(ymn[i]*re330[i]-re331[i])+(xmn[i]*re332[i]-re333[i])-(zmn[i]*re330[i]-re334[i]);
	}
	re330.clear();
	re331.clear();
	re332.clear();
	re333.clear();
	re334.clear();

	//re4
	Vec_Com re340,re341,re342,re343;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*xmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,re340);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re341);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re342);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re343);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = -dfymn[i]*(xmn[i]*re340[i]-re341[i])+(xmn[i]*re342[i]-re343[i]);
	}
	re340.clear();
	re341.clear();
	re342.clear();
	re343.clear();

	//re5
	Vec_Com re350,re351,re352;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*zmn[i];
		Tvec[1][i] = vec[5][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,re350);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re351);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re352);
	for ( i=0; i<Num; i++ )
	{
		Ti[5][i] = (1./permi)*(dfymn[i]*(zmn[i]*re350[i]-re351[i])+(ymn[i]*re350[i]-re352[i]));
	}
	re350.clear();
	re351.clear();
	re352.clear();
	
	for ( i=0; i<Num; i++ )
	{
		re[3][i] = Ti[0][i] + Ti[1][i] + Ti[3][i] + Ti[4][i] + Ti[5][i];
	}

	//第五个积分方程
	//re0
	Vec_Com re400,re401;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,vec[0],k2,dxy,nxy,rd_num,re400);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re401);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = i_*(k1*wave_impe1)*(-re400[i]-dfxmn[i]*re401[i]);
	}
	re400.clear();
	re401.clear();

	//re1
	Vec_Com re410;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_g_b_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re410);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = i_*(k1*wave_impe1)*(-dfxmn[i]*re410[i]);
	}
	re410.clear();

	//re3
	Vec_Com re430,re431,re432,re433;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
		Tvec[1][i] = vec[3][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[3][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,re430);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re431);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re432);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re433);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = (ymn[i]*re431[i]-re432[i])-dfxmn[i]*(ymn[i]*re430[i]-re433[i]);
	}
	re430.clear();
	re431.clear();
	re432.clear();
	re433.clear();

	//re4
	Vec_Com re440,re441,re442,re443,re444;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*zmn[i];
		Tvec[1][i] = vec[4][i]*dfymn[i];
		Tvec[2][i] = vec[4][i]*dfymn[i]*ymn[i];
		Tvec[3][i] = vec[4][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,re440);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re441);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re442);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re443);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,re444);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = -(zmn[i]*re440[i]-re441[i])+(ymn[i]*re442[i]-re443[i])+dfxmn[i]*(xmn[i]*re440[i]-re444[i]);
	}
	re440.clear();
	re441.clear();
	re442.clear();
	re443.clear();
	re444.clear();

	//re5
	Vec_Com re450,re451,re452;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[5][i]*zmn[i];
		Tvec[1][i] = vec[5][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[5],k2,dxy,nxy,rd_num,re450);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re451);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re452);
	for ( i=0; i<Num; i++ )
	{
		Ti[5][i] = (1./permi)*(-dfxmn[i]*(zmn[i]*re450[i]-re451[i])+(xmn[i]*re450[i]-re452[i]));
	}
	re450.clear();
	re451.clear();
	re452.clear();

	for ( i=0; i<Num; i++ )
	{
		re[4][i] = Ti[0][i] + Ti[1][i] + Ti[3][i] + Ti[4][i] + Ti[5][i];
	}

	//第六个积分方程
	//re0
	Vec_Com re500,re501,re502,re503,re504,re505;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[0][i]*dfxmn[i];
		Tvec[1][i] = vec[0][i]*dfxmn[i]*ymn[i];
		Tvec[2][i] = vec[0][i]*dfxmn[i]*xmn[i];
		Tvec[3][i] = vec[0][i]*zmn[i];
		Tvec[4][i] = vec[0][i]*ymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[0],k2,dxy,nxy,rd_num,re500);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re501);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re502);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re503);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,re504);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,re505);
	for ( i=0; i<Num; i++ )
	{
		Ti[0][i] = -dfxmn[i]*(ymn[i]*re501[i]-re502[i])+dfymn[i]*((xmn[i]*re501[i]-re503[i])-(zmn[i]*re500[i]-re504[i]))-(ymn[i]*re500[i]-re505[i]);
	}
	re500.clear();
	re501.clear();
	re502.clear();
	re503.clear();
	re504.clear();
	re505.clear();

	//re1
	Vec_Com re510,re511,re512,re513,re514,re515;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[1][i]*zmn[i];
		Tvec[1][i] = vec[1][i]*dfymn[i];
		Tvec[2][i] = vec[1][i]*dfymn[i]*xmn[i];
		Tvec[3][i] = vec[1][i]*dfymn[i]*ymn[i];
		Tvec[4][i] = vec[1][i]*xmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[1],k2,dxy,nxy,rd_num,re510);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re511);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re512);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re513);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[3],k2,dxy,nxy,rd_num,re514);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[4],k2,dxy,nxy,rd_num,re515);
	for ( i=0; i<Num; i++ )
	{
		Ti[1][i] = dfxmn[i]*(zmn[i]*re510[i]-re511[i])+dfymn[i]*(xmn[i]*re512[i]-re513[i])-dfxmn[i]*(ymn[i]*re512[i]-re514[i])+(xmn[i]*re510[i]-re515[i]);
	}
	re510.clear();
	re511.clear();
	re512.clear();
	re513.clear();
	re514.clear();
	re515.clear();

	//re2
	Vec_Com re520,re521,re522,re523;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[2][i]*xmn[i];
		Tvec[1][i] = vec[2][i]*ymn[i];
		Tvec[2][i] = vec[2][i]*zmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[2],k2,dxy,nxy,rd_num,re520);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re521);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[1],k2,dxy,nxy,rd_num,re522);
	Two_D_Diel_Multi_Zw_Vec_G_a_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[2],k2,dxy,nxy,rd_num,re523);
	for ( i=0; i<Num; i++ )
	{
		Ti[2][i] = dfxmn[i]*(xmn[i]*re520[i]-re521[i])+dfymn[i]*(ymn[i]*re520[i]-re522[i])-(zmn[i]*re520[i]-re523[i]);
	}
	re520.clear();
	re521.clear();
	re522.clear();
	re523.clear();

	//re3
	Vec_Com re530,re531;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[3][i]*dfxmn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[3],k2,dxy,nxy,rd_num,re530);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re531);
	for ( i=0; i<Num; i++ )
	{
		Ti[3][i] = -i_*(k1*wave_impe1)*(1./permi)*(dfxmn[i]*re530[i]-re531[i]);
	}
	re530.clear();
	re531.clear();

	//re4
	Vec_Com re540,re541;
	for ( i=0; i<Num; i++ )
	{
		Tvec[0][i] = vec[4][i]*dfymn[i];
	}
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,vec[4],k2,dxy,nxy,rd_num,re540);
	Two_D_Diel_Multi_Zw_Vec_g_b_Child(' ',Para_G_a_Num,xmn,ymn,zmn,Tvec[0],k2,dxy,nxy,rd_num,re541);
	for ( i=0; i<Num; i++ )
	{
		Ti[4][i] = -i_*(k1*wave_impe1)*(1./permi)*(dfymn[i]*re540[i]-re541[i]);
	}
	
	re540.clear();
	re541.clear();

	for ( i=0; i<Num; i++ )
	{
		re[5][i] = Ti[0][i] + Ti[1][i] + Ti[2][i] + Ti[3][i] + Ti[4][i];
	}

	//////////////////////////////////////////////////////////////////////////
	
	for ( i=0; i<5; i++ )
	{
		Tvec[i].clear();
	}
	for ( i=0; i<6; i++ )
	{
		Ti[i].clear();
	}

}
void resolve_matrix::Multi_vec_vec(Vec_Com &va, Vec_Com &vb, Vec_Com &re )
{
	int num = va.size();
	re.assign(num);
	int i;
	for ( i=0; i<num; i++ )
	{
		re[i] = va[i] * vb[i];
	}
}
void resolve_matrix::Multi_vec_num(Vec_Com &va, Complex nb, Vec_Com &re )
{
	int num = va.size();
	re.assign(num);
	int i;
	for ( i=0; i<num; i++ )
	{
		re[i] = va[i] * nb;
	}
}
void resolve_matrix::Two_D_Diel_Multi_Zw_Vec_G_a_Child(char a, int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex kk,double dxy,int nxy,int rd_num,Vec_Com &re)
{
	int i;
	int Num = nxy*nxy;//总元素个数
	re.assign(Num);

	switch( Para_Num )
	{
	case 0:
		{//只有平面矩阵，在其它中，都没有加平面矩阵
			Vec_Com Bccb0;
			Bccb0.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					if ( a == ' ' )
					{
						Bccb0[i] = Para_Two_D_Dielectric_G_a(0,rou,kk);
					}
					else
					{
						Bccb0[i] = conj(Para_Two_D_Dielectric_G_a(0,rou,kk));
					}					
				}
			}
			
			//Zw1
			Multi_BTTB_Vec(Bccb0,vec,nxy,re);
			for ( i=0; i<Num; i++ )
			{
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb0.clear();
			
			break;
		}
	case 1://以下都没有加平面矩阵
		{
			Vec_Com Bccb1;
			Bccb1.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
				}
			}
			
			Vec_Com fu1,fu2;
			fu1.assign(Num);
			fu2.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i];
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			
			fu1.clear();
			fu2.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			
			break;
		}
	case 2:
		{
			Vec_Com Bccb1,Bccb2;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_G_a(2,rou,kk)/pow(rou,4);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i];

				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();

			break;
		}
	case 3:
		{
			Vec_Com Bccb1,Bccb2,Bccb3;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_G_a(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_G_a(3,rou,kk)/pow(rou,6);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i];
				
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();

			break;
		}
	case 4:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_G_a(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_G_a(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_G_a(4,rou,kk)/pow(rou,8);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
					fu8[i]=fu7[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
				
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
				
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);	
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
				
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
						
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i];

				re[i] = dxy*dxy* re[i];
			}
				
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
				
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
				
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			
			break;
		}
	case 5:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			double rou;
					
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_G_a(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_G_a(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_G_a(4,rou,kk)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Dielectric_G_a(5,rou,kk)/pow(rou,10);
				}
			}
					
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
			}
					
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
					
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
					
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
					
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
					
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
					
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
							
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
							
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i];
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
					
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
					
			re11.clear();
			re12.clear();
			re13.clear();	
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();

			break;
		}
	case 6:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5,Bccb6;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			Bccb6.assign(Num);
			double rou;
						
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_G_a(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_G_a(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_G_a(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_G_a(4,rou,kk)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Dielectric_G_a(5,rou,kk)/pow(rou,10);
					Bccb6[i] = Para_Two_D_Dielectric_G_a(6,rou,kk)/pow(rou,12);
				}
			}
						
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10,fu11,fu12;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
			fu11.assign(Num);
			fu12.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
				fu11[i]=fu10[i]*fmn[i];
				fu12[i]=fu11[i]*fmn[i];
			}
						
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511,
				re61,re62,re63,re64,re65,re66,re67,re68,re69,re610,re611,re612,re613;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
						
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
						
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
						
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
						
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
						
			//Zw6
			Multi_BTTB_Vec(Bccb6,vec,nxy,re61);
			Multi_BTTB_Vec(Bccb6,fu1,nxy,re62);
			Multi_BTTB_Vec(Bccb6,fu2,nxy,re63);
			Multi_BTTB_Vec(Bccb6,fu3,nxy,re64);
			Multi_BTTB_Vec(Bccb6,fu4,nxy,re65);
			Multi_BTTB_Vec(Bccb6,fu5,nxy,re66);
			Multi_BTTB_Vec(Bccb6,fu6,nxy,re67);
			Multi_BTTB_Vec(Bccb6,fu7,nxy,re68);
			Multi_BTTB_Vec(Bccb6,fu8,nxy,re69);
			Multi_BTTB_Vec(Bccb6,fu9,nxy,re610);
			Multi_BTTB_Vec(Bccb6,fu10,nxy,re611);
			Multi_BTTB_Vec(Bccb6,fu11,nxy,re612);
			Multi_BTTB_Vec(Bccb6,fu12,nxy,re613);
						
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
								
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
								
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
								
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
						
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i]
							
					+pow(fmn[i],12)*re61[i]-12.0*pow(fmn[i],11)*re62[i]+66.0*pow(fmn[i],10)*re63[i]
					-220.0*pow(fmn[i],9)*re64[i]+495.0*pow(fmn[i],8)*re65[i]-792.0*pow(fmn[i],7)*re66[i]
					+924.0*pow(fmn[i],6)*re67[i]-790.0*pow(fmn[i],5)*re68[i]+487.0*pow(fmn[i],4)*re69[i]
					-208.0*pow(fmn[i],3)*re610[i]+58.0*pow(fmn[i],2)*re611[i]-8.0*fmn[i]*re612[i]+re613[i];

					re[i] = dxy*dxy* re[i];
			}
						
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
			Bccb6.clear();
						
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
			fu11.clear();
			fu12.clear();
						
			re11.clear();
			re12.clear();
			re13.clear();
						
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
						
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
					
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();
						
			re61.clear();
			re62.clear();
			re63.clear();
			re64.clear();
			re65.clear();
			re66.clear();
			re67.clear();
			re68.clear();
			re69.clear();
			re610.clear();
			re611.clear();
			re612.clear();
			re613.clear();
						
			break;
		}
	}
}
void resolve_matrix::Two_D_Diel_Multi_Zw_Vec_g_b_Child(char a, int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex kk,double dxy,int nxy,int rd_num,Vec_Com &re)
{
	int i;
	int Num = nxy*nxy;//总元素个数
	re.assign(Num);

	switch( Para_Num ) 
	{
	case 0:
		{//只有平面矩阵，在其它中，都没有加平面矩阵
			Vec_Com Bccb0;
			Bccb0.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					if ( a == ' ' )
					{
						Bccb0[i] = Para_Two_D_Dielectric_g_b(0,rou,kk);
					}
					else
					{
						Bccb0[i] = conj(Para_Two_D_Dielectric_g_b(0,rou,kk));
					}					
				}
			}	
			
			//Zw1
			Multi_BTTB_Vec(Bccb0,vec,nxy,re);
			for ( i=0; i<Num; i++ )
			{
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb0.clear();
			
			break;
		}
	case 1:
		{
			Vec_Com Bccb1;
			Bccb1.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
				}
			}
			
			Vec_Com fu1,fu2;
			fu1.assign(Num);
			fu2.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i];
				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			
			fu1.clear();
			fu2.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			
			break;
		}
	case 2:
		{
			Vec_Com Bccb1,Bccb2;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_g_b(2,rou,kk)/pow(rou,4);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i];

				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();

			break;
		}
	case 3:
		{
			Vec_Com Bccb1,Bccb2,Bccb3;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_g_b(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_g_b(3,rou,kk)/pow(rou,6);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i];

				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();

			break;
		}
	case 4:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			double rou;
				
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_g_b(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_g_b(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_g_b(4,rou,kk)/pow(rou,8);
				}
			}
				
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
				
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
			}
				
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49;
				
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
				
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
				
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);	
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
				
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
						
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i];

				re[i]=dxy*dxy*re[i];
			}
				
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
				
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
				
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();

			break;
		}
	case 5:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			double rou;
					
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_g_b(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_g_b(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_g_b(4,rou,kk)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Dielectric_g_b(5,rou,kk)/pow(rou,10);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
			}
					
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
					
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
					
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
					
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
					
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
					
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
							
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
							
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i];

					re[i]=dxy*dxy*re[i];
			}
					
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
					
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
					
			re11.clear();
			re12.clear();
			re13.clear();	
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();

			break;
		}
	case 6:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5,Bccb6;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			Bccb6.assign(Num);
			double rou;
						
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Dielectric_g_b(1,rou,kk)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Dielectric_g_b(2,rou,kk)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Dielectric_g_b(3,rou,kk)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Dielectric_g_b(4,rou,kk)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Dielectric_g_b(5,rou,kk)/pow(rou,10);
					Bccb6[i] = Para_Two_D_Dielectric_g_b(6,rou,kk)/pow(rou,12);
				}
			}
						
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10,fu11,fu12;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
			fu11.assign(Num);
			fu12.assign(Num);
						
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
				fu11[i]=fu10[i]*fmn[i];
				fu12[i]=fu11[i]*fmn[i];
			}
						
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511,
				re61,re62,re63,re64,re65,re66,re67,re68,re69,re610,re611,re612,re613;
						
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
						
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
						
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
						
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
						
			//Zw6
			Multi_BTTB_Vec(Bccb6,vec,nxy,re61);
			Multi_BTTB_Vec(Bccb6,fu1,nxy,re62);
			Multi_BTTB_Vec(Bccb6,fu2,nxy,re63);
			Multi_BTTB_Vec(Bccb6,fu3,nxy,re64);
			Multi_BTTB_Vec(Bccb6,fu4,nxy,re65);
			Multi_BTTB_Vec(Bccb6,fu5,nxy,re66);
			Multi_BTTB_Vec(Bccb6,fu6,nxy,re67);
			Multi_BTTB_Vec(Bccb6,fu7,nxy,re68);
			Multi_BTTB_Vec(Bccb6,fu8,nxy,re69);
			Multi_BTTB_Vec(Bccb6,fu9,nxy,re610);
			Multi_BTTB_Vec(Bccb6,fu10,nxy,re611);
			Multi_BTTB_Vec(Bccb6,fu11,nxy,re612);
			Multi_BTTB_Vec(Bccb6,fu12,nxy,re613);
						
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
								
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
								
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
								
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i]
								
					+pow(fmn[i],12)*re61[i]-12.0*pow(fmn[i],11)*re62[i]+66.0*pow(fmn[i],10)*re63[i]
					-220.0*pow(fmn[i],9)*re64[i]+495.0*pow(fmn[i],8)*re65[i]-792.0*pow(fmn[i],7)*re66[i]
					+924.0*pow(fmn[i],6)*re67[i]-790.0*pow(fmn[i],5)*re68[i]+487.0*pow(fmn[i],4)*re69[i]
					-208.0*pow(fmn[i],3)*re610[i]+58.0*pow(fmn[i],2)*re611[i]-8.0*fmn[i]*re612[i]+re613[i];

  				re[i]=dxy*dxy*re[i];
			}
						
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
			Bccb6.clear();
						
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
			fu11.clear();
			fu12.clear();
						
			re11.clear();
			re12.clear();
			re13.clear();
				
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
						
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
						
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
						
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();
						
			re61.clear();
			re62.clear();
			re63.clear();
			re64.clear();
			re65.clear();
			re66.clear();
			re67.clear();
			re68.clear();
			re69.clear();
			re610.clear();
			re611.clear();
			re612.clear();
			re613.clear();
	
			break;
		}
	}
}
Complex resolve_matrix::Para_Two_D_Dielectric_G_a(int Bool,double R,Complex k)
{
	return ( - Para_Two_D_Conductor_Vector(Bool,R,k) );
}
Complex resolve_matrix::Para_Two_D_Dielectric_g_b(int Bool,double R,Complex k)
{
	return ( Para_Two_D_Conductor_Scalar(Bool,R,k) );
}
//////////////////////////////////////////////////////////////////////////
//
Complex resolve_matrix::G_0_fs_R(double R,Complex k)
{
	Complex re;
	re=exp(i_*k*R)/(4.0*PI*R);
	return re;
}
//标量波入射二维导体粗糙面时的弱矩阵向量积
void resolve_matrix::Two_D_Conductor_Scalar_Multi_Zw_Vec(int Para_Num, Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re)
{
	int i;
	int Num=nxy*nxy;//总元素个数
	re.assign(Num);

	switch(Para_Num) 
	{
	case 1:
		{
			Vec_Com Bccb1;
			Bccb1.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
				}
			}
			
			Vec_Com fu1,fu2;
			fu1.assign(Num);
			fu2.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i];
				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			
			fu1.clear();
			fu2.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			
			break;
		}
	case 2:
		{
			Vec_Com Bccb1,Bccb2;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Scalar(2,rou,k)/pow(rou,4);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i];

				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();

			break;
		}
	case 3:
		{
			Vec_Com Bccb1,Bccb2,Bccb3;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Scalar(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Scalar(3,rou,k)/pow(rou,6);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i];

				re[i]=dxy*dxy*re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();

			break;
		}
	case 4:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			double rou;
				
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Scalar(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Scalar(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Scalar(4,rou,k)/pow(rou,8);
				}
			}
				
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
				
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
			}
				
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49;
				
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
				
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
				
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);	
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
				
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
						
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i];

				re[i]=dxy*dxy*re[i];
			}
				
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
				
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
				
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();

			break;
		}
	case 5:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			double rou;
					
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Scalar(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Scalar(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Scalar(4,rou,k)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Conductor_Scalar(5,rou,k)/pow(rou,10);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
			}
					
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
					
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
					
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
					
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
					
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
					
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
							
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
							
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i];

					re[i]=dxy*dxy*re[i];
			}
					
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
					
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
					
			re11.clear();
			re12.clear();
			re13.clear();	
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();

			break;
		}
	case 6:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5,Bccb6;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			Bccb6.assign(Num);
			double rou;
						
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Scalar(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Scalar(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Scalar(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Scalar(4,rou,k)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Conductor_Scalar(5,rou,k)/pow(rou,10);
					Bccb6[i] = Para_Two_D_Conductor_Scalar(6,rou,k)/pow(rou,12);
				}
			}
						
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10,fu11,fu12;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
			fu11.assign(Num);
			fu12.assign(Num);
						
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
				fu11[i]=fu10[i]*fmn[i];
				fu12[i]=fu11[i]*fmn[i];
			}
						
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511,
				re61,re62,re63,re64,re65,re66,re67,re68,re69,re610,re611,re612,re613;
						
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
						
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
						
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
						
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
						
			//Zw6
			Multi_BTTB_Vec(Bccb6,vec,nxy,re61);
			Multi_BTTB_Vec(Bccb6,fu1,nxy,re62);
			Multi_BTTB_Vec(Bccb6,fu2,nxy,re63);
			Multi_BTTB_Vec(Bccb6,fu3,nxy,re64);
			Multi_BTTB_Vec(Bccb6,fu4,nxy,re65);
			Multi_BTTB_Vec(Bccb6,fu5,nxy,re66);
			Multi_BTTB_Vec(Bccb6,fu6,nxy,re67);
			Multi_BTTB_Vec(Bccb6,fu7,nxy,re68);
			Multi_BTTB_Vec(Bccb6,fu8,nxy,re69);
			Multi_BTTB_Vec(Bccb6,fu9,nxy,re610);
			Multi_BTTB_Vec(Bccb6,fu10,nxy,re611);
			Multi_BTTB_Vec(Bccb6,fu11,nxy,re612);
			Multi_BTTB_Vec(Bccb6,fu12,nxy,re613);
						
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
								
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
								
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
								
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i]
								
					+pow(fmn[i],12)*re61[i]-12.0*pow(fmn[i],11)*re62[i]+66.0*pow(fmn[i],10)*re63[i]
					-220.0*pow(fmn[i],9)*re64[i]+495.0*pow(fmn[i],8)*re65[i]-792.0*pow(fmn[i],7)*re66[i]
					+924.0*pow(fmn[i],6)*re67[i]-790.0*pow(fmn[i],5)*re68[i]+487.0*pow(fmn[i],4)*re69[i]
					-208.0*pow(fmn[i],3)*re610[i]+58.0*pow(fmn[i],2)*re611[i]-8.0*fmn[i]*re612[i]+re613[i];

  				re[i]=dxy*dxy*re[i];
			}
						
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
			Bccb6.clear();
						
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
			fu11.clear();
			fu12.clear();
						
			re11.clear();
			re12.clear();
			re13.clear();
				
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
						
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
						
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
						
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();
						
			re61.clear();
			re62.clear();
			re63.clear();
			re64.clear();
			re65.clear();
			re66.clear();
			re67.clear();
			re68.clear();
			re69.clear();
			re610.clear();
			re611.clear();
			re612.clear();
			re613.clear();
	
			break;
		}
	}
}

//二维标量波入射时格林函数g1,2(R)的泰勒展开系数
Complex resolve_matrix::Para_Two_D_Conductor_Scalar(int Bool,double R,Complex k)
{
	Complex re;
	Complex t0=i_*k*R;
	switch(Bool)
	{
	case 0:
		{
			re = G_0_fs_R(R,k);
			break;
		}
	case 1:
		{
			re = (t0-1.0)*G_0_fs_R(R,k)/2.0;
			break;
		}
	case 2:
		{
			re = (pow(t0,2)-3.0*t0+3.0)*G_0_fs_R(R,k)/8.0;
			break;
		}
	case 3:
		{
			re = (pow(t0,3)-6.0*pow(t0,2)+15.0*t0-15.0)*G_0_fs_R(R,k)/(16.0*6.0);
			break;
		}
	case 4:
		{
			re = (pow(t0,4)-10.0*pow(t0,3)+45.0*pow(t0,2)-105.0*t0+105.0)*G_0_fs_R(R,k)/(32.0*6.0*24.0);
			break;
		}
	case 5:
		{
			re = (pow(t0,5)-15.0*pow(t0,4)+105.0*pow(t0,3)-420.0*pow(t0,2)+945.0*t0-945.0)*G_0_fs_R(R,k)/(64.0*6.0*24.0*24.0*5.0);
			break;
		}
	case 6:
		{
			re = (pow(t0,6)-21.0*pow(t0,5)+210.0*pow(t0,4)-1260.0*pow(t0,3)+945.0*5.0*pow(t0,2)-945.0*11.0*t0+945.0*11.0)*G_0_fs_R(R,k)/(128.0*6.0*24.0*24.0*5.0*24.0*5.0*6.0);
			break;
		}
	}
	return re;
}
//////////////////////////////////////////////////////////////////////////
//矢量波入射二维导体粗糙面时的弱矩阵向量积
void resolve_matrix::Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(int Para_Num, Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn,Vec_Dou &dfymn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re)
{
	int i;
	int Num=nxy*nxy;//总元素个数
	re.assign(2.0*Num);
	Vec_Com Re_LU,Re_RU,Re_LD,Re_RD,Vec_UP,Vec_Down;
	Re_LU.assign(Num);
	Re_RU.assign(Num);
	Re_LD.assign(Num);
	Re_RD.assign(Num);
	Vec_UP.assign(Num);
	Vec_Down.assign(Num);

	for (i=0; i<Num; i++)
	{
		Vec_UP[i] = vec[i];
		Vec_Down[i] = vec[i+Num];
	}
	
	//左上
	Vec_Com vec11,vec12,vec13,vec14;
	vec11.assign(Num);
	vec12.assign(Num);
	vec13.assign(Num);
	vec14.assign(Num);
	for (i=0; i<Num; i++)
	{
		vec11[i] = Vec_UP[i]*ymn[i];
		vec12[i] = Vec_UP[i]*dfxmn[i];
		vec13[i] = Vec_UP[i]*xmn[i]*dfxmn[i];
		vec14[i] = Vec_UP[i]*fmn[i];
	}
	Vec_Com re10,re11,re12,re13,re14;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,Vec_UP,k,dxy,nxy,rd_num,re10);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec11,k,dxy,nxy,rd_num,re11);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec12,k,dxy,nxy,rd_num,re12);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec13,k,dxy,nxy,rd_num,re13);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec14,k,dxy,nxy,rd_num,re14);
	for (i=0; i<Num; i++)
	{
		Re_LU[i] = -dfymn[i]*(ymn[i]*re10[i]-re11[i])-(xmn[i]*re12[i]-re13[i])+(fmn[i]*re10[i]-re14[i]);
	}

	//右上
	Vec_Com vec21,vec22,vec23;
	vec21.assign(Num);
	vec22.assign(Num);
	vec23.assign(Num);
	for (i=0; i<Num; i++)
	{
		vec21[i] = Vec_Down[i]*xmn[i];
		vec22[i] = Vec_Down[i]*dfymn[i];
		vec23[i] = Vec_Down[i]*xmn[i]*dfymn[i];
	}
	Vec_Com re20,re21,re22,re23;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,Vec_Down,k,dxy,nxy,rd_num,re20);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec21,k,dxy,nxy,rd_num,re21);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec22,k,dxy,nxy,rd_num,re22);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec23,k,dxy,nxy,rd_num,re23);
	for (i=0; i<Num; i++)
	{
		Re_RU[i] = dfymn[i]*(xmn[i]*re20[i]-re21[i])-(xmn[i]*re22[i]-re23[i]);
	}

	//左下
	Vec_Com vec31,vec32,vec33;
	vec31.assign(Num);
	vec32.assign(Num);
	vec33.assign(Num);
	for (i=0; i<Num; i++)
	{
		vec31[i] = Vec_UP[i]*ymn[i];
		vec32[i] = Vec_UP[i]*dfxmn[i];
		vec33[i] = Vec_UP[i]*ymn[i]*dfxmn[i];
	}
	Vec_Com re30,re31,re32,re33;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,Vec_UP,k,dxy,nxy,rd_num,re30);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec31,k,dxy,nxy,rd_num,re31);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec32,k,dxy,nxy,rd_num,re32);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec33,k,dxy,nxy,rd_num,re33);
	for (i=0; i<Num; i++)
	{
		Re_LD[i] = dfxmn[i]*(ymn[i]*re30[i]-re31[i])-(ymn[i]*re32[i]-re33[i]);
	}

	//右下
	Vec_Com vec41,vec42,vec43,vec44;
	vec41.assign(Num);
	vec42.assign(Num);
	vec43.assign(Num);
	vec44.assign(Num);
	for (i=0; i<Num; i++)
	{
		vec41[i] = Vec_Down[i]*xmn[i];
		vec42[i] = Vec_Down[i]*dfymn[i];
		vec43[i] = Vec_Down[i]*ymn[i]*dfymn[i];
		vec44[i] = Vec_Down[i]*fmn[i];
	}
	Vec_Com re40,re41,re42,re43,re44;
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,Vec_Down,k,dxy,nxy,rd_num,re40);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec41,k,dxy,nxy,rd_num,re41);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec42,k,dxy,nxy,rd_num,re42);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec43,k,dxy,nxy,rd_num,re43);
	Two_D_Conductor_Vector_Multi_Zw_Vec_Child(Para_Num,xmn,ymn,fmn,vec44,k,dxy,nxy,rd_num,re44);
	for (i=0; i<Num; i++)
	{
		Re_RD[i] = -dfxmn[i]*(xmn[i]*re40[i]-re41[i])-(ymn[i]*re42[i]-re43[i])+(fmn[i]*re40[i]-re44[i]);
	}

	//总的列向量
	for (i=0; i<Num; i++)
	{
		re[i] = Re_LU[i] + Re_RU[i];
		re[i+Num] = Re_LD[i] + Re_RD[i];
	}

	Re_LU.clear();
	Re_RU.clear();
	Re_LD.clear();
	Re_RD.clear();
	Vec_UP.clear();
	Vec_Down.clear();
	vec11.clear();
	vec12.clear();
	vec13.clear();
	vec14.clear();
	re10.clear();
	re11.clear();
	re12.clear();
	re13.clear();
	re14.clear();
	vec21.clear();
	vec22.clear();
	vec23.clear();
	re20.clear();
	re21.clear();
	re22.clear();
	re23.clear();
	vec31.clear();
	vec32.clear();
	vec33.clear();
	re30.clear();
	re31.clear();
	re32.clear();
	re33.clear();
	vec41.clear();
	vec42.clear();
	vec43.clear();
	vec44.clear();
	re40.clear();
	re41.clear();
	re42.clear();
	re43.clear();
	re44.clear();

}
//矢量波入射导体粗糙面时弱矩阵向量积的子函数：弱矩阵与任意向量的乘积
void resolve_matrix::Two_D_Conductor_Vector_Multi_Zw_Vec_Child(int Para_Num,Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Com &vec,Complex k,double dxy,int nxy,int rd_num,Vec_Com &re)
{
    int i;
	int Num=nxy*nxy;//总元素个数
	re.assign(Num);

	switch(Para_Num) 
	{
	case 0:
		{//只有平面矩阵，在其它中，都没有加平面矩阵
			Vec_Com Bccb0;
			Bccb0.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb0[i] = Para_Two_D_Conductor_Vector(0,rou,k);
				}
			}	
			
			//Zw1
			Multi_BTTB_Vec(Bccb0,vec,nxy,re);
			for ( i=0; i<Num; i++ )
			{
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb0.clear();
			
			break;
		}
	case 1://以下都没有加平面矩阵
		{
			Vec_Com Bccb1;
			Bccb1.assign(Num);
			double rou;
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				//cout<<rou<<endl;
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
				}
			}
			
			Vec_Com fu1,fu2;
			fu1.assign(Num);
			fu2.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i];
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			
			fu1.clear();
			fu2.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			
			break;
		}
	case 2:
		{
			Vec_Com Bccb1,Bccb2;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Vector(2,rou,k)/pow(rou,4);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i];

				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();

			break;
		}
	case 3:
		{
			Vec_Com Bccb1,Bccb2,Bccb3;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Vector(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Vector(3,rou,k)/pow(rou,6);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
			
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
			
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
			
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
					
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i];
				
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();

			break;
		}
	case 4:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			double rou;
			
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Vector(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Vector(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Vector(4,rou,k)/pow(rou,8);
				}
			}
			
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
					fu8[i]=fu7[i]*fmn[i];
			}
			
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49;
			
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
				
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
				
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
				
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);	
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
				
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
						
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i];

				re[i] = dxy*dxy* re[i];
			}
				
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
				
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
				
			re11.clear();
			re12.clear();
			re13.clear();
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			
			break;
		}
	case 5:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			double rou;
					
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Vector(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Vector(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Vector(4,rou,k)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Conductor_Vector(5,rou,k)/pow(rou,10);
				}
			}
					
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
			}
					
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
					
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
					
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
					
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
					
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
					
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
							
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
						
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
							
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i];
				re[i] = dxy*dxy* re[i];
			}
			
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
					
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
					
			re11.clear();
			re12.clear();
			re13.clear();	
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();

			break;
		}
	case 6:
		{
			Vec_Com Bccb1,Bccb2,Bccb3,Bccb4,Bccb5,Bccb6;
			Bccb1.assign(Num);
			Bccb2.assign(Num);
			Bccb3.assign(Num);
			Bccb4.assign(Num);
			Bccb5.assign(Num);
			Bccb6.assign(Num);
			double rou;
						
			for (i=0; i<Num; i++)
			{
				rou=sqrt(pow(xmn[0]-xmn[i],2.0)+pow(ymn[0]-ymn[i],2.0));
				rou=rou/dxy;
				if(rou <= rd_num || fabs(rou-rd_num) < 1.0e-12){}
				else
				{
					rou = rou * dxy;
					Bccb1[i] = Para_Two_D_Conductor_Vector(1,rou,k)/pow(rou,2);
					Bccb2[i] = Para_Two_D_Conductor_Vector(2,rou,k)/pow(rou,4);
					Bccb3[i] = Para_Two_D_Conductor_Vector(3,rou,k)/pow(rou,6);
					Bccb4[i] = Para_Two_D_Conductor_Vector(4,rou,k)/pow(rou,8);
					Bccb5[i] = Para_Two_D_Conductor_Vector(5,rou,k)/pow(rou,10);
					Bccb6[i] = Para_Two_D_Conductor_Vector(6,rou,k)/pow(rou,12);
				}
			}
						
			Vec_Com fu1,fu2,fu3,fu4,fu5,fu6,fu7,fu8,fu9,fu10,fu11,fu12;
			fu1.assign(Num);
			fu2.assign(Num);
			fu3.assign(Num);
			fu4.assign(Num);
			fu5.assign(Num);
			fu6.assign(Num);
			fu7.assign(Num);
			fu8.assign(Num);
			fu9.assign(Num);
			fu10.assign(Num);
			fu11.assign(Num);
			fu12.assign(Num);
					
			for (i=0; i<Num; i++)
			{
				fu1[i]=fmn[i]*vec[i];
				fu2[i]=fu1[i]*fmn[i];
				fu3[i]=fu2[i]*fmn[i];
				fu4[i]=fu3[i]*fmn[i];
				fu5[i]=fu4[i]*fmn[i];
				fu6[i]=fu5[i]*fmn[i];
				fu7[i]=fu6[i]*fmn[i];
				fu8[i]=fu7[i]*fmn[i];
				fu9[i]=fu8[i]*fmn[i];
				fu10[i]=fu9[i]*fmn[i];
				fu11[i]=fu10[i]*fmn[i];
				fu12[i]=fu11[i]*fmn[i];
			}
						
			Vec_Com re11,re12,re13,
				re21,re22,re23,re24,re25,
				re31,re32,re33,re34,re35,re36,re37,
				re41,re42,re43,re44,re45,re46,re47,re48,re49,
				re51,re52,re53,re54,re55,re56,re57,re58,re59,re510,re511,
				re61,re62,re63,re64,re65,re66,re67,re68,re69,re610,re611,re612,re613;
					
			//Zw1
			Multi_BTTB_Vec(Bccb1,vec,nxy,re11);
			Multi_BTTB_Vec(Bccb1,fu1,nxy,re12);
			Multi_BTTB_Vec(Bccb1,fu2,nxy,re13);
						
			//Zw2
			Multi_BTTB_Vec(Bccb2,vec,nxy,re21);
			Multi_BTTB_Vec(Bccb2,fu1,nxy,re22);
			Multi_BTTB_Vec(Bccb2,fu2,nxy,re23);
			Multi_BTTB_Vec(Bccb2,fu3,nxy,re24);
			Multi_BTTB_Vec(Bccb2,fu4,nxy,re25);
						
			//Zw3
			Multi_BTTB_Vec(Bccb3,vec,nxy,re31);
			Multi_BTTB_Vec(Bccb3,fu1,nxy,re32);
			Multi_BTTB_Vec(Bccb3,fu2,nxy,re33);
			Multi_BTTB_Vec(Bccb3,fu3,nxy,re34);
			Multi_BTTB_Vec(Bccb3,fu4,nxy,re35);
			Multi_BTTB_Vec(Bccb3,fu5,nxy,re36);
			Multi_BTTB_Vec(Bccb3,fu6,nxy,re37);
						
			//Zw4
			Multi_BTTB_Vec(Bccb4,vec,nxy,re41);
			Multi_BTTB_Vec(Bccb4,fu1,nxy,re42);
			Multi_BTTB_Vec(Bccb4,fu2,nxy,re43);
			Multi_BTTB_Vec(Bccb4,fu3,nxy,re44);
			Multi_BTTB_Vec(Bccb4,fu4,nxy,re45);
			Multi_BTTB_Vec(Bccb4,fu5,nxy,re46);
			Multi_BTTB_Vec(Bccb4,fu6,nxy,re47);
			Multi_BTTB_Vec(Bccb4,fu7,nxy,re48);
			Multi_BTTB_Vec(Bccb4,fu8,nxy,re49);
						
			//Zw5
			Multi_BTTB_Vec(Bccb5,vec,nxy,re51);
			Multi_BTTB_Vec(Bccb5,fu1,nxy,re52);
			Multi_BTTB_Vec(Bccb5,fu2,nxy,re53);
			Multi_BTTB_Vec(Bccb5,fu3,nxy,re54);
			Multi_BTTB_Vec(Bccb5,fu4,nxy,re55);
			Multi_BTTB_Vec(Bccb5,fu5,nxy,re56);
			Multi_BTTB_Vec(Bccb5,fu6,nxy,re57);
			Multi_BTTB_Vec(Bccb5,fu7,nxy,re58);
			Multi_BTTB_Vec(Bccb5,fu8,nxy,re59);
			Multi_BTTB_Vec(Bccb5,fu9,nxy,re510);
			Multi_BTTB_Vec(Bccb5,fu10,nxy,re511);
						
			//Zw6
			Multi_BTTB_Vec(Bccb6,vec,nxy,re61);
			Multi_BTTB_Vec(Bccb6,fu1,nxy,re62);
			Multi_BTTB_Vec(Bccb6,fu2,nxy,re63);
			Multi_BTTB_Vec(Bccb6,fu3,nxy,re64);
			Multi_BTTB_Vec(Bccb6,fu4,nxy,re65);
			Multi_BTTB_Vec(Bccb6,fu5,nxy,re66);
			Multi_BTTB_Vec(Bccb6,fu6,nxy,re67);
			Multi_BTTB_Vec(Bccb6,fu7,nxy,re68);
			Multi_BTTB_Vec(Bccb6,fu8,nxy,re69);
			Multi_BTTB_Vec(Bccb6,fu9,nxy,re610);
			Multi_BTTB_Vec(Bccb6,fu10,nxy,re611);
			Multi_BTTB_Vec(Bccb6,fu11,nxy,re612);
			Multi_BTTB_Vec(Bccb6,fu12,nxy,re613);
						
			for (i=0; i<Num; i++)
			{
				re[i]=pow(fmn[i],2.0)*re11[i]-2.0*fmn[i]*re12[i]+re13[i]
								
					+pow(fmn[i],4.0)*re21[i]-4.0*pow(fmn[i],3.0)*re22[i]+6.0*pow(fmn[i],2.0)*re23[i]-4.0*fmn[i]*re24[i]+re25[i]
								
					+pow(fmn[i],6.0)*re31[i]-6.0*pow(fmn[i],5.0)*re32[i]+15.0*pow(fmn[i],4.0)*re33[i]
					-20.0*pow(fmn[i],3.0)*re34[i]+15.0*pow(fmn[i],2.0)*re35[i]-6.0*fmn[i]*re36[i]+re37[i]
								
					+pow(fmn[i],8.0)*re41[i]-8.0*pow(fmn[i],7.0)*re42[i]+28.0*pow(fmn[i],6.0)*re43[i]
					-56.0*pow(fmn[i],5.0)*re44[i]+70.0*pow(fmn[i],4.0)*re45[i]-56.0*pow(fmn[i],3.0)*re46[i]
					+28.0*pow(fmn[i],2.0)*re47[i]-6.0*fmn[i]*re48[i]+re49[i]
						
					+pow(fmn[i],10)*re51[i]-10.0*pow(fmn[i],9)*re52[i]+45.0*pow(fmn[i],8)*re53[i]
					-120.0*pow(fmn[i],7)*re54[i]+210.0*pow(fmn[i],6)*re55[i]-252.0*pow(fmn[i],5)*re56[i]
					+210.0*pow(fmn[i],4)*re57[i]-118.0*pow(fmn[i],3)*re58[i]+41.0*pow(fmn[i],2)*re59[i]
					-8.0*fmn[i]*re510[i]+re511[i]
							
					+pow(fmn[i],12)*re61[i]-12.0*pow(fmn[i],11)*re62[i]+66.0*pow(fmn[i],10)*re63[i]
					-220.0*pow(fmn[i],9)*re64[i]+495.0*pow(fmn[i],8)*re65[i]-792.0*pow(fmn[i],7)*re66[i]
					+924.0*pow(fmn[i],6)*re67[i]-790.0*pow(fmn[i],5)*re68[i]+487.0*pow(fmn[i],4)*re69[i]
					-208.0*pow(fmn[i],3)*re610[i]+58.0*pow(fmn[i],2)*re611[i]-8.0*fmn[i]*re612[i]+re613[i];

					re[i] = dxy*dxy* re[i];
			}
						
			Bccb1.clear();
			Bccb2.clear();
			Bccb3.clear();
			Bccb4.clear();
			Bccb5.clear();
			Bccb6.clear();
						
			fu1.clear();
			fu2.clear();
			fu3.clear();
			fu4.clear();
			fu5.clear();
			fu6.clear();
			fu7.clear();
			fu8.clear();
			fu9.clear();
			fu10.clear();
			fu11.clear();
			fu12.clear();
						
			re11.clear();
			re12.clear();
			re13.clear();
						
			re21.clear();
			re22.clear();
			re23.clear();
			re24.clear();
			re25.clear();
						
			re31.clear();
			re32.clear();
			re33.clear();
			re34.clear();
			re35.clear();
			re36.clear();
			re37.clear();
			
			re41.clear();
			re42.clear();
			re43.clear();
			re44.clear();
			re45.clear();
			re46.clear();
			re47.clear();
			re48.clear();
			re49.clear();
					
			re51.clear();
			re52.clear();
			re53.clear();
			re54.clear();
			re55.clear();
			re56.clear();
			re57.clear();
			re58.clear();
			re59.clear();
			re510.clear();
			re511.clear();
						
			re61.clear();
			re62.clear();
			re63.clear();
			re64.clear();
			re65.clear();
			re66.clear();
			re67.clear();
			re68.clear();
			re69.clear();
			re610.clear();
			re611.clear();
			re612.clear();
			re613.clear();
						
			break;
		}
	}
	//////////////////////////////////////////////////////////////////////////
}
//二维矢量波入射时G1,2(R)的泰勒展开系数
Complex resolve_matrix::Para_Two_D_Conductor_Vector(int Bool,double R,Complex k)
{
	Complex re;
	Complex t0 = exp(i_*k*R) / (4.0*PI*R);
	switch(Bool)
	{
	case 0:
		{
			re = t0*(i_*k*R-1.0)/pow(R,2.0);
			break;
		}
	case 1:
		{
			re = t0*(-pow(k,2.0) - 3.0*i_*k/R + 3.0/pow(R,2.0))/2.0;
			break;
		}
	case 2:
		{
			re = (1.0/factor(2))*t0*(-i_*pow(k,3.0)*R + 6.0*pow(k,2.0) + 15.0*i_*k/R - 15.0/pow(R,2.0))/4.0;
			break;
		}
	case 3:
		{
			re = (pow(k,4.0)*pow(R,2.0) + 10.0*i_*pow(k,3.0)*R - 45.0*pow(k,2.0) - 105.0*i_*k/R + 105.0/pow(R,2.0))/8.0;
			break;
		}
	case 4:
		{
			re = (1.0/(factor(4)*factor(3)*factor(2)))*t0*(i_*pow(k,5.0)*pow(R,3.0) - 15.0*pow(k,4.0)*pow(R,2.0) - 105.0*i_*pow(k,3.0)*R + 420.0*pow(k,2.0) + 945.0*i_*k/R - 945.0/pow(R,2.0))/16.0;
			break;
		}
	case 5:
		{
			re = (1.0/(factor(5)*factor(4)*factor(3)*factor(2)))*t0*(-pow(k,6.0)*pow(R,4.0) - 21.0*i_*pow(k,5.0)*pow(R,3.0) + 210.0*pow(k,4.0)*pow(R,2.0) + 1260.0*i_*pow(k,3.0)*R - 4725.0*pow(k,2.0) - 10395.0*i_*k/R + 10395.0/pow(R,2.0))/32.0;
			break;
		}
	case 6:
		{
			re = (1.0/(factor(6)*factor(5)*factor(4)*factor(3)*factor(2)))*t0*(-i_*pow(k,7.0)*pow(R,5.0) + 28.0*pow(k,6.0)*pow(R,4.0) + 378.0*i_*pow(k,5.0)*pow(R,3.0) - 3150.0*pow(k,4.0)*pow(R,2.0) - 17325.0*i_*pow(k,3.0)*R + 62370.0*pow(k,2.0) + 135135.0*i_*k/R - 135135.0/pow(R,2.0))/64.0;
			break;
		}
	}
	return re;
}
//Bi-CGSTAB algorithm with preconditioning :P138
void resolve_matrix::Two_D_Cond_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs,Vec_Com &zfs,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &rex)
{
	int i;
	int Num=nxy*nxy;
	rex.assign(Num);

	Vec_Com vector_r,vector_p,vector_pt,vector_r0p,vector_s,vector_st;
	vector_r.assign(Num);
	vector_p.assign(Num);
	vector_pt.assign(Num);
	vector_r0p.assign(Num);
	vector_s.assign(Num);
	vector_st.assign(Num);

	Vec_Com vector_t,vector_v,t1,t2;
	vector_t.assign(Num);
	vector_v.assign(Num);
		
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	double norm_0=sqrt(L2_Norm(vmn));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(Num);
	for (i=0; i<Num; i++)
	{
		Diag_Matrix_Inverse[i] = 1.0 / zs_diag[i];
		//cout<<vector_r0p[i]<<"   "<<vector_r[i]<<endl;
	}
	
	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r);
		//cout<<rou;
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<Num; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<Num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		NearMultiplyVector(' ',zs,vector_pt,t1);
		//Multi_Zs_Vec(zs,vector_pt,rd_num,nxy,t1);
		Multi_BTTB_Vec(zfs,vector_pt,nxy,t2);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = t1[i] + t2[i];
		}
		temp1=Multi_vt_v(vector_r0p,vector_v);
		alfa=rou/temp1;
		
		for (i=0; i<Num; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
			//cout<<vector_v[i]<<endl;
		}

		Norm=sqrt(L2_Norm(vector_s));
		
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<Num; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<Num; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}
		NearMultiplyVector(' ',zs,vector_st,t1);
		//Multi_Zs_Vec(zs,vector_st,rd_num,nxy,t1);
		Multi_BTTB_Vec(zfs,vector_st,nxy,t2);
		for (i=0; i<Num; i++)
		{
			vector_t[i] = t1[i] + t2[i];
		}
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for (i=0; i<Num; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eeps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

}
void resolve_matrix::Two_D_Cond_Vector_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs_lu,vector <Zmn_Cell> &zs_ru,vector <Zmn_Cell> &zs_ld,vector <Zmn_Cell> &zs_rd,
													 Vec_Com &zfs_lu,Vec_Com &zfs_ru,Vec_Com &zfs_ld,Vec_Com &zfs_rd,
													 Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d,
													 int nxy,int rd_num,double eeps,Vec_Com &rex)
{
	int i;
	int Num=nxy*nxy;
	rex.assign(2*Num);

	Vec_Com vector_r,vector_p,vector_pt,vector_r0p,vector_s,vector_st;
	vector_r.assign(2*Num);
	vector_p.assign(2*Num);
	vector_pt.assign(2*Num);
	vector_r0p.assign(2*Num);
	vector_s.assign(2*Num);
	vector_st.assign(2*Num);

	Vec_Com vector_t,vector_v,t1,t2,t3,t4;
	vector_t.assign(2*Num);
	vector_v.assign(2*Num);
		
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	Vec_Com vmn;
	vmn.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		vmn[i] = vmn_u[i];
		vmn[i+Num] = vmn_d[i];
	}
	double norm_0=sqrt(L2_Norm(vmn));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		Diag_Matrix_Inverse[i] = 1.0 / zs_d_lu[i];
		Diag_Matrix_Inverse[i+Num] = 1.0 / zs_d_rd[i];
		//cout<<vector_r0p[i]<<"   "<<vector_r[i]<<endl;
	}
	Vec_Com vector_pt_Up,vector_pt_Down,vector_st_Up,vector_st_Down;
	vector_pt_Up.assign(Num);
	vector_pt_Down.assign(Num);
	vector_st_Up.assign(Num);
	vector_st_Down.assign(Num);
	
	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r);
		//cout<<rou;
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<2*Num; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<2*Num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_pt_Up[i] = vector_pt[i];
			vector_pt_Down[i] = vector_pt[i+Num];
		}
		Two_D_Vector_NearMultiplyVector(' ',zs_lu,vector_pt_Up,t1);
		Two_D_Vector_NearMultiplyVector(' ',zs_ru,vector_pt_Down,t2);
        Two_D_Vector_NearMultiplyVector(' ',zs_ld,vector_pt_Up,t3);
		Two_D_Vector_NearMultiplyVector(' ',zs_rd,vector_pt_Down,t4);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = t1[i] + t2[i];
			vector_v[i+Num] = t3[i] + t4[i];
		}
		Multi_BTTB_Vec(zfs_lu,vector_pt_Up,nxy,t1);
		Multi_BTTB_Vec(zfs_ru,vector_pt_Down,nxy,t2);
		Multi_BTTB_Vec(zfs_ld,vector_pt_Up,nxy,t3);
		Multi_BTTB_Vec(zfs_rd,vector_pt_Down,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = vector_v[i] + t1[i] + t2[i];
			vector_v[i+Num] = vector_v[i+Num] + t3[i] + t4[i];
		}
		temp1=Multi_vt_v(vector_r0p,vector_v);
		alfa=rou/temp1;
		
		for (i=0; i<2*Num; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
			//cout<<vector_v[i]<<endl;
		}

		Norm=sqrt(L2_Norm(vector_s));
		
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<2*Num; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<2*Num; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_st_Up[i] = vector_st[i];
			vector_st_Down[i] = vector_st[i+Num];
		}
		Two_D_Vector_NearMultiplyVector(' ',zs_lu,vector_st_Up,t1);
		Two_D_Vector_NearMultiplyVector(' ',zs_ru,vector_st_Down,t2);
        Two_D_Vector_NearMultiplyVector(' ',zs_ld,vector_st_Up,t3);
		Two_D_Vector_NearMultiplyVector(' ',zs_rd,vector_st_Down,t4);
		
		for (i=0; i<Num; i++)
		{
			vector_t[i] = t1[i] + t2[i];
			vector_t[i+Num] = t3[i] + t4[i];
		}
		Multi_BTTB_Vec(zfs_lu,vector_st_Up,nxy,t1);
		Multi_BTTB_Vec(zfs_ru,vector_st_Down,nxy,t2);
		Multi_BTTB_Vec(zfs_ld,vector_st_Up,nxy,t3);
		Multi_BTTB_Vec(zfs_rd,vector_st_Down,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_t[i] = vector_t[i] + t1[i] + t2[i];
			vector_t[i+Num] = vector_t[i+Num] + t3[i] + t4[i];
		}
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for (i=0; i<2*Num; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eeps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	vector_r.clear();
	vector_p.clear();
	vector_pt.clear();
	vector_r0p.clear();
	vector_s.clear();
	vector_st.clear();
	vector_t.clear();
	vector_v.clear();
	t1.clear();
	t2.clear();
	t3.clear();
	t4.clear();
	vmn.clear();
	Diag_Matrix_Inverse.clear();
	vector_pt_Up.clear();
	vector_pt_Down.clear();
	vector_st_Up.clear();
	vector_st_Down.clear();

}
void resolve_matrix::Two_D_Cond_Vector_Bi_CGSTAB_Pre(Vec_Com &zs_lu,Vec_Com &zs_ru,Vec_Com &zs_ld,Vec_Com &zs_rd,
													 Vec_Com &zfs_lu,Vec_Com &zfs_ru,Vec_Com &zfs_ld,Vec_Com &zfs_rd,
													 Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d,
													 int nxy,int rd_num,double eeps,Vec_Com &rex)
{
	int i;
	int Num=nxy*nxy;
	rex.assign(2*Num);

	Vec_Com vector_r,vector_p,vector_pt,vector_r0p,vector_s,vector_st;
	vector_r.assign(2*Num);
	vector_p.assign(2*Num);
	vector_pt.assign(2*Num);
	vector_r0p.assign(2*Num);
	vector_s.assign(2*Num);
	vector_st.assign(2*Num);

	Vec_Com vector_t,vector_v,t1,t2,t3,t4;
	vector_t.assign(2*Num);
	vector_v.assign(2*Num);
		
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	Vec_Com vmn;
	vmn.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		vmn[i] = vmn_u[i];
		vmn[i+Num] = vmn_d[i];
	}
	double norm_0=sqrt(L2_Norm(vmn));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		Diag_Matrix_Inverse[i] = 1.0 / zs_d_lu[i];
		Diag_Matrix_Inverse[i+Num] = 1.0 / zs_d_rd[i];
		//cout<<vector_r0p[i]<<"   "<<vector_r[i]<<endl;
	}
	Vec_Com vector_pt_Up,vector_pt_Down,vector_st_Up,vector_st_Down;
	vector_pt_Up.assign(Num);
	vector_pt_Down.assign(Num);
	vector_st_Up.assign(Num);
	vector_st_Down.assign(Num);
	
	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r);
		//cout<<rou;
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<2*Num; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<2*Num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_pt_Up[i] = vector_pt[i];
			vector_pt_Down[i] = vector_pt[i+Num];
		}
		Multi_Zs_Vec(zs_lu,vector_pt_Up,rd_num,nxy,t1);
		Multi_Zs_Vec(zs_ru,vector_pt_Down,rd_num,nxy,t2);
		Multi_Zs_Vec(zs_ld,vector_pt_Up,rd_num,nxy,t3);
		Multi_Zs_Vec(zs_rd,vector_pt_Down,rd_num,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = t1[i] + t2[i];
			vector_v[i+Num] = t3[i] + t4[i];
		}
		Multi_BTTB_Vec(zfs_lu,vector_pt_Up,nxy,t1);
		Multi_BTTB_Vec(zfs_ru,vector_pt_Down,nxy,t2);
		Multi_BTTB_Vec(zfs_ld,vector_pt_Up,nxy,t3);
		Multi_BTTB_Vec(zfs_rd,vector_pt_Down,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = vector_v[i] + t1[i] + t2[i];
			vector_v[i+Num] = vector_v[i+Num] + t3[i] + t4[i];
		}
		temp1=Multi_vt_v(vector_r0p,vector_v);
		alfa=rou/temp1;
		
		for (i=0; i<2*Num; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
			//cout<<vector_v[i]<<endl;
		}

		Norm=sqrt(L2_Norm(vector_s));
		
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<2*Num; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<2*Num; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_st_Up[i] = vector_st[i];
			vector_st_Down[i] = vector_st[i+Num];
		}
		
		Multi_Zs_Vec(zs_lu,vector_st_Up,rd_num,nxy,t1);
		Multi_Zs_Vec(zs_ru,vector_st_Down,rd_num,nxy,t2);
		Multi_Zs_Vec(zs_ld,vector_st_Up,rd_num,nxy,t3);
		Multi_Zs_Vec(zs_rd,vector_st_Down,rd_num,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_t[i] = t1[i] + t2[i];
			vector_t[i+Num] = t3[i] + t4[i];
		}
		Multi_BTTB_Vec(zfs_lu,vector_st_Up,nxy,t1);
		Multi_BTTB_Vec(zfs_ru,vector_st_Down,nxy,t2);
		Multi_BTTB_Vec(zfs_ld,vector_st_Up,nxy,t3);
		Multi_BTTB_Vec(zfs_rd,vector_st_Down,nxy,t4);
		for (i=0; i<Num; i++)
		{
			vector_t[i] = vector_t[i] + t1[i] + t2[i];
			vector_t[i+Num] = vector_t[i+Num] + t3[i] + t4[i];
		}
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for (i=0; i<2*Num; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eeps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	vector_r.clear();
	vector_p.clear();
	vector_pt.clear();
	vector_r0p.clear();
	vector_s.clear();
	vector_st.clear();
	vector_t.clear();
	vector_v.clear();
	t1.clear();
	t2.clear();
	t3.clear();
	t4.clear();
	vmn.clear();
	Diag_Matrix_Inverse.clear();
	vector_pt_Up.clear();
	vector_pt_Down.clear();
	vector_st_Up.clear();
	vector_st_Down.clear();

}
//SMFSIA计算二维矢量波入射导体粗糙面散射中的矩阵方程解法
void resolve_matrix::Two_D_Cond_Vector_Bi_CGSTAB_Pre(vector <Zmn_Cell> &zs_lu,vector <Zmn_Cell> &zs_ru,vector <Zmn_Cell> &zs_ld,vector <Zmn_Cell> &zs_rd,
													 Vec_Com &zs_d_lu,Vec_Com &zs_d_rd,Vec_Com &vmn_u,Vec_Com &vmn_d,
													 Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &fmn,Vec_Dou &dfxmn,
													 Vec_Dou &dfymn,Complex kw,double dxy,int nxy,int rd_num,double eeps,
													 Vec_Com &rex)
{
	int i;
	int Num=nxy*nxy;
	rex.assign(2*Num);

	Vec_Com vector_r,vector_p,vector_pt,vector_r0p,vector_s,vector_st;
	vector_r.assign(2*Num);
	vector_p.assign(2*Num);
	vector_pt.assign(2*Num);
	vector_r0p.assign(2*Num);
	vector_s.assign(2*Num);
	vector_st.assign(2*Num);

	Vec_Com vector_t,vector_v,t1,t2,t3,t4;
	vector_t.assign(2*Num);
	vector_v.assign(2*Num);
		
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	Vec_Com vmn;
	vmn.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		vmn[i] = vmn_u[i];
		vmn[i+Num] = vmn_d[i];
	}
	double norm_0=sqrt(L2_Norm(vmn));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(2*Num);
	for (i=0; i<Num; i++)
	{
		Diag_Matrix_Inverse[i] = 1.0 / zs_d_lu[i];
		Diag_Matrix_Inverse[i+Num] = 1.0 / zs_d_rd[i];
		//cout<<vector_r0p[i]<<"   "<<vector_r[i]<<endl;
	}
	Vec_Com vector_pt_Up,vector_pt_Down,vector_st_Up,vector_st_Down;
	vector_pt_Up.assign(Num);
	vector_pt_Down.assign(Num);
	vector_st_Up.assign(Num);
	vector_st_Down.assign(Num);
	
	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r);
		//cout<<rou;
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<2*Num; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<2*Num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_pt_Up[i] = vector_pt[i];
			vector_pt_Down[i] = vector_pt[i+Num];
		}
		Two_D_Vector_NearMultiplyVector(' ',zs_lu,vector_pt_Up,t1);
		Two_D_Vector_NearMultiplyVector(' ',zs_ru,vector_pt_Down,t2);
        Two_D_Vector_NearMultiplyVector(' ',zs_ld,vector_pt_Up,t3);
		Two_D_Vector_NearMultiplyVector(' ',zs_rd,vector_pt_Down,t4);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = t1[i] + t2[i];
			vector_v[i+Num] = t3[i] + t4[i];
		}
		Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(0,xmn,ymn,fmn,dfxmn,dfymn,vector_pt,kw,dxy,nxy,rd_num,t1);
		for (i=0; i<Num; i++)
		{
			vector_v[i] = vector_v[i] + t1[i];
			vector_v[i+Num] = vector_v[i+Num] + t1[i+Num];
		}
		temp1=Multi_vt_v(vector_r0p,vector_v);
		alfa=rou/temp1;
		
		for (i=0; i<2*Num; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
			//cout<<vector_v[i]<<endl;
		}

		Norm=sqrt(L2_Norm(vector_s));
		
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<2*Num; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<2*Num; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}
		for (i=0; i<Num; i++)
		{
			vector_st_Up[i] = vector_st[i];
			vector_st_Down[i] = vector_st[i+Num];
		}
		Two_D_Vector_NearMultiplyVector(' ',zs_lu,vector_st_Up,t1);
		Two_D_Vector_NearMultiplyVector(' ',zs_ru,vector_st_Down,t2);
        Two_D_Vector_NearMultiplyVector(' ',zs_ld,vector_st_Up,t3);
		Two_D_Vector_NearMultiplyVector(' ',zs_rd,vector_st_Down,t4);
		
		for (i=0; i<Num; i++)
		{
			vector_t[i] = t1[i] + t2[i];
			vector_t[i+Num] = t3[i] + t4[i];
		}
		Two_D_Conductor_Vector_Multi_Zw_Vec_Vec(0,xmn,ymn,fmn,dfxmn,dfymn,vector_st,kw,dxy,nxy,rd_num,t1);
		for (i=0; i<Num; i++)
		{
			vector_t[i] = vector_t[i] + t1[i];
			vector_t[i+Num] = vector_t[i+Num] + t1[i+Num];
		}
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for (i=0; i<2*Num; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eeps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	vector_r.clear();
	vector_p.clear();
	vector_pt.clear();
	vector_r0p.clear();
	vector_s.clear();
	vector_st.clear();
	vector_t.clear();
	vector_v.clear();
	t1.clear();
	t2.clear();
	t3.clear();
	t4.clear();
	vmn.clear();
	Diag_Matrix_Inverse.clear();
	vector_pt_Up.clear();
	vector_pt_Down.clear();
	vector_st_Up.clear();
	vector_st_Down.clear();

}
//SMCG计算二维介质粗糙面散射中的矩阵方程解法
void resolve_matrix::Two_D_Diel_Bi_CGSTAB_Pre(Vec_Zmn_Cell zs[6][6],Vec_Com zs_diag[6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,
											  Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,int nxy,int rd_num,double eeps,Vec_Com rex[6])
{
	int i,k;
	int Num=nxy*nxy;
	for ( i=0; i<6; i++ )
	{
		rex[i].assign(Num);
	}

	Vec_Com vector_r[6],vector_p[6],vector_r0p[6],vector_s[6],vector_pt[6],vector_st[6],vector_t[6],vector_v[6];
	for ( k=0; k<6; k++ )
	{
		vector_r[k].assign(Num);
		vector_p[k].assign(Num);
		vector_r0p[k].assign(Num);
		vector_s[k].assign(Num);
		vector_pt[k].assign(Num);
		vector_st[k].assign(Num);
		vector_t[k].assign(Num);
		vector_v[k].assign(Num);
	}
	
	Vec_Com tre[6];
	
		
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;

	double norm_0=sqrt(L2_Norm(vmn,6));

	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0p);
	
	//预处理矩阵为对角阵
/*	Vec_Com Diag_Matrix_Inverse[6];
	for ( k=0; k<6; k++ )
	{
		Diag_Matrix_Inverse[k].assign(Num);
	}
	
	for ( k=0; k<6; k++ )
	{
		for (i=0; i<Num; i++)
		{
			Diag_Matrix_Inverse[k][i] = 1.0 / zs_diag[k][i];
		}
	}*/
	
	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou = Multi_vt_v(vector_r0p,vector_r);
		//cout<<rou;
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p);
		}
		else
		{
			beta = (rou/rou0)*(alfa/value_w);
			for ( k=0; k<6; k++ )
			{
				for (i=0; i<Num; i++)
				{
					vector_p[k][i] = vector_r[k][i]+beta*(vector_p[k][i]-value_w*vector_v[k][i]);
				}
			}
			
		}
		Two_D_Diel_CGM(zs,vector_p,nxy,eeps,vector_pt);
/*		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				vector_pt[k][i] = Diag_Matrix_Inverse[k][i]*vector_p[k][i];
			}
		}*/
		
		Two_D_Diel_Zs_Vector_Multiply(' ',zs,vector_pt,vector_v);
		Two_D_Diel_Multi_Zfs_Vec(' ',xmn,ymn,zmn,dfxmn,dfymn,vector_pt,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				vector_v[k][i] += tre[k][i];
			}
		}
		temp1=Multi_vt_v(vector_r0p,vector_v);
		alfa=rou/temp1;
		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				vector_s[k][i] = vector_r[k][i]-alfa*vector_v[k][i];
			}
		}		

		Norm=sqrt(L2_Norm(vector_s,6));
		
		if ( Norm < 1.0e-6 )
		{
			for ( k=0; k<6; k++ )
			{
				for (i=0; i<Num; i++)
				{
					rex[k][i] += alfa*vector_pt[k][i];
				}
			}
			break;
		}
		Two_D_Diel_CGM(zs,vector_s,nxy,eeps,vector_st);
/*		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				vector_st[k][i] = Diag_Matrix_Inverse[k][i]*vector_s[k][i];
			}
		}*/
		Two_D_Diel_Zs_Vector_Multiply(' ', zs,vector_st,vector_t);
		Two_D_Diel_Multi_Zfs_Vec(' ',xmn,ymn,zmn,dfxmn,dfymn,vector_st,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				vector_t[k][i] += tre[k][i];
			}
		}
		temp1=Multi_vt_v(vector_t,vector_s);
		temp2=Multi_vt_v(vector_t,vector_t);
		value_w=temp1/temp2;
		for ( k=0; k<6; k++ )
		{
			for (i=0; i<Num; i++)
			{
				rex[k][i] += alfa*vector_pt[k][i]+value_w*vector_st[k][i];
				vector_r[k][i] = vector_s[k][i]-value_w*vector_t[k][i];
			}
		}
		
		
		Norm=sqrt(L2_Norm(vector_r,6))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eeps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	for ( k=0; k<6; k++ )
	{
		vector_r[k].clear();
		vector_p[k].clear();
		vector_r0p[k].clear();
		vector_s[k].clear();
		vector_pt[k].clear();
		vector_st[k].clear();
		vector_t[k].clear();
		vector_v[k].clear();
		tre[k].clear();
//		Diag_Matrix_Inverse[k].clear();
	}

}
void resolve_matrix::Two_D_Diel_GPBI_CG(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,
										Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,
										int nxy,int rd_num,double eeps,Vec_Com rex[6])
{
	int k,i;
	int Num = nxy * nxy;
	
	for ( i=0; i<6; i++ )
	{
		rex[i].assign(Num);
	}
	Vec_Com tre[6];
	Vec_Com vector_r[6],vector_r0[6],vector_u[6],vector_z[6],vector_p[6],vector_t[6],vector_y[6],vector_w[6],vector_s[6];
	Vec_Com vector_v[6],vector_q[6];
	for ( i=0; i<6; i++ )
	{
		tre[i].assign(Num);
		vector_r[i].assign(Num);
		vector_r0[i].assign(Num);
		vector_u[i].assign(Num);
		vector_z[i].assign(Num);
		vector_p[i].assign(Num);
		vector_t[i].assign(Num);
		vector_y[i].assign(Num);
		vector_w[i].assign(Num);
		vector_s[i].assign(Num);
	}
	Complex rou,rou0,alfa,beta,value_u1,value_u2,value_u3,value_u4,
		value_u5,value_gama,value_eta,value_tao,temp;
	double Norm;

	double norm_0=sqrt(L2_Norm(vmn,6));
	Assign_Vecor(vmn,vector_r);
	Assign_Vecor(vmn,vector_r0);

	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,6)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0,vector_r);
		if ( iter == 1 )
		{
			Assign_Vecor(vector_r,vector_p);
			Two_D_Diel_Zs_Vector_Multiply(' ',zs,vector_p,vector_q);
			Two_D_Diel_Multi_Zw_Vec(0,0,xmn,ymn,zmn,dfxmn,dfymn,vector_p,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_q[k][i] += tre[k][i];
				}
			}
//			vector_q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp=Multi_vt_v(vector_r0,vector_q);
			alfa=rou/temp;
			for ( k=0; k<6; k++ )
			{
				for (i=0; i<Num; i++)
				{
					vector_t[k][i] = vector_r[k][i] - alfa * vector_q[k][i];
				}
			}
			
			Two_D_Diel_Zs_Vector_Multiply(' ',zs,vector_t,vector_v);
			Two_D_Diel_Multi_Zw_Vec(0,0,xmn,ymn,zmn,dfxmn,dfymn,vector_t,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_v[k][i] += tre[k][i];
				}
			}
//			vector_v=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_y[k][i] = alfa * vector_q[k][i] - vector_r[k][i];
				}
			}
			
			value_u2 = Multi_vt_v(vector_v,vector_t);
			value_u5 = Multi_vt_v(vector_v,vector_v);
			value_gama = value_u2 / value_u5;
			value_eta = 0.0;
		}
		else
		{
			beta = (rou / rou0) * (alfa / value_gama);
			for ( k=0; k<6; k++ )
			{
				for (i=0; i<Num; i++)
				{
					vector_w[k][i] = vector_v[k][i] + beta * vector_q[k][i];
					vector_p[k][i] = vector_r[k][i] + beta * (vector_p[k][i] - vector_u[k][i]);
				}
			}
			
			Two_D_Diel_Zs_Vector_Multiply(' ',zs,vector_p,vector_q);
			Two_D_Diel_Multi_Zw_Vec(0,0,xmn,ymn,zmn,dfxmn,dfymn,vector_p,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
		
/*			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					cout<<tre[k][i]<<"  "<<k<<"  "<<i<<endl;
				}
			}*/
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_q[k][i] += tre[k][i];
				}
			}
			
//			vector_q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp = Multi_vt_v(vector_r0,vector_q);
			alfa = rou / temp;
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_s[k][i] = vector_t[k][i] - vector_r[k][i];
					vector_t[k][i] = vector_r[k][i] - alfa * vector_q[k][i];
				}
			}

			Two_D_Diel_Zs_Vector_Multiply(' ',zs,vector_t,vector_v);
			Two_D_Diel_Multi_Zw_Vec(0,0,xmn,ymn,zmn,dfxmn,dfymn,vector_t,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,tre);
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_v[k][i] += tre[k][i];
				}
			}
//			vector_v=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for ( k=0; k<6; k++ )
			{
				for ( i=0; i<Num; i++ )
				{
					vector_y[k][i] = vector_s[k][i] - alfa * (vector_w[k][i] - vector_q[k][i]);
				}
			}
			value_u1 = Multi_vt_v(vector_y,vector_y);
			value_u2 = Multi_vt_v(vector_v,vector_t);
			value_u3 = Multi_vt_v(vector_y,vector_t);
			value_u4 = Multi_vt_v(vector_v,vector_y);
			value_u5 = Multi_vt_v(vector_v,vector_v);
			value_tao = value_u5 * value_u1 - conj(value_u4) * value_u4;
			value_gama = (value_u1 * value_u2 - value_u3 * value_u4) / value_tao;
			value_eta = (value_u5 * value_u3 - conj(value_u4) * value_u2) / value_tao;
		}
		
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				vector_u[k][i] = value_gama * vector_q[k][i] + value_eta * (vector_s[k][i] + beta * vector_u[k][i]);
				vector_z[k][i] = value_gama * vector_r[k][i] + value_eta * vector_z[k][i] - alfa * vector_u[k][i];
				rex[k][i] += alfa * vector_p[k][i] + vector_z[k][i];
				vector_r[k][i] = vector_t[k][i] - value_eta * vector_y[k][i] - value_gama * vector_v[k][i];
			}
		}

		Norm=sqrt(L2_Norm(vector_r,6))/norm_0;
		cout<<"GPBi-CG steps "<<iter<<",  res = "<<Norm<<endl;
		
		if ( Norm < eps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total GPBi-CG Steps = "<<iter<<", to achieve "<<Norm<<endl;

	for ( i=0; i<6; i++ )
	{
		tre[i].clear();
		vector_r[i].clear();
		vector_r0[i].clear();
		vector_u[i].clear();
		vector_z[i].clear();
		vector_p[i].clear();
		vector_t[i].clear();
		vector_y[i].clear();
		vector_w[i].clear();
		vector_s[i].clear();
		vector_v[i].clear();
		vector_q[i].clear();
	}
	
}
void resolve_matrix::Two_D_Diel_CGM(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,
									Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,
									int nxy,int rd_num,double eps,Vec_Com rex[6])
{  
	int k,i;
	int Num = nxy * nxy;
	
	for ( i=0; i<6; i++ )
	{
		rex[i].assign(Num);
	}
	double ebs;
	double excitationAbs;
	double alpha, beta;
	
	Vec_Com p[6],r[6];
	for ( i=0; i<6; i++ )
	{
		p[i].assign(Num);
		r[i].assign(Num);
	}

	Vec_Com t0[6],t1[6];
	
	int step=0;
	
	excitationAbs=sqrt(L2_Norm(vmn,6));   //平方和开根号
	Assign_Vecor(vmn,r);
	Two_D_Diel_Zs_Vector_Multiply('*',zs,vmn,t0);
	Two_D_Diel_Multi_Zfs_Vec('*',xmn,ymn,zmn,dfxmn,dfymn,vmn,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			t0[k][i] += t1[k][i];
		}
	}
//	t0=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,Bm,CX,DX,fbw,lbw,dbw,ndg,ncg);
	beta=1.0/L2_Norm(t0,6);   //平方和
	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			p[k][i] = beta * t0[k][i]; 
		}
	}
//	for(i=0;i<2*ndg;i++)
//		p[i]=beta*t0[i];
	
	while ( true )
	{
		out_resual<<step<<"   "<<sqrt(L2_Norm(r,6)) / excitationAbs<<endl;
		step++;
		Two_D_Diel_Zs_Vector_Multiply(' ',zs,p,t0);
		Two_D_Diel_Multi_Zfs_Vec(' ',xmn,ymn,zmn,dfxmn,dfymn,p,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				t1[k][i] += t0[k][i];
			}
		}
//		t1=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
		alpha=1.0/L2_Norm(t1,6);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				rex[k][i] += alpha * p[k][i];
				r[k][i] -= alpha * t1[k][i];
			}
		}
//		for(i=0;i<2*ndg;i++)
//		{
//			rex[i]+=alpha*p[i];
//			r[i]=r[i]-alpha*t1[i];
//		}
		ebs=sqrt(L2_Norm(r,6))/excitationAbs;

		cout << "steps= " <<step << ", ebs= " <<ebs << endl;
		
		if(ebs < eps )//停止准则:0.01即可得到精确解
		{
			break;     //终止整个循环
		}
		
		Two_D_Diel_Zs_Vector_Multiply('*',zs,r,t0);
		Two_D_Diel_Multi_Zfs_Vec('*',xmn,ymn,zmn,dfxmn,dfymn,r,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				t0[k][i] += t1[k][i];
			}
		}
//		t2=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
		beta=1.0/L2_Norm(t0,6);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				p[k][i] += beta * t0[k][i];
			}
		}
//		for(i=0;i<2*ndg;i++)
//			p[i]+=beta*t2[i];
	}

	for ( i=0; i<6; i++ )
	{
		p[i].clear();
		r[i].clear();
		t0[i].clear();
		t1[i].clear();
	}
}
void resolve_matrix::Two_D_Diel_CGM(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],int nxy,double eps,Vec_Com rex[6])
{  
	int k,i;
	int Num = nxy * nxy;
	
	for ( i=0; i<6; i++ )
	{
		rex[i].assign(Num);
	}
	double ebs;
	double excitationAbs;
	double alpha, beta;
	
	Vec_Com p[6],r[6];
	for ( i=0; i<6; i++ )
	{
		p[i].assign(Num);
		r[i].assign(Num);
	}

	Vec_Com t0[6];
	
	int step=0;
	
	excitationAbs=sqrt(L2_Norm(vmn,6));   //平方和开根号
	Assign_Vecor(vmn,r);
	Two_D_Diel_Zs_Vector_Multiply('*',zs,vmn,t0);
	beta=1.0/L2_Norm(t0,6);   //平方和
	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			p[k][i] = beta * t0[k][i]; 
		}
	}
	
	while ( true )
	{
//		out_resual<<step<<"   "<<sqrt(L2_Norm(r,6)) / excitationAbs<<endl;
		step++;
		Two_D_Diel_Zs_Vector_Multiply(' ',zs,p,t0);
		alpha=1.0/L2_Norm(t0,6);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				rex[k][i] += alpha * p[k][i];
				r[k][i] -= alpha * t0[k][i];
			}
		}
		ebs=sqrt(L2_Norm(r,6))/excitationAbs;

		//cout << "CGM steps= " <<step << ", ebs= " <<ebs << endl;

		
		if(ebs < eps )//停止准则:0.01即可得到精确解
		{
			break;     //终止整个循环
		}
		
		Two_D_Diel_Zs_Vector_Multiply('*',zs,r,t0);
		beta=1.0/L2_Norm(t0,6);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				p[k][i] += beta * t0[k][i];
			}
		}
	}

	for ( i=0; i<6; i++ )
	{
		p[i].clear();
		r[i].clear();
		t0[i].clear();
	}
}
void resolve_matrix::Two_D_Diel_CGLS(Vec_Zmn_Cell zs[6][6],Vec_Com vmn[6],Vec_Dou &xmn,Vec_Dou &ymn,Vec_Dou &zmn,Vec_Dou &dfxmn,
									 Vec_Dou &dfymn,double wave_num1,Complex wave_num2,double wave_Im1,Complex permi,double dxy,
									 int nxy,int rd_num,double eps,Vec_Com rex[6])
{
	int k,i;
	int Num = nxy * nxy;
	
	for ( i=0; i<6; i++ )
	{
		rex[i].assign(Num);
	}
	double beta,alpha,te;
	
	Vec_Com p[6],r[6];
	for ( i=0; i<6; i++ )
	{
		p[i].assign(Num);
		r[i].assign(Num);
	}
	
	Vec_Com t0[6],t1[6];

	//
	Assign_Vecor(vmn,r);
	Two_D_Diel_Zs_Vector_Multiply('*',zs,r,t0);
	Two_D_Diel_Multi_Zfs_Vec('*',xmn,ymn,zmn,dfxmn,dfymn,r,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
	for ( k=0; k<6; k++ )
	{
		for ( i=0; i<Num; i++ )
		{
			t0[k][i] += t1[k][i];
		}
	}
//  s=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
	Assign_Vecor(t0,p);
	double norm_0=sqrt(L2_Norm(vmn,6));
	int Steps=0;
	while( sqrt(L2_Norm(r,6))/norm_0 > eps )
	{
		//out_resual<<Steps<<"   "<<sqrt(L2_Norm(r,6)) / norm_0<<endl;
		Steps++;
		cout<<"Steps="<<Steps<<" ,R_Norm="<<sqrt(L2_Norm(r,6))/norm_0<<endl;
		
		te=L2_Norm(t0,6);
		Two_D_Diel_Zs_Vector_Multiply(' ',zs,p,t0);
		Two_D_Diel_Multi_Zfs_Vec(' ',xmn,ymn,zmn,dfxmn,dfymn,p,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				t1[k][i] += t0[k][i];
			}
		}
//		q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
		alpha=te/L2_Norm(t1,6);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				rex[k][i] += alpha * p[k][i];
				r[k][i] -= alpha * t1[k][i];
			}
		}
		Two_D_Diel_Zs_Vector_Multiply('*',zs,r,t0);
		Two_D_Diel_Multi_Zfs_Vec('*',xmn,ymn,zmn,dfxmn,dfymn,r,wave_num1,wave_num2,wave_Im1,permi,dxy,nxy,rd_num,t1);
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				t0[k][i] += t1[k][i];
			}
		}
//		s=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
		beta=L2_Norm(t0,6)/te;
		for ( k=0; k<6; k++ )
		{
			for ( i=0; i<Num; i++ )
			{
				p[k][i] = t0[k][i] + beta * p[k][i];
			}
		}			
	}
	cout<<"Steps="<<Steps<<" ,R_Norm="<<sqrt(L2_Norm(r,6))/norm_0<<endl;

	for ( i=0; i<6; i++ )
	{
		p[i].clear();
		r[i].clear();
		t0[i].clear();
		t1[i].clear();
	}
}
void resolve_matrix::Two_D_Cond_CGM(Vec_Com &zs,Vec_Com &zfs,Vec_Com zs_inver,Vec_Com zfs_inver,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &resx)
{
	int i;
	double rfa,beta;
	int num=nxy*nxy;
	Vec_Com r,d;
	r.assign(num);
	d.assign(num);
	resx.assign(num);

	//矩阵的共轭转置
	
	//第0步
	Vec_Com t1,t2;
	Multi_Zs_Vec(zs_inver,vmn,rd_num,nxy,t1);
	Multi_BTTB_Vec(zfs_inver,vmn,nxy,t2);
	for (i=0; i<num; i++)
	{
		t1[i] = t1[i] + t2[i];
		//cout<<t1[i]<<endl;
	}
	Assign_Vecor(vmn,r);//0步
	Assign_Vecor(t1,d);//1步

	double norm_0=sqrt(L2_Norm(vmn));

	double sum_rr,sum_r,sum_m;
	sum_rr=L2_Norm(t1);
	double Norm=sqrt(L2_Norm(r)) / norm_0;
	//第一步到第n步
	int iter=0;
	while( Norm > eeps)
	{
		iter++;
		
		cout<<"CGM steps "<<iter<<",  res = "<<Norm<<endl;
	
		Multi_Zs_Vec(zs,d,rd_num,nxy,t1);
		Multi_BTTB_Vec(zfs,d,nxy,t2);
		for (i=0; i<num; i++)
		{
			t1[i] = t1[i] + t2[i];
			//cout<<t1[i]<<endl;
		}
		sum_m=L2_Norm(t1);
		rfa=sum_rr/sum_m;
		
	    for(i=0; i<num; i++)
		{
			resx[i]=resx[i]+rfa*d[i];
			r[i]=r[i]-rfa*t1[i];
		}
		//
		Multi_Zs_Vec(zs_inver,r,rd_num,nxy,t1);
    	Multi_BTTB_Vec(zfs_inver,r,nxy,t2);
    	for (i=0; i<num; i++)
		{
    		t1[i] = t1[i] + t2[i];
		}
        sum_r=L2_Norm(t1);
		beta=sum_r/sum_rr;
		
		for(i=0; i<num; i++)
		{
			d[i]=t1[i]+beta*d[i];//2
		}
		sum_rr=sum_r;
		Norm = sqrt(L2_Norm(r)) / norm_0;
	}

	r.clear();
	d.clear();
	t1.clear();
	t2.clear();
	
}
void resolve_matrix::Two_D_Cond_CGS_Pre(Vec_Com &zs,Vec_Com &zfs,Vec_Com &zs_diag,Vec_Com &vmn,int nxy,int rd_num,double eeps,Vec_Com &resx)
{
	int i;
	int num = nxy*nxy;
	resx.assign(num);
	Vec_Com vector_r0,vector_r,vector_u,vector_p,vector_q,vector_pt,vector_ut,vector_uq,
		vector_vp,t1,t2;
	vector_r0.assign(num);
	vector_r.assign(num);
	vector_u.assign(num);
	vector_p.assign(num);
	vector_q.assign(num);
	vector_pt.assign(num);
	vector_ut.assign(num);
	vector_uq.assign(num);
	vector_vp.assign(num);
	
	Complex alfa,beta,rou0,rou,t3;

	double Norm;
	double norm_0=sqrt(L2_Norm(vmn));
	
	//预处理矩阵为对角阵
	Vec_Com Diag_Matrix_Inverse;
	Diag_Matrix_Inverse.assign(num);
	for (i=0; i<num; i++)
	{
		Diag_Matrix_Inverse[i] = 1.0 / zs_diag[i];
	}
	
	Assign_Vecor(vmn,vector_r0);
	Assign_Vecor(vmn,vector_r);

	int iter=0;
	while ( true )
	{
		//out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0,vector_r);
		if ( iter == 1)
		{
			Assign_Vecor(vector_r0,vector_u);
			Assign_Vecor(vector_u,vector_p);
		}
		else
		{
			beta=rou/rou0;
			for (i=0; i<num; i++)
			{
				vector_u[i]=vector_r[i]+beta*vector_q[i];
				vector_p[i]=vector_u[i]+beta*(vector_q[i]+beta*vector_p[i]);
			}
		}
		
		for (i=0; i<num; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		
		Multi_Zs_Vec(zs,vector_pt,rd_num,nxy,t1);
		Multi_BTTB_Vec(zfs,vector_pt,nxy,t2);
		for (i=0; i<num; i++)
		{
			vector_vp[i] = t1[i] + t2[i];
		}

		t3=Multi_vt_v(vector_r0,vector_vp);
		alfa=rou/t3;
		for (i=0; i<num; i++)
		{
	    	vector_q[i]=vector_u[i]-alfa*vector_vp[i];
	    	vector_uq[i]=vector_u[i]+vector_q[i];
	        vector_ut[i]=Diag_Matrix_Inverse[i]*vector_uq[i];
		}

		Multi_Zs_Vec(zs,vector_ut,rd_num,nxy,t1);
		Multi_BTTB_Vec(zfs,vector_ut,nxy,t2);
		for (i=0; i<num; i++)
		{
			t1[i] = t1[i] + t2[i];
		}

    	for (i=0; i<num; i++)
		{
    		resx[i] += alfa*vector_ut[i];
    		vector_r[i] -= alfa*t1[i];
		}

		Norm=sqrt(L2_Norm(vector_r))/norm_0;
		cout<<"CGS"<<" Steps = "<<iter<<",  res = "<<Norm<<endl;
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r))<<endl;
		if (Norm < eps)
		{
			break;
		}
		
		rou0=rou;		
	}
	cout<<"Total CGS Steps = "<<iter<<", to achieve "<<Norm<<endl;
	
	vector_r0.clear();
	vector_r.clear();
	vector_u.clear();
	vector_p.clear();
	vector_q.clear();
	vector_pt.clear();
	vector_ut.clear();
	vector_uq.clear();
	vector_vp.clear();
	t1.clear();
	t2.clear();
	Diag_Matrix_Inverse.clear();
}
void resolve_matrix::Multi_Zs_Vec(Vec_Com &Zs,Vec_Com &vec,int RD_num,int nxy,Vec_Com &result)
{
	int i,j;
	result.assign(nxy*nxy);
	int Band_Ser,Band_Num;//带状矩阵的序号和子带矩阵的数目
	int Zs_bw;
	int Index_up=0;
	for(Band_Ser=0;Band_Ser < RD_num;Band_Ser++)
	{
		for(Band_Num=0;Band_Num < nxy-Band_Ser;Band_Num++)
		{
			Zs_bw=int(sqrt(RD_num*RD_num-Band_Ser*Band_Ser));
			for(i=Band_Num*nxy;i < nxy+Band_Num*nxy;i++)
			{
				for(j=Band_Num*nxy+Band_Ser*nxy;j < nxy+Band_Num*nxy+Band_Ser*nxy;j++)
				{
					if(abs((i-Band_Num*nxy)-(j-Band_Num*nxy-Band_Ser*nxy)) <= Zs_bw)
					{
						result[i]+=Zs[Index_up]*vec[j];//右半部分(包括对角线)与列向量相乘
						if(Band_Ser != 0)              //左半部分与列向量相乘
						{
							result[j]+=Zs[Index_up]*vec[i];
						}
						Index_up++;
					}
				}
			}
		}
	}
	for(i=0;i < nxy*nxy-RD_num*nxy;i++)
	{
		result[i]+=Zs[Index_up]*vec[i+RD_num*nxy];//右半部分
		result[i+RD_num*nxy]+=Zs[Index_up]*vec[i];//左半部分
		Index_up++;
	}
}
void resolve_matrix::NearMultiplyVector(char a, vector <Zmn_Cell> &zs, vector<Complex> &vec, vector<Complex> &result)
{//zs只需存储上一半
	int Z_num;
	int nxy = vec.size();
	result.assign(nxy,0.0);

	if(a == ' ')
	{
		for(Z_num=0; Z_num<zs.size(); Z_num++)
		{
			if ( zs[Z_num].Zmn_i == zs[Z_num].Zmn_j )
			{
				result[zs[Z_num].Zmn_i] += zs[Z_num].Zmn_z*vec[zs[Z_num].Zmn_j];
			}
			else
			{
				result[zs[Z_num].Zmn_i] += zs[Z_num].Zmn_z*vec[zs[Z_num].Zmn_j];
				result[zs[Z_num].Zmn_j] += zs[Z_num].Zmn_z*vec[zs[Z_num].Zmn_i];
			}
			//outresult<<zs[Z_num].Zmn_i<<"   "<<zs[Z_num].Zmn_j<<"   "<<zs[Z_num].Zmn_z<<endl;	
		}
	}
	else
	{//矩阵转置共轭的矩阵向量积
		for(Z_num=0; Z_num<zs.size(); Z_num++)
		{
			if ( zs[Z_num].Zmn_i == zs[Z_num].Zmn_j )
			{
				result[zs[Z_num].Zmn_i] += conj(zs[Z_num].Zmn_z)*vec[zs[Z_num].Zmn_j];
			}
			else
			{
				result[zs[Z_num].Zmn_i] += conj(zs[Z_num].Zmn_z)*vec[zs[Z_num].Zmn_j];
				result[zs[Z_num].Zmn_j] += conj(zs[Z_num].Zmn_z)*vec[zs[Z_num].Zmn_i];
			}
		}
	}
}
void resolve_matrix::Two_D_Vector_NearMultiplyVector(char a, vector <Zmn_Cell> &zs, vector<Complex> &vec, vector<Complex> &result)
{//zs要求全部存储
	int Z_num;
	int nxy = vec.size();
	result.assign(nxy,0.0);

	if(a == ' ')
	{
		for(Z_num=0; Z_num<zs.size(); Z_num++)
		{
			result[zs[Z_num].Zmn_i] += zs[Z_num].Zmn_z*vec[zs[Z_num].Zmn_j];
		}
	}
	else
	{//矩阵转置共轭的矩阵向量积
		for(Z_num=0; Z_num<zs.size(); Z_num++)
		{
			result[zs[Z_num].Zmn_j] += conj(zs[Z_num].Zmn_z)*vec[zs[Z_num].Zmn_i];
		}
	}
}
void resolve_matrix::Multiply_Full_Matrix_Vector(char a, Vec_Com &zm, Vec_Com &vec, Vec_Com &result )
{
	int i,j,Num;
	Num = vec.size();
	result.assign(Num,0.0);

	if(a == ' ')
	{
		for ( i=0; i<Num; i++ )
		{
			for ( j=0; j<Num; j++ )
			{
				result[i] += zm[i*Num+j]*vec[j];
			}
		}
	}
	else
	{
		for ( i=0; i<Num; i++ )
		{
			for ( j=0; j<Num; j++ )
			{
				result[i] += conj(zm[j*Num+i])*vec[j];
			}
		}
	}
}
//The GPBi-CG algorithm with preconditioning :P143
void resolve_matrix::GPBI_CG_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	
	Complex *vector_r=new Complex[2*ndg];
	Complex *vector_r0=new Complex[2*ndg];
	Complex *vector_u=new Complex[2*ndg];
	Complex *vector_z=new Complex[2*ndg];
	Complex *vector_p=new Complex[2*ndg];
	Complex *vector_t=new Complex[2*ndg];
	Complex *vector_y=new Complex[2*ndg];
	Complex *vector_w=new Complex[2*ndg];
	Complex *vector_s=new Complex[2*ndg];
	Complex *vector_v,*vector_q;
		
	Complex rou,rou0,alfa,beta,value_u1,value_u2,value_u3,value_u4,
		value_u5,value_gama,value_eta,value_tao,temp;
	double Norm;

	//预处理矩阵: inv(K)*A, inv(K)*b
	//预处理矩阵为对角阵
	Complex *Diag_Matrix_Inverse=Preconditoning_Matrix_Diag_Inverse(asmn,b1mn,ndg,lbw,fbw);
	Complex *Diag_M1=new Complex[ndg];
	Complex *Diag_M2=new Complex[ndg];
	for (i=0; i<ndg; i++)
	{
		Diag_M1[i] = Diag_Matrix_Inverse[i];
		Diag_M2[i] = Diag_Matrix_Inverse[i+ndg];
	}
	Complex *asmn_M = Multi_DiagM_BandM(Diag_M1,asmn,ndg,fbw);
	Complex *bsmn_M = Multi_DiagM_BandM(Diag_M1,bsmn,ndg,fbw);
	Complex *a1mn_M = Multi_DiagM_BandM(Diag_M2,a1mn,ndg,lbw);
	Complex *b1mn_M = Multi_DiagM_BandM(Diag_M2,b1mn,ndg,lbw);
	Complex *aabmn_M = Multi_DiagM_ABM(Diag_M1,aabmn,ndg,ncg,dbw,fbw);
	Complex *babmn_M = Multi_DiagM_ABM(Diag_M1,babmn,ndg,ncg,dbw,fbw);
	Complex *bm_M =new Complex[ndg];
	for (i=0; i<ndg; i++)
	{
		bm_M[i] = bm[i] * Diag_M1[i];
	}


	double norm_0=sqrt(L2_Norm(bm_M,ndg));
	Assign_Vecor(bm_M,vector_r,ndg);
	Assign_Vecor(bm_M,vector_r0,ndg);

	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,2*ndg)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0,vector_r,2*ndg);
		if ( iter == 1 )
		{
			Assign_Vecor(vector_r,vector_p,2*ndg);
			vector_q=Multi_M_V(asmn_M,bsmn_M,a1mn_M,b1mn_M,aabmn_M,babmn_M,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp=Multi_vt_v(vector_r0,vector_q,2*ndg);
			alfa=rou/temp;
			for (i=0; i<2*ndg; i++)
			{
				vector_t[i] = vector_r[i] - alfa * vector_q[i];
			}
			vector_v=Multi_M_V(asmn_M,bsmn_M,a1mn_M,b1mn_M,aabmn_M,babmn_M,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for (i=0; i<2*ndg; i++)
			{
				vector_y[i] = alfa * vector_q[i] - vector_r[i];
			}
			value_u2 = Multi_vt_v(vector_v,vector_t,2*ndg);
			value_u5 = Multi_vt_v(vector_v,vector_v,2*ndg);
			value_gama = value_u2 / value_u5;
			value_eta = 0.0;
		}
		else
		{
			beta = (rou / rou0) * (alfa / value_gama);
			for (i=0; i<2*ndg; i++)
			{
				vector_w[i] = vector_v[i] + beta * vector_q[i];
				vector_p[i] = vector_r[i] + beta * (vector_p[i] - vector_u[i]);
			}
			
			vector_q=Multi_M_V(asmn_M,bsmn_M,a1mn_M,b1mn_M,aabmn_M,babmn_M,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp = Multi_vt_v(vector_r0,vector_q,2*ndg);
			alfa = rou / temp;
			for (i=0; i<2*ndg; i++)
			{
				vector_s[i] = vector_t[i] - vector_r[i];
				vector_t[i] = vector_r[i] - alfa * vector_q[i];
			}
			vector_v=Multi_M_V(asmn_M,bsmn_M,a1mn_M,b1mn_M,aabmn_M,babmn_M,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for (i=0; i<2*ndg; i++)
			{
				vector_y[i] = vector_s[i] - alfa * (vector_w[i] - vector_q[i]);
			}
			value_u1 = Multi_vt_v(vector_y,vector_y,2*ndg);
			value_u2 = Multi_vt_v(vector_v,vector_t,2*ndg);
			value_u3 = Multi_vt_v(vector_y,vector_t,2*ndg);
			value_u4 = Multi_vt_v(vector_v,vector_y,2*ndg);
			value_u5 = Multi_vt_v(vector_v,vector_v,2*ndg);
			value_tao = value_u5 * value_u1 - conj(value_u4) * value_u4;
			value_gama = (value_u1 * value_u2 - value_u3 * value_u4) / value_tao;
			value_eta = (value_u5 * value_u3 - conj(value_u4) * value_u2) / value_tao;
		}
		
		for (i=0; i<2*ndg; i++)
		{
			vector_u[i] = value_gama * vector_q[i] + value_eta * (vector_s[i] + beta * vector_u[i]);
			vector_z[i] = value_gama * vector_r[i] + value_eta * vector_z[i] - alfa * vector_u[i];
			rex[i] += alfa * vector_p[i] + vector_z[i];
			vector_r[i] = vector_t[i] - value_eta * vector_y[i] - value_gama * vector_v[i];
		}

		Norm=sqrt(L2_Norm(vector_r,2*ndg))/norm_0;
		cout<<"GPBi-CG steps "<<iter<<",  res = "<<Norm<<endl;
		
		if ( Norm < eps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total GPBi-CG Steps = "<<iter<<", to achieve "<<Norm<<endl;

	delete [] vector_r;
	delete [] vector_r0;
	delete [] vector_u;
	delete [] vector_z;
	delete [] vector_p;
	delete [] vector_t;
	delete [] vector_y;
	delete [] vector_w;
	delete [] vector_s;
	delete [] vector_v;
	delete [] vector_q;
	delete [] Diag_Matrix_Inverse;
	delete [] Diag_M1;
	delete [] Diag_M2;
	delete [] asmn_M;
	delete [] bsmn_M;
	delete [] a1mn_M;
	delete [] b1mn_M;
	delete [] aabmn_M;
	delete [] babmn_M;
}
//The GPBi-CG algorithm without preconditioning :P143
void resolve_matrix::GPBI_CG(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	
	Complex *vector_r=new Complex[2*ndg];
	Complex *vector_r0=new Complex[2*ndg];
	Complex *vector_u=new Complex[2*ndg];
	Complex *vector_z=new Complex[2*ndg];
	Complex *vector_p=new Complex[2*ndg];
	Complex *vector_t=new Complex[2*ndg];
	Complex *vector_y=new Complex[2*ndg];
	Complex *vector_w=new Complex[2*ndg];
	Complex *vector_s=new Complex[2*ndg];
	Complex *vector_v,*vector_q;
		
	Complex rou,rou0,alfa,beta,value_u1,value_u2,value_u3,value_u4,
		value_u5,value_gama,value_eta,value_tao,temp;
	double Norm;

	double norm_0=sqrt(L2_Norm(bm,ndg));
	Assign_Vecor(bm,vector_r,ndg);
	Assign_Vecor(bm,vector_r0,ndg);

	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,2*ndg)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0,vector_r,2*ndg);
		if ( iter == 1 )
		{
			Assign_Vecor(vector_r,vector_p,2*ndg);
			vector_q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp=Multi_vt_v(vector_r0,vector_q,2*ndg);
			alfa=rou/temp;
			for (i=0; i<2*ndg; i++)
			{
				vector_t[i] = vector_r[i] - alfa * vector_q[i];
			}
			vector_v=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for (i=0; i<2*ndg; i++)
			{
				vector_y[i] = alfa * vector_q[i] - vector_r[i];
			}
			value_u2 = Multi_vt_v(vector_v,vector_t,2*ndg);
			value_u5 = Multi_vt_v(vector_v,vector_v,2*ndg);
			value_gama = value_u2 / value_u5;
			value_eta = 0.0;
		}
		else
		{
			beta = (rou / rou0) * (alfa / value_gama);
			for (i=0; i<2*ndg; i++)
			{
				vector_w[i] = vector_v[i] + beta * vector_q[i];
				vector_p[i] = vector_r[i] + beta * (vector_p[i] - vector_u[i]);
			}
			
			vector_q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_p,CX,DX,fbw,lbw,dbw,ndg,ncg);
			temp = Multi_vt_v(vector_r0,vector_q,2*ndg);
			alfa = rou / temp;
			for (i=0; i<2*ndg; i++)
			{
				vector_s[i] = vector_t[i] - vector_r[i];
				vector_t[i] = vector_r[i] - alfa * vector_q[i];
			}
			vector_v=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_t,CX,DX,fbw,lbw,dbw,ndg,ncg);
			for (i=0; i<2*ndg; i++)
			{
				vector_y[i] = vector_s[i] - alfa * (vector_w[i] - vector_q[i]);
			}
			value_u1 = Multi_vt_v(vector_y,vector_y,2*ndg);
			value_u2 = Multi_vt_v(vector_v,vector_t,2*ndg);
			value_u3 = Multi_vt_v(vector_y,vector_t,2*ndg);
			value_u4 = Multi_vt_v(vector_v,vector_y,2*ndg);
			value_u5 = Multi_vt_v(vector_v,vector_v,2*ndg);
			value_tao = value_u5 * value_u1 - conj(value_u4) * value_u4;
			value_gama = (value_u1 * value_u2 - value_u3 * value_u4) / value_tao;
			value_eta = (value_u5 * value_u3 - conj(value_u4) * value_u2) / value_tao;
		}
		
		for (i=0; i<2*ndg; i++)
		{
			vector_u[i] = value_gama * vector_q[i] + value_eta * (vector_s[i] + beta * vector_u[i]);
			vector_z[i] = value_gama * vector_r[i] + value_eta * vector_z[i] - alfa * vector_u[i];
			rex[i] += alfa * vector_p[i] + vector_z[i];
			vector_r[i] = vector_t[i] - value_eta * vector_y[i] - value_gama * vector_v[i];
		}

		Norm=sqrt(L2_Norm(vector_r,2*ndg))/norm_0;
		cout<<"GPBi-CG steps "<<iter<<",  res = "<<Norm<<endl;
		
		if ( Norm < eps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total GPBi-CG Steps = "<<iter<<", to achieve "<<Norm<<endl;

	delete [] vector_r;
	delete [] vector_r0;
	delete [] vector_u;
	delete [] vector_z;
	delete [] vector_p;
	delete [] vector_t;
	delete [] vector_y;
	delete [] vector_w;
	delete [] vector_s;
	delete [] vector_v;
	delete [] vector_q;
}
//Bi-CGSTAB algorithm with preconditioning :P138
void resolve_matrix::BI_CGSTAB_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	Complex *vector_r=new Complex[2*ndg];
	Complex *vector_p=new Complex[2*ndg];
	Complex *vector_pt=new Complex[2*ndg];
	Complex *vector_r0p=new Complex[2*ndg];
	Complex *vector_s=new Complex[2*ndg];
	Complex *vector_st=new Complex[2*ndg];
	Complex *vector_t;
	Complex *vector_v;
	
	Complex rou,rou0,alfa,beta,value_w,temp1,temp2;
	double Norm;
	double norm_0=sqrt(L2_Norm(bm,ndg));

	Assign_Vecor(bm,vector_r,ndg);
	Assign_Vecor(bm,vector_r0p,ndg);
	
	//预处理矩阵为对角阵
	Complex *Diag_Matrix_Inverse=Preconditoning_Matrix_Diag_Inverse(asmn,b1mn,ndg,lbw,fbw);
	
	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,2*ndg)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0p,vector_r,2*ndg);
		if ( iter == 1)
		{
			Assign_Vecor(vector_r,vector_p,2*ndg);
		}
		else
		{
			beta=(rou/rou0)*(alfa/value_w);
			for (i=0; i<2*ndg; i++)
			{
				vector_p[i]=vector_r[i]+beta*(vector_p[i]-value_w*vector_v[i]);
			}
		}
		for (i=0; i<2*ndg; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		vector_v=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_pt,CX,DX,fbw,lbw,dbw,ndg,ncg);
		temp1=Multi_vt_v(vector_r0p,vector_v,2*ndg);
		alfa=rou/temp1;
		for (i=0; i<2*ndg; i++)
		{
			vector_s[i]=vector_r[i]-alfa*vector_v[i];
		}

		Norm=sqrt(L2_Norm(vector_s,2*ndg));
		if ( Norm < 1.0e-6 )
		{
			for (i=0; i<2*ndg; i++)
			{
				rex[i]=rex[i]+alfa*vector_pt[i];
			}
			break;
		}
		
		for (i=0; i<2*ndg; i++)
		{
			vector_st[i]=Diag_Matrix_Inverse[i]*vector_s[i];
		}
		vector_t=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_st,CX,DX,fbw,lbw,dbw,ndg,ncg);
		temp1=Multi_vt_v(vector_t,vector_s,2*ndg);
		temp2=Multi_vt_v(vector_t,vector_t,2*ndg);
		value_w=temp1/temp2;
		for (i=0; i<2*ndg; i++)
		{
			rex[i] += alfa*vector_pt[i]+value_w*vector_st[i];
			vector_r[i]=vector_s[i]-value_w*vector_t[i];
		}
		
		Norm=sqrt(L2_Norm(vector_r,2*ndg))/norm_0;
		cout<<"Bi-CGSTAB_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		if ( Norm < eps )
		{
			break;
		}

		rou0=rou;
	}
	cout<<"Total Bi-CGSTAB_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	delete [] Diag_Matrix_Inverse;
	delete [] vector_r;
	delete [] vector_p;
	delete [] vector_pt;
	delete [] vector_r0p;
	delete [] vector_s;
	delete [] vector_st;
	delete [] vector_t;
	delete [] vector_v;
}
//CGS:P104
void resolve_matrix::CGS(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	Complex *vector_r0=new Complex[2*ndg];
	Complex *vector_r=new Complex[2*ndg];
	Complex *vector_u=new Complex[2*ndg];
	Complex *vector_p=new Complex[2*ndg];
	Complex *vector_q=new Complex[2*ndg];
	Complex *vector_pt=new Complex[2*ndg];
	Complex *vector_ut=new Complex[2*ndg];
	Complex *vector_uq=new Complex[2*ndg];
	
	Complex *vector_vp,*temp1;
	Complex alfa,beta,rou0,rou,temp2;

	double Norm;
	double norm_0=sqrt(L2_Norm(bm,ndg));
	//double eps=1.0e-2;
	
	//预处理矩阵为对角阵
	Complex *Diag_Matrix_Inverse=Preconditoning_Matrix_Diag_Inverse(asmn,b1mn,ndg,lbw,fbw);
	
	Assign_Vecor(bm,vector_r0,ndg);
	Assign_Vecor(bm,vector_r,ndg);

	int iter=0;
	while ( true )
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,2*ndg)) / norm_0<<endl;
		iter++;
		rou=Multi_vt_v(vector_r0,vector_r,2*ndg);
		if ( iter == 1)
		{
			Assign_Vecor(vector_r0,vector_u,2*ndg);
			Assign_Vecor(vector_u,vector_p,2*ndg);
		}
		else
		{
			beta=rou/rou0;
			for (i=0; i<2*ndg; i++)
			{
				vector_u[i]=vector_r[i]+beta*vector_q[i];
				vector_p[i]=vector_u[i]+beta*(vector_q[i]+beta*vector_p[i]);
			}
		}
		
		for (i=0; i<2*ndg; i++)
		{
			vector_pt[i]=Diag_Matrix_Inverse[i]*vector_p[i];
		}
		vector_vp=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_pt,CX,DX,fbw,lbw,dbw,ndg,ncg);
		temp2=Multi_vt_v(vector_r0,vector_vp,2*ndg);
		alfa=rou/temp2;
		for (i=0; i<2*ndg; i++)
		{
	    	vector_q[i]=vector_u[i]-alfa*vector_vp[i];
	    	vector_uq[i]=vector_u[i]+vector_q[i];
	        vector_ut[i]=Diag_Matrix_Inverse[i]*vector_uq[i];
		}
		temp1=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,vector_ut,CX,DX,fbw,lbw,dbw,ndg,ncg);
    	for (i=0; i<2*ndg; i++)
		{
    		rex[i] += alfa*vector_ut[i];
    		vector_r[i] -= alfa*temp1[i];
		}

		Norm=sqrt(L2_Norm(vector_r,2*ndg))/norm_0;
		cout<<"CGS"<<" Steps = "<<iter<<",  res = "<<Norm<<endl;
		out_resual<<iter<<"  "<<sqrt(L2_Norm(vector_r,2*ndg))<<endl;
		if (Norm < eps)
		{
			break;
		}
		
		rou0=rou;		
	}
	cout<<"Total CGS Steps = "<<iter<<", to achieve "<<Norm<<endl;
	
	delete [] Diag_Matrix_Inverse;
	delete [] vector_r0;
	delete [] vector_r;
	delete [] vector_u;
	delete [] vector_p;
	delete [] vector_q;
	delete [] vector_pt;
	delete [] vector_ut;
	delete [] vector_uq;
	delete [] vector_vp;
	delete [] temp1;
}
//预处理的双共轭梯度法：不收敛
void resolve_matrix::Bi_CG_Pre(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
			Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
			Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
			Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,
			int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	double Norm;
	Complex rho_1, rho_2, alpha, beta, temp;
	
	Complex *z=new Complex[2*ndg];
	Complex *ztilde=new Complex[2*ndg];
	Complex *p=new Complex[2*ndg];
	Complex *ptilde=new Complex[2*ndg];
	Complex *r=new Complex[2*ndg];
	Complex *rtilde=new Complex[2*ndg];
	Complex *q,*qtilde;
	
	double normb = sqrt(L2_Norm(bm,ndg));
	Assign_Vecor(bm,r,ndg);
	Assign_Vecor(bm,rtilde,ndg);
	//预处理矩阵为对角阵
	Complex *Diag_Matrix_Inverse = Preconditoning_Matrix_Diag_Inverse(asmn,b1mn,ndg,lbw,fbw);
	Complex *Diag_Tran_Matrix_Inverse = new Complex[2*ndg];
	for (i=0; i<2*ndg; i++)
	{
		Diag_Tran_Matrix_Inverse[i] = 1.0 / conj(1.0 / Diag_Matrix_Inverse[i]);
	}

	int iter = 0;
	while ( true ) 
	{
		out_resual<<iter<<"  "<<sqrt(L2_Norm(r,2*ndg)) / normb<<endl;
		iter++;
		for (i=0; i<2*ndg; i++)
		{
			z[i] = Diag_Matrix_Inverse[i] * r[i];
			ztilde[i] = Diag_Tran_Matrix_Inverse[i] * rtilde[i];
		}
		rho_1 = Multi_vt_v(rtilde,z,2*ndg);//原文有问题，原文中为 multi_vec_vec(ztilde,z,rho_1);
		if (iter == 1)
		{
			Assign_Vecor(z,p,2*ndg);
			Assign_Vecor(ztilde,ptilde,2*ndg);		
		}
		else 
		{
			beta = rho_1 / rho_2;
			for (i=0; i<2*ndg; i++)
			{
				p[i] = z[i] + beta * p[i];
				ptilde[i] = ztilde[i] + conj(beta) * ptilde[i];
			}
		}
		q = Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
		qtilde = Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,ptilde,CX,DX,fbw,lbw,dbw,ndg,ncg);
		temp = Multi_vt_v(ptilde,q,2*ndg);
		alpha = rho_1 / temp;
		for (i=0; i<2*ndg; i++)
		{
			rex[i] += alpha * p[i];
			r[i] -= alpha * q[i];
			rtilde[i] -= conj(alpha) * qtilde[i];
		}
				
		Norm = sqrt(L2_Norm(r,2*ndg)) / normb;
		cout<<"Bi_CG_Pre steps "<<iter<<",  res = "<<Norm<<endl;
		
		if (Norm < eps ) 
		{
			break;
		}
		rho_2 = rho_1;
	}
	
	cout<<"Total Bi-CG_Preconditioning Steps = "<<iter<<", to achieve "<<Norm<<endl;

	delete [] z;
	delete [] ztilde;
	delete [] p;
	delete [] ptilde;
	delete [] r;
	delete [] rtilde;
	delete [] q;
	delete [] qtilde;
	delete [] Diag_Matrix_Inverse;
	delete [] Diag_Tran_Matrix_Inverse;
}
void resolve_matrix::Preconditioning_CGM(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
						 Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
						 Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
						 Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,
						 int lbw,int fbw,int dbw)
{	
}
Complex *resolve_matrix::Preconditoning_Matrix_Diag_Inverse(Complex *asmn,Complex *b1mn,int ndg,int lbw,int fbw)
{
	Complex *Diag_Matrix=new Complex[2*ndg];
	int i;
	for (i=0;i < ndg;i++)
	{
		if (i <= fbw)
			Diag_Matrix[i]=1.0/asmn[i*(2*fbw+1)+i];
		else
			Diag_Matrix[i]=1.0/asmn[i*(2*fbw+1)+fbw];
		if (i <= lbw)
			Diag_Matrix[i+ndg]=1.0/b1mn[i*(2*lbw+1)+i];
		else
			Diag_Matrix[i+ndg]=1.0/b1mn[i*(2*lbw+1)+lbw];
	}
	
	return Diag_Matrix;
}
//根据不完全 Cholesky 分解所得的下三角和上三角矩阵解方程
Complex *resolve_matrix::Resolve_Cholesky_LU(Complex *C_Left_Up,Complex *C_Left_Down,Complex *C_Light_Down,Complex *Bm,int ndg,int lbw,int fbw)
{
	Complex *result=new Complex[2*ndg];
	Complex *re=new Complex[2*ndg];
	int i,j;
	//下三角矩阵部分 回代过程
	//上半部分
	result[0]=Bm[0]/C_Left_Up[0];
	for (i=1;i < ndg;i++)
	{
		if (i <= fbw)
		{
			for (j=0;j < i;j++)
			{
				result[i]+=C_Left_Up[i*(fbw+1)+j]*result[j];
			}
			result[i]=(Bm[i]-result[i])/C_Left_Up[i*(fbw+1)+i];
		}
		else
		{
			for (j=i-fbw;j < i;j++)
			{
				result[i]+=C_Left_Up[i*(fbw+1)+j-(i-fbw)]*result[j];
			}
			result[i]=(Bm[i]-result[i])/C_Left_Up[i*(fbw+1)+fbw];
		}		
	}
	//下半部分
	for (i=0;i < ndg;i++)
	{
		if (i <= lbw)
		{
			for (j=0;j <= i+lbw;j++)
			{
				result[i+ndg]+=C_Left_Down[i*(2*lbw+1)+j]*result[j];
			}
			for (j=0;j < i;j++)
			{
				result[i+ndg]+=C_Light_Down[i*(lbw+1)+j]*result[j+ndg];
			}
			result[i+ndg]=(Bm[i+ndg]-result[i+ndg])/C_Light_Down[i*(lbw+1)+i];
		}
		else
		{
			for (j=i-lbw;j <= i+lbw && j < ndg;j++)
			{
				result[i+ndg]+=C_Left_Down[i*(2*lbw+1)+j-(i-lbw)]*result[j];
			}
			for (j=i-lbw;j < i;j++)
			{
				result[i+ndg]+=C_Light_Down[i*(lbw+1)+j-(i-lbw)]*result[j+ndg];
			}
			result[i+ndg]=(Bm[i+ndg]-result[i+ndg])/C_Light_Down[i*(lbw+1)+lbw];
		}
	}
	//上三角矩阵部分(即为下三角矩阵的转置共轭) 回代过程
	//下半部分
	re[2*ndg-1]=result[2*ndg-1]/conj(C_Light_Down[(ndg-1)*(lbw+1)+lbw]);
	for (i=ndg-2;i >= 0;i--)
	{
		if (i < ndg-lbw)
		{
			for (j=i+1;j <= i+lbw;j++)
			{
				if (j <= lbw)
				{
					re[i+ndg]+=conj(C_Light_Down[j*(lbw+1)+i])*re[j+ndg];
				}
				else
				{
					re[i+ndg]+=conj(C_Light_Down[j*(lbw+1)+i-(j-lbw)])*re[j+ndg];
				}
			}			
		}
		else
		{
			for (j=i+1;j < ndg;j++)
			{
				if (j <= lbw)
				{
					re[i+ndg]+=conj(C_Light_Down[j*(lbw+1)+i])*re[j+ndg];
				}
				else
				{
					re[i+ndg]+=conj(C_Light_Down[j*(lbw+1)+i-(j-lbw)])*re[j+ndg];
				}
			}
		}
		if (i <= lbw)
		{
			re[i+ndg]=(result[i+ndg]-re[i+ndg])/conj(C_Light_Down[i*(lbw+1)+i]);
		}
		else
		{
			re[i+ndg]=(result[i+ndg]-re[i+ndg])/conj(C_Light_Down[i*(lbw+1)+lbw]);
		}	
	}
	//上半部分
	for (i=ndg-1;i >= 0;i--)
	{
		if (i <= lbw)
		{
			for (j=0;j <= i+lbw;j++)
			{
				if (j <= lbw)
				{
					re[i]+=conj(C_Left_Down[j*(2*lbw+1)+i])*re[j+ndg];
				}
				else
				{
					re[i]+=conj(C_Left_Down[j*(2*lbw+1)+i-(j-lbw)])*re[j+ndg];
				}
			}
		}
		else
		{
			for (j=i-lbw;j <= i+lbw && j < ndg;j++)
			{
				if (j <= lbw)
				{
					re[i]+=conj(C_Left_Down[j*(2*lbw+1)+i])*re[j+ndg];
				}
				else
				{
					re[i]+=conj(C_Left_Down[j*(2*lbw+1)+i-(j-lbw)])*re[j+ndg];
				}
			}
		}
		if (i < ndg-fbw)
		{
			for (j=i+1;j <= i+fbw;j++)
			{
				if (j <= fbw)
				{
					re[i]+=conj(C_Left_Up[j*(fbw+1)+i])*re[j];
				}
				else
				{
					re[i]+=conj(C_Left_Up[j*(fbw+1)+i-(j-fbw)])*re[j];
				}
			}
		}
		else
		{
			for (j=i+1;j < ndg;j++)
			{
				if (j <= fbw)
				{
					re[i]+=conj(C_Left_Up[j*(fbw+1)+i])*re[j];
				}
				else
				{
					re[i]+=conj(C_Left_Up[j*(fbw+1)+i-(j-fbw)])*re[j];
				}
			}
		}
		if (i <= fbw)
		{
			re[i]=(result[i]-re[i])/conj(C_Left_Up[i*(fbw+1)+i]);
		}
		else
		{
			re[i]=(result[i]-re[i])/conj(C_Left_Up[i*(fbw+1)+fbw]);
		}
	}

	delete [] result;

	return re;
}
void resolve_matrix::F_GMRES_K(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw,double eps_m,double eps_p)
{
	int i,j,L;
	rex=new Complex[2*ndg];
	int Gmres_K=2*ndg/10;
	int maxit=2*ndg;
	double r2,c;
	Complex temp,R_temp,s,rhs;
	Complex *x0,*r,*q,*w,*y,*rex_temp,*z;
	Vec_Com Q,Zm,R,RHS,S;
	Vec_Dou C;
	q=new Complex[2*ndg];
	x0=new Complex[2*ndg];
	z=new Complex[2*ndg];
	double norm_0=sqrt(L2_Norm(bm,ndg));
	//预处理矩阵为对角阵
	Complex *Diag_Matrix_Inverse=Preconditoning_Matrix_Diag_Inverse(asmn,b1mn,ndg,lbw,fbw);
	//
	////////预处理矩阵为不完全Cholesky分解矩阵
	//不完全Cholesky分解，生成下三角矩阵
/*	Complex *Diag=Creat_Diag(asmn,a1mn,b1mn,ndg,lbw,fbw);
	Complex *D_Up=new Complex[ndg];
	Complex *D_Down=new Complex[ndg];
    for(i=0;i<ndg;i++)
	{
		D_Up[i]=Diag[i];
		D_Down[i]=Diag[i+ndg];
	}
	Complex *U_Tria_asmn=Creat_triang(asmn,D_Up,ndg,fbw);
	Complex *U_Tria_b1mn=Creat_triang(b1mn,D_Down,ndg,lbw);
	Diag=Creat_Inv_Diag(Diag,ndg);
	for(i=0;i<ndg;i++)
	{
		D_Up[i]=Diag[i];
		D_Down[i]=Diag[i+ndg];
	}
	//C
	Complex *C_Left_Up=Multi_HalfBM_Diag(U_Tria_asmn,D_Up,ndg,fbw);
	Complex *C_Light_Down=Multi_HalfBM_Diag(U_Tria_b1mn,D_Down,ndg,lbw);
	Complex *C_Left_Down=Multi_BM_Diag(a1mn,D_Up,ndg,lbw);
	/////////////////////////////////////////////////*/
	int iter=0;
	int ok=0;

	while(iter < maxit)
	{
		
		//Step 1
		for(j=0;j<2*ndg;j++)
			x0[j]=rex[j];
		r=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,x0,CX,DX,fbw,lbw,dbw,ndg,ncg);
		//multiply(am,x0,r);
		for (j=0;j<ndg;j++)
		{
			r[j]=bm[j]-r[j];
			r[j+ndg]=-r[j+ndg];
		}
		r2=sqrt(L2_Norm(r,2*ndg));
		//cout<<r2;
		//Step 2
		for(j=0;j<2*ndg;j++)
		{
			q[j]=r[j]/r2;
			Q.push_back(q[j]);
		}
		rhs=r2;
		RHS.push_back(rhs);
		for(i=0;i<Gmres_K;i++)
		{
			//Step 3
			for(j=0;j<2*ndg;j++)
			{
				q[j]=Q[i*2*ndg+j];
			}
			//
			for (j=0;j < 2*ndg;j++)
			{
				z[j]=Diag_Matrix_Inverse[j]*q[j];
			}
			//不完全 Cholesky 分解预处理阵
			//z=Resolve_Cholesky_LU(C_Left_Up,C_Left_Down,C_Light_Down,q,ndg,lbw,fbw);
			//z=GMRES_K(a1mn,b1mn,asmn,bsmn,aabmn,babmn,CX,DX,q,ndg,ncg,lbw,fbw,dbw,eps_m);
			for(j=0;j<2*ndg;j++)
			{
				Zm.push_back(z[j]);
				//z[j]=Zm[i*2*ndg+j];
			}
			w=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,z,CX,DX,fbw,lbw,dbw,ndg,ncg);
			//multiply(am,q,w);
			
			//Step 4
			for(L=0;L<=i;L++)//
			{
				temp=0.0;
				for(j=0;j<2*ndg;j++)
					temp+=conj(w[j])*Q[L*2*ndg+j];
				R.push_back(temp);
				for(j=0;j<2*ndg;j++)
				{
					w[j]-=temp*Q[L*2*ndg+j];
					//cout<<w[j]<<endl;
				}
			}
			//Step 5
			R_temp=sqrt(L2_Norm(w,2*ndg));
			//cout<<R_temp;
			//R.push_back(temp);
			for(j=0;j<2*ndg;j++)
			{
				q[j]=w[j]/R_temp;
				Q.push_back(q[j]);
			}
			//Step 6
			for(L=0;L<=i-1;L++)
			{
				temp=C[L]*R[((i+1)*i)/2+L]+conj(S[L])*R[((i+1)*i)/2+L+1];
				R[((i+1)*i)/2+L+1]=S[L]*R[((i+1)*i)/2+L]-C[L]*R[((i+1)*i)/2+L+1];
				R[((i+1)*i)/2+L]=temp;
			}
			//Step 7
			temp=R[((i+1)*i)/2+i];
			if(sqrt(temp.real()*temp.real()+temp.imag()*temp.imag()) <= 1.0e-16)
			{
				c=0.0;
				C.push_back(c);
				s=1.0;
				S.push_back(s);
			}
			else
			{
				temp=R_temp/R[((i+1)*i)/2+i];
				c=1.0/sqrt(1.0+(temp.real()*temp.real()+temp.imag()*temp.imag()));
				C.push_back(c);
				s=C[i]*temp;
				S.push_back(s);
			}
			
			//Step 8
			temp=C[i]*RHS[i];
			rhs=S[i]*RHS[i];
			RHS.push_back(rhs);
			RHS[i]=temp;
			R[((i+1)*i)/2+i]=C[i]*R[((i+1)*i)/2+i]+conj(S[i])*R_temp;
			R_temp=0.0;
			
			//Step 9
			
			
			y=Resolve_Tri_Matrix_GMRES(R,RHS,i+1);
						
			//Step 10
			
			rex_temp=multiply_mn_GMRES(Zm,y,2*ndg,i+1);
			
			for(j=0;j<2*ndg;j++)
				rex[j]=x0[j]+rex_temp[j];
			//for (j=0;j<m;j++)
			//{
			//	cout<<rex[j]<<"  "<<i<<endl;
			//}
			//Step 11
			rex_temp=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,rex,CX,DX,fbw,lbw,dbw,ndg,ncg);
			//multiply(am,rex,rex_temp);
			
			for (j=0;j<ndg;j++)
			{
				r[j]=bm[j]-rex_temp[j];
				r[j+ndg]=-rex_temp[j+ndg];
			}
	    	//for(j=0;j<m;j++)
		    //	r[j]=bm[j]-rex_temp[j];
			r2=sqrt(L2_Norm(r,2*ndg))/norm_0;
			cout<<"FGMRES"<<" Steps="<<iter*Gmres_K+i+1<<", iter = "<<iter<<",  res="<<r2<<endl;
			out_resual<<iter*Gmres_K+i+1<<"   "<<r2<<endl;
			if(r2 <= eps_p)
			{
				ok=1;
				break;
			}
			//Step 12
		}
		Zm.clear();
		Q.clear();
		R.clear();
		RHS.clear();
		C.clear();
		S.clear();
		if(ok == 1)
			break;
		iter++;
		
	}
	if(iter >= maxit)
		iter-=1;
	cout<<"Total FGMRES("<<Gmres_K<<") Steps = "<<iter*Gmres_K+i+1<<", iter = "<<iter<<", to achieve "<<r2<<endl;

	delete [] x0;
	delete [] r;
	delete [] q;
	delete [] w;
	delete [] y;
	delete [] rex_temp;
	delete [] z;
	delete [] Diag_Matrix_Inverse;
}
Complex *resolve_matrix::GMRES_K(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,		
			Complex *aabmn,Complex *babmn,Vec_Dou &CX,Vec_Dou &DX,
			Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw,double tol)
{
	int i,j,L;
	Complex *re=new Complex[2*ndg];
	int Gmres_K=2*ndg/10;
	//double tol=1.0e-2;
	int maxit=2*ndg;
	double r2,c;
	Complex temp,R_temp,s,rhs;
	Complex *x0,*r,*q,*w,*y,*rex_temp;
	Vec_Com Q,R,RHS,S;
	Vec_Dou C;
	q=new Complex[2*ndg];
	x0=new Complex[2*ndg];
	double norm_0=sqrt(L2_Norm(bm,ndg));
	int iter=0;
	int ok=0;

	while(iter < maxit)
	{
		//Step 1
		for(j=0;j<2*ndg;j++)
			x0[j]=re[j];
		r=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,x0,CX,DX,fbw,lbw,dbw,ndg,ncg);
		//multiply(am,x0,r);
		for (j=0;j<ndg;j++)
		{
			r[j]=bm[j]-r[j];
			r[j+ndg]=-r[j+ndg];
		}
				
		r2=sqrt(L2_Norm(r,2*ndg));
		//cout<<r2;
		//Step 2
		for(j=0;j<2*ndg;j++)
		{
			q[j]=r[j]/r2;
			Q.push_back(q[j]);
		}
		rhs=r2;
		RHS.push_back(rhs);
		for(i=0;i<Gmres_K;i++)
		{
			//Step 3
			for(j=0;j<2*ndg;j++)
			{
				q[j]=Q[i*2*ndg+j];
			}
			w=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,q,CX,DX,fbw,lbw,dbw,ndg,ncg);
			//multiply(am,q,w);
			
			//Step 4
			for(L=0;L<=i;L++)//
			{
				temp=0.0;
				for(j=0;j<2*ndg;j++)
					temp+=conj(w[j])*Q[L*2*ndg+j];
				R.push_back(temp);
				for(j=0;j<2*ndg;j++)
				{
					w[j]-=temp*Q[L*2*ndg+j];
					//cout<<w[j]<<endl;
				}
			}
			//Step 5
			R_temp=sqrt(L2_Norm(w,2*ndg));
			//cout<<R_temp;
			//R.push_back(temp);
			for(j=0;j<2*ndg;j++)
			{
				q[j]=w[j]/R_temp;
				Q.push_back(q[j]);
			}
			//Step 6
			for(L=0;L<=i-1;L++)
			{
				temp=C[L]*R[((i+1)*i)/2+L]+conj(S[L])*R[((i+1)*i)/2+L+1];
				R[((i+1)*i)/2+L+1]=S[L]*R[((i+1)*i)/2+L]-C[L]*R[((i+1)*i)/2+L+1];
				R[((i+1)*i)/2+L]=temp;
			}
			//Step 7
			temp=R[((i+1)*i)/2+i];
			if(sqrt(temp.real()*temp.real()+temp.imag()*temp.imag()) < 1.0e-16)
			{
				c=0.0;
				C.push_back(c);
				s=1.0;
				S.push_back(s);
			}
			else
			{
				temp=R_temp/R[((i+1)*i)/2+i];
				c=1.0/sqrt(1.0+(temp.real()*temp.real()+temp.imag()*temp.imag()));
				C.push_back(c);
				s=C[i]*temp;
				S.push_back(s);
			}
			
			//Step 8
			temp=C[i]*RHS[i];
			rhs=S[i]*RHS[i];
			RHS.push_back(rhs);
			RHS[i]=temp;
			R[((i+1)*i)/2+i]=C[i]*R[((i+1)*i)/2+i]+conj(S[i])*R_temp;
			R_temp=0.0;
			
			//Step 9
			
			
			y=Resolve_Tri_Matrix_GMRES(R,RHS,i+1);
						
			//Step 10
			
			rex_temp=multiply_mn_GMRES(Q,y,2*ndg,i+1);
			
			for(j=0;j<2*ndg;j++)
				re[j]=x0[j]+rex_temp[j];
			//for (j=0;j<m;j++)
			//{
			//	cout<<rex[j]<<"  "<<i<<endl;
			//}
			//Step 11
			rex_temp=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,re,CX,DX,fbw,lbw,dbw,ndg,ncg);
			//multiply(am,rex,rex_temp);
			
			for (j=0;j<ndg;j++)
			{
				r[j]=bm[j]-rex_temp[j];
				r[j+ndg]=-rex_temp[j+ndg];
			}
			r2=sqrt(L2_Norm(r,2*ndg))/norm_0;
			cout<<"GMRES"<<" Steps = "<<iter*Gmres_K+i+1<<", iter = "<<iter<<",  res = "<<r2<<endl;
			out_resual<<iter*Gmres_K+i+1<<"   "<<r2<<endl;
			if(r2 <= tol)
			{
				ok=1;
				break;
			}
			//Step 12
		}
		Q.clear();
		R.clear();
		RHS.clear();
		C.clear();
		S.clear();
		if(ok == 1)
			break;
		iter++;
	}
	if(iter >= maxit)
		iter-=1;
	cout<<"Total GMRES("<<Gmres_K<<") Steps = "<<iter*Gmres_K+i+1<<", iter = "<<iter<<", to achieve "<<r2<<endl;

	delete [] x0;
	delete [] r;
	delete [] q;
	delete [] w;
	delete [] y;
	delete [] rex_temp;

	return re;
}
Complex *resolve_matrix::Resolve_Tri_Matrix_GMRES(Vec_Com &am,Vec_Com &bm,int m)
{
	int i,j;
	Complex *result=new Complex[m];
	result[m-1]=bm[m-1]/am[((m-1+1)*(m-1))/2+m-1];
	for(i=m-2;i>=0;i--)
	{
		for(j=i+1;j<m;j++)
			result[i]=result[i]+am[((j+1)*j)/2+i]*result[j];
		result[i]=(bm[i]-result[i])/am[((i+1)*i)/2+i];
	}
	return result;
}
Complex *resolve_matrix::multiply_mn_GMRES(Vec_Com &am,Complex *bm,int m,int n)
{
	int i,j;
	Complex *result=new Complex[m];
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			result[i]+=am[j*m+i]*bm[j];
		}
	}
	return result;
}

void resolve_matrix::CGLS(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
			Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
			Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
			Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,
			int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	double beta,te,alpha;
	Complex *r=new Complex[2*ndg];
	Complex *p=new Complex[2*ndg];
	Complex *s,*q;
	Vec_Com Resual;
	Resual.assign(2*ndg);

	//
	for(i=0;i<ndg;i++)
		r[i]=bm[i];
	
    s=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
	for(i=0;i<2*ndg;i++)
		p[i]=s[i];
	double norm_0=sqrt(L2_Norm(bm,ndg));
    double eps=1e-2;
	int Steps=0;
	while( sqrt(L2_Norm(r,2*ndg))/norm_0 > eps && Steps < 2*ndg )
	{
		out_resual<<Steps<<"   "<<sqrt(L2_Norm(r,2*ndg)) / norm_0<<endl;
		Steps++;
		cout<<"Steps="<<Steps<<" ,R_Norm="<<sqrt(L2_Norm(r,2*ndg))/norm_0<<endl;
		
		te=L2_Norm(s,2*ndg);
		q=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
		alpha=te/L2_Norm(q,2*ndg);
		for(i=0;i<2*ndg;i++)
		{
			rex[i]+=alpha*p[i];
			r[i]-=alpha*q[i];
		}
		s=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
		beta=L2_Norm(s,2*ndg)/te;
		for(i=0;i<2*ndg;i++)
			p[i]=s[i]+beta*p[i];
				
	}
	cout<<"Steps="<<Steps<<" ,R_Norm="<<sqrt(L2_Norm(r,2*ndg))/norm_0<<endl;
}
void resolve_matrix::LSQR(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
			Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
			Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
			Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{
	int i;
	rex=new Complex[2*ndg];
	
	double rss;
	double alpha,beta,pha,rou,phap,roup,Ci,Si,Cita;
	Complex *Ui=new Complex[2*ndg];
	Complex *Vi=new Complex[2*ndg];
	Complex *temp;
	Vec_Com r,U,V,W;
	r.assign(2*ndg);
	U.assign(2*ndg);
	V.assign(2*ndg);
	W.assign(2*ndg);
		
	//
	beta=sqrt(L2_Norm(bm,ndg));
	for(i=0;i<ndg;i++)
		Ui[i]=bm[i]/beta;
	temp=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,Ui,CX,DX,fbw,lbw,dbw,ndg,ncg);
	for(i=0;i<2*ndg;i++)
		V[i]=temp[i];
	//multiply(amt,Ui,V);
	alpha=sqrt(L2_Norm(V));
    for(i=0;i<2*ndg;i++)
	{
		Vi[i]=V[i]/alpha;
		W[i]=Vi[i];
	}
    phap=beta;
	roup=alpha;


	int iterNum=0;
	double Norm_b=sqrt(L2_Norm(bm,ndg));
	double eps=1.0e-2;
	//cout<<Norm_b<<endl;
	while ( true )
	{
		temp=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,rex,CX,DX,fbw,lbw,dbw,ndg,ncg);
		//multiply(am,x,temp);
		for(i=0;i<ndg;i++)
		{
			r[i]=bm[i]-temp[i];
			r[i+ndg]=-temp[i+ndg];
		}
			
		rss=sqrt(L2_Norm(r))/Norm_b;//余量
		iterNum++;
		cout<<"iterNum= "<<iterNum<<"  , rss= "<<rss<<endl;
		
		if(rss < eps || iterNum > 2*ndg)
		{
			break;     //终止整个循环
		}

		temp=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,Vi,CX,DX,fbw,lbw,dbw,ndg,ncg);
		//multiply(am,Vi,temp);
		for(i=0;i<2*ndg;i++)
			U[i]=temp[i]-alpha*Ui[i];
		beta=sqrt(L2_Norm(U));
		for(i=0;i<2*ndg;i++)
			Ui[i]=U[i]/beta;
		temp=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,Ui,CX,DX,fbw,lbw,dbw,ndg,ncg);
		//multiply(amt,Ui,temp);
		for(i=0;i<2*ndg;i++)
			V[i]=temp[i]-beta*Vi[i];
		alpha=sqrt(L2_Norm(V));
		for(i=0;i<2*ndg;i++)
			Vi[i]=V[i]/alpha;
		rou=sqrt(roup*roup+beta*beta);
		Ci=roup/rou;
		Si=beta/rou;
		Cita=Si*alpha;
		roup=-Ci*alpha;
		pha=Ci*phap;
		phap=Si*phap;
		for(i=0;i<2*ndg;i++)
		{
			rex[i]+=(pha/rou)*W[i];
			W[i]=Vi[i]-(Cita/rou)*W[i];
		}
	}
}

void resolve_matrix::resolve_bm(Vec_Com &zm,Vec_Com &vm,int n,int bw)
{//b为左端带状矩阵，d为右端列向量，n为矩阵维数，bw为带状矩阵半带宽
	int i,j,k;
	Complex T;
	int mbw=2*bw+1;
	rex=new Complex[n];
	Complex *zmm=new Complex[n*mbw];
	//将Z矩阵和V矩阵赋到另外的矩阵中，保证不改变其值
	for(i=0;i<n;i++)
	{
		for(j=0;j<mbw;j++)
		{
			zmm[i*mbw+j]=zm[i*mbw+j];
		}
		rex[i]=vm[i];
	}
	//消元
	for(i=0;i<n-1;i++)
	{
		if(i < bw)
		{
			for(j=i+1;j<i+bw+1 && j<n;j++)
			{
				if(j <= bw)
				{
					T=zmm[j*mbw+i]/zmm[i*mbw+i];
					for(k=i+1;k<i+bw+1;k++)
						zmm[j*mbw+k]=zmm[j*mbw+k]-zmm[i*mbw+k]*T;
				}
				else
				{
					T=zmm[j*mbw+i-(j-bw)]/zmm[i*mbw+i];
					for(k=i+1;k<i+bw+1;k++)
						zmm[j*mbw+k-(j-bw)]=zmm[j*mbw+k-(j-bw)]-zmm[i*mbw+k]*T;
				}
				rex[j]=rex[j]-rex[i]*T;
			}
		}
		else
		{
			for(j=i+1;j<i+bw+1 && j<n;j++)
			{
				T=zmm[j*mbw+i-(j-bw)]/zmm[i*mbw+bw];
				for(k=i+1;k<i+bw+1;k++)
					zmm[j*mbw+k-(j-bw)]=zmm[j*mbw+k-(j-bw)]-zmm[i*mbw+k-(i-bw)]*T;
				rex[j]=rex[j]-rex[i]*T;
			}
		}
	}
	//回代过程
	rex[n-1]=rex[n-1]/zmm[(n-1)*mbw+bw];
	for(i=n-2;i>=0;i--)
	{
		T=0.0;
		if(i>bw)
		{
			for(j=i+1;j<i+bw+1 && j<n;j++)
				T=T+zmm[i*mbw+j-(i-bw)]*rex[j];
			rex[i]=(rex[i]-T)/zmm[i*mbw+bw];
		}
		else
		{                            
			for(j=i+1;j<i+bw+1 && j<n;j++)
				T=T+zmm[i*mbw+j]*rex[j];
			rex[i]=(rex[i]-T)/zmm[i*mbw+i];
		}
	}
	delete [] zmm;
}
void resolve_matrix::resolve_m(Vec_Com &a,Vec_Com &b,Vec_Com &re)
{
	int l,k,i,j,is,u,v;
    double p,d;
	Complex t0;
	int m = b.size();
	re.assign(m);
	int *js=new int[m];
    for (k=0;k<=m-2;k++)
	{
		d=0.0;
		for (i=k;i<=m-1;i++)
        for (j=k;j<=m-1;j++)
        {
			u=i*m+j;
            p=norm(a[u]);
            if (p>d) {d=p;js[k]=j;is=i;}
        }
        if (d+1.0==1.0)
        {
			cout<<"err**fail\n";
        }
        if (is!=k)
        {
			for (j=k;j<=m-1;j++)
            {
				u=k*m+j; v=is*m+j;
                t0=a[u]; a[u]=a[v]; a[v]=t0;
            }
			t0=b[k]; b[k]=b[is]; b[is]=t0;
        }
        if (js[k]!=k)
			for (i=0;i<=m-1;i++)
            {
				u=i*m+k; v=i*m+js[k];
				t0=a[u]; a[u]=a[v]; a[v]=t0;
            }
        v=k*m+k;
        for (j=k+1;j<=m-1;j++)
        {
			u=k*m+j;
			a[u]=a[u]/a[v];
		}
		b[k]=b[k]/a[v];
        for (i=k+1;i<=m-1;i++)
        {
			u=i*m+k;
            for (j=k+1;j<=m-1;j++)
            {
				v=k*m+j; l=i*m+j;
				a[l]=a[l]-a[u]*a[v];
            }
			b[i]=b[i]-a[u]*b[k];
        }
	}
    u=(m-1)*m+m-1;
    d=norm(a[u]);
    if (d+1.0==1.0)
    {
		cout<<"err**fail\n";
    }
	b[m-1]=b[m-1]/a[u];
    for (i=m-2;i>=0;i--)
    for (j=i+1;j<=m-1;j++)
    {
		u=i*m+j;
		b[i]=b[i]-a[u]*b[j];
    }
    js[m-1]=m-1;
    for (k=m-1;k>=0;k--)
      if (js[k]!=k)
        {
		  t0=b[k];b[k]=b[js[k]];b[js[k]]=t0;
        }
	for(i=0;i<m;i++)
	{
		re[i]=b[i];
	}
}
void resolve_matrix::resolve_m(Vec_Com &a,Vec_Com &b,int m)
{
	int l,k,i,j,is,u,v;
    double p,d;
	Complex t0;
	rex=new Complex[m];
	int *js=new int[m];
    for (k=0;k<=m-2;k++)
	{
		d=0.0;
		for (i=k;i<=m-1;i++)
        for (j=k;j<=m-1;j++)
        {
			u=i*m+j;
            p=norm(a[u]);
            if (p>d) {d=p;js[k]=j;is=i;}
        }
        if (d+1.0==1.0)
        {
			cout<<"err**fail\n";
            //return(0);
        }
        if (is!=k)
        {
			for (j=k;j<=m-1;j++)
            {
				u=k*m+j; v=is*m+j;
                t0=a[u]; a[u]=a[v]; a[v]=t0;
                //p=a[u].Imag; a[u].Imag=a[v].Imag; a[v].Imag=p;
            }
			t0=b[k]; b[k]=b[is]; b[is]=t0;
            //p=b[k].Real; b[k].Real=b[is].Real; b[is].Real=p;
            //p=b[k].Imag; b[k].Imag=b[is].Imag; b[is].Imag=p;
        }
        if (js[k]!=k)
			for (i=0;i<=m-1;i++)
            {
				u=i*m+k; v=i*m+js[k];
				t0=a[u]; a[u]=a[v]; a[v]=t0;
                //p=a[u].Real; a[u].Real=a[v].Real; a[v].Real=p;
                //p=a[u].Imag; a[u].Imag=a[v].Imag; a[v].Imag=p;
            }
        v=k*m+k;
        for (j=k+1;j<=m-1;j++)
        {
			u=k*m+j;
            //p=a[u].Real*a[v].Real; q=-a[u].Imag*a[v].Imag;
            //s=(a[v].Real-a[v].Imag)*(a[u].Real+a[u].Imag);
            //a[u].Real=(p-q)/d; a[u].Imag=(s-p-q)/d;
			a[u]=a[u]/a[v];
		}
        //p=b[k].Real*a[v].Real; q=-b[k].Imag*a[v].Imag;
        //s=(a[v].Real-a[v].Imag)*(b[k].Real+b[k].Imag);
        //b[k].Real=(p-q)/d; b[k].Imag=(s-p-q)/d;
		b[k]=b[k]/a[v];
        for (i=k+1;i<=m-1;i++)
        {
			u=i*m+k;
            for (j=k+1;j<=m-1;j++)
            {
				v=k*m+j; l=i*m+j;
                //p=a[u].Real*a[v].Real; q=a[u].Imag*a[v].Imag;
                //s=(a[u].Real+a[u].Imag)*(a[v].Real+a[v].Imag);
                //a[l].Real=a[l].Real-p+q;
                //a[l].Imag=a[l].Imag-s+p+q;
				a[l]=a[l]-a[u]*a[v];
            }
            //p=a[u].Real*b[k].Real; q=a[u].Imag*b[k].Imag;
            //s=(a[u].Real+a[u].Imag)*(b[k].Real+b[k].Imag);
            //b[i].Real=b[i].Real-p+q; b[i].Imag=b[i].Imag-s+p+q;
			b[i]=b[i]-a[u]*b[k];
        }
	}
    u=(m-1)*m+m-1;
    d=norm(a[u]);
    if (d+1.0==1.0)
    {
		cout<<"err**fail\n";
        //return(0);
    }
    //p=a[u].Real*b[m-1].Real; q=-a[u].Imag*b[m-1].Imag;
    //s=(a[u].Real-a[u].Imag)*(b[m-1].Real+b[m-1].Imag);
    //b[m-1].Real=(p-q)/d; b[m-1].Imag=(s-p-q)/d;
	b[m-1]=b[m-1]/a[u];
    for (i=m-2;i>=0;i--)
    for (j=i+1;j<=m-1;j++)
    {
		u=i*m+j;
        //p=a[u].Real*b[j].Real; q=a[u].Imag*b[j].Imag;
        //s=(a[u].Real+a[u].Imag)*(b[j].Real+b[j].Imag);
        //b[i].Real=b[i].Real-p+q;
        //b[i].Imag=b[i].Imag-s+p+q;
		b[i]=b[i]-a[u]*b[j];
    }
    js[m-1]=m-1;
    for (k=m-1;k>=0;k--)
      if (js[k]!=k)
        {
		  t0=b[k];b[k]=b[js[k]];b[js[k]]=t0;
		  //p=b[k].Real; b[k].Real=b[js[k]].Real; b[js[k]].Real=p;
          //p=b[k].Imag; b[k].Imag=b[js[k]].Imag; b[js[k]].Imag=p;
        }
	for(i=0;i<m;i++)
	{
		rex[i]=b[i];
		//rex[i].Imag=b[i].Imag;
		//rex[i].display();
	}
    //return(1);
}
void resolve_matrix::resolve_cgm(Complex *am,Complex *bm,int n,int bw)
{
	int i;
	double rfa,beta;
	//int mbw=2*bw+1;
	Complex *r=new Complex[n];//余量
	Complex *d=new Complex[n];
	Complex *amt;
	rex=new Complex[n];
	//矩阵的共轭转置
	amt=Transpose_BandMatrix(am,n,bw);
	
	//第0步
	Complex *t0;
	t0=BandMatrix_multiply(amt,bm,n,bw);
	for(i=0;i<n;i++)
	{
		//rex[i]=0.0;//0步
		r[i]=bm[i];//0步
		d[i]=t0[i];//1步
	}
	double sum_rr,sum_r,sum_m;
	sum_rr=L2_Norm(t0,n);
	//第一步到第n步
	Complex *t1,*t2;
	int t=0;
	while(sqrt(L2_Norm(r,n)) > eps)
	{
		t++;
		cout<<sqrt(L2_Norm(r,n))<<endl;
	
		t1=BandMatrix_multiply(am,d,n,bw);
		sum_m=L2_Norm(t1,n);
		rfa=sum_rr/sum_m;
		
	    for(i=0;i<n;i++)
		{
			rex[i]=rex[i]+rfa*d[i];
			r[i]=r[i]-rfa*t1[i];
		}
		t2=BandMatrix_multiply(amt,r,n,bw);
        sum_r=L2_Norm(t2,n);
		beta=sum_r/sum_rr;
		
		for(i=0;i<n;i++)
		{
			d[i]=t2[i]+beta*d[i];//2
		}
		sum_rr=sum_r;
	}
	delete [] r;
	delete [] d;
	delete [] amt;
	delete [] t0;
	delete [] t1;
	delete [] t2;
}


void resolve_matrix::resolve_pbtg_cgm(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
									  Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
									  Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
									  Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{//bm应该是ndg维
	rex=new Complex[2*ndg];//方程的解
	int i;
	Complex *t1,*t2;
	int n1=ndg/ncg;
	double rfa,beta;
	
	Complex *r=new Complex[2*ndg];//余量
	Complex *d;
	Vec_Com Resual;
	Resual.assign(2*ndg);

	//第0步
	//d=BandMatrix_multiply(amt,bm,n,b);
	Complex *Bm=new Complex[2*ndg];
	for(i=0;i<ndg;i++)
		Bm[i]=bm[i];
	d=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,Bm,CX,DX,fbw,lbw,dbw,ndg,ncg);
	
	for(i=0;i<ndg;i++)
	{
		//0步
		r[i]=bm[i];//0步下半部分为零
	}
	double sum_rr,sum_r,sum_m;
	sum_rr=L2_Norm(d,2*ndg);
	double norm_0=sqrt(L2_Norm(bm,ndg));
	//double Con_N=L2_Norm(bm,ndg);
	//cout<<sum_rr<<endl;
	//第一步到第n步
	int step=0;
	while(sqrt(L2_Norm(r,2*ndg))/norm_0 > eps && step<2*ndg)
	{
		out_resual<<step<<"   "<<sqrt(L2_Norm(r,2*ndg)) / norm_0<<endl;
		step++;
		//out_test<<L2_Norm(r,2*ndg)<<endl;
		cout<<sqrt(L2_Norm(r,2*ndg))/norm_0<<"     "<<step<<endl;
		

		//t1=multiply(am,d,n,bw);
		t1=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,d,CX,DX,fbw,lbw,dbw,ndg,ncg);
	    
		sum_m=L2_Norm(t1,2*ndg);
		rfa=sum_rr/sum_m;
		//cout<<rfa<<endl;

		for(i=0;i<2*ndg;i++)
		{
			rex[i]=rex[i]+rfa*d[i];
			r[i]=r[i]-rfa*t1[i];
			//cout<<r[i]<<endl;
		}
		
		
		//t2=BandMatrix_multiply(amt,r,n,bw);
		t2=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
         
		
        sum_r=L2_Norm(t2,2*ndg);
		beta=sum_r/sum_rr;
		
		for(i=0;i<2*ndg;i++)
		{
			d[i]=t2[i]+beta*d[i];//2
			//cout<<rex[i]<<"        "<<i<<endl;
		}
		sum_rr=sum_r;
	}
	delete [] r;
	delete [] Bm;
	delete [] d;
	delete [] t1;
	delete [] t2;

}
void resolve_matrix::CG(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
						Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
						Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
						Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{  
	rex=new Complex[2*ndg];
	int i;
	double ebs;
	double excitationAbs;
	double alpha, beta;
	
	Complex *p=new Complex[2*ndg];
	Complex *r=new Complex[2*ndg];

	Complex *t0,*t1,*t2;
	
	int step=0;
	
	excitationAbs=sqrt(L2_Norm(bm,ndg));   //平方和开根号
	for(i=0;i<ndg;i++)
		r[i]=bm[i];
	Complex *Bm=new Complex[2*ndg];
	for(i=0;i<ndg;i++)
		Bm[i]=bm[i];
	t0=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,Bm,CX,DX,fbw,lbw,dbw,ndg,ncg);
	beta=1.0/L2_Norm(t0,2*ndg);   //平方和
	for(i=0;i<2*ndg;i++)
		p[i]=beta*t0[i];
	
	while ( true )
	{
		out_resual<<step<<"   "<<sqrt(L2_Norm(r,2*ndg)) / excitationAbs<<endl;
		step++;
		t1=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
		alpha=1.0/L2_Norm(t1,2*ndg);
		for(i=0;i<2*ndg;i++)
		{
			rex[i]+=alpha*p[i];
			r[i]=r[i]-alpha*t1[i];
		}
		ebs=sqrt(L2_Norm(r,2*ndg))/excitationAbs;

		cout << "steps= " <<step << ", ebs= " <<ebs << endl;
		
		if(ebs < eps || step > 2*ndg)//停止准则:0.01即可得到精确解
		{
			break;     //终止整个循环
		}
		
		t2=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
		beta=1.0/L2_Norm(t2,2*ndg);
		for(i=0;i<2*ndg;i++)
			p[i]+=beta*t2[i];
	}
	delete [] p;
	delete [] r;
	delete [] Bm;
	delete [] t0;
	delete [] t1;
	delete [] t2;
	
}
Complex *resolve_matrix::Multi_Tran_M_V(Complex *asmnt,Complex *bsmnt,Complex *a1mnt,Complex *b1mnt,Complex *aabmnt,Complex *babmnt,Complex *VR,Vec_Dou &CX,Vec_Dou &DX,int fbw,int lbw,int dbw,int ndg,int ncg)
{
	Complex *re=new Complex[2*ndg];
	Complex *d1=new Complex[ndg];
	Complex *d2=new Complex[ndg];
	Complex *d3,*t0,*t1,*t2,*t3;
	int i;
	int n1=ndg/ncg;
	for(i=0;i<ndg;i++)
	{
		d1[i]=VR[i];
		d2[i]=VR[i+ndg];
	}
	d3=Average(d1,n1,ncg);
	//上
	t0=BandMatrix_multiply(asmnt,d1,ndg,fbw);//asmnt与r的 上 半部分d1相乘，维数为ndg
	t1=Matrix_multiply(aabmnt,d3,fbw,dbw,ncg,n1);//aabmnt与r的 上 半部分d3相乘，维数为ncg
	//t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//插值到ndg维，上
	t1=Lg3_Interpolate(t1,CX,DX,ncg,ndg);
	t2=BandMatrix_multiply(a1mnt,d2,ndg,lbw);//a1mnt与r的 下 半部分d2相乘，维数为ndg
	t0=Plus_Matrix(t0,t1,ndg);//上
	t0=Plus_Matrix(t0,t2,ndg);//上

    //下
	t1=BandMatrix_multiply(bsmnt,d1,ndg,fbw);//bsmnt与r的 上 半部分d1相乘，维数为ndg 下 
	t2=Matrix_multiply(babmnt,d3,fbw,dbw,ncg,n1);//babmnt与r的 上 半部分d3相乘，维数为ncg 下 
	//t2=Linear_Interpolate(t2,CX,DX,ncg,ndg);//插值到ndg维， 下
	t2=Lg3_Interpolate(t2,CX,DX,ncg,ndg);
	t3=BandMatrix_multiply(b1mnt,d2,ndg,lbw);//b1mnt与r的 下 半部分d2相乘，维数为ndg 下 
	t1=Plus_Matrix(t1,t2,ndg);//相加， 下 
	t1=Plus_Matrix(t1,t3,ndg);//相加， 下 
		

	re=Unite_Matrix(t0,t1,ndg);//上下合并

	delete [] d1;
	delete [] d2;
	delete [] d3;
	delete [] t0;
	delete [] t1;
	delete [] t2;
	delete [] t3;
	
	return re;
}
Complex *resolve_matrix::Multi_M_V(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *aabmn,Complex *babmn,Complex *VR,Vec_Dou &CX,Vec_Dou &DX,int fbw,int lbw,int dbw,int ndg,int ncg)
{
	Complex *re;
	Complex *d1=new Complex[ndg];
	Complex *d2=new Complex[ndg];
	Complex *d3,*d4,*t0,*t1,*t2,*t3;
	int i;
	int n1=ndg/ncg;
	for(i=0;i<ndg;i++)
	{
		d1[i]=VR[i];
		d2[i]=VR[i+ndg];
	}
	d3=Average(d1,n1,ncg);
	d4=Average(d2,n1,ncg);
	//上
	t0=BandMatrix_multiply(asmn,d1,ndg,fbw);//asmn与d的 上 半部分d1相乘，维数为ndg 上 
	t1=Matrix_multiply(aabmn,d3,fbw,dbw,ncg,n1);//aabmn与d的 上 半部分d3相乘，维数为ncg 上 
	//t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//插值到ndg维，上
	t1=Lg3_Interpolate(t1,CX,DX,ncg,ndg);
	t2=BandMatrix_multiply(bsmn,d2,ndg,fbw);//bsmn与d的 下 半部分d2相乘，维数为ndg 上 
	t3=Matrix_multiply(babmn,d4,fbw,dbw,ncg,n1);//babmn与d的 下 半部分d4相乘，维数为ncg 上 
	//t3=Linear_Interpolate(t3,CX,DX,ncg,ndg);//插值到ndg维， 上
	t3=Lg3_Interpolate(t3,CX,DX,ncg,ndg);
	t0=Plus_Matrix(t0,t1,ndg);//上
	t0=Plus_Matrix(t0,t2,ndg);//上
	t0=Plus_Matrix(t0,t3,ndg);//上
	
    //下
	t1=BandMatrix_multiply(a1mn,d1,ndg,lbw);//a1mn与d的 上 半部分d1相乘，维数为ndg下
	t2=BandMatrix_multiply(b1mn,d2,ndg,lbw);//b1mn与d的 下 半部分d2相乘，维数为ndg下
	t3=Plus_Matrix(t1,t2,ndg);//下

	re=Unite_Matrix(t0,t3,ndg);//上下合并

	
	delete [] d1;
	delete [] d2;
	delete [] d3;
	delete [] d4;
	delete [] t0;
	delete [] t1;
	delete [] t2;
	delete [] t3;
	

	return re;
}
Complex *resolve_matrix::Multi_A_R(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *VR,int ndg,int bw)
{
	Complex *re=new Complex[2*ndg];
	Complex *d1=new Complex[ndg];
	Complex *d2=new Complex[ndg];
	Complex *t1,*t2,*t3,*t4;
	int i;
	for(i=0;i<ndg;i++)
	{
		d1[i]=VR[i];
		d2[i]=VR[i+ndg];
	}
	/////////////////////////////////////
	//上
	t1=BandMatrix_multiply(asmn,d1,ndg,bw);//asmn与d的 上 半部分d1相乘，维数为ndg 上 
	t2=BandMatrix_multiply(bsmn,d2,ndg,bw);//bsmn与d的 下 半部分d2相乘，维数为ndg 上 
	t1=Plus_Matrix(t1,t2,ndg);//上
	//下
	t3=BandMatrix_multiply(a1mn,d1,ndg,bw);//a1mn与d的 上 半部分d1相乘，维数为ndg下
	t4=BandMatrix_multiply(b1mn,d2,ndg,bw);//b1mn与d的 下 半部分d2相乘，维数为ndg下
	t3=Plus_Matrix(t3,t4,ndg);//下
	re=Unite_Matrix(t1,t3,ndg);//上下合并

	delete [] d1;
	delete [] d2;
	delete [] t1;
	delete [] t2;
	delete [] t3;
	delete [] t4;

	return re;
}
void resolve_matrix::resolve_pbtg_HJG_cgm(Complex *a1mn,Complex *b1mn,Complex *asmn,Complex *bsmn,
										  Complex *a1mnt,Complex *b1mnt,Complex *asmnt,Complex *bsmnt,
										  Complex *aabmn,Complex *babmn,Complex *aabmnt,Complex *babmnt,
										  Vec_Dou &CX,Vec_Dou &DX,Complex *bm,int ndg,int ncg,int lbw,int fbw,int dbw)
{//方程组右端列向量bm应该是2*ndg维
	rex=new Complex[2*ndg];//方程的解
	int i;
	int n1=ndg/ncg;
	double rfa,beta;
			
	Complex *r=new Complex[2*ndg];//余量
	Complex *p=new Complex[2*ndg];

	
	//第0步
	Complex *t1,*t2;
	/*Complex *d1=new Complex[ndg];
	Complex *d2=new Complex[ndg];
	Complex *d3,*d4;
	
	/////////////////////////////////////
	d3=Average(bm,n1,ncg);
	
	/////////////////////////////////////
	//t0=BandMatrix_multiply(amt,bm,n,bw);
	//上
	t0=BandMatrix_multiply(asmnt,bm,ndg,fbw);//asmnt乘以bm，上
	t1=Matrix_multiply(aabmnt,d3,fbw,dbw,ncg,n1);//aabmnt与bm的上半部分d3相乘，维数为ncg，上
	t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//将所得的维数为ncg的列向量线形插值到维数为ndg，上
	d1=Plus_Matrix(t0,t1,ndg);//所得的上半部分的两个矩阵相加，上
	

	//下
	t0=BandMatrix_multiply(bsmnt,bm,ndg,fbw);//bsmnt乘以bm，下
	t1=Matrix_multiply(babmnt,d3,fbw,dbw,ncg,n1);//babmnt与bm的上半部分d3相乘，维数为ncg，下
	t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//将所得的维数为ncg的列向量线形插值到维数为ndg，下
	d2=Plus_Matrix(t0,t1,ndg);//列向量相加，下
	
	p=Unite_Matrix(d1,d2,ndg);//将两个ndg维的列向量合并成一个2*ndg维的列向量//1步下半部分为零*/
	p=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,bm,CX,DX,fbw,lbw,dbw,ndg,ncg);
	
	for(i=0;i<ndg;i++)
	{
		r[i]=bm[i];//0步下半部分为零
	}
	double sum_rr,sum_r,sum_m;
	sum_rr=L2_Norm(r,2*ndg);
	//double Con_N=L2_Norm(bm,ndg);
	//cout<<Con_N<<endl;
	//第一步到第n步
	int t=0;
	while(sqrt(sum_rr) > eps && t<2*ndg)
	{
		t++;
		cout<<sum_rr<<endl;

		sum_m=L2_Norm(p,2*ndg);
		rfa=sum_rr/sum_m;
		//t1=multiply(am,p,n,bw);
		/*for(i=0;i<ndg;i++)
		{
			d1[i]=p[i];
			d2[i]=p[i+ndg];
			
		}
		/////////////////////////////////////
		d3=Average(d1,n1,ncg);
		d4=Average(d2,n1,ncg);
		//for(i=0;i<ncg;i++)
		//	outresult<<d3[i]<<"          "<<d4[i]<<"        "<<i<<endl;
		/////////////////////////////////////
		//上
		t0=BandMatrix_multiply(asmn,d1,ndg,fbw);//asmn与d的 上 半部分d1相乘，维数为ndg 上 
		t1=Matrix_multiply(aabmn,d3,fbw,dbw,ncg,n1);//aabmn与d的 上 半部分d3相乘，维数为ncg 上 
		t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//插值到ndg维，上
		t2=BandMatrix_multiply(bsmn,d2,ndg,fbw);//bsmn与d的 下 半部分d2相乘，维数为ndg 上 
		t3=Matrix_multiply(babmn,d4,fbw,dbw,ncg,n1);//babmn与d的 下 半部分d4相乘，维数为ncg 上 
		t3=Linear_Interpolate(t3,CX,DX,ncg,ndg);//插值到ndg维， 上 
		t0=Plus_Matrix(t0,t1,ndg);//上
		t0=Plus_Matrix(t0,t2,ndg);//上
		t0=Plus_Matrix(t0,t3,ndg);//上
		

    	//下
		t1=BandMatrix_multiply(a1mn,d1,ndg,lbw);//a1mn与d的 上 半部分d1相乘，维数为ndg下
		t2=BandMatrix_multiply(b1mn,d2,ndg,lbw);//b1mn与d的 下 半部分d2相乘，维数为ndg下
		t3=Plus_Matrix(t1,t2,ndg);//下

		t1=Unite_Matrix(t0,t3,ndg);//上下合并*/
		t1=Multi_M_V(asmn,bsmn,a1mn,b1mn,aabmn,babmn,p,CX,DX,fbw,lbw,dbw,ndg,ncg);
	
		
		//cout<<sum_rr<<endl;

		for(i=0;i<2*ndg;i++)
		{
			rex[i]=rex[i]+rfa*p[i];
			r[i]=r[i]-rfa*t1[i];
			//cout<<r[i]<<endl;
		}
		sum_r=L2_Norm(r,2*ndg);
		beta=sum_r/sum_rr;
		//t2=BandMatrix_multiply(amt,r,n,bw);
		/*for(i=0;i<ndg;i++)
		{
			d1[i]=r[i];
			d2[i]=r[i+ndg];
		}
		/////////////////////////////////////
		d3=Average(d1,n1,ncg);
		/////////////////////////////////////
		
		//上
		t0=BandMatrix_multiply(asmnt,d1,ndg,fbw);//asmnt与r的 上 半部分d1相乘，维数为ndg
		t1=Matrix_multiply(aabmnt,d3,fbw,dbw,ncg,n1);//aabmnt与r的 上 半部分d3相乘，维数为ncg
		t1=Linear_Interpolate(t1,CX,DX,ncg,ndg);//插值到ndg维，上
		t2=BandMatrix_multiply(a1mnt,d2,ndg,lbw);//a1mnt与r的 下 半部分d2相乘，维数为ndg
		t0=Plus_Matrix(t0,t1,ndg);//上
		t0=Plus_Matrix(t0,t2,ndg);//上

    	//下
		t1=BandMatrix_multiply(bsmnt,d1,ndg,fbw);//bsmnt与r的 上 半部分d1相乘，维数为ndg 下 
		t2=Matrix_multiply(babmnt,d3,fbw,dbw,ncg,n1);//babmnt与r的 上 半部分d3相乘，维数为ncg 下 
	    t2=Linear_Interpolate(t2,CX,DX,ncg,ndg);//插值到ndg维， 下 
		t3=BandMatrix_multiply(b1mnt,d2,ndg,lbw);//b1mnt与r的 下 半部分d2相乘，维数为ndg 下 
		t1=Plus_Matrix(t1,t2,ndg);//相加， 下 
		t1=Plus_Matrix(t1,t3,ndg);//相加， 下 
		

		t2=Unite_Matrix(t0,t1,ndg);//上下合并*/
		t2=Multi_Tran_M_V(asmnt,bsmnt,a1mnt,b1mnt,aabmnt,babmnt,r,CX,DX,fbw,lbw,dbw,ndg,ncg);
	        
		for(i=0;i<2*ndg;i++)
		{
			p[i]=t2[i]+beta*p[i];//2
			//cout<<rex[i]<<"        "<<i<<endl;
		}
		sum_rr=sum_r;
	}
	delete [] r;
	delete [] p;
	delete [] t1;
	delete [] t2;

}
Complex *resolve_matrix::Transpose_BandMatrix(Complex *am,int n,int bw)
{
	//带状矩阵的共轭转置
	int i,j,u,v;
	int mbw=2*bw+1;
	Complex *amt=new Complex[n*mbw];
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<mbw;j++)
		{
			u=i*mbw+j;
			if(i <= bw)
			{
				if(j <= bw)
				{
					v=j*mbw+i;
					amt[u]=conj(am[v]);//
				}
				else
				{
					v=j*mbw+i-(j-bw);
					amt[u]=conj(am[v]);
				}
			}
			else
			{
				if(j+i-bw <= bw)
				{
					v=(j+i-bw)*mbw+i;
		    		amt[u]=conj(am[v]);
				}
				else if(j+i-bw > bw && j+i-bw < n)
				{
					v=(j+i-bw)*mbw+i-((j+i-bw)-bw);
		    		amt[u]=conj(am[v]);
					//cout<<u<<"             "<<v<<endl;
				}
			}
		}
	}
	return amt;
}

Complex *resolve_matrix::Multi_DiagM_BandM(Complex *DiagM,Complex *BandM,int m,int bw)
{
	int i,j,u;
	int mbw = 2 * bw + 1;
	Complex *re=new Complex[m*mbw];
	for(i=0; i<m; i++)
	{
		for(j=0; j<mbw; j++)
		{
			u = i * mbw + j;
			re[u] = BandM[u] * DiagM[i];
		}
	}
	return re;
}
Complex *resolve_matrix::Multi_DiagM_ABM(Complex *DiagM,Complex *ABM,int md,int mc,int rd_bw,int rf_bw)
{
	int i,j,u;
	int nn1 = md / mc;
	Complex *re=new Complex[mc*mc];
	for(i=0;i<mc;i++)
	{
		for(j=0;j<mc;j++)
		{
			u=i*mc+j;
			if(abs(i-j) < rd_bw/nn1 && abs(i-j) > rf_bw/nn1)
			{
				re[u] = ABM[u] * DiagM[i];
			}
		}
	}
	return re;
}
Complex *resolve_matrix::BandMatrix_multiply(Complex *a,Complex *b,int n,int bw)
{
	int i,l;
	int mbw=2*bw+1;
	Complex *re=new Complex[n];
	for(i=0;i<n;i++)
	{
		if(i <= bw)
		{
			re[i]=0.0;
		   	for(l=0;l<=i+bw;l++)
			{
				re[i]=re[i]+a[i*mbw+l]*b[l];
			}
		}
		else if(i > bw && i < n-bw)
		{
			re[i]=0.0;
		   	for(l=0;l<mbw;l++)
			{
				re[i]=re[i]+a[i*mbw+l]*b[l+i-bw];
			}
		}
		else
		{
			re[i]=0.0;
		   	for(l=0;l<n-(i-bw);l++)//n减去每行零的个数(i-bw)
			{
				re[i]=re[i]+a[i*mbw+l]*b[l+i-bw];
			}
		}
	}
	return re;
}
Complex *resolve_matrix::Matrix_multiply(Complex *a,Complex *b,int fbw,int dbw,int n,int nn1)
{
	int i,j;
	Complex *re=new Complex[n];
	for(i=0;i<n;i++)
	{
		//re[i]=0.0;
		for(j=0;j<n;j++)
		{
			if(abs(i-j) > fbw/nn1 && abs(i-j) < dbw/nn1)
			{
				re[i]=re[i]+a[i*n+j]*b[j];
			}
		}
	}
	return re;
}


Complex *resolve_matrix::Creat_Diag(Complex *asmn,Complex *a1mn,Complex *b1mn,int ndg,int lbw,int fbw)
{
	int j,k;
	Complex *re=new Complex[2*ndg];
	
	re[0]=asmn[0];
	for(j=1;j<ndg;j++)
	{
		if(j <= fbw)
		{	
			for(k=0;k<j;k++)
			{//注意：只有有k的元素是累加的
				re[j]=re[j]+pow(asmn[j*(2*fbw+1)+k],2)/re[k];
			}
			re[j]=asmn[j*(2*fbw+1)+j]-re[j];
		}
		else
		{
			for(k=j-fbw;k<j;k++)
			{
				re[j]=re[j]+pow(asmn[j*(2*fbw+1)+k-(j-fbw)],2)/re[k];
			}
			re[j]=asmn[j*(2*fbw+1)+fbw]-re[j];
		}
	}
	//下半对角阵
	for(k=0;k <= lbw;k++)
	{
		re[ndg]=re[ndg]+pow(a1mn[k],2)/re[k];
	}
	re[ndg]=b1mn[0]-re[ndg];
	for(j=ndg+1;j<2*ndg;j++)
	{
		if(j-ndg <= lbw)
		{	
			for(k=0;k <= lbw+j-ndg;k++)
			{
				re[j]=re[j]+pow(a1mn[(j-ndg)*(2*lbw+1)+k],2)/re[k];//
			}
			for(k=0;k < j-ndg;k++)
			{
				re[j]=re[j]+pow(b1mn[(j-ndg)*(2*lbw+1)+k],2)/re[k+ndg];//
			}
			re[j]=b1mn[(j-ndg)*(2*lbw+1)+(j-ndg)]-re[j];
		}
		else
		{
			if(j-ndg < ndg-lbw)
			{
				for(k=0;k < 2*lbw+1;k++)
				{
					re[j]=re[j]+pow(a1mn[(j-ndg)*(2*lbw+1)+k],2)/re[k+j-ndg-lbw];//
				}
			}
			else
			{
				for(k=j-ndg-lbw;k < ndg;k++)
				{
					re[j]=re[j]+pow(a1mn[(j-ndg)*(2*lbw+1)+k-(j-ndg-lbw)],2)/re[k];//
				}
			}
			for(k=j-ndg-lbw;k < j-ndg;k++)
			{
				re[j]=re[j]+pow(b1mn[(j-ndg)*(2*lbw+1)+k-(j-ndg-lbw)],2)/re[k+ndg];//
			}
			re[j]=b1mn[(j-ndg)*(2*lbw+1)+lbw]-re[j];
		}
	}
	return re;
}
Complex *resolve_matrix::Creat_triang(Complex *asmn,Complex *Dm,int ndg,int fbw)
{
	Complex *re=new Complex[ndg*(fbw+1)];
	int i,j,t;
	int mfbw=2*fbw+1;
	for(i=0;i<ndg;i++)
	{
		if(i <= fbw)
		{
			for(j=0;j <= i;j++)
			{
				t=i*(fbw+1)+j;
				if(i == j)
					re[t]=Dm[i];
				else
					re[t]=asmn[i*mfbw+j];
				//cout<<re[t]<<"    ";
			}
		}
		else
		{
			for(j=0;j <= fbw;j++)
			{
				t=i*(fbw+1)+j;
				if(j == fbw)
					re[t]=Dm[i];
				else
					re[t]=asmn[i*mfbw+j];
				//cout<<re[t]<<"    ";
			}
		}
		//cout<<endl;
	}
	return re;
}


Complex *resolve_matrix::Creat_Inv_Diag(Complex *Dm,int ndg)
{
	Complex *re=new Complex[2*ndg];
	for(int i=0;i<2*ndg;i++)
		re[i]=sqrt(1.0/Dm[i]);
	return re;
}
//下三角的带状阵与对角阵相乘
Complex *resolve_matrix::Multi_HalfBM_Diag(Complex *Half_BM,Complex *Dm_U_D,int ndg,int bw)
{
	Complex *re=new Complex[ndg*(bw+1)];
	int u;
	for(int i=0;i<ndg;i++)
	{
		for(int j=0;j<bw+1;j++)
		{
			u=i*(bw+1)+j;
			if(i<=bw)
				re[u]=Half_BM[u]*Dm_U_D[j];
			else
				re[u]=Half_BM[u]*Dm_U_D[j+i-bw];
			//cout<<re[u]<<"    ";
		}
		//cout<<endl;
	}
	return re;
}
//带状矩阵与对角阵相乘
Complex *resolve_matrix::Multi_BM_Diag(Complex *BM,Complex *Dm_U_D,int ndg,int bw)
{
	Complex *re=new Complex[ndg*(2*bw+1)];
	int u;
	for(int i=0;i<ndg;i++)
	{
		for(int j=0;j<2*bw+1;j++)
		{
			u=i*(2*bw+1)+j;
			if(i<=bw)
				re[u]=BM[u]*Dm_U_D[j];
			else
				re[u]=BM[u]*Dm_U_D[j+i-bw];
			//cout<<re[u]<<"    ";
		}
		//cout<<endl;
	}
	return re;
}
//下三角带状矩阵的求逆，生成满的下三角阵
Complex *resolve_matrix::Tria_Decom(Complex *TM,int ndg,int bw)
{
	int i,j,k;
	Complex *re=new Complex[((ndg+1)*ndg)/2];
	for(i=0;i<ndg;i++)
	{
		if(i<=bw)
			re[((i+1)*i)/2+i]=1.0/TM[i*(bw+1)+i];
		else
			re[((i+1)*i)/2+i]=1.0/TM[i*(bw+1)+bw];
	}
	//注意：k是行，i是列
	for(i=0;i<ndg;i++)
	{
		for(k=i+1;k<ndg;k++)
		{
			if(k<=bw)
			{
				for(j=i;j<k;j++)
				{
					re[(k*(k+1))/2+i]=re[(k*(k+1))/2+i]-TM[k*(bw+1)+j]*re[(j*(j+1))/2+i]/TM[k*(bw+1)+k];
				}
			}
			else
			{
				if(i<k-bw)
				{
					for(j=0;j<bw;j++)
					{
			    		re[(k*(k+1))/2+i]=re[(k*(k+1))/2+i]-TM[k*(bw+1)+j]*re[(j+(k-bw))*((j+(k-bw))+1)/2+i]/TM[k*(bw+1)+bw];
					}
				}
				else
				{//k是行
					for(j=i-(k-bw);j<bw;j++)//?
					{
			    		re[(k*(k+1))/2+i]=re[(k*(k+1))/2+i]-TM[k*(bw+1)+j]*re[((j+(k-bw))*((j+(k-bw))+1))/2+i]/TM[k*(bw+1)+bw];
					}
				}			
			}
		}
	}
	return re;
}
//带状矩阵与下三角阵相乘
Complex *resolve_matrix::Multi_BM_TM(Complex *BM,Complex *TM,int ndg,int bw)
{
	Complex *re=new Complex[ndg*ndg-((ndg-bw)*(ndg-bw-1))/2];
	int i,j,k,t;
	int mbw=2*bw+1;
	for(i=0;i<ndg;i++)
	{
		if(i<ndg-bw-1)
		{
			for(j=0;j<=i+bw;j++)
			{
				t=((2*bw+i+1)*i)/2+j;
				if(i<=bw)
				{
					for(k=j;k<=i+bw;k++)
					{
						re[t]=re[t]+BM[i*mbw+k]*TM[((k+1)*k)/2+j];
					}
				}
				else
				{
					if(i-bw>=j)
					{
			    		for(k=i-bw;k<=i+bw;k++)
						{
				    		re[t]=re[t]+BM[i*mbw+k-(i-bw)]*TM[((k+1)*k)/2+j];
						}
					}
			    	else
					{
			    		for(k=j;k<=i+bw;k++)
						{
				    		re[t]=re[t]+BM[i*mbw+k-(i-bw)]*TM[((k+1)*k)/2+j];
						}
					}
				}
			}
		}
		else
		{
			for(j=0;j<ndg;j++)
			{
				t=i*ndg+j-((ndg-bw)*(ndg-bw-1))/2;
				if(i-bw>=j)
				{
					for(k=i-bw;k<ndg;k++)
					{
						re[t]=re[t]+BM[i*mbw+k-(i-bw)]*TM[((k+1)*k)/2+j];
					}
				}
				else
				{
					for(k=j;k<ndg;k++)
					{
						re[t]=re[t]+BM[i*mbw+k-(i-bw)]*TM[((k+1)*k)/2+j];
					}
				}
			}
		}
	}
	return re;
}
//下三角阵与上三角阵为带状结构的矩阵相乘
Complex *resolve_matrix::Multi_TM_LM(Complex *TM,Complex *LM,int ndg,int bw)
{
	Complex *re=new Complex[ndg*ndg-((ndg-bw)*(ndg-bw-1))/2];
	int i,j,k,t;
	for(i=0;i<ndg;i++)
	{
		if(i<ndg-bw-1)
		{
			for(j=0;j<=i+bw;j++)
			{
				t=((2*bw+i+1)*i)/2+j;
				if(j<=bw)
				{
					for(k=0;k<=i;k++)
					{//k此时既是行，又是列，要分两种情况进行计数
						if(k<ndg-bw-1)
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[((2*bw+1+k)*k)/2+j];
						else
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[k*ndg+j-((ndg-bw)*(ndg-bw-1))/2];
					}
				}
				else
				{
					for(k=j-bw;k<=i;k++)
					{
						if(k<ndg-bw-1)
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[((2*bw+1+k)*k)/2+j];
						else
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[k*ndg+j-((ndg-bw)*(ndg-bw-1))/2];
					}
				}
				re[t]=0.0-re[t];
			}
		}
		else
		{
			for(j=0;j<ndg;j++)
			{
				t=i*ndg+j-((ndg-bw)*(ndg-bw-1))/2;
				if(j<=bw)
				{
					for(k=0;k<=i;k++)
					{
						if(k<ndg-bw-1)
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[((2*bw+1+k)*k)/2+j];
						else
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[k*ndg+j-((ndg-bw)*(ndg-bw-1))/2];
					}
				}
				else
				{
					for(k=j-bw;k<=i;k++)
					{
						if(k<ndg-bw-1)
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[((2*bw+1+k)*k)/2+j];
						else
							re[t]=re[t]+TM[((i+1)*i)/2+k]*LM[k*ndg+j-((ndg-bw)*(ndg-bw-1))/2];
					}
				}
				re[t]=0.0-re[t];
			}
		}
	}
	return re;
}
//乔列斯基分解的下三角矩阵的逆与列向量相乘，得到列向量
Complex *resolve_matrix::Multi_CN_Vector(Complex *Inv_Left_up_Tria,Complex *Inv_Light_Down_Tria,Complex *Inv_Left_Down,Complex *Vector,int ndg,int bw)
{
	Complex *re=new Complex[2*ndg];
	int i,j;
	for(i=0;i<ndg;i++)
	{
		for(j=0;j<=i;j++)
		{
			re[i]=re[i]+Inv_Left_up_Tria[((i+1)*i)/2+j]*Vector[j];
			re[i+ndg]=re[i+ndg]+Inv_Light_Down_Tria[((i+1)*i)/2+j]*Vector[j+ndg];
		}
		if(i<ndg-bw-1)
		{
			for(j=0;j<=bw+i;j++)
			{
				re[i+ndg]=re[i+ndg]+Inv_Left_Down[((2*bw+i+1)*i)/2+j]*Vector[j];
			}
		}
		else
		{
			for(j=0;j<ndg;j++)
			{
				re[i+ndg]=re[i+ndg]+Inv_Left_Down[i*ndg+j-((ndg-bw)*(ndg-bw-1))/2]*Vector[j];
			}
		}
	}
	return re;
}
//乔列斯基分解的下三角矩阵的转置与列向量相乘，得到列向量
Complex *resolve_matrix::Multi_C_Trans_Vector(Complex *Left_up_Tria,Complex *Light_Down_Tria,Complex *Light_Up,Complex *Vector,int ndg,int bw)
{
	Complex *re=new Complex[2*ndg];
	int i,j;
	for(i=0;i<ndg;i++)
	{
		for(j=i;j<=bw+i && j<ndg;j++)
		{
			if(j<=bw)
				re[i]=re[i]+conj(Left_up_Tria[j*(bw+1)+i])*Vector[j];
			else
				re[i]=re[i]+conj(Left_up_Tria[j*(bw+1)+i-(j-bw)])*Vector[j];
		}
		if(i<=bw)
		{
			for(j=0;j<=bw+i;j++)
			{
				if(j<=bw)
			    	re[i]=re[i]+conj(Light_Up[j*(bw+1)+i])*Vector[j+ndg];
		    	else
		    		re[i]=re[i]+conj(Light_Up[j*(bw+1)+i-(j-bw)])*Vector[j+ndg];
			}
		}
		else
		{
			for(j=i-bw;j<=bw+i && j<ndg;j++)
			{
				if(j<=bw)
			    	re[i]=re[i]+conj(Light_Up[j*(bw+1)+i])*Vector[j+ndg];
		    	else
		    		re[i]=re[i]+conj(Light_Up[j*(bw+1)+i-(j-bw)])*Vector[j+ndg];
			}
		}
		//re的下半部分
		for(j=i;j<=bw+i && j<ndg;j++)
		{
			if(j<=bw)
		    	re[i+ndg]=re[i+ndg]+conj(Light_Down_Tria[j*(bw+1)+i])*Vector[j+ndg];
	    	else
		    	re[i+ndg]=re[i+ndg]+conj(Light_Down_Tria[j*(bw+1)+i-(j-bw)])*Vector[j+ndg];
		}
	}
	return re;
}
//乔列斯基分解的下三角矩阵的逆的转置与列向量相乘，得到列向量
Complex *resolve_matrix::Multi_CN_Trans_Vector(Complex *Left_up_Tria,Complex *Light_Down_Tria,Complex *Light_Up,Complex *Vector,int ndg,int bw)
{
	Complex *re=new Complex[2*ndg];
	int i,j;
	for(i=0;i<ndg;i++)
	{
		for(j=i;j<ndg;j++)
		{
			re[i]=re[i]+conj(Left_up_Tria[((j+1)*j)/2+i])*Vector[j];
		}
		if(i<=bw)
		{
			for(j=0;j<ndg;j++)
			{
				if(j<ndg-bw-1)
					re[i]=re[i]+conj(Light_Up[((2*bw+1+j)*j)/2+i])*Vector[j+ndg];
				else
					re[i]=re[i]+conj(Light_Up[j*ndg+i-((ndg-bw)*(ndg-bw-1))/2])*Vector[j+ndg];
			}
		}
		else
		{
			for(j=i-bw;j<ndg;j++)
			{
				if(j<ndg-bw-1)
					re[i]=re[i]+conj(Light_Up[((2*bw+1+j)*j)/2+i])*Vector[j+ndg];
				else
					re[i]=re[i]+conj(Light_Up[j*ndg+i-((ndg-bw)*(ndg-bw-1))/2])*Vector[j+ndg];
			}
		}
		//re的下半部分
		for(j=i;j<ndg;j++)
		{
			re[i+ndg]=re[i+ndg]+conj(Light_Down_Tria[((j+1)*j)/2+i])*Vector[j+ndg];
		}
	}
	return re;
}
void resolve_matrix::resolve_bmia_cag_ilucg(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *bm,int ndg,int bw,double ebs)
{
	rex=new Complex[2*ndg];//方程的解
	int i;
	double rfa,beta;
	int mbw=2*bw+1;//
	double sum_rr,sum_r,sum_m;
	//////////////////////////////////////////
	//不完全Cholesky分解，生成下三角矩阵
	Complex *Diag=Creat_Diag(asmn,a1mn,b1mn,ndg,bw,bw);
	Complex *D_Up=new Complex[ndg];
	Complex *D_Down=new Complex[ndg];
    for(i=0;i<ndg;i++)
	{
		D_Up[i]=Diag[i];
		D_Down[i]=Diag[i+ndg];
	}
	Complex *U_Tria_asmn=Creat_triang(asmn,D_Up,ndg,bw);
	Complex *U_Tria_b1mn=Creat_triang(b1mn,D_Down,ndg,bw);
	Diag=Creat_Inv_Diag(Diag,ndg);
	for(i=0;i<ndg;i++)
	{
		D_Up[i]=Diag[i];
		D_Down[i]=Diag[i+ndg];
	}
	//C
	Complex *C_Left_Up=Multi_HalfBM_Diag(U_Tria_asmn,D_Up,ndg,bw);
	Complex *C_Light_Down=Multi_HalfBM_Diag(U_Tria_b1mn,D_Down,ndg,bw);
	Complex *C_Left_Down=Multi_BM_Diag(a1mn,D_Up,ndg,bw);
	//C逆
	Complex *CN_Left_Up=Tria_Decom(C_Left_Up,ndg,bw);
	Complex *CN_Light_Down=Tria_Decom(C_Light_Down,ndg,bw);
	Complex *CN_Left_Down=Multi_BM_TM(C_Left_Down,CN_Left_Up,ndg,bw);
	CN_Left_Down=Multi_TM_LM(CN_Light_Down,CN_Left_Down,ndg,bw);
	//////////////////////////////////////////////////////////////////////////////
	//ILUCG
	Complex *r=new Complex[2*ndg];//余量
	Complex *p=new Complex[2*ndg];

	
	//第0步
	Complex *t1;
	Complex *ap;

	/////////////////////////////////////
	//t0=BandMatrix_multiply(amt,bm,n,bw);
		
	for(i=0;i<2*ndg;i++)
	{
		r[i]=bm[i];//0步下半部分为零
	}
	p=Multi_UALR(asmn,bsmn,a1mn,b1mn,CN_Left_Up,CN_Light_Down,CN_Left_Down,r,ndg,bw);
	t1=Multi_CN_Vector(CN_Left_Up,CN_Light_Down,CN_Left_Down,r,ndg,bw);
	sum_rr=L2_Norm(t1,2*ndg);

	
	//double Con_N=L2_Norm(bm,ndg);
	//cout<<Con_N<<endl;
	//第一步到第n步
	int t=0;
	while(L2_Norm(r,2*ndg) > ebs)
	{
		t++;
		cout<<L2_Norm(r,2*ndg)<<"     "<<t<<endl;
		//
		
		t1=Multi_C_Trans_Vector(C_Left_Up,C_Light_Down,C_Left_Down,p,ndg,bw);
		sum_m=L2_Norm(t1,2*ndg);
		rfa=sum_rr/sum_m;
		//cout<<sum_rr<<endl;
		//multiply(am,d,n,bw);
		ap=Multi_A_R(asmn,bsmn,a1mn,b1mn,p,ndg,bw);

		for(i=0;i<2*ndg;i++)
		{
			rex[i]=rex[i]+rfa*p[i];
			r[i]=r[i]-rfa*ap[i];
			//cout<<r[i]<<endl;
		}
		
		t1=Multi_CN_Vector(CN_Left_Up,CN_Light_Down,CN_Left_Down,r,ndg,bw);
		sum_r=L2_Norm(t1,2*ndg);
		beta=sum_r/sum_rr;
		t1=Multi_UALR(asmn,bsmn,a1mn,b1mn,CN_Left_Up,CN_Light_Down,CN_Left_Down,r,ndg,bw);
		for(i=0;i<2*ndg;i++)
		{
			p[i]=t1[i]+beta*p[i];//2
			//cout<<rex[i]<<"        "<<i<<endl;
		}
		sum_rr=sum_r;
	}
	delete [] r;
	delete [] p;
	delete t1;
	delete ap;
	delete Diag;
	delete [] D_Up;
	delete [] D_Down;
	delete U_Tria_asmn;
	delete U_Tria_b1mn;
	delete C_Left_Up;
	delete C_Light_Down;
	delete C_Left_Down;
	delete CN_Left_Up;
	delete CN_Light_Down;
	delete CN_Left_Down;
}
Complex *resolve_matrix::Multi_UALR(Complex *asmn,Complex *bsmn,Complex *a1mn,Complex *b1mn,Complex *CN_L_U,Complex *CN_Li_D,Complex *CN_L_D,Complex *VR,int ndg,int bw)
{
	Complex *re=new Complex[2*ndg];
	int i;
	Complex *a1mnt,*b1mnt,*asmnt,*bsmnt;
	Complex *t1=new Complex[ndg];
	Complex *t2=new Complex[ndg];
	Complex *t3,*t4,*t5,*t6;
		
	//矩阵的共轭转置
	
	asmnt=Transpose_BandMatrix(asmn,ndg,bw);
	bsmnt=Transpose_BandMatrix(bsmn,ndg,bw);
	a1mnt=Transpose_BandMatrix(a1mn,ndg,bw);
	b1mnt=Transpose_BandMatrix(b1mn,ndg,bw);

	re=Multi_CN_Vector(CN_L_U,CN_Li_D,CN_L_D,VR,ndg,bw);
	re=Multi_CN_Trans_Vector(CN_L_U,CN_Li_D,CN_L_D,re,ndg,bw);
	for(i=0;i<ndg;i++)
	{
		t1[i]=re[i];
		t2[i]=re[i+ndg];
	}
	//上
	t3=BandMatrix_multiply(asmnt,t1,ndg,bw);//asmnt与r的 上 半部分d1相乘，维数为ndg
	t4=BandMatrix_multiply(a1mnt,t2,ndg,bw);//a1mnt与r的 下 半部分d2相乘，维数为ndg
	t3=Plus_Matrix(t3,t4,ndg);//上
   	//下
	t5=BandMatrix_multiply(bsmnt,t1,ndg,bw);//bsmnt与r的 上 半部分d1相乘，维数为ndg 下 
	t6=BandMatrix_multiply(b1mnt,t2,ndg,bw);//b1mnt与r的 下 半部分d2相乘，维数为ndg 下 
	t5=Plus_Matrix(t5,t6,ndg);//相加， 下 
	re=Unite_Matrix(t3,t5,ndg);//上下合并
	re=Multi_CN_Vector(CN_L_U,CN_Li_D,CN_L_D,re,ndg,bw);
	re=Multi_CN_Trans_Vector(CN_L_U,CN_Li_D,CN_L_D,re,ndg,bw);

	delete a1mnt;
	delete b1mnt;
	delete asmnt;
	delete bsmnt;
	delete [] t1;
	delete [] t2;
	delete t3;
	delete t4;
	delete t5;
	delete t6;

	return re;
}