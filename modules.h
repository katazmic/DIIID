#ifndef MD_HD/Users/katyghantous/Desktop/420/TOYM1.5 2/modules.h
#define Md_HD





/********************************************/
/* the plasma profiles are set from transp  */
/* while the machine parameters are         */
/* specified in the main file.              */
/********************************************/


typedef struct{

  double *T_r;  
  double *Te_r;
  double *beta_pc_r;
  double *beta_elec;
  double *beta_beam;
  double *del_beta_beam;
  double *neut_ini; // check since it's without *a in 
  double *q;
  double *del_q;
  double *shear;
  double *ni;
  double *r;
  double sgm;
  double R;
  double a;
  double B; 
} discharge_;


typedef struct{

  double *vA;
  double *Di;
  double *Do;
  double *wtae;
 
  double *crit;
  double *gamma_b; 
  double *gamma_e;
  double *gamma_i;
  double *gamma_iLT;
  double *gamma_rad; 
  double *gamma_extra; 
  double *sgm_gauss; 
  double *gamma_EP_prm; 
} TAE_;


 
typedef struct{ 
  double *Db;
  double *orb_wd;
  double *w_dia;
  double *wd_hs;
  double gamm;
  double E_EP;
  double mEP;
  double Chi0;
  double Z;
  double mu; /* mass ion/mass proton = 2 for deutrium */
  double v_EP_0;
  double w_c;
  double rho_EP;
  int edge;
} EP_;

typedef struct{
  double *profile;
  double *neutron_profile;
  double EP_loss;
  double neutron_loss;
} relaxed_;

typedef struct{
  double r_a;
  double gamma_b;
  double gamma_e;
  double gamma_i;
  double gamma_iLT;
  double gamma_rad;
  double gamma_extra;
} point_;

typedef struct{
  int r;
  int q;
  int gamma_b; int gamma_e;  int gamma_i;  int  gamma_iLT;int  gamma_rad; int gamma_extra;
  
  int beta_alpha_ini;  int beta_alpha_rel; 
  int beta_alpha_crit;int del_beta_alpha;

  int beta_beam_ini;int beta_beam_rel;
  int beta_beam_crit;int del_beta_beam;
} visualization_;



/***********************************************************************/
/* set_discharge function uses trasnp file to interpolate for NUM long */
/* arrays for the profiles.NUM id givrn in main                        */
/***********************************************************************/
void set_discharge(discharge_ *discharge, char *fname,int Ntransp,  double sgm, int NUM);

void malloc_discharge(discharge_ **discharge, int NUM);

void get_data(FILE *ftransp, const char * name_data,  int Ntransp, double *ary);
double get_data_NAME(FILE *someFile, const char * nameData);
int get_data_NAME_array(double *vlu, int Npts,  FILE *someFile, const char * nameData);
int get_Ntransp(FILE *someFile);

void NormalizeArray(TAE_ *TAE, discharge_ discharge, int Npt, point_ *point,int NUM);

void NormalizeArray_Fac(TAE_ *TAE, discharge_ discharge, int Np, point_ *point,int NUM);
void set_EP( EP_ *EP, char *fname, discharge_ discharge,  double mEP, double Z, double E_EP , double gamm, double bndFac, int NUM);
void malloc_EP(EP_ **EP, int NUM);


void set_TAE(TAE_ *TAE, EP_ EP, discharge_ discharge,int NUM); 
void malloc_TAE(TAE_ **TAE, int NUM);

void malloc_relaxed(relaxed_ **relaxed, int NUM);
int compute_relaxed( relaxed_ * relaxed, discharge_ discharge, TAE_ TAE, EP_ EP, int NUM);
void compute_neut_rlxd( relaxed_ *relaxed,discharge_ discharge,  int NUM);
void compute_loss( relaxed_ * relaxed, discharge_ discharge, int NUM);


void heaviside_smoothen(relaxed_ *relaxed, EP_ EP,int NUM);
void smoothen_gauss_rates(TAE_ * TAE, discharge_ discharge, int NUM);


void read_input_visualization(visualization_ *visualization, char *fname);


void malloc_visualization(visualization_ **visualization);


void read_input_points(point_ *point1, point_ *point2, char *fname, int *shot, int *discharge);
void read_input_points_array(point_ *point, int N, char *fname, int *shot, int *discharge, int * NUM, double * boundary);
void set_point(point_ *point, double r, double gamma_b, double gamma_e, double gamma_i, double gamma_iLT, double gamma_rad,double gamma_extra);
void malloc_point(point_ ** point);
void malloc_point_array(point_ ** point,int N);

void Normalize(TAE_ *TAE, discharge_ discharge, point_ point1, point_ point2,int NUM);


void free_discharge(discharge_ *discharge);
void free_TAE(TAE_ *TAE);
void free_EP(EP_ *EP);
void free_relaxed(relaxed_ *relaxed);
void free_point(point_ *point);
void free_point_array(point_ *point);
void free_visualization(visualization_ *visualization);

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);


#endif
