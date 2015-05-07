#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<errno.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_errno.h>
#include "modules.h"
#define PI 3.1415
#define NR_END 1
#define FREE_ARG char*
     


/****************************************************************************/
/*                  Reading from the INPUT file                             */
/****************************************************************************/
void read_input_visualization(visualization_ *visualization, char *fname)
{
  FILE *ftransp;
  int r, gamma_b, gamma_e,  gamma_i,  gamma_iLT,  gamma_rad,gamma_extra, beta_alpha_ini, beta_alpha_crit, beta_alpha_rel,  beta_beam_ini, beta_beam_crit, beta_beam_rel,q,  del_beta_beam,  del_beta_alpha;
 

  ftransp = fopen(fname,"r");   
  if(ftransp==NULL){
    printf("\n Provide INPUT file! Refer to README for more information.\n \n");
    abort();
  }


  r = get_data_NAME(ftransp, "r");  
  fclose(ftransp);
  

  ftransp = fopen(fname,"r");   
  gamma_b = get_data_NAME(ftransp, "P_gmma_b");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  gamma_e = get_data_NAME(ftransp, "P_gmma_e");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_i = get_data_NAME(ftransp, "P_gmma_i");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_iLT = get_data_NAME(ftransp, "P_gmma_iLT");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_rad = get_data_NAME(ftransp, "P_gmma_rad");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_extra = get_data_NAME(ftransp, "P_gmma_extra");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  beta_beam_ini = get_data_NAME(ftransp, "beta_beam_initial");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  beta_alpha_ini = get_data_NAME(ftransp, "beta_alpha_initial");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  beta_beam_rel = get_data_NAME(ftransp, "beta_beam_relaxed");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  beta_alpha_rel = get_data_NAME(ftransp, "beta_alpha_relaxed");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  beta_beam_crit = get_data_NAME(ftransp, "beta_beam_crit_gradient");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  del_beta_beam = get_data_NAME(ftransp, "del_beta_beam");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  beta_alpha_crit = get_data_NAME(ftransp, "beta_alpha_crit_gradient");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  del_beta_alpha = get_data_NAME(ftransp, "del_beta_alpha");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  q = get_data_NAME(ftransp, "safety_factor");
  fclose(ftransp);



  visualization->r = r;

  visualization->gamma_b = gamma_b;
  visualization->gamma_e = gamma_e;
  visualization->gamma_i = gamma_i;
  visualization->gamma_iLT = gamma_iLT;
  visualization->gamma_rad = gamma_rad;
  visualization->gamma_extra = gamma_extra;
 
  visualization->q = q;
  visualization->beta_beam_crit = beta_beam_crit;
  visualization->del_beta_beam = del_beta_beam;

  visualization->beta_beam_ini = beta_beam_ini;
  visualization->beta_beam_rel = beta_beam_rel;


  visualization->beta_alpha_crit = beta_alpha_crit;
  visualization->del_beta_alpha = del_beta_alpha;

  visualization->beta_alpha_ini = beta_alpha_ini;
  visualization->beta_alpha_rel = beta_alpha_rel;

}
/****************************************************************************/
/*               allocating memory for vizualization                        */
/****************************************************************************/

void malloc_visualization(visualization_ ** visualization)
{
  *visualization = malloc(sizeof(visualization_));
}

/****************************************************************************/
/*                  Reading from the INPUT file                             */
/****************************************************************************/
void read_input_points(point_ *point1, point_ *point2, char *fname, int *shot, int *discharge)
{
  FILE *ftransp;
  double r, gamma_b, gamma_e,  gamma_i,  gamma_iLT,  gamma_rad;
 

  ftransp = fopen(fname,"r");   
  if(ftransp==NULL){
    printf("\n Provide INPUT file! Refer to README for more information.\n \n");
    abort();
  }
  r = get_data_NAME(ftransp, "r_a1");  
  fclose(ftransp);
  

   ftransp = fopen(fname,"r");   
   gamma_b = get_data_NAME(ftransp, "gamma_b1");
   fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_e = get_data_NAME(ftransp, "gamma_e1");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_i = get_data_NAME(ftransp, "gamma_i1");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_iLT = get_data_NAME(ftransp, "gamma_iLT1");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_rad = get_data_NAME(ftransp, "gamma_rad1");
  fclose(ftransp);


  point1->r_a = r;
  point1->gamma_b = gamma_b;
  point1->gamma_e = gamma_e;
  point1->gamma_i = gamma_i;
  point1->gamma_iLT = gamma_iLT;
  point1->gamma_rad = gamma_rad;

  ftransp = fopen(fname,"r");   
  r = get_data_NAME(ftransp, "r_a2");
  fclose(ftransp);

  ftransp = fopen(fname,"r");     
  gamma_b = get_data_NAME(ftransp, "gamma_b2");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_e = get_data_NAME(ftransp, "gamma_e2");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  gamma_i = get_data_NAME(ftransp, "gamma_i2");
  fclose(ftransp);
  
  ftransp = fopen(fname,"r");   
  gamma_iLT = get_data_NAME(ftransp, "gamma_iLT2");
  fclose(ftransp);

  
  ftransp = fopen(fname,"r");   
  gamma_rad = get_data_NAME(ftransp, "gamma_rad2");
  fclose(ftransp);


 
  point2->r_a = r;
  point2->gamma_b = gamma_b;
  point2->gamma_e = gamma_e;
  point2->gamma_i = gamma_i;
  point2->gamma_iLT = gamma_iLT;
  point2->gamma_rad = gamma_rad;



  ftransp = fopen(fname,"r");   
  *shot = (int) get_data_NAME(ftransp, "time");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  *discharge = (int) get_data_NAME(ftransp, "discharge");
  fclose(ftransp);


}

/****************************************************************************/
/*                  Reading from the INPUT file                             */
/****************************************************************************/
void read_input_points_array(point_ *point, int N, char *fname, int *shot, int *discharge, int * NUM, double * boundary)
{

  int i;
  FILE *ftransp;
  double r[N], gamma_b[N], gamma_e[N],  gamma_i[N],  gamma_iLT[N],  gamma_rad[N], gamma_extra[N];
 

  ftransp = fopen(fname,"r");   
  if(ftransp==NULL){
    printf("\n Provide INPUT file! Refer to README for more information.\n \n");
    abort();
  }
  get_data_NAME_array(r,N,ftransp, "radius");  
  fclose(ftransp);
  

   ftransp = fopen(fname,"r");   
   get_data_NAME_array(gamma_b,N,ftransp, "gamma_b");
   fclose(ftransp);

  ftransp = fopen(fname,"r");   
  get_data_NAME_array(gamma_e,N,ftransp, "gamma_e_coll");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  get_data_NAME_array(gamma_i,N,ftransp, "gamma_iD");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  get_data_NAME_array(gamma_iLT,N,ftransp, "gamma_iLT");
  fclose(ftransp);

  ftransp = fopen(fname,"r");   
  get_data_NAME_array(gamma_rad,N,ftransp, "gamma_rad");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  get_data_NAME_array(gamma_extra,N,ftransp, "gamma_extra");
  fclose(ftransp);


  for(i=0;i<N;i++){

    (point+i)->r_a = r[i];
    (point+i)->gamma_b = gamma_b[i];
    (point+i)->gamma_e = gamma_e[i];
    (point+i)->gamma_i = gamma_i[i];
    (point+i)->gamma_iLT = gamma_iLT[i];
    (point+i)->gamma_rad = gamma_rad[i];
    (point+i)->gamma_extra = gamma_extra[i];
  }


  ftransp = fopen(fname,"r");   
  *shot = (int) get_data_NAME(ftransp, "time");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  *discharge = (int) get_data_NAME(ftransp, "discharge");
  fclose(ftransp);
  
  ftransp = fopen(fname,"r");   
   *NUM = (int) get_data_NAME(ftransp, "NUM of GRIDS");
  fclose(ftransp);


  ftransp = fopen(fname,"r");   
  *boundary = (double) get_data_NAME(ftransp, "BOUNDARY");
  fclose(ftransp);

 
}
/****************************************************************************/
/*  Setting up the various components for finding the critical gradient for */ 
/*  getting marginal stability                                              */
/****************************************************************************/
 


/****************************************************************************/
/* Set_discharge function reads the profiles from TRANSP file and fits the  */
/* points NUM long arrays to work with                                      */
/****************************************************************************/
void set_discharge(discharge_ *discharge, char *fname, int Ntransp,  double sgm, int NUM)
{

  int i,ird;


  FILE *ftransp;      
  
  ////////////////////////
  double ary[Ntransp],rtransp[Ntransp], y2[Ntransp], r[NUM],tempary[Ntransp];
  double temp,y1,yn,R,a,B,BR;
  ////////////////////////


  //////for del_beta//////////



      gsl_interp *workspace;
      gsl_interp_accel *accel;
      int EXT=1;
      int EXT2 = 3;
      int Ntemp = 5;  // or is it 25?
      int DIM = EXT+EXT2+Ntemp;

      double rTemp[DIM], aryTemp[DIM];

      int NUMEXT;
  /////////////
  ftransp = fopen(fname,"r");   
  BR = get_data_NAME(ftransp, "BZXR");
  fclose(ftransp);
 
      
  ftransp = fopen(fname,"r"); 
  // for some reason it runs if this printf is activated??
  printf(" in set_discharge Ntransp=%d\n", Ntransp);
  get_data(ftransp, " RMJMP ", Ntransp,tempary);
  fclose(ftransp);
  R = (tempary[0] +tempary[Ntransp-1])/2;
  
  discharge->R = R;      
 
  
  B = BR/R;
  discharge->B = B;
 
  
  ftransp = fopen(fname,"r"); 
  get_data(ftransp, " RMNMP ", Ntransp,tempary);
  fclose(ftransp);
  
  a = tempary[Ntransp-1];
  discharge->a = a;

  discharge->sgm = sgm;

  /****************************************************************************/
      /*-------------------
	getting radius
	-------------------*/
      ftransp = fopen(fname,"r");   

 
      get_data(ftransp, " X ", Ntransp,rtransp);
 
     fclose(ftransp);
      for(i=0;i<Ntransp;i++)
	rtransp[i] = rtransp[i]*a;

      for(i=0;i<NUM;i++)
	{
	  r[i] = (0.01 + (((float) i)/( (float) NUM))*0.99)*a;
	  discharge->r[i] = r[i];
	}
      

      /*----------------------------------------
	polynomial coef for the ion temperature
	-----------------------------------------*/   
      
      ftransp = fopen(fname,"r");   
      get_data(ftransp, "ION TEMPERATURE", Ntransp,ary);
      fclose(ftransp);
      
      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);
     
      spline(rtransp, ary, Ntransp,  y1, yn, y2);      
      for(i=0;i<NUM;i++)
	{temp = 0;
	  splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->T_r[i] = temp;
	}


     /*----------------------------------------
     polynomial coef for the electron temperature
     ----------------------------------------*/
     ftransp = fopen(fname,"r");   
     get_data(ftransp, "ELECTRON TEMPERATURE", Ntransp,ary);
     fclose(ftransp);

     y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
     yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);

     spline(rtransp, ary, Ntransp,  y1, yn, y2);
   
     for(i=0;i<NUM;i++)
       {
	 temp = 0;
	 splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	 discharge->Te_r[i] = temp;
       }

      /*-------------------------------------
      polynomial coeficients for plasma beta
      ---------------------------------------*/

      ftransp = fopen(fname,"r");   
      get_data(ftransp, "PLASMA BETA TOROIDAL", Ntransp,ary);
      fclose(ftransp);

      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);

      spline(rtransp, ary, Ntransp,  y1, yn, y2);
    
      for(i=0;i<NUM;i++)
	{
	  temp = 0;
	  splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->beta_pc_r[i] = temp;
	}
 
      /*------------------------------------
      polynomial for beta of electron     
      -------------------------------------*/
      ftransp = fopen(fname,"r");   
      get_data(ftransp, "ELECTRON BETA TOROIDAL", Ntransp,ary);
      fclose(ftransp);

      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);
      
      spline(rtransp, ary, Ntransp,  y1, yn, y2);
    

      for(i=0;i<NUM;i++)
	{temp = 0;
	  splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->beta_elec[i] = temp;
	}

      /*------------------------------
      poly for ion density 
      -------------------------------------  */   
      ftransp = fopen(fname,"r");   
      get_data(ftransp, "TOTAL ION DENSITY", Ntransp,ary);
      fclose(ftransp);

      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);

      spline(rtransp, ary, Ntransp,  y1, yn, y2);
    

      for(i=0;i<NUM;i++)
	{temp = 0;
	  splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->ni[i] = temp;
	}

      /*-------------------------------
      poly soof to fit beta_beam 
      -----------------------------------*/

      workspace = gsl_interp_alloc(gsl_interp_polynomial, DIM);
      accel = gsl_interp_accel_alloc();  

      ftransp = fopen(fname,"r");   
      get_data(ftransp, "BEAM BETA TOROIDAL", Ntransp,ary);
      fclose(ftransp);
   
  
      for(ird = EXT;ird<(DIM-EXT2);ird++)
      	{
	  aryTemp[ird] = ary[(Ntransp/Ntemp)*(ird-EXT)];
	  rTemp[ird] = r[(NUM/Ntemp)*(ird-EXT)];
	}
     
      for(ird = EXT-1; ird>-1;ird--)
	{
	  aryTemp[ird] = ary[0];
	  rTemp[ird] = (ird-EXT)*(rTemp[EXT+1]-rTemp[EXT]);
	}
      
      for(ird = DIM-EXT2; ird<DIM;ird++)
	{
	  aryTemp[ird] = ary[Ntransp-1];
	  rTemp[ird] = (ird-EXT)*(rTemp[EXT+1]-rTemp[1]);
	}


      ////////////////
         gsl_interp_init(workspace, rTemp, aryTemp, DIM);
      ////////////// 
      NUMEXT = NUM+ (int) ((rTemp[EXT]-rTemp[0])/(r[2]-r[1])) + (int) ((rTemp[DIM-1]-rTemp[Ntemp+EXT-1])/(r[2]-r[1]));
   
      float rEXT[NUMEXT], beta_beamEXT[NUMEXT];
      
      for(i = 0; i<NUMEXT;i++)
	rEXT[i]= rTemp[0] + i*(r[2]-r[1]);


      for(i=0;i<NUMEXT;i++)
	{
	  temp = 0;
	  temp = gsl_interp_eval(workspace, rTemp, aryTemp, rEXT[i], NULL);
	  beta_beamEXT[i] = temp;
	}
      for(i=0;i<NUM;i++)
	discharge->beta_beam[i] = beta_beamEXT[i + ((int) ((rTemp[EXT]-rTemp[0])/(r[2]-r[1])))];
      

    gsl_interp_free(workspace);

    gsl_interp_accel_free(accel);
    
      /**************************
            del_beta_beam 
      ****************************/ 
      for(i= 1;i<NUM;i++)
	discharge->del_beta_beam[i] = (discharge->beta_beam[i]-discharge->beta_beam[i-1])/(r[i]-r[i-1]);

      discharge->del_beta_beam[0] = discharge->del_beta_beam[1];


      /* --------------------------
      Poly soof for the q profile
      ---------------------------*/
 
      ftransp = fopen(fname,"r");   
      get_data(ftransp, "Q ", Ntransp,ary);
      fclose(ftransp);

      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);

      spline(rtransp, ary, Ntransp,  y1, yn, y2);
      

      for(i=0;i<NUM;i++)
	{temp = 0;
	splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->q[i] = temp;
	  if(i !=0)
	   {
	     discharge->del_q[i] = (discharge->q[i]-discharge->q[i-1])/(r[i]-r[i-1]);
	     discharge->shear[i] = r[i]*discharge->del_q[i]/discharge->q[i];	  
	   }
	}
      discharge->del_q[0] = discharge->del_q[1];
      discharge->shear[0] = discharge->shear[1];

    /*-----------------------------
     poly fit for initial neutron distribution 
     ------------------------------*/

      ftransp = fopen(fname,"r");   
      get_data(ftransp, "TOTAL NEUTRONS ", Ntransp,ary);
      fclose(ftransp);


      y1 = (ary[1]-ary[0])/(rtransp[1]-rtransp[0]);
      yn = (ary[Ntransp-1] - ary[Ntransp-2])/(rtransp[Ntransp-1]-rtransp[Ntransp-2]);

      spline(rtransp, ary, Ntransp,  y1, yn, y2);
      

      for(i=0;i<NUM;i++)
	{
	  temp = 0;
	  //  discharge->r[i] = (0.01 + (i/( (float) NUM))*0.99);
	  splint(rtransp, ary, y2,Ntransp, r[i], &temp);	
	  discharge->neut_ini[i] = temp;
	}

}


void malloc_discharge(discharge_ **discharge, int NUM)
{
  *discharge = malloc(sizeof(discharge_));
  (*discharge)->T_r = malloc(NUM*sizeof(double));
  (*discharge)->Te_r = malloc(NUM*sizeof(double));
  (*discharge)->beta_pc_r = malloc(NUM*sizeof(double));
  (*discharge)->beta_elec = malloc(NUM*sizeof(double));
  (*discharge)->beta_beam = malloc(NUM*sizeof(double));
  (*discharge)->del_beta_beam = malloc(NUM*sizeof(double));
  (*discharge)->neut_ini = malloc(NUM*sizeof(double));
  (*discharge)->q = malloc(NUM*sizeof(double));
  (*discharge)->del_q = malloc(NUM*sizeof(double));
  (*discharge)->shear = malloc(NUM*sizeof(double));
  (*discharge)->ni = malloc(NUM*sizeof(double));
  (* discharge)->r = malloc(NUM*sizeof(double));

}


/****************************************************************************/
/* get_data reads the values from TRANSP                                    */
/****************************************************************************/
void get_data(FILE *ftransp, const char * name_data, int Ntransp, double *ary)
{
  int dataExists;
  int m,i,j, k;
  char dat;
  char data[Ntransp*Ntransp];
  char temp[Ntransp];
  int dg;
  dg = 0;
  m = 0;
  dataExists = 0;
  dat = fgetc(ftransp);


 
 	   dat = fgetc(ftransp);
	  
 	   while(dat != EOF && dataExists == 0)
	     {
	       i=0;     
	       dat = fgetc(ftransp);

	       if(dat == name_data[i])
		 {
		   
		   while(dat == name_data[i])
		     {
		       dat = fgetc(ftransp);
		       i++;
		       
		     }
		   
		   if(name_data[i] == '\0')
		     { 
		       dataExists = 1;		  
		     }
		 }
	     }
	   
	   while(dat !='\n')
	     {
	     dat = fgetc(ftransp);
	     }

	   while(dat == ' '|| dat == '\n')
	     {
	     dat = fgetc(ftransp);
	     }
	  
	   
	   while(m<Ntransp)
	     {
	       
	       while((dat != ' ') && (dat != '\n'))
		 {
		   data[dg] = dat;
		   dat = fgetc(ftransp);		 
		   dg++;
		 }
	       
	        data[dg] = '\t';
		dg++;
	       
	       while((dat == ' ') || (dat == '\n'))
		 dat = fgetc(ftransp);	
	       
	       m++;
	       
	    }
	   data[dg+1] = '\0';   
	   i=0;
	   j=0;
	   
	   while(data[i] !='\0')
	     {k=0;
	     while(data[i]!='\t')
	       {
		 temp[k] = data[i];
		 i++;
		 k++;
	       }
	     temp[k] = '\0';
	     
	     ary[j]= atof(temp);
	     i++;
	     j++;
	     }


}




/****************************************************************************/
/* set_TAE builds TAE stucture specific to the discarge  and computes the   */ 
/* various damping and growth rates associated as well the critical gradient*/
/****************************************************************************/
void set_TAE( TAE_ *TAE, EP_ EP, discharge_ discharge,int NUM)
{
  int i;
  double r[NUM];

  double  Vsnd , rhosnd2,   v_iD,   v_iT, v_e,nq,m, x_iD, x_iT, x_b, x_e, n_e, nu_omg, Zeff, I1, I2, fac, chi0sq , v_b0sq, v_Asq, D_i1,  D_i2, a, zeta, p, Gfow, B,gamm,epsln;
  
  epsln=0.000001;
 
  B = discharge.B;
  gamm =  EP.gamm;

  for(i=0;i<NUM;i++)
    {
      r[i] = discharge.r[i];


      Vsnd = 9.79*pow(10,5)*sqrt(gamm*discharge.Te_r[i]/EP.mEP);
      rhosnd2 = Vsnd*Vsnd/(EP.w_c*EP.w_c);
      
      v_iD = 9.79*pow(10,5)*sqrt(2*discharge.T_r[i]/2);
      v_iT = 9.79*pow(10,5)*sqrt(2*discharge.T_r[i]/3);

      v_e = 4.19*pow(10,7)*sqrt(discharge.Te_r[i]);
     
      TAE->vA[i] = sqrt(1/(EP.mEP*discharge.ni[i]))*2.18*pow(10,11)*B*10000;


      //      nq = r[i]*EP.w_c/EP.v_EP_0;
      nq =  discharge.R*EP.w_c/(EP.v_EP_0*discharge.q[i]); // r[i]*EP.w_c/EP.v_EP_0;
      m = fabs(nq-0.5) ; // instead of nq as initial code..
    
      x_iD = TAE->vA[i]/(3*v_iD);
      x_iT = TAE->vA[i]/(3*v_iT);
      x_b =  TAE->vA[i]/EP.v_EP_0;
      x_e =  TAE->vA[i]/v_e;

   
      n_e = pow(10,11)*B*10000*B*10000*discharge.beta_elec[i]/(4.01*discharge.Te_r[i]);

      if(discharge.beta_elec[i]<0)
	n_e = 0;


      nu_omg = 2.91*pow(10,-6)*16*pow(discharge.Te_r[i],-1.5)*n_e*2*discharge.q[i]*(discharge.R)/TAE->vA[i];
      
      Zeff = 1.3; /* in case of D3-D the dominant impurity is carbon*/

      I1=  0.43*Zeff +1.06;

      I2 = 1.03*Zeff + 2.3;

      fac = I1*pow((8*discharge.shear[i]*nq*sqrt(rhosnd2)/(5*r[i]*r[i]/discharge.R)),2) + I2*discharge.q[i]*discharge.q[i]*(8*discharge.beta_pc_r[i]/(1+discharge.sgm)); 
 

      TAE->gamma_e[i] =   -(sqrt(PI/2)/4*(fac)*sqrt(nu_omg))/(sqrt(fabs(log(16*sqrt(r[i]/(discharge.R*nu_omg)))))*sqrt(fabs(log(16*sqrt(r[i]/(discharge.R*nu_omg)))))*sqrt(fabs(log(16*sqrt(r[i]/(discharge.R*nu_omg)))))) - epsln;



//  TAE->gamma_e[i] = -(sqrt(PI/2)/4*(fac)*sqrt(nu_omg))*pow((log(16*sqrt(r[i]/(discharge.R*nu_omg)))),(-1.5));
       
      TAE->gamma_i[i] = -discharge.q[i]*discharge.q[i]*discharge.sgm*sqrt(PI)*pow(x_iD,3)*exp(-x_iD*x_iD)/(18*(1+discharge.sgm/4)) -epsln ;


   
      TAE->gamma_iLT[i] = -discharge.q[i]*discharge.q[i]*discharge.sgm*sqrt(PI)*pow(x_iT,3)*exp(-x_iT*x_iT)/(18*(1+discharge.sgm/4)) -epsln ;
      
      TAE->gamma_rad[i] = -3*pow( fabs(sqrt(rhosnd2)/r[i]*discharge.shear[i]*nq*(nq+1)/(2*nq+1)/1.4142 ),0.67) +epsln;
 
      TAE->gamma_extra[i] = -1 ;
  
      ////////////////// for growth rates //////////////////////    
      // this is only vA/3 resonance and only for the above the plato point which is 80c in Fulop's paper PPCF 96

      chi0sq = EP.Chi0*EP.Chi0;
      v_b0sq = EP.v_EP_0*EP.v_EP_0;
      v_Asq = (TAE->vA[i])*(TAE->vA[i]);

      TAE->wtae[i] = (TAE->vA[i])/(2*discharge.q[i]*discharge.R);
      // this is the width of the TAE mode
      D_i1 = 2.5*PI*r[i]*r[i]/(8*m*discharge.R);
      D_i2 = pow(r[i],1.5)/(sqrt(fabs(discharge.shear[i])*discharge.R)*m);
 

      if(D_i1>D_i2)
	TAE->Di[i] = D_i1;
      else
	TAE->Di[i] = D_i2;

      TAE->Do[i] = EP.v_EP_0/(EP.w_c);

      a = (EP.w_dia[i])/(TAE->wtae[i]);
  
      zeta = (1+chi0sq)*(discharge.q[i])*(EP.v_EP_0)/((EP.Chi0)*(EP.w_c));

      p = 2.0; 
      Gfow = pow((1+pow(EP.Db[i]/TAE->Di[i],p)),(-1/p));
     
      TAE->gamma_b[i] = - EP.Db[i]*Gfow*discharge.q[i]*discharge.q[i]*discharge.beta_beam[i]*a*(1+chi0sq)*(1+ chi0sq)*24*(TAE->Do[i])*(TAE->Do[i])/pow(zeta,3)*((TAE->vA[i])/(3*EP.Chi0*EP.v_EP_0)) +epsln;

      TAE->sgm_gauss[i] = TAE->Di[i]/2;/*averging the grwth rate over mode width */
      TAE->gamma_EP_prm[i] = TAE->gamma_b[i]/discharge.del_beta_beam[i];

      TAE->crit[i] =  (TAE->gamma_e[i]+ TAE->gamma_rad[i] +TAE->gamma_i[i]+TAE->gamma_iLT[i])*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);

         
    }
 
 
}


void malloc_TAE(TAE_ **TAE, int NUM)
{
  *TAE = malloc(sizeof(TAE_));
  
  (*TAE)->vA = malloc(NUM*sizeof(double));
  (*TAE)->Di = malloc(NUM*sizeof(double));
  (*TAE)->Do = malloc(NUM*sizeof(double));
  (*TAE)->wtae = malloc(NUM*sizeof(double));


  (*TAE)->gamma_b = malloc(NUM*sizeof(double));
  (*TAE)->gamma_e = malloc(NUM*sizeof(double));
  (*TAE)->gamma_i = malloc(NUM*sizeof(double));
  (*TAE)->gamma_iLT = malloc(NUM*sizeof(double));
  (*TAE)->gamma_rad = malloc(NUM*sizeof(double));
  (*TAE)->gamma_extra = malloc(NUM*sizeof(double));
  (*TAE)->crit = malloc(NUM*sizeof(double));
  (*TAE)->sgm_gauss = malloc(NUM*sizeof(double));
  (*TAE)->gamma_EP_prm = malloc(NUM*sizeof(double));  


}

/****************************************************************************/
/* set_EP builds the EP structure depending on the species                  */
/****************************************************************************/
void set_EP( EP_ *EP, char *fname, discharge_ discharge, double mEP, double Z, double E_EP, double gamm,double bndFac, int NUM)
{
  int i;
  FILE *ftransp;
  double Chi0;
  ftransp = fopen(fname,"r");   
  Chi0 = get_data_NAME(ftransp, "COFRC_D");
  fclose(ftransp);


  EP->E_EP = E_EP;
  EP->mEP = mEP;
  EP->Chi0 = Chi0;
  EP->Z = Z;
  EP->gamm = gamm;
  
  EP->w_c =  9.58*pow(10,3)*(Z*(discharge.B)*10000/mEP) ;/* deutrium cyclotron freq nu = 2. p.28 formulary*/

  EP->v_EP_0 = 9.79*pow(10,5)*sqrt(E_EP*pow(10,6)/mEP);

  EP->rho_EP = 9.79*pow(10,5)*sqrt(E_EP*pow(10,6)/mEP)/( 9.58*pow(10,3)*(discharge.B*10000/mEP));

  for( i=0;i<NUM;i++)
    {   
      EP->Db[i] = (EP->v_EP_0/EP->w_c)*discharge.q[i];

      EP->w_dia[i] =  (EP->v_EP_0)*discharge.del_beta_beam[i]/discharge.beta_beam[i];

      EP->orb_wd[i] = discharge.q[i]*(EP->rho_EP);
      EP->wd_hs[i] = bndFac*discharge.q[i]*(EP->rho_EP)/(discharge.a) + 2*(EP->rho_EP)/(discharge.a);
     
    }
 
  i=NUM-1;
  while(0.5*EP->wd_hs[i]>((discharge.a-discharge.r[i])/discharge.a))
    {
      i--;
    }
  if(bndFac == 0)
    i=NUM-1;
  
  EP->edge = i;
 

}

void malloc_EP(EP_ **EP, int NUM)
{
  *EP = malloc(sizeof(EP_));
  (*EP)->orb_wd = malloc(NUM*sizeof(double));
  (*EP)->Db = malloc(NUM*sizeof(double));
  (*EP)->w_dia = malloc(NUM*sizeof(double));
  (*EP)->wd_hs = malloc(NUM*sizeof(double));
}

/****************************************************************************/
/*  Setting up the various components for finding the pressure profiles of  */ 
/*  EP at marginal stability                                                */
/****************************************************************************/


/****************************************************************************/
/*  Setting up the normalizations points                                    */
/****************************************************************************/
void set_point(point_ *point, double r, double gamma_b, double gamma_e, double gamma_i, double gamma_iLT, double gamma_rad, double gamma_extra)
{
  point->r_a = r;
  point->gamma_b = gamma_b;
  point->gamma_e = gamma_e;
  point->gamma_i = gamma_i;
  point->gamma_iLT = gamma_iLT;
  point->gamma_rad = gamma_rad;
  point->gamma_extra = gamma_extra;
}


void malloc_point(point_ ** point)
{
  *point = malloc(sizeof(point_));
}

void malloc_point_array(point_ ** point,int N)
{
  *point = calloc(sizeof(point_),N);
}

/****************************************************************************/
/*       Normalizing the rates according to NOVA computed rates using       */
/*                            an array of points.                            */
/****************************************************************************/
void NormalizeArray(TAE_ *TAE, discharge_ discharge, int Np, point_ *point,int NUM)
{
  int i,p;
  //  int Npt;
  //  Npt = sizeof(*point)/sizeof(point_);
  int pt[Np];
  double r,rp[Np];
  double gamma_b[NUM], gamma_e[NUM], gamma_i[NUM],gamma_iLT[NUM], gamma_rad[NUM], gamma_extra[NUM];
  
  
  for(p=0;p<Np;p++)
    pt[p] = (int) ((point[p].r_a)*((double) NUM));
 
  

  for(i=0;i<NUM;i++)
    {
      gamma_b[i] = TAE->gamma_b[i];
      gamma_e[i] = TAE->gamma_e[i];
      gamma_i[i] = TAE->gamma_i[i];
      gamma_iLT[i] = TAE->gamma_iLT[i];
      gamma_rad[i] = TAE->gamma_rad[i];
      gamma_extra[i] = TAE->gamma_extra[i];
    }
  printf(" The following numbers are factors for normalization with NOVA at the given radii:\n");
  for(i=0;i<Np;i++)
    printf("%d: r/a= %lf\t gammabFAC = %lf \t gamma_eFAC = %lf \t gamma_iFAC = %lf \t gamma_radFAC = %lf \n", i+1, point[i].r_a,point[i].gamma_b/TAE->gamma_b[pt[i]],-1*point[i].gamma_e/TAE->gamma_e[pt[i]], -1*point[i].gamma_i/TAE->gamma_i[pt[i]], -1*point[i].gamma_rad/TAE->gamma_rad[pt[i]]);



  p=0;

  for(i = 0; i<pt[p]; i++)
    {
      if(point[p].gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b/TAE->gamma_b[pt[p]]);
      if(point[p].gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point[p].gamma_e/TAE->gamma_e[pt[p]]);
      if(point[p].gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point[p].gamma_i/TAE->gamma_i[pt[p]]);
      if(point[p].gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point[p].gamma_iLT/TAE->gamma_iLT[pt[p]]);
       if(point[p].gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point[p].gamma_rad/TAE->gamma_rad[pt[p]]);
       if(point[p].gamma_extra != -1)
	gamma_extra[i] = TAE->gamma_extra[i]*fabs(point[p].gamma_extra/TAE->gamma_extra[pt[p]]);

    }
  while(p<Np-1){
  
    while(i < pt[p+1])
      {
	r = discharge.r[i];
	rp[p] = point[p].r_a*discharge.a;
	rp[p+1] = point[p+1].r_a*discharge.a;
	
      if(point[p].gamma_b != -1)
        gamma_b[i] =  point[p].gamma_b + (r -rp[p])*(point[p+1].gamma_b - point[p].gamma_b)/(rp[p+1]-rp[p]);
      //	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b/TAE->gamma_b[pt[p]] + (r -rp[p])*(point[p+1].gamma_b/TAE->gamma_b[pt[p+1]] -point[p].gamma_b/TAE->gamma_b[pt[p]])/(rp[p+1]-rp[p]));
      //      printf("hello...i=%d, p=%d pt_p=%d r=%g ggi=%g pointp=%g\n",i,p,pt[p],r,gamma_b[i],point[p].gamma_b);

      if(point[p].gamma_e != -1)
        gamma_e[i] =  point[p].gamma_e + (r -rp[p])*(point[p+1].gamma_e - point[p].gamma_e)/(rp[p+1]-rp[p]);
      //       gamma_e[i] =  TAE->gamma_e[i]*(fabs(point[p].gamma_e/TAE->gamma_e[pt[p]]) + (r -rp[p])*(fabs(point[p+1].gamma_e/TAE->gamma_e[pt[p+1]]) - fabs(point[p].gamma_e/TAE->gamma_e[pt[p]]))/(rp[p+1]-rp[p]));
    
//nng ion Landau damping 

//      printf("hello...i=%d, p=%d r=%g ggi=%g\n",i,p,r,TAE->gamma_i[i]);
      if(point[p].gamma_i != -1)
        gamma_i[i] =  point[p].gamma_i + (r -rp[p])*(point[p+1].gamma_i - point[p].gamma_i)/(rp[p+1]-rp[p]);

       if(point[p].gamma_iLT != -1)
	 gamma_iLT[i] =  point[p].gamma_iLT + (r -rp[p])*(point[p+1].gamma_iLT - point[p].gamma_iLT)/(rp[p+1]-rp[p]);
	 //	 gamma_iLT[i] = TAE->gamma_iLT[i]*(fabs(point[p].gamma_iLT/TAE->gamma_iLT[pt[p]]) + (r -rp[p])*(fabs(point[p+1].gamma_iLT/TAE->gamma_iLT[pt[p+1]]) - fabs(point[p].gamma_iLT/TAE->gamma_iLT[pt[p]]))/(rp[p+1]-rp[p]));


       if(point[p].gamma_rad != -1)
        gamma_rad[i] =  point[p].gamma_rad + (r -rp[p])*(point[p+1].gamma_rad - point[p].gamma_rad)/(rp[p+1]-rp[p]);
	 //	 gamma_rad[i] = TAE->gamma_rad[i]*(fabs(point[p].gamma_rad/TAE->gamma_rad[pt[p]]) + (r -rp[p])*(fabs(point[p+1].gamma_rad/TAE->gamma_rad[pt[p+1]]) - fabs(point[p].gamma_rad/TAE->gamma_rad[pt[p]]))/(rp[p+1]-rp[p]));

       if(point[p].gamma_extra != -1)
        gamma_extra[i] =  point[p].gamma_extra + (r -rp[p])*(point[p+1].gamma_extra - point[p].gamma_extra)/(rp[p+1]-rp[p]);
	 //	 gamma_extra[i] = TAE->gamma_extra[i]*(fabs(point[p].gamma_extra/TAE->gamma_extra[pt[p]]) + (r -rp[p])*(fabs(point[p+1].gamma_extra/TAE->gamma_extra[pt[p+1]]) - fabs(point[p].gamma_extra/TAE->gamma_extra[pt[p]]))/(rp[p+1]-rp[p]));


       i++;
    }
    p++;



  }

  for(i = pt[p]; i<NUM; i++)
    { 
      if(point[p].gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b/TAE->gamma_b[pt[p]]);

      if(point[p].gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point[p].gamma_e/TAE->gamma_e[pt[p]]);
      if(point[p].gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point[p].gamma_i/TAE->gamma_i[pt[p]]);
       if(point[p].gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point[p].gamma_iLT/TAE->gamma_iLT[pt[p]]);
       if(point[p].gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point[p].gamma_rad/TAE->gamma_rad[pt[p]]);

       if(point[p].gamma_extra != -1)
	gamma_extra[i] = TAE->gamma_extra[i]*fabs(point[p].gamma_extra/TAE->gamma_extra[pt[p]]);

    }
 

  for(i=0;i<NUM;i++)
    {
      for(p=0;p<Np;p++)    
	if(point[p].gamma_e == -100)
	  gamma_e[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_i == -100)
	gamma_i[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_iLT == -100)
	  gamma_iLT[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_rad == -100)
	  gamma_rad[i] = 0;

      for(p=0;p<Np;p++)    
	if(point[p].gamma_extra == -100)
	  gamma_extra[i] = 0;
// filling back the arrays with the damping/growth rates
      TAE->gamma_b[i] = gamma_b[i];
      TAE->gamma_e[i] = gamma_e[i];
      TAE->gamma_i[i] = gamma_i[i];
      TAE->gamma_iLT[i] = gamma_iLT[i];
      TAE->gamma_rad[i] = gamma_rad[i];
      TAE->gamma_extra[i] = gamma_extra[i];

      TAE->crit[i] =  (TAE->gamma_e[i]+ TAE->gamma_rad[i] + TAE->gamma_i[i] +TAE->gamma_extra[i]  )*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);

    }
 
 

};



/****************************************************************************/
/*       Normalizing the rates according to NOVA computed rates using       */
/*                            an array of points.                            */
/****************************************************************************/
void NormalizeArray_Fac(TAE_ *TAE, discharge_ discharge, int Np, point_ *point,int NUM)
{
  int i,p;
  //  int Npt;
  //  Npt = sizeof(*point)/sizeof(point_);
  int pt[Np];
  double r,rp[Np];
  double gamma_b[NUM], gamma_e[NUM], gamma_i[NUM],gamma_iLT[NUM], gamma_rad[NUM], gamma_extra[NUM];
  
  
  for(p=0;p<Np;p++)
    pt[p] = (int) ((point[p].r_a)*((double) NUM));
 
  

  for(i=0;i<NUM;i++)
    {
      gamma_b[i] = TAE->gamma_b[i];
      gamma_e[i] = TAE->gamma_e[i];
      gamma_i[i] = TAE->gamma_i[i];
      gamma_iLT[i] = TAE->gamma_iLT[i];
      gamma_rad[i] = TAE->gamma_rad[i];
      gamma_extra[i] = TAE->gamma_extra[i];
    }
 
  p=0;

  for(i = 0; i<pt[p]; i++)
    {
      if(point[p].gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b);
      if(point[p].gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point[p].gamma_e);
      if(point[p].gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point[p].gamma_i);
      if(point[p].gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point[p].gamma_iLT);
       if(point[p].gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point[p].gamma_rad);
       if(point[p].gamma_extra != -1)
	gamma_extra[i] = TAE->gamma_extra[i]*fabs(point[p].gamma_extra);

    }
  while(p<Np-1){
  
    while(i < pt[p+1])
      {
	r = discharge.r[i];
	rp[p] = point[p].r_a*discharge.a;
	rp[p+1] = point[p+1].r_a*discharge.a;
	
      if(point[p].gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b + (r -rp[p])*(point[p+1].gamma_b -point[p].gamma_b)/(rp[p+1]-rp[p]));

      if(point[p].gamma_e != -1)
       gamma_e[i] =  TAE->gamma_e[i]*(fabs(point[p].gamma_e) + (r -rp[p])*(fabs(point[p+1].gamma_e) - fabs(point[p].gamma_e))/(rp[p+1]-rp[p]));
    

      if(point[p].gamma_i != -1)
       gamma_i[i] =  TAE->gamma_i[i]*(fabs(point[p].gamma_i) + (r -rp[p])*(fabs(point[p+1].gamma_i) - fabs(point[p].gamma_i))/(rp[p+1]-rp[p]));


       if(point[p].gamma_iLT != -1)
	 gamma_iLT[i] = TAE->gamma_iLT[i]*(fabs(point[p].gamma_iLT) + (r -rp[p])*(fabs(point[p+1].gamma_iLT) - fabs(point[p].gamma_iLT))/(rp[p+1]-rp[p]));


       if(point[p].gamma_rad != -1)
	 gamma_rad[i] = TAE->gamma_rad[i]*(fabs(point[p].gamma_rad) + (r -rp[p])*(fabs(point[p+1].gamma_rad) - fabs(point[p].gamma_rad))/(rp[p+1]-rp[p]));

       if(point[p].gamma_extra != -1)
	 gamma_extra[i] = TAE->gamma_extra[i]*(fabs(point[p].gamma_extra) + (r -rp[p])*(fabs(point[p+1].gamma_extra) - fabs(point[p].gamma_extra))/(rp[p+1]-rp[p]));


       i++;
    }
    p++;



  }

  for(i = pt[p]; i<NUM; i++)
    { 
      if(point[p].gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point[p].gamma_b);

      if(point[p].gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point[p].gamma_e);
      if(point[p].gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point[p].gamma_i);
       if(point[p].gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point[p].gamma_iLT);
       if(point[p].gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point[p].gamma_rad);

       if(point[p].gamma_extra != -1)
	gamma_extra[i] = TAE->gamma_extra[i]*fabs(point[p].gamma_extra);

    }
 

  for(i=0;i<NUM;i++)
    {
      for(p=0;p<Np;p++)    
	if(point[p].gamma_e == -100)
	  gamma_e[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_i == -100)
	gamma_i[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_iLT == -100)
	  gamma_iLT[i] = 0;
      for(p=0;p<Np;p++)    
	if(point[p].gamma_rad == -100)
	  gamma_rad[i] = 0;

      for(p=0;p<Np;p++)    
	if(point[p].gamma_extra == -100)
	  gamma_extra[i] = 0;

      TAE->gamma_b[i] = gamma_b[i];
      TAE->gamma_e[i] = gamma_e[i];
      TAE->gamma_i[i] = gamma_i[i];
      TAE->gamma_iLT[i] = gamma_iLT[i];
      TAE->gamma_rad[i] = gamma_rad[i];
      TAE->gamma_extra[i] = gamma_extra[i];

      TAE->crit[i] =  (TAE->gamma_e[i]+ TAE->gamma_rad[i] + TAE->gamma_i[i] +TAE->gamma_extra[i]  )*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);

    }
 
 

};


/****************************************************************************/
/* Normalizing the rates according to NOVA computed rates.                  */
/****************************************************************************/
void Normalize(TAE_ * TAE, discharge_ discharge, point_ point1, point_ point2,int NUM)
{
  int i, pt1, pt2;
  double r,r1,r2;
  double gamma_b[NUM], gamma_e[NUM], gamma_i[NUM],gamma_iLT[NUM], gamma_rad[NUM];

  pt1 = (int) ((point1.r_a)*((double) NUM));
  pt2 = (int) ((point2.r_a)*((double) NUM));  

  for(i=0;i<NUM;i++)
    {
      gamma_b[i] = TAE->gamma_b[i];
      gamma_e[i] = TAE->gamma_e[i];
      gamma_i[i] = TAE->gamma_i[i];
      gamma_iLT[i] = TAE->gamma_iLT[i];
      gamma_rad[i] = TAE->gamma_rad[i];
    }

  printf("gammabFAC = %lf \n gamma_eFAC = %lf \n gamma_radFAC = %lf \n  ",point1.gamma_b/TAE->gamma_b[pt1], point1.gamma_e/TAE->gamma_e[pt1], point1.gamma_rad/TAE->gamma_rad[pt1]);

  printf("gammabFAC2 = %lf \n gamma_eFAC2 = %lf \n gamma_radFAC2 = %lf \n  ",point1.gamma_b/TAE->gamma_b[pt2], point1.gamma_e/TAE->gamma_e[pt2], point1.gamma_rad/TAE->gamma_rad[pt2]);

  for(i = 0; i<pt1; i++)
    {
      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point1.gamma_b/TAE->gamma_b[pt1]);
      if(point1.gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point1.gamma_e/TAE->gamma_e[pt1]);
      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point1.gamma_i/TAE->gamma_i[pt1]);
      if(point1.gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point1.gamma_iLT/TAE->gamma_iLT[pt1]);
       if(point1.gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point1.gamma_rad/TAE->gamma_rad[pt1]);

    }

  for(i = pt1;i<pt2;i++)
    {
      r = discharge.r[i];
      r1 = point1.r_a*discharge.a;
      r2 = point2.r_a*discharge.a;

      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point1.gamma_b/TAE->gamma_b[pt1] + (r -r1)*(point2.gamma_b/TAE->gamma_b[pt2] -point1.gamma_b/TAE->gamma_b[pt1])/(r2-r1));

      if(point1.gamma_e != -1)
       gamma_e[i] =  TAE->gamma_e[i]*(fabs(point1.gamma_e/TAE->gamma_e[pt1]) + (r -r1)*(fabs(point2.gamma_e/TAE->gamma_e[pt2]) - fabs(point1.gamma_e/TAE->gamma_e[pt1]))/(r2-r1));
    

      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*(fabs(point1.gamma_i/TAE->gamma_i[pt1]) + (r -r1)*(fabs(point2.gamma_i/TAE->gamma_i[pt2]) - fabs(point1.gamma_i/TAE->gamma_i[pt1]))/(r2-r1));


       if(point1.gamma_iLT != -1)
	 gamma_iLT[i] = TAE->gamma_iLT[i]*(fabs(point1.gamma_iLT/TAE->gamma_iLT[pt1]) + (r -r1)*(fabs(point2.gamma_iLT/TAE->gamma_iLT[pt2]) - fabs(point1.gamma_iLT/TAE->gamma_iLT[pt1]))/(r2-r1));


       if(point1.gamma_rad != -1)
	 gamma_rad[i] = TAE->gamma_rad[i]*(fabs(point1.gamma_rad/TAE->gamma_rad[pt1]) + (r -r1)*(fabs(point2.gamma_rad/TAE->gamma_rad[pt2]) - fabs(point1.gamma_rad/TAE->gamma_rad[pt1]))/(r2-r1));


    }

  for(i = pt2; i<NUM; i++)
    { 
      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point2.gamma_b/TAE->gamma_b[pt2]);

      if(point1.gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point2.gamma_e/TAE->gamma_e[pt2]);
      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point2.gamma_i/TAE->gamma_i[pt2]);
       if(point1.gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point2.gamma_iLT/TAE->gamma_iLT[pt2]);
       if(point1.gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point2.gamma_rad/TAE->gamma_rad[pt2]);
    }

  for(i=0;i<NUM;i++)
    {

      if(point1.gamma_e == -100)
	gamma_e[i] = 0;

      if(point1.gamma_i == -100)
	gamma_i[i] = 0;

      if(point1.gamma_iLT == -100)
	gamma_iLT[i] = 0;

      if(point1.gamma_rad == -100)
	gamma_rad[i] = 0;

      TAE->gamma_b[i] = gamma_b[i];
      TAE->gamma_e[i] = gamma_e[i];
      TAE->gamma_i[i] = gamma_i[i];
      TAE->gamma_iLT[i] = gamma_iLT[i];
      TAE->gamma_rad[i] = gamma_rad[i];

      TAE->crit[i] =  (TAE->gamma_e[i]+ TAE->gamma_rad[i] + TAE->gamma_i[i]  )*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);

    
    }

}

/****************************************************************************/
/* Normalizing the rates according to a given normalization factor          */
/****************************************************************************/
void Normalize_Fac(TAE_ *TAE, discharge_ discharge, point_ point1, point_ point2,int NUM)
{
  int i, pt1, pt2;
  double r,r1,r2;
  double gamma_b[NUM], gamma_e[NUM], gamma_i[NUM],gamma_iLT[NUM], gamma_rad[NUM];

  pt1 = (int) ((point1.r_a)*((double) NUM));
  pt2 = (int) ((point2.r_a)*((double) NUM));
  

  for(i=0;i<NUM;i++)
    {
      gamma_b[i] = TAE->gamma_b[i];
      gamma_e[i] = TAE->gamma_e[i];
      gamma_i[i] = TAE->gamma_i[i];
      gamma_iLT[i] = TAE->gamma_iLT[i];
      gamma_rad[i] = TAE->gamma_rad[i];
    }
  for(i = 0; i<pt1; i++)
    {
      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point1.gamma_b);
      if(point1.gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point1.gamma_e);
      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point1.gamma_i);
      if(point1.gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point1.gamma_iLT);
       if(point1.gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point1.gamma_rad);

    }

  for(i = pt1;i<pt2;i++)
    {
      r = discharge.r[i];
      r1 = point1.r_a*discharge.a;
      r2 = point2.r_a*discharge.a;

      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point1.gamma_b + (r -r1)*(point2.gamma_b -point1.gamma_b))/(r2-r1);

      if(point1.gamma_e != -1)
	gamma_e[i] =  TAE->gamma_e[i]*(fabs(point1.gamma_e) + (r -r1)*(fabs(point2.gamma_e) - fabs(point1.gamma_e))/(r2-r1));
    

      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*(fabs(point1.gamma_i) + (r -r1)*(fabs(point2.gamma_i) - fabs(point1.gamma_i))/(r2-r1));


       if(point1.gamma_iLT != -1)
	 gamma_iLT[i] = TAE->gamma_iLT[i]*(fabs(point1.gamma_iLT) + (r -r1)*(fabs(point2.gamma_iLT) - fabs(point1.gamma_iLT))/(r2-r1));


       if(point1.gamma_rad != -1)
	 gamma_rad[i] = TAE->gamma_rad[i]*(fabs(point1.gamma_rad) + (r -r1)*(fabs(point2.gamma_rad) - fabs(point1.gamma_rad))/(r2-r1));


    }
  for(i = pt2; i<NUM; i++)
    { 
      if(point1.gamma_b != -1)
	gamma_b[i] = TAE->gamma_b[i]*(point2.gamma_b);

      if(point1.gamma_e != -1)
	gamma_e[i] = TAE->gamma_e[i]*fabs(point2.gamma_e);
      if(point1.gamma_i != -1)
	gamma_i[i] = TAE->gamma_i[i]*fabs(point2.gamma_i);
       if(point1.gamma_iLT != -1)
	gamma_iLT[i] = TAE->gamma_iLT[i]*fabs(point2.gamma_iLT);
       if(point1.gamma_rad != -1)
	gamma_rad[i] = TAE->gamma_rad[i]*fabs(point2.gamma_rad);
    }

  for(i=0;i<NUM;i++)
    {

      if(point1.gamma_e == -100)
	gamma_e[i] = 0;

      if(point1.gamma_i == -100)
	gamma_i[i] = 0;

      if(point1.gamma_iLT == -100)
	gamma_iLT[i] = 0;

      if(point1.gamma_rad == -100)
	gamma_rad[i] = 0;

      TAE->gamma_b[i] = gamma_b[i];
      TAE->gamma_e[i] = gamma_e[i];
      TAE->gamma_i[i] = gamma_i[i];
      TAE->gamma_iLT[i] = gamma_iLT[i];
      TAE->gamma_rad[i] = gamma_rad[i];

      TAE->crit[i] =  (TAE->gamma_e[i]+ TAE->gamma_rad[i] + TAE->gamma_i[i]  )*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);

    
    }

}


/****************************************************************************/
/* Computes the relaxed pressure profiles given the critical gradient       */
/****************************************************************************/
int compute_relaxed( relaxed_ * relaxed, discharge_ discharge, TAE_ TAE, EP_ EP, int NUM)
{

  int i,edge;
  double rlxd[NUM], crit[NUM], beta[NUM], del_beta[NUM],r[NUM];


  int Conserved_Particles = 1 ;
  int k = 0 ;
  int bndry1 = 0; 
  int bndry2;  // bndry2ini is not used!? why!?!
  int bndryLOSS;
 


  int continous = 1;
  int var;
  float Prtcl_Rel = 0;
  float Prtcl_Beta = 0;
  float Prtcl_Diff_Ini;


  for(i =0;i<NUM;i++)
    {
      crit[i] = TAE.crit[i];
      beta[i] = discharge.beta_beam[i];
      del_beta[i] = discharge.del_beta_beam[i];
      r[i] = discharge.r[i];
      rlxd[i] = 0;
    }

  edge = EP.edge;
   
  bndryLOSS=  edge;


  /***********************************************************************/
  /* dependign on region of instability, the beta profile either RELAXED */
  /* without any losses, relaxes but there is LOSS of particles to teh   */
  /* edge,  or is STABLE everywhere already.                             */
  /***********************************************************************/
  enum Statetype {STABLE,LOSS,RELAXED} state;
  
  state = RELAXED;
  /************************************************************************/
  /* bndry1 and bndry2 are the boundaries of the unstable region which    */
  /* relaxes quasilinearly, dependign on their existance intially and as  */
  /* the relaxation proceeds, the profile will either be stable, or relax */
  /* quasilearly. Dependign on whehter bndry2 reaches in the loss region  */
  /* or not, there will be loss of fast ions or not. The following two    */
  /*loops specify the boundaries intially.                                */
  /************************************************************************/

  while((crit[bndry1]) > -1*(del_beta[bndry1]))
    {
      bndry1++;
    
      
      if(bndry1 > bndryLOSS)
 	{
 	  state = STABLE;
 	  break;
 	}    
    }
 
  if(state == STABLE)
    {
    
      for(var=0;var<NUM;var++)
	rlxd[var] = beta[var];


  for(i=0;i<NUM;i++)
    relaxed->profile[i] = rlxd[i];      

      return 1;
    }
  bndry2 = bndry1+1;

 
  if(state != STABLE)
    {
      bndry2 = bndryLOSS;
      while(fabs(crit[bndry2]) > fabs(del_beta[bndry2]))
	bndry2--;

      bndry1 = bndry2 -1;
          
      while(fabs(crit[bndry1]) < fabs(del_beta[bndry1]))
	bndry1--;
      
    }
  /***********************************************************************/
  /***********************************************************************/
  /*              double checking on the biggest unstable region.        */
  /***********************************************************************/
  /***********************************************************************/
  int bndry1t,bndry2t;
  
  bndry2t = bndry1-1; 
  while((fabs(crit[bndry2t]) > fabs(del_beta[bndry2t])) & bndry2t>2)
	bndry2t--;
      
  bndry1t = bndry2t - 1;
   
  while((fabs(crit[bndry1t]) < fabs(del_beta[bndry1t])) & bndry1t>1)
    bndry1t--;
  if((bndry2t-bndry1t)>(bndry2-bndry1))
    {
      bndry2 = bndry2t;
      bndry1 = bndry1t;

    }
  if(bndry2 == bndryLOSS)
    {
      state = LOSS;
    }
  
   
  //  printf("bndry1 is %d\t bndry2 is%d\n",bndry1,bndry2);
  /***********************************************************************/
  /***********************************************************************/
 
  /************************************************************************/
  /* if region is considered loss then integration from bndryLoss to r=0  */
  /************************************************************************/
  if(state == LOSS)
    { 
      var = bndryLOSS;

      while( (var>0))
	{
	  rlxd[var] = rlxd[var+1] + fabs(crit[var])*(r[var] - r[var-1]);
	  
	  if(rlxd[var]>beta[var])
	    {
	      var++;
	      break;
	    }
	  
	  var--;
	} 

      while(var>0)
	{
	  rlxd[var]=beta[var];
	  var--;
	}

    }
  /************************************************************************/
  /* if the region of instability lies completely within the reactor, the */
  /* profile relaxes pushing the boundaries of relaxed region. If in the  */
  /*  process bndry2 is pushed to loss region, state is rendered LOSS.    */
  /************************************************************************/
  if(state == RELAXED)
    {
  
      /************************/
 
      /************************/


      while(Conserved_Particles !=0 & state == RELAXED)
 	{
	  continous = 1;
       	  
 	  bndry2 = bndry2 + 1;
	  
	  if(bndry2 > bndryLOSS)
	    {
	      state = LOSS;
	      break;
	    }

	  for(var = bndryLOSS; var> bndry2;var--)
	    rlxd[var] = beta[var];

	 	  
	  while(var>(bndry1))
	    {
	      rlxd[var] = rlxd[var+1] + fabs(crit[var])*(r[var] - r[var-1]);
	      var--;
	     
	    }

	  while( (var>0))
	    {
	      rlxd[var] = rlxd[var+1] + fabs(crit[var])*(r[var] - r[var-1]);
	   
	      if(rlxd[var]>beta[var])
		break;

	      var--;
	    } 		 
 
	  var++;
	  bndry1 = var;
	 
	  while(var>0)
	    {
	      rlxd[var] = beta[var];
	      var--;
	    }
	  
	  if(k == 0)
	    {
	      for(var = bndry1;var<bndry2;var++)
		{
		  Prtcl_Rel = Prtcl_Rel + rlxd[var]*r[var]*(r[var]-r[var-1]);
		  Prtcl_Beta = Prtcl_Beta +  beta[var]*r[var]*(r[var]-r[var-1]);
		  var++;
		}
	      Prtcl_Diff_Ini = Prtcl_Rel - Prtcl_Beta;
	    }

 
	  /* code to compute losses and compare it to intial */
	
        if(k > 0)
	    {
	   
	      Prtcl_Rel = 0;
	      Prtcl_Beta = 0;
	      
	      for(var = bndry1;var<bndry2;var++)
		{
		  Prtcl_Rel = Prtcl_Rel + rlxd[var]*r[var]*(r[var]-r[var-1]);
		  Prtcl_Beta = Prtcl_Beta + beta[var]*r[var]*(r[var]-r[var-1]);
		  var++;
		}
	      
	      if(  ((Prtcl_Rel - Prtcl_Beta)*Prtcl_Diff_Ini) < 0)
		Conserved_Particles = 0; 
	      	    
	    }
	k++;
	}
    }      
 
  for(i=0;i<NUM;i++)
    relaxed->profile[i] = 0.25*rlxd[i]+0.75*beta[i];
  //rlxd[i]; 

  printf("Test 4 Valpha %g %g \n",TAE.vA[1],EP.v_EP_0);
  //  vA[1],alpha.v_EP_0 EP.v_EP_0
  return 1;
}


/****************************************************************************/
/*              Computes the relaxed profiles of the neutrons               */
/****************************************************************************/
void compute_neut_rlxd( relaxed_ *relaxed,discharge_ discharge,  int NUM)
{
  int i;
  for( i= 0; i<NUM;i++)
   relaxed->neutron_profile[i] = discharge.neut_ini[i]*relaxed->profile[i]/discharge.beta_beam[i];


  
}


void malloc_relaxed(relaxed_ **relaxed, int NUM)
{ 
  *relaxed = malloc(sizeof(relaxed_));
  (*relaxed)->profile = malloc(NUM*sizeof(double));
  (*relaxed)->neutron_profile = malloc(NUM*sizeof(double));
}


/****************************************************************************/
/*              Code to compute the value of the loss                       */
/****************************************************************************/
void compute_loss( relaxed_ * relaxed,discharge_ discharge, int NUM)
{
  int i;
  double r[NUM], relxd[NUM],beta_bm[NUM], neut_rlxd[NUM], neut[NUM];
  double loss;

  double ini=0;
  double fin = 0;
  for(i = 0; i<NUM; i++)
    {
      r[i] = discharge.r[i];
      relxd[i] = relaxed->profile[i];
      neut_rlxd[i] = relaxed->neutron_profile[i];
      beta_bm[i] = discharge.beta_beam[i]; 
      neut[i] = discharge.neut_ini[i]; 
    }

  for (i = 0; i<NUM-1; i++)
    { 
      ini = ini+ beta_bm[i]*r[i]*(r[i+1]-r[i]);
      fin = fin+relxd[i]*r[i]*(r[i+1]-r[i]);
    }
  loss = (ini-fin)/ini*100;

  relaxed->EP_loss = loss;
 
for (i = 0; i<NUM-1; i++)
    { 
      ini = ini+ neut[i]*r[i]*(r[i+1]-r[i]);
      fin = fin+neut_rlxd[i]*r[i]*(r[i+1]-r[i]);
    }
  loss = (ini-fin)/ini*100;

  relaxed->neutron_loss = loss;

}
/****************************************************************************/
/*                       Smoothening functions                              */
/****************************************************************************/

void heaviside_smoothen(relaxed_ *relaxed, EP_ EP, int NUM)
{
  double integN, integP,norm,wd;
  int i,j;
  double funcN[NUM], funcEP[NUM], width_frac[NUM];
  int edge = EP.edge;

  for(i = 0;i<NUM;i++)
    {
      funcEP[i] = relaxed->profile[i];
      funcN[i] = relaxed->neutron_profile[i];
      width_frac[i] = EP.wd_hs[i];
    }
  for (i = edge;i>-1; i--)
    {
      wd = width_frac[i]*NUM;
      integN = 0;
      integP = 0;
      norm = 0;
      for(j = (int) -wd/2;j<(int) wd/2;j++)
	{
	  if((i+j)>0 && (i+j)<NUM)
	    {
	      integN = integN + funcN[i+j];
	      integP = integP + funcEP[i+j];
	      norm = norm + 1;
	    }
        }
      relaxed->neutron_profile[i] = integN/norm;
      relaxed->profile[i] = integP/norm;
    }

}


void FLR_effect(relaxed_* relaxed,EP_ EP, discharge_ discharge, int NUM)
{
  int i,j;
  double rlxd[NUM];
  double a = discharge.a;

  double width  = 4*EP.rho_EP*((double) NUM)/a;
  int len = (int) width;
  
  for(i = 0;i<NUM - len; i++)
    {
      for(j = 0; j<width; j++)
	{
	  rlxd[i] = rlxd[i]+rlxd[i+j];
	}
      
      rlxd[i] = rlxd[i]/(width);
    }
  for (i = NUM-len; i<NUM; i++)
    { 
      
      for(j = 0; j<(width - (i - NUM+len )); j++)
	{
	  rlxd[i] = rlxd[i]+rlxd[i+j];
	}
      
      relaxed->profile[i] = rlxd[i]/(width);
    }

}



  /**************************************************************/
  /*       Gaussian averaging the rates over mode width         */
  /**************************************************************/
void smoothen_gauss_rates(TAE_ * TAE, discharge_ discharge, int NUM)
{

  double integb,intege,integi,integiLT,integrad,norm,pwr,wd;
  
  double *func_temp, *gamma_b, *gamma_e, *gamma_rad, *gamma_i, *gamma_iLT,*r;

  gamma_b = calloc(NUM,sizeof(double));
  gamma_e = calloc(NUM,sizeof(double));
  gamma_i = calloc(NUM,sizeof(double));
  gamma_iLT = calloc(NUM,sizeof(double));
  gamma_rad = calloc(NUM,sizeof(double));
  r = calloc(NUM,sizeof(double));
  func_temp = calloc(NUM,sizeof(double));

  int i,j;
  for (i = 0;i<NUM; i++)
    {
      r[i] = discharge.r[i];
      wd = TAE->sgm_gauss[i];
   
      integb = 0; intege=0; integi = 0; integiLT = 0; integrad = 0;
      norm = 0;

      if(wd > 5*(r[2]-r[1]))
	{
	  for(j = 0;j<NUM;j++)
	    {
	    
	      pwr = -0.5*(r[i] - r[j])*(r[i] - r[j])/(wd*wd);
	   
	      integb = integb + TAE->gamma_b[j]*exp(pwr);
	      intege = intege +TAE->gamma_e[j]*exp(pwr);
	      integi = integi +TAE->gamma_i[j]*exp(pwr);
	      integiLT = integiLT +TAE->gamma_iLT[j]*exp(pwr);
	      integrad = integrad +TAE->gamma_rad[j]*exp(pwr);
	      norm = norm + exp(pwr);
	      
	    }
	}
      else
	{
	  gamma_b[i] = TAE->gamma_b[i];
	  gamma_e[i] = TAE->gamma_e[i];
	  gamma_i[i] = TAE->gamma_i[i]; 
	  gamma_iLT[i] = TAE->gamma_iLT[i];
	  gamma_rad[i] = TAE->gamma_rad[i];
	}
      gamma_b[i] = integb/norm;
      gamma_e[i] = intege/norm;
      gamma_i[i] = integi/norm;
      gamma_iLT[i] = integiLT/norm;
      gamma_rad[i] = integrad/norm;
    }
  
  for(i=0;i<NUM;i++)
    {
      TAE->gamma_b[i] = gamma_b[i];
      TAE->gamma_e[i] = gamma_e[i];
      TAE->gamma_i[i] = gamma_i[i];
      TAE->gamma_iLT[i] = gamma_iLT[i];
      TAE->gamma_rad[i] = gamma_rad[i];
      TAE->crit[i] = (TAE->gamma_e[i]+ TAE->gamma_rad[i])*(discharge.del_beta_beam[i])/(TAE->gamma_b[i]);
    }

  free(func_temp);
  free(gamma_b);
  free(gamma_e);
  free(gamma_rad);
  free(gamma_i);
  free(gamma_iLT);
  free(r);


}

/****************************************************************************/
/*              Functions to free the allocated memories                    */
/****************************************************************************/
void free_discharge(discharge_ *discharge)
{

  free(discharge->T_r);
  free(discharge->Te_r);
  free(discharge->beta_pc_r);
  free(discharge->beta_elec);
  free(discharge->beta_beam);
  free(discharge->del_beta_beam);
  free(discharge->neut_ini);
  free(discharge->q);
  free(discharge->del_q);
  free(discharge->shear);
  free(discharge->ni);
  free(discharge->r);
  free(discharge);

  
}

void free_TAE(TAE_ *TAE)
{

  free(TAE->vA);
  free(TAE->Di);
  free(TAE->Do);
  free(TAE->wtae);
  free(TAE->crit);
  free(TAE->gamma_b);
  free(TAE->gamma_e);
  free(TAE->gamma_rad);
  free(TAE->gamma_i);
  free(TAE->gamma_iLT);
  free(TAE->sgm_gauss);
  free(TAE->gamma_EP_prm);  
  free(TAE);
}

void free_EP(EP_ *EP)
{

  free(EP->Db);
  free(EP->orb_wd);
  free(EP->w_dia);
  free(EP->wd_hs);
  free(EP);
}

void free_relaxed(relaxed_ *relaxed)
{
  free(relaxed->neutron_profile);
  free(relaxed->profile);
  free(relaxed);
}

void free_point(point_ *point)
{
  free(point);
}

void free_visualization(visualization_ * visualization)
{
  free(visualization);
}

/****************************************************************************/
/* functions from NR (numerical recipes) for interpolation                  */
/****************************************************************************/
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u=vector(1,n-1);
  if (yp1 > 0.99e30) 
    y2[1]=u[1]=0.0;
  else
    {
      y2[1] = -0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
  
  for (i=2;i<=n-1;i++) 
    { 
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
  if (ypn > 0.99e30) 
    qn=un=0.0; 
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--) 
    y2[k]=y2[k]*y2[k+1]+u[k]; 
  free_vector(u,1,n-1);
}


void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
  //	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
	  
//////////////////////////////  vectors ////////////////////////////////
	  
double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) printf("allocation failure in vector()");
	return v-nl+NR_END;
}



void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}




/****************************************************************************/
/* get_data_NAME returns a double of the value of what follows name_data if */
/* it exists in file someFile,                                              */
/****************************************************************************/
int get_Ntransp(FILE *someFile)
{
  
  char Ntc[2];
  
  char a,b,c,d;
 
  int Ntransp;


  a = fgetc(someFile);
  b = fgetc(someFile);
  c = fgetc(someFile);
  d = fgetc(someFile);
 

  while(d != EOF)
    {      
      if(c == ' ' && d == 'X')
	break;
     
      a = b;
      b = c;
      c = d;
      d = fgetc(someFile);
    }

  Ntc[0] = a;
  Ntc[1] = b;
  Ntransp =  (int) atof(Ntc);

  return Ntransp; 

}
/****************************************************************************/
/* get_data_NAME returns a double of the value of what follows name_data if */
/* it exists in file someFile,                                              */
/****************************************************************************/
double get_data_NAME(FILE *someFile, const char * nameData)
{
  double vlu;
  int dataExists;
  char dat;
  int dg,i;
  char data[20]; // = malloc(20*sizeof(char));
 
  dataExists = 0;
  dat = fgetc(someFile);
 
  while(dat != EOF && dataExists == 0)
    {
      i=0;     
      dat = fgetc(someFile);
      
      if(dat == nameData[i])
	{
	  while(dat == nameData[i])
	    {	    
	      dat = fgetc(someFile);
	      i++;
	    }
	  
	  if(nameData[i] == '\0')
	    { 
	      dataExists = 1;		  
	    }
	}
    }
  
  dg =0;  
  dat = fgetc(someFile);

  while((dat == '\t') || (dat == ' '))
    {
      dat = fgetc(someFile);		 
    }
  while((dat != '\t') && (dat != ' ') && (dat != EOF) && (dat != '\0')  )
    {
      data[dg] = dat;
      dat = fgetc(someFile);		 
      dg++;
    }

  vlu =(double) atof(data);

  if(dataExists ==0){
    printf("\n You goofed! Make sure the correct spelling of %s exists in the INPUT file. \n \n",nameData);
    abort();
  }
 
  return vlu;
 
    
}
/****************************************************************************/
/* get_data_NAME returns a double of the value of what follows name_data if */
/* it exists in file someFile,                                              */
/****************************************************************************/
int get_data_NAME_array(double *vlu, int Npts,  FILE *someFile, const char * nameData)
{
 
  int dataExists;
  char dat;
  int dg,i,j,k;
  char data[20]; // = malloc(20*sizeof(char));
  char cmp[50];
 
  dataExists = 0;
  dat = fgetc(someFile);
 
  while(dat != EOF && dataExists == 0)
    {
      i=0;     
      dat = fgetc(someFile);
      
      if(dat == nameData[i])
	{
	  while(dat == nameData[i])
	    {	    
	      cmp[i]=dat;
	      dat = fgetc(someFile);
	      i++;
	    }
	  
	  if(nameData[i] == '\0')
	    {  
	     
	      dataExists = 1;		  
	    }
	}
    }
  
  k=0;
  for(j=0;j<50;j++)
    {
      if(cmp[j]!='\0')
	k++;
    }

  if(strncmp(cmp,nameData,k)!=0){
   
     printf("\n You goofed! Make sure the correct spelling of %s exists in the INPUT file. \n \n",nameData);
    abort();
  }
 
 
  
  dg =0;  
  dat = fgetc(someFile);

  while((dat == '\t') || (dat == ' '))
    {
      dat = fgetc(someFile);		 
    }

  for(i=0;i<Npts;i++){
    for(dg=0;dg<20;dg++)
      data[dg]=0;
    dg=0;
    while((dat != '\t'&& dat!='\n') )
      {
	data[dg] = dat;
	dat = fgetc(someFile);		 
	dg++;
      }
    vlu[i] = (double) atof(data);
   
   
    if(dat=='\n' & i!=Npts-1){
      printf("\n You goofed! Make sure you have the %d values for teh rates in the INPUT file. Seperate them by tabs. \n \n",Npts);
    abort();
    }
    dat = fgetc(someFile);
 
  }
  printf("\n");
  

  return 0;
 
    
}
