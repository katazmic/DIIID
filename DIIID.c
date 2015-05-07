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



int main(int argc, char **argv)
{
  
  int Npt;
  int fac;
  discharge_ *diiid211;
  EP_ *beam;
  TAE_ *TAE;
  relaxed_ *relaxed;
  point_ *point1;
  point_ *point2;
  visualization_ *visualization;

  
  point_ *point;

  int i;


  char fname[200];
  int Ntransp;
  int shot, discharge;
  int NUM;
  double sgm = 0.8;
  double mb,Eb,gamm,Z; 
  double boundary;

  FILE *prnt;
  prnt = fopen("Results1_5D.txt","w");


  FILE *pgnu;
  pgnu = fopen("plot.gnu","w");


  char  *INname = "INPUT";
  gamm = 5/3;
 

  FILE *ftransp;
  FILE *fInp;
  double r, gamma_b, gamma_e,  gamma_i,  gamma_iLT,  gamma_rad;
 


  mb = 2;/* ratio of m_H to beam mass*/
  Eb = 0.08; /*in mega electron volts*/
  Z = 1; /*charge of fast ion in units of charge of proton. */
  ////////////////////////////////  NORMALIZATION /////////////////////////

 
  fInp = fopen(INname,"r");   
  Npt =  (int) get_data_NAME(fInp,"Npoints");
  fclose(fInp);

  malloc_point_array(&point,Npt);
  
  read_input_points_array(point, Npt, INname,&shot, &discharge, &NUM, &boundary);
  
  ///////////////////////////////// VIZUALIZATION /////////////////////////

 
  malloc_visualization(&visualization);
  read_input_visualization(visualization, INname);

  /////////////////////////////////////////////////////////////////////////


  if(shot == 360 || shot==780 ||shot==120 || shot == 310 || shot ==410 || shot == 460 || shot == 510 || shot==560)
    sprintf(fname, "122117_transp_t0%d.dat",shot);
  else
    sprintf(fname, "transp_t0%d.dat",shot);

  ftransp = fopen(fname, "r");
 
  // IF SPECIAL FILES DONT EXIST, OPENS XXX_TRANSP_YYY.DAT
  if (ftransp == NULL){ 
    sprintf(fname, "%d_transp_%d.dat",discharge,shot);
    ftransp = fopen(fname, "r");
  }


  // IF THE FILE NAME IS NOT FOUND, IT OPENS TRANSP.DAT
  if (ftransp == NULL){ 
    sprintf(fname, "transp.dat");
    ftransp = fopen(fname, "r");
  }

  // IF FILE STILL NOT FOUND, ABORT!
  if (ftransp == NULL){ 
	printf("\n FILE NOT FOUND! \n --> Check INPUT file and specify the integers for discharge and time that mathces the transp file name is (discharge)_transp_(shot).dat\n Refer to README for details. \n \n");
	  abort();
    
  }
  
  Ntransp = get_Ntransp(ftransp);
  
  malloc_discharge(&diiid211,NUM);
  set_discharge(diiid211, fname,  Ntransp,  sgm, NUM);


  malloc_EP(&beam,NUM);
  set_EP(beam, fname, *diiid211,  mb, Z,  Eb , gamm,boundary, NUM);

 
  malloc_TAE(&TAE,NUM);
  set_TAE(TAE,  *beam, *diiid211, NUM);
 

  // smoothen_gauss_rates(TAE, *diiid211, NUM); (not working)
 

  fInp = fopen(INname,"r");   
  
  fac = get_data_NAME(fInp, "Normalization_Scheme");
  
  fclose(fInp);
  if(fac>2 || fac<1){
    printf("\n Provide correct parameter for normalization scheme (1 or 2) in INPUT file.\n Refer to README for more details. \n\n");
    abort();}

  if(fac == 2 )
    NormalizeArray_Fac(TAE, *diiid211, Npt,point, NUM);

  if(fac==1)
    NormalizeArray(TAE, *diiid211, Npt,point, NUM);


  malloc_relaxed(&relaxed,NUM);
  compute_relaxed(relaxed, *diiid211, *TAE, *beam,  NUM);
  compute_neut_rlxd(relaxed, *diiid211, NUM);

  heaviside_smoothen(relaxed, *beam, NUM);
   
 
  compute_loss(relaxed, *diiid211,NUM);


  /***************************************************************/  
  
  printf("%lf is the fast ion loss \n",relaxed->EP_loss);
  printf("%lf is the neutron loss \n",relaxed->neutron_loss);
  int rt,sf,crta,crtb,prfa,prfb;
  char rateNm[400], prfNm[100], crtNm[100];
  strcpy(rateNm, "");  strcpy(prfNm, "");strcpy(crtNm, "");
  int nm;

  rt =0;
  sf=0;
  crta=0;
  prfa=0;
  crtb=0;
  prfb=0;
      

  for (i = 2; i<NUM; i++)
    {
      fprintf( prnt,"%lf \t",  diiid211->r[i] );


      //////////////// RATES //////////////////////////
      if( visualization->gamma_b ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "growth rate X ");
	}
	fprintf( prnt,"%lf \t", TAE->gamma_b[i]  );}

      if( visualization->gamma_e ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "electron collision damping rate X ");
	}
	fprintf( prnt,"%lf \t", TAE->gamma_e[i]  );}


      if(visualization->gamma_i ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "ion landau damping rate (D) X ");
	}
	fprintf( prnt,"%lf \t",  TAE->gamma_i[i] );}


      if(visualization->gamma_iLT ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "ion landau damping rate (T) X ");
	}
	fprintf( prnt,"%lf \t",  TAE->gamma_iLT[i] );}


      if(visualization->gamma_rad ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "Radiative damping rate X ");
	}
	fprintf( prnt,"%lf \t",  TAE->gamma_rad[i] );}


      if(visualization->gamma_extra ==1){
	if(i==2){
	  rt++;
	  strcat(rateNm, "Extra damping rate X ");
	}
	fprintf( prnt,"%lf \t",  TAE->gamma_extra[i] );}
 

      ///////// /////  CRITICCALL      ////// ////// //////
      if(visualization->beta_beam_crit ==1){
	if(i==2){
	  crtb++;
	  strcat(crtNm, "critical beam beta gradient X ");
	}
	fprintf( prnt,"%lf \t", TAE->crit[i]  );}

      if(visualization->del_beta_beam ==1){
	if(i==2){
	  crtb++;
	  strcat(crtNm, "intial beam beta gradient X ");
	}
	fprintf( prnt,"%lf \t", -1*diiid211->del_beta_beam[i] );}
      
      //////////////////// PROFILES ////////////////////
      if(visualization->beta_beam_ini ==1){
        if(i==2){
	  prfb++;
	  strcat(prfNm, "intial beam beta profile X ");
	}
	fprintf( prnt,"%lf \t", diiid211->beta_beam[i]  );}

      if(visualization->beta_beam_rel ==1){
	if(i==2){
	  prfb++;
	  strcat(prfNm, "relaxed beam beta profile X ");
	}
	fprintf( prnt,"%lf \t", relaxed->profile[i]  );}


      ///////////
      if(visualization->beta_alpha_crit ==1){
	crta++;
      	fprintf( prnt,"%lf \t", TAE->crit[i]  );
	printf("FIX ALPHA FOR PLOTTING!!!! \n");
      }

      if(visualization->del_beta_alpha ==1){
	crta++;
	fprintf( prnt,"%lf \t",  -1*diiid211->del_beta_beam[i] );
	printf("FIX ALPHA FOR PLOTTING!!!! \n");
      } 
      
      ////////
      if(visualization->beta_alpha_ini ==1){
	prfa++;
	fprintf( prnt,"%lf \t",  diiid211->beta_beam[i] );
	printf("FIX ALPHA FOR PLOTTING!!!! \n");
      }

      if(visualization->beta_alpha_rel ==1){
	prfa++;
	fprintf( prnt,"%lf \t", relaxed->profile[i]  );
	printf("FIX ALPHA FOR PLOTTING!!!! \n");
      }

      ////// ////// ////// SAFETY FACRTOR  ////// ////// //////
      if(visualization->q ==1){

	if(i==2){
	  sf++;
	}
	fprintf( prnt,"%lf \t",  diiid211->q[i] );}


      fprintf( prnt,"\n"  );
     
    }

  int wnd=0;
  
  fprintf(pgnu,"set terminal x11 %d \n" ,wnd);
  i=0;
  /////////////////  RATES   /////
  if(rt>0){
    fprintf(pgnu,"set yrange[-0.2:0.2]\n");
    fprintf(pgnu,"plot ");
  }
  while(i<rt-1){
    i++;

    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" ",i+1 );
    while(rateNm[nm]!='X'){
      fprintf(pgnu,"%c", rateNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines,");

  }
  if(rt>0){
    i++;

    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" ",i+1 );
    while(rateNm[nm]!='X'){
      fprintf(pgnu,"%c", rateNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines;");
    
  }
 
  

  //////////////////////  CRITICAL  //////

  nm=0;
  if(i>0&crtb>0){
    wnd++;
    fprintf(pgnu,"\n set terminal x11 %d \n",wnd );
  }
  if(crtb>0){
    fprintf(pgnu,"set yrange[-0.00:0.0005]\n");
    fprintf(pgnu,"plot ");
  }
  i=0;
  while(i<crtb-1){
    i++;
    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" " ,i+1+rt);
    while(crtNm[nm]!='X'){
      fprintf(pgnu,"%c", crtNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines,");


  }
  if(crtb>0){
    i++;
    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" " ,i+1+rt );
    while(crtNm[nm]!='X'){
      fprintf(pgnu,"%c", crtNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines;");

  }
  nm=0;

  ////////////////////  PROFILES  //////////
  if(i>0&prfb>0){
    wnd++;
    fprintf(pgnu,"\n set terminal x11 %d \n",wnd );
 
  }
  if(prfb>0){
    fprintf(pgnu,"set yrange[0:0.015]\n");
    fprintf(pgnu,"plot ");
  }
  i=0;
  while(i<prfb-1){
    i++;
    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" ",i+1+rt+crtb );
    while(prfNm[nm]!='X'){
      fprintf(pgnu,"%c", prfNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines,");
  }
if(prfb>0){
    i++;
    fprintf(pgnu,"\"Results1_5D.txt\" using 1:%d title \" ",i+1+rt+crtb );
    while(prfNm[nm]!='X'){
      fprintf(pgnu,"%c", prfNm[nm]);
      nm++;}
    nm++;
    fprintf(pgnu," \" with lines;");
 }

  /////////////////// SAFETY FACTOR /////////////////
  i=1;
  if(sf>0){
    wnd++;
    fprintf(pgnu,"set terminal x11 %d \n" ,wnd);
    fprintf(pgnu,"set yrange[4:10]\n");
    fprintf(pgnu,"plot \"Results1_5D.txt\" using 1:%d title \" Safety Factor \"  with lines;\n",i+1+rt+crtb+prfb);
  }
  
////////////////////////


  fclose(prnt);
  fclose(pgnu);
  free_point(point);
  
  free_discharge(diiid211);
  free_EP(beam);
  free_TAE(TAE);
  free_relaxed(relaxed);
  //free_point(point1);
  //free_point(point2);
  free_visualization(visualization);
  return 1;
}
