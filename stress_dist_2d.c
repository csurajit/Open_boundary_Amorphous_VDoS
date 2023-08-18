#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>


#define TP (400)
double double_TP=TP*1.e0;
#define dim (2)
//...LJ parameters...
#define sig_AA (1.e0)  
#define eps_AA (1.e0)
#define sig_AB (0.8e0*sig_AA)  
#define eps_AB (1.5e0*eps_AA)  
#define sig_BB (0.88e0*sig_AA)  
#define eps_BB (0.5e0*eps_AA)  
#define r_cut_AA (2.5e0*sig_AA)
#define r_cut_AB (2.5e0*sig_AB)
#define r_cut_BB (2.5e0*sig_BB)
#define c_0 (0.040490237952e0)
#define c_2 (-0.00970155098112e0)
#define c_4 (0.0006201261686784e0)
//double  pi=acos(-1.e0);

int main(){
    int i, j, k, n, m, l, shit;
    double total_sig_xx_0 = 0.e0, total_sig_yy_0 = 0.e0, total_sig_xy_0 = 0.e0,  total_sig_xx = 0.e0, total_sig_yy = 0.e0, total_sig_xy = 0.e0; //total stresses
    double pi = acos(-1.e0);    
    //---------- for force and stress calculation ----
    int tag_i, tag_j;
    double xi,yi,r2i,r6i,r;
    double dx,dy,r2,s2,fr,store_x,store_y,pe, dist_x, dist_y;
    double cmx, cmy;
    double ft;
    //==================== BINNING ===============================
    int no_of_x_bin = 200;
    double x_max, x_min, delta_x;
    x_max = 30.e0;
    x_min = -30.e0;
    delta_x = (x_max-x_min)/(no_of_x_bin*1.e0);
    //printf("%E\n",delta_x);
    
    int no_of_y_bin = 200;
    double y_max, y_min, delta_y;
    y_max = 30.e0;
    y_min = -30.e0;
    delta_y = (y_max-y_min)/(no_of_y_bin*1.e0);
    //printf("%E\n",delta_y);
    
    int no_of_sq_bin = no_of_x_bin*no_of_y_bin;
    
    int no_of_rad_bin = 100;
    double delta_r, dr;
    delta_r = (30.e0)/(no_of_rad_bin*1.e0);
    //printf("%E\n",delta_r);
    
    
    //=== stress in each square grid ===   
    double *sxx_sq_0, *syy_sq_0, *sxy_sq_0;
    int *num_sq;
    double *sxx_rad_0, *syy_rad_0, *sxy_rad_0;
    int *num_rad;
   
    sxx_sq_0 = malloc(no_of_sq_bin * sizeof(*sxx_sq_0));
    syy_sq_0 = malloc(no_of_sq_bin * sizeof(*syy_sq_0));
    sxy_sq_0 = malloc(no_of_sq_bin * sizeof(*sxy_sq_0));
    
    num_sq = malloc(no_of_sq_bin * sizeof(*num_sq));
    //=========================
    for(n=0;n<no_of_sq_bin;n++){
        sxx_sq_0[n]=0.e0;
        syy_sq_0[n]=0.e0;
        sxy_sq_0[n]=0.e0;
            
        num_sq[n]=0;   
    }
    //=================================================== radial-distribution ========================================================================   
    sxx_rad_0 = malloc(no_of_rad_bin * sizeof(*sxx_rad_0));
    syy_rad_0 = malloc(no_of_rad_bin * sizeof(*syy_rad_0));
    sxy_rad_0 = malloc(no_of_rad_bin * sizeof(*sxy_rad_0));
    
    num_rad = malloc(no_of_rad_bin * sizeof(*num_rad));
    //=========================
    for(n=0;n<no_of_rad_bin;n++){
        sxx_rad_0[n]=0.e0;
        syy_rad_0[n]=0.e0;
        sxy_rad_0[n]=0.e0;
    
        num_rad[n]=0;
    }    
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    unsigned int ens_id=0,no_config=100000, config, load_config=0, dummy;
    int *tag_0;
    double *x_0, *y_0, *fx_0, *fy_0, *sig_xx_0, *sig_yy_0, *sig_xy_0;

    load_config=0,ens_id=0;
    double dummy2, dummy3, dummy4;
    for(dummy=100000;dummy<100000+no_config;dummy++){
    	printf("%06u\n",dummy);
    	

        char PositionsFileName[128];
        FILE *positionFile;
        sprintf(PositionsFileName,"../min_configs_%d/%06u_configuration.txt",TP,dummy);
        positionFile=fopen(PositionsFileName,"rb");
        
        
        char eigenvaluesFileName[128];
       	FILE *eigenvaluesFile;
       	sprintf(eigenvaluesFileName,"../eigenvalue_%d/eigenvalues_of_%06u_ens_id.txt",TP,dummy);
       	eigenvaluesFile=fopen(eigenvaluesFileName,"rb");
       	if(positionFile==NULL || eigenvaluesFile==NULL) {
       	        		printf("ERROR No File at :'%s'\n",PositionsFileName);
        }

        else{
        	
        	for (i=0; i < 3; i++){
        		fscanf(eigenvaluesFile, "%*[^\n]\n");
    		}
    		fscanf(eigenvaluesFile,"%lf %lf %lf",&(dummy2),&(dummy3),&(dummy4));
            	fclose(eigenvaluesFile);
        
            x_0= malloc(TP * sizeof(*x_0));
            y_0= malloc(TP * sizeof(*y_0));
            tag_0= malloc(TP * sizeof(*tag_0));
            for(k=0;k<TP;k++){
                fscanf(positionFile,"%d %lf %lf\n",&(tag_0[k]),&(x_0[k]),&(y_0[k]));
            }
            fclose(positionFile);
            
            
            if(dummy2>1e-5){
            cmx=0.e0, cmy=0.e0;
            for(i=0;i<TP;i++){
                cmx+=x_0[i];
                cmy+=y_0[i];
            }
            cmx/=double_TP;
            cmy/=double_TP;
            //printf("%E\t%E\n",cmx,cmy);
            
            fx_0= malloc(TP * sizeof(*fx_0));
            fy_0= malloc(TP * sizeof(*fy_0));
            sig_xx_0= malloc(TP * sizeof(*sig_xx_0));
            sig_yy_0= malloc(TP * sizeof(*sig_yy_0));
            sig_xy_0= malloc(TP * sizeof(*sig_xy_0));
            
            
            for(i=0;i<TP;i++){
                fx_0[i]=0.e0;
                fy_0[i]=0.e0;
                sig_xx_0[i]=0.e0;
                sig_yy_0[i]=0.e0;
                sig_xy_0[i]=0.e0;
            }
            pe = 0.e0;
            
            for(i=0;i<TP-1;i++){
                xi=x_0[i];
                yi=y_0[i];
                tag_i= tag_0[i];
                for(j=i+1;j<TP;j++){
                    tag_j = tag_0[j];
                    dx=x_0[j]-xi;
                    dy=y_0[j]-yi;
                    r2=(dx*dx)+(dy*dy);
                    if(tag_i==1 && tag_j==1 && r2<(r_cut_AA*r_cut_AA)){
                        s2= sig_AA*sig_AA;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_AA*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_AA*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx_0[j] += (store_x);  //
                        fx_0[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy_0[j] += (store_y);  //Final force calculation componantwise
                        fy_0[i] -= (store_y);  //
                        //-----------------------------------
                        sig_xx_0[i]+=store_x*dx;
                        sig_yy_0[i]+=store_y*dy;
                        sig_xy_0[i]+=store_x*dy;
                        sig_xx_0[j]+=store_x*dx;
                        sig_yy_0[j]+=store_y*dy;
                        sig_xy_0[j]+=store_x*dy;
                        
                    }
                    if (tag_i==0 && tag_j==0 && r2<(r_cut_BB*r_cut_BB)){
                        s2= sig_BB*sig_BB;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_BB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_BB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx_0[j] += (store_x);  //
                        fx_0[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy_0[j] += (store_y);  //Final force calculation componantwise
                        fy_0[i] -= (store_y);  //

                        //-----------------------------------
                        sig_xx_0[i]+=store_x*dx;
                        sig_yy_0[i]+=store_y*dy;
                        sig_xy_0[i]+=store_x*dy;
                        sig_xx_0[j]+=store_x*dx;
                        sig_yy_0[j]+=store_y*dy;
                        sig_xy_0[j]+=store_x*dy;
                        
                    }
                    if (tag_i!=tag_j && r2<(r_cut_AB*r_cut_AB)){
                        s2= sig_AB*sig_AB;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_AB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_AB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx_0[j] += (store_x);  //
                        fx_0[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy_0[j] += (store_y);  //Final force calculation componantwise
                        fy_0[i] -= (store_y);  //
                        //-----------------------------------
                        sig_xx_0[i]+=store_x*dx;
                        sig_yy_0[i]+=store_y*dy;
                        sig_xy_0[i]+=store_x*dy;
                        sig_xx_0[j]+=store_x*dx;
                        sig_yy_0[j]+=store_y*dy;
                        sig_xy_0[j]+=store_x*dy;
                        
                    }
                }
            }
            for(i=0;i<TP;i++){ 
                    dx = x_0[i]-cmx;
                    dy = y_0[i]-cmy;
                    dr = sqrt((dx*dx)+(dy*dy));
                    j = (int) ((dx-x_min)/delta_x);
                    k = (int) ((dy-y_min)/delta_y);
                    n=j+(no_of_y_bin*k);
                    m=(int) (dr/delta_r);
                    
                    sxx_sq_0[n]+= sig_xx_0[i];
                    syy_sq_0[n]+= sig_yy_0[i];
                    sxy_sq_0[n]+= sig_xy_0[i];
                    num_sq[n]++;
                    
                    sxx_rad_0[m]+= sig_xx_0[i];
                    syy_rad_0[m]+= sig_yy_0[i];
                    sxy_rad_0[m]+= sig_xy_0[i];
                    
                    num_rad[m]++;
            }
            

            free(fx_0);
            free(fy_0);
            free(sig_xx_0);
            free(sig_xy_0);
            free(sig_yy_0);
            
             ens_id++;
            }
            free(x_0);
            free(y_0);
            free(tag_0);
            
      	}      
            
	}
    printf("%06u\t%E\t%E\t%E\t%d\n",ens_id,delta_x,delta_y,delta_r,no_of_sq_bin);
    for(n=0;n<no_of_sq_bin;n++){
    	//printf("%E\t%E\t%E\n",sxx_sq_0[n],syy_sq_0[n],sxy_sq_0[n]);
        sxx_sq_0[n]/=(delta_x*delta_y*ens_id*1.e0);
        syy_sq_0[n]/=(delta_x*delta_y*ens_id*1.e0);
        sxy_sq_0[n]/=(delta_x*delta_y*ens_id*1.e0);
        
    }
    char sq_stress_FileName[128];
    FILE *sq_stress;
    sprintf(sq_stress_FileName,"%d_sq_stress_of_particles_%.2lf_bin.txt",TP,delta_x);
    sq_stress=fopen(sq_stress_FileName,"w");
    for(i=0;i<no_of_x_bin;i++){
        for(j=0;j<no_of_y_bin;j++){
            m = i+(j*no_of_y_bin);
            if(num_sq[m]!=0) fprintf(sq_stress,"%E\t%E\t%E\t%E\t%E\t%E\t%E\t%d\n",((i+0.5e0)*delta_x)+x_min,((j+0.5e0)*delta_y)+y_min,sxx_sq_0[m],syy_sq_0[m],sxy_sq_0[m],(sxx_sq_0[m]+syy_sq_0[m])/2.e0,(sxx_sq_0[m]-syy_sq_0[m]),num_sq[m]);
        }
    }
    
    
    for(n=0;n<no_of_rad_bin;n++){
        sxx_rad_0[n]/=(pi*(((n+1)*(n+1))-(n*n))*delta_r*delta_r*ens_id*1.e0);
        syy_rad_0[n]/=(pi*(((n+1)*(n+1))-(n*n))*delta_r*delta_r*ens_id*1.e0);
        sxy_rad_0[n]/=(pi*(((n+1)*(n+1))-(n*n))*delta_r*delta_r*ens_id*1.e0);
        
    }

    char rad_stress_FileName[128];
    FILE *rad_stress;
    sprintf(rad_stress_FileName,"%d_rad_stress_of_particles_%.2lf_bin.txt",TP,delta_r);
    rad_stress=fopen(rad_stress_FileName,"w");
    for(m=0;m<no_of_rad_bin;m++){
        if(num_rad[m]!=0) fprintf(rad_stress,"%E\t%E\t%E\t%E\t%E\t%E\t%d\n",((m+0.5e0)*delta_r),sxx_rad_0[m],syy_rad_0[m],sxy_rad_0[m],(sxx_rad_0[m]+syy_rad_0[m])/2.e0,(sxx_rad_0[m]-syy_rad_0[m]),num_rad[m]);
    }
    

    


    fclose(sq_stress);
    fclose(rad_stress);
    
    free(sxx_rad_0);
    free(syy_rad_0);
    free(sxy_rad_0);
    
	free(num_sq);
    free(num_rad);
}
