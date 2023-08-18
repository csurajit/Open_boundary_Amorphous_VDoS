#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>


#define TP (4096)
double double_TP=TP*1.e0;
#define dim (3)
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
    double total_sig_xx = 0.e0, total_sig_yy = 0.e0, total_sig_zz = 0.e0,  total_sig_xy = 0.e0, total_sig_yz = 0.e0, total_sig_xz = 0.e0; //total stresses
    double pi = acos(-1.e0);    
    //---------- for force and stress calculation ----
    int tag_i, tag_j;
    double xi,yi,zi,r2i,r6i,r;
    double dx,dy,dz,r2,s2,fr,store_x,store_y,store_z, pe, dist_x, dist_y, dist_z;
    double cmx, cmy, cmz;
    double ft;
    
    //==================== BINNING ===============================
    int no_of_rad_bin = 100;
    double delta_r, dr;
    delta_r = (28.e0)/(no_of_rad_bin*1.e0);
    //printf("%E\n",delta_r);
    
    //=== stress in each square grid ===
    double *sxx_rad, *syy_rad, *szz_rad;
    int *num_rad;
   
    //=================================================== radial-distribution ========================================================================   
    sxx_rad = malloc(no_of_rad_bin * sizeof(*sxx_rad));
    syy_rad = malloc(no_of_rad_bin * sizeof(*syy_rad));
    szz_rad = malloc(no_of_rad_bin * sizeof(*szz_rad));
    
    num_rad = malloc(no_of_rad_bin * sizeof(*num_rad));
    //=========================
    for(n=0;n<no_of_rad_bin;n++){
        sxx_rad[n]=0.e0;
        syy_rad[n]=0.e0;
        szz_rad[n]=0.e0;
    
        num_rad[n]=0;
    }    
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    unsigned int ens_id=0,no_config=100000, config, load_config=0, dummy;
    int *tag;
    double *x, *y, *z, *fx, *fy, *fz, *sig_xx, *sig_yy, *sig_zz;
	double Tsig_xx, Tsig_yy, Tsig_zz;
    load_config=0,ens_id=0;
    double dummy2, dummy3, dummy4;
    for(dummy=0;dummy<no_config;dummy++){
    	printf("%06u\n",dummy);
    	
        char PositionsFileName[128];
        FILE *positionFile;
        sprintf(PositionsFileName,"first_non_zero_mode_of_%d_particles_%06u_ens_id.txt",TP,dummy);
        positionFile=fopen(PositionsFileName,"rb");
        
        char eigenvaluesFileName[128];
           FILE *eigenvaluesFile;
           sprintf(eigenvaluesFileName,"eigenvalues_of_%d_particles_%06u_ens_id.txt",TP,dummy);
           eigenvaluesFile=fopen(eigenvaluesFileName,"rb");
           if(positionFile==NULL || eigenvaluesFile==NULL) {
                           printf("ERROR No File at :'%s'\n",PositionsFileName);
        }

        else{
        	
            x= malloc(TP * sizeof(*x));
            y= malloc(TP * sizeof(*y));
            z= malloc(TP * sizeof(*z));
            tag= malloc(TP * sizeof(*tag));
            for(k=0;k<TP;k++){
                fscanf(positionFile,"%d %lf %lf %lf %lf %lf %lf\n",&(tag[k]),&(x[k]),&(y[k]),&(z[k]),&dummy2,&dummy3,&dummy4);
            }
            fclose(positionFile);
            for (i=0; i < 6; i++){
                fscanf(eigenvaluesFile, "%*[^\n]\n");
            }
            fscanf(eigenvaluesFile,"%lf %lf %lf",&(dummy2),&(dummy3),&(dummy4));
            fclose(eigenvaluesFile);
            
            if(dummy2>1e-5){
            cmx=0.e0, cmy=0.e0, cmz=0.e0;
            for(i=0;i<TP;i++){
                cmx+=x[i];
                cmy+=y[i];
                cmz+=z[i];
            }
            cmx/=double_TP;
            cmy/=double_TP;
            cmz/=double_TP;
            //printf("%E\t%E\n",cmx,cmy);
            
            fx= malloc(TP * sizeof(*fx));
            fy= malloc(TP * sizeof(*fy));
            fz= malloc(TP * sizeof(*fz));
            sig_xx= malloc(TP * sizeof(*sig_xx));
            sig_yy= malloc(TP * sizeof(*sig_yy));
            sig_zz= malloc(TP * sizeof(*sig_zz));
            
            
            for(i=0;i<TP;i++){
                fx[i]=0.e0;
                fy[i]=0.e0;
                fz[i]=0.e0;
                sig_xx[i]=0.e0;
                sig_yy[i]=0.e0;
                sig_zz[i]=0.e0;
            }
            pe = 0.e0;
            Tsig_xx=0.e0;
            Tsig_yy=0.e0;
            Tsig_zz=0.e0;
            for(i=0;i<TP-1;i++){
                xi=x[i];
                yi=y[i];
                zi=z[i];
                tag_i= tag[i];
                for(j=i+1;j<TP;j++){
                    tag_j = tag[j];
                    dx=x[j]-xi;
                    dy=y[j]-yi;
                    dz=z[j]-zi;
                    r2=(dx*dx)+(dy*dy)+(dz*dz);
                    if(tag_i==1 && tag_j==1 && r2<(r_cut_AA*r_cut_AA)){
                        s2= sig_AA*sig_AA;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_AA*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_AA*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx[j] += (store_x);  //
                        fx[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy[j] += (store_y);  //
                        fy[i] -= (store_y);  //
                        store_z=fr*dz;
                        fz[j] += (store_z);  //
                        fz[i] -= (store_z);  //
                        //-----------------------------------
                        sig_xx[i]+=store_x*dx;
                        sig_yy[i]+=store_y*dy;
                        sig_zz[i]+=store_z*dz;
                        sig_xx[j]+=store_x*dx;
                        sig_yy[j]+=store_y*dy;
                        sig_zz[j]+=store_z*dz;
                        Tsig_xx+=store_x*dx;
                        Tsig_yy+=store_y*dy;
                        Tsig_zz+=store_z*dz;
                    }
                    if (tag_i==0 && tag_j==0 && r2<(r_cut_BB*r_cut_BB)){
                        s2= sig_BB*sig_BB;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_BB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_BB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx[j] += (store_x);  //
                        fx[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy[j] += (store_y);  //
                        fy[i] -= (store_y);  //
                        store_z=fr*dz;
                        fz[j] += (store_z);  //
                        fz[i] -= (store_z);  //
                        //-----------------------------------
                        sig_xx[i]+=store_x*dx;
                        sig_yy[i]+=store_y*dy;
                        sig_zz[i]+=store_z*dz;
                        sig_xx[j]+=store_x*dx;
                        sig_yy[j]+=store_y*dy;
                        sig_zz[j]+=store_z*dz;
                        Tsig_xx+=store_x*dx;
                        Tsig_yy+=store_y*dy;
                        Tsig_zz+=store_z*dz;
                    }
                    if (tag_i!=tag_j && r2<(r_cut_AB*r_cut_AB)){
                        s2= sig_AB*sig_AB;
                        r2i = 1.e0/r2;
                        r6i = r2i*r2i*r2i*s2*s2*s2;
                        pe+=4.e0*eps_AB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                        fr =4.e0*eps_AB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                            
                        store_x=fr*dx;
                        fx[j] += (store_x);  //
                        fx[i] -= (store_x);  //
                        store_y=fr*dy;
                        fy[j] += (store_y);  //
                        fy[i] -= (store_y);  //
                        store_z=fr*dz;
                        fz[j] += (store_z);  //
                        fz[i] -= (store_z);  //
                        //-----------------------------------
                        sig_xx[i]+=store_x*dx;
                        sig_yy[i]+=store_y*dy;
                        sig_zz[i]+=store_z*dz;
                        sig_xx[j]+=store_x*dx;
                        sig_yy[j]+=store_y*dy;
                        sig_zz[j]+=store_z*dz;
                        Tsig_xx+=store_x*dx;
                        Tsig_yy+=store_y*dy;
                        Tsig_zz+=store_z*dz;
                    }
                }
            }
	   
            for(i=0;i<TP;i++){ 
                    dx = x[i]-cmx;
                    dy = y[i]-cmy;
                    dz = z[i]-cmz;
                    dr = sqrt((dx*dx)+(dy*dy)+(dz*dz));
                    
                    m=(int) (dr/delta_r);
                    
                    sxx_rad[m]+= sig_xx[i];
                    syy_rad[m]+= sig_yy[i];
                    szz_rad[m]+= sig_zz[i];
                    
                    num_rad[m]++;
            }
		
            free(fx);
            free(fy);
            free(fz);
            free(sig_xx);
            free(sig_zz);
            free(sig_yy);
            
             ens_id++;
            }
            free(x);
            free(y);
            free(z);
            free(tag);
	}
    }
    
    for(n=0;n<no_of_rad_bin;n++){
        sxx_rad[n]/=((4.e0/3.e0)*pi*(((n+1)*(n+1)*(n+1))-(n*n*n))*delta_r*delta_r*delta_r*ens_id*1.e0);
        syy_rad[n]/=((4.e0/3.e0)*pi*(((n+1)*(n+1)*(n+1))-(n*n*n))*delta_r*delta_r*delta_r*ens_id*1.e0);
        szz_rad[n]/=((4.e0/3.e0)*pi*(((n+1)*(n+1)*(n+1))-(n*n*n))*delta_r*delta_r*delta_r*ens_id*1.e0);
        
    }

    char rad_stress_FileName[128];
    FILE *rad_stress;
    sprintf(rad_stress_FileName,"%d_rad_stress_of_particles_%.2lf_bin.txt",TP,delta_r);
    rad_stress=fopen(rad_stress_FileName,"w");
    for(m=0;m<no_of_rad_bin;m++){
        fprintf(rad_stress,"%E\t%E\t%E\t%E\t%E\t%d\n",((m+0.5e0)*delta_r),sxx_rad[m],syy_rad[m],szz_rad[m],(sxx_rad[m]+syy_rad[m]+szz_rad[m])/3.e0,num_rad[m]);
    }
    


    fclose(rad_stress);
    
    free(sxx_rad);
    free(syy_rad);
    free(szz_rad);
    
	free(num_rad);
}
