#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>

#define TP (4096)
#define double_TP (4096.e0)
#define dim (3)
//..................
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
double r_cut_ij, sig_ij, eps_ij;

int main(){
    unsigned int no_config = 250000, ens_id, load_config=0, config;
    double *x, *y, *z, *delta_x, *delta_y, *delta_z;
    int *tag;
    double *omega_sq, *participation, *omega;
    int dimN = TP*dim;
    int i,j,k;

    double *B2, *PR, *sur_contr;
    B2 = malloc(no_config * sizeof(*B2));
    PR = malloc(no_config * sizeof(PR));
    sur_contr = malloc(no_config * sizeof(sur_contr));
    for(ens_id=0;ens_id<no_config;ens_id++){
        B2[ens_id]=0.e0;
	PR[ens_id]=0.e0;
        sur_contr[ens_id]=0.e0;
    }
    double eigenvalue, participation_ratio,cmx,cmy,cmz;
    double dx, dy, dz, r2, r, phi, phi_r, phi_rr, phi_rrr, phi_rrrr, s2, s4, s6, s12, ri2, ri6, ri8, ri14, ri10, ri12, ri16;

    int cor1, cor2, cor3, cor4;
    double si_ij_1,si_ij_2,si_ij_3,si_ij_4,delta_r_ij_1,delta_r_ij_2,delta_r_ij_3,delta_r_ij_4;
    char NQLM_FileName[128];
    FILE *NQLM_File;    
    sprintf(NQLM_FileName,"parti_surf_of_%06u_ens_id_%d_particles_cut_out_dynamics_damping.txt",no_config,TP);
    NQLM_File=fopen(NQLM_FileName,"w");
    ens_id =0;
    for(config=0;config<no_config;config++){
        
        char eigenvalue_FileName[128];
        FILE *eigenvalue_File; 
        sprintf(eigenvalue_FileName,"eigenvalues_of_%d_particles_%06u_ens_id.txt",TP,config);
        eigenvalue_File = fopen(eigenvalue_FileName,"r");
        

        char eigenvector_FileName[128];
        FILE *eigenvector_File;        
        sprintf(eigenvector_FileName,"first_non_zero_mode_of_%d_particles_%06u_ens_id.txt",TP,config);
        eigenvector_File=fopen(eigenvector_FileName,"r");

       
        if(eigenvalue_File == NULL || eigenvector_File == NULL){
            printf("ERROR: No file at: '%d'\n",config);
            if (eigenvalue_File != NULL){
                fclose(eigenvalue_File);
            }
            if (eigenvector_File != NULL){
                fclose(eigenvector_File);
            }
         
        }
        else
        {
    
            omega_sq = malloc(dimN * sizeof(*omega_sq));
            omega = malloc(dimN * sizeof(*omega));
            participation = malloc(dimN * sizeof(*participation));
            
            x = malloc(TP * sizeof(*x));
            y = malloc(TP * sizeof(*y));
            z = malloc(TP * sizeof(*z));
            delta_x = malloc(TP * sizeof(*delta_x));
            delta_y = malloc(TP * sizeof(*delta_y));
            delta_z = malloc(TP * sizeof(*delta_z));
            tag = malloc(TP * sizeof(*tag));
            k=0;
            while(  fscanf(eigenvalue_File,"%lf %lf %lf",&(omega_sq[k]),&(omega[k]),&(participation[k])) !=EOF){
                  k++;
            }
            fclose(eigenvalue_File);
        
            k=0;
            while ( fscanf(eigenvector_File,"%d %lf %lf %lf %lf %lf %lf",&(tag[k]),&(x[k]),&(y[k]),&(z[k]),&(delta_x[k]),&(delta_y[k]),&(delta_z[k])) != EOF){
                k++;
            }
            fclose(eigenvector_File);
        
        
            // ------------- to distingush the correct configs --------------
            for(i=0;i<dimN;i++){
                if(omega_sq[i]<1e-5 && omega_sq[i]<(-1e-10)){
                    break;
                }
                if(omega_sq[i]>1e-5){
                    eigenvalue = omega_sq[i];
                    participation_ratio = participation[i];
                    break;
                }                   
            }
            //----------------------------------------------------------
            //----------------- surface modes --------------------------
            cmx = 0.e0, cmy =0.e0, cmz=0.e0;
	     

		   for(i=0;i<TP;i++){
		   		cmx+=x[i];
				cmy+=y[i];
				cmz+=z[i];
		    }
			cmx/=(TP*1.e0);
			cmy/=(TP*1.e0);
			cmz/=(TP*1.e0);   
	    	double e2, e4, deno = 0.e0, nume = 0.e0, nume_sur = 0.e0; 
            for(k=0;k<TP;k++){
                dx=x[k]-cmx;
                dy=y[k]-cmy;
		    	dz=z[k]-cmz;
                r2=((dx*dx)+(dy*dy)+(dz*dz));               
		        e2 = (delta_x[k]*delta_x[k])+(delta_y[k]*delta_y[k])+(delta_z[k]*delta_z[k]);
		    	e4 = e2*e2;
		    	deno+=e4;
                nume+=e2;
				if(sqrt(r2)>6.e0) nume_sur+=e2;
		    }
		    PR[ens_id]=(1.e0/(double_TP*deno));
            sur_contr[ens_id]=nume_sur;

        	for(i=0;i<TP-1;i++){
            	for(j=i+1;j<TP;j++){
                	dx = x[j]-x[i];
                	dy = y[j]-y[i];
					dz = z[j]-z[i];
                	if(tag[i]==1 && tag[j] == 1){
                       	r_cut_ij = r_cut_AA;
                       	sig_ij = sig_AA;
                       	eps_ij = eps_AA;
                     }
		             if(tag[i]==0 && tag[j] == 0){
		             	r_cut_ij = r_cut_BB;
		             	sig_ij = sig_BB;
		             	eps_ij = eps_BB;
		             }
                 	if(tag[i]!=tag[j]){
                       	r_cut_ij = r_cut_AB;
                       	sig_ij = sig_AB;
                       	eps_ij = eps_AB;
                     }
                	r2 = (dx*dx)+(dy*dy)*(dz*dz);
                if(r2<r_cut_ij*r_cut_ij){
                    r = sqrt(r2);
                    s2= sig_ij*sig_ij;
                    ri2 = 1.e0/r2;
                    s4 = s2*s2;
                    s6 = s4*s2;
                    s12= s6*s6;
                    ri6 = ri2*ri2*ri2;
                    ri8 = ri6*ri2;
                    ri10 = ri8*ri2;
                    ri12 = ri10*ri2;
                    ri14 = ri12*ri2;
                    ri16 = ri14*ri2;
                    phi = 4.e0*eps_ij*((s12*ri12)-(s6*ri6)+c_0+(c_2*r2/s2)+(c_4*r2*r2/s4));
                    phi_r = 4.e0*eps_ij*(-(12.e0*s12*ri12/r)+(6.e0*s6*ri6/r)+(2.e0*c_2*r/s2)+(4.e0*c_4*r2*r/s4));
                    phi_rr = 4.e0*eps_ij*((12.e0*13.e0*s12*ri14)-(6.e0*7.e0*s6*ri8)+(2.e0*c_2/s2)+(12.e0*c_4*r2/s4));
                    for(cor1=0;cor1<3;cor1++){
                    	for(cor2=0;cor2<3;cor2++){
                        	if(cor1==0)	{
                                delta_r_ij_1=dx;
                                si_ij_1 = delta_x[j]-delta_x[i];
                            }
                            if(cor1==1){
                                delta_r_ij_1=dy;
                                si_ij_1 = delta_y[j]-delta_y[i];
                            }	
							if(cor1==2){
                                delta_r_ij_1=dz;
                                si_ij_1 = delta_z[j]-delta_z[i];
                            }
                            
                            if(cor2==0){
                                delta_r_ij_2=dx;
                                si_ij_2 = delta_x[j]-delta_x[i];			
                            }	
                            if(cor2==1){
                                delta_r_ij_2=dy;
                                si_ij_2 = delta_y[j]-delta_y[i];
                            }
                            if(cor2==2){
                                delta_r_ij_2=dz;
                                si_ij_2 = delta_z[j]-delta_z[i];
                            }
                            
                            B2[ens_id]+=(((phi_rr*ri2)-(phi_r*ri2/r))*delta_r_ij_1*si_ij_1*delta_r_ij_2*si_ij_2);
                            if(cor1==cor2) B2[ens_id]+=((phi_r/r)*si_ij_1*si_ij_2);
                            
                  	}	
       	           }
                    
            	}
       		 }
	}
            if(((fabs(B2[ens_id]-eigenvalue))<1e-5)){
                fprintf(NQLM_File,"%E\t%E\t%E\t%E\t%06u\n",eigenvalue,PR[ens_id],B2[ens_id],sur_contr[ens_id],config);
                load_config++;
            } 
            
			ens_id++;

            free(x);
            free(y);
            free(tag);
            free(delta_x);
            free(delta_y);
            free(omega_sq);
            free(omega);
            free(participation);
        }
    }
    
    free(B2);
    free(PR);
    free(sur_contr);
    fclose(NQLM_File);
 
}

