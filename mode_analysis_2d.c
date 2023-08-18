#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>

#define TP (400)
#define double_TP (400.e0)
#define dim (2)
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

int *nn;// *sur_par;
#define neigh_dist (1.4e0)
int main(){
int f_up;
    unsigned int no_config = 606000, ens_id, load_config=0, config;
    double *x, *y, *delta_x, *delta_y;
    int *tag;
    double *omega_sq, *participation, *omega;
    int dimN = TP*dim;
    int i,j,k;

    double *B2, *B3, *B4, *PR, *sur_contr;
    B2 = malloc(no_config * sizeof(*B2));
    B3 = malloc(no_config * sizeof(*B3));
    B4 = malloc(no_config * sizeof(*B4));
    PR = malloc(no_config * sizeof(PR));
    sur_contr = malloc(no_config * sizeof(sur_contr));
    for(ens_id=0;ens_id<no_config;ens_id++){
        B2[ens_id]=0.e0;
        B3[ens_id]=0.e0;
        B4[ens_id]=0.e0;
	PR[ens_id]=0.e0;
        sur_contr[ens_id]=0.e0;
    }
    double eigenvalue, participation_ratio;
    double dx, dy, r2, r, phi, phi_r, phi_rr, phi_rrr, phi_rrrr, s2, s4, s6, s12, ri2, ri6, ri8, ri14, ri10, ri12, ri16;
    double Fou_deri_1, Fou_deri_2, Fou_deri_3, Tri_deri_1, Tri_deri_2;

    int cor1, cor2, cor3, cor4;
    double si_ij_1,si_ij_2,si_ij_3,si_ij_4,delta_r_ij_1,delta_r_ij_2,delta_r_ij_3,delta_r_ij_4;
    char NQLM_FileName[128];
    FILE *NQLM_File;    
    sprintf(NQLM_FileName,"NQLM_of_%06u_ens_id_%d_particles_cut_out_dynamics_damping.txt",no_config,TP);
    NQLM_File=fopen(NQLM_FileName,"w");


	ens_id =0;
    for(config=0;config<no_config;config++){
        
        char eigenvalue_FileName[128];
        FILE *eigenvalue_File; 
        sprintf(eigenvalue_FileName,"eigenvalues_of_%d_particles_%06u_ens_id.txt",TP,config);
        eigenvalue_File = fopen(eigenvalue_FileName,"rb");
        

        char eigenvector_FileName[128];
        FILE *eigenvector_File;        
        sprintf(eigenvector_FileName,"first_non_zero_mode_of_%d_particles_%06u_ens_id.txt",TP,config);
        eigenvector_File=fopen(eigenvector_FileName,"rb");

        char configuration_FileName[128];
        FILE *configuration_File;
        sprintf(configuration_FileName,"min_configs_%d_TP_%06u_ens_id.txt",TP,config);
        configuration_File=fopen(configuration_FileName,"rb");    

        if(eigenvalue_File == NULL || eigenvector_File == NULL || configuration_File == NULL){
            printf("ERROR: No file at: '%d'\n",config);
            if (eigenvalue_File != NULL){
                fclose(eigenvalue_File);
            }
            if (eigenvector_File != NULL){
                fclose(eigenvector_File);
            }
            if (configuration_File != NULL){
                fclose(configuration_File);
            }
        }
        else
        {
    
            omega_sq = malloc(dimN * sizeof(*omega_sq));
            omega = malloc(dimN * sizeof(*omega));
            participation = malloc(dimN * sizeof(*participation));
            
            x = malloc(TP * sizeof(*x));
            y = malloc(TP * sizeof(*y));
            delta_x = malloc(TP * sizeof(*delta_x));
            delta_y = malloc(TP * sizeof(*delta_y));
            tag = malloc(TP * sizeof(*tag));
            nn= malloc(TP * sizeof(*nn));
            k=0;
            while(  fscanf(eigenvalue_File,"%lf %lf %lf",&(omega_sq[k]),&(omega[k]),&(participation[k])) !=EOF){
                  k++;
            }
            fclose(eigenvalue_File);
        
            k=0;
            while ( fscanf(eigenvector_File,"%lf %lf %lf %lf",&(x[k]),&(y[k]),&(delta_x[k]),&(delta_y[k])) != EOF){
                k++;
            }
            fclose(eigenvector_File);
        
            k=0;
            while ( fscanf(configuration_File,"%d %lf %lf",&(tag[k]),&(x[k]),&(y[k])) != EOF){
                k++;
            }
            fclose(configuration_File);

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
            for(i=0;i<TP;i++){
                nn[i]=0;
            }

            for(k=0;k<(TP-1);k++){
                for(j=k+1;j<TP;j++){
                    dx=x[j]-x[k];
                    dy=y[j]-y[k];
                    r2=((dx*dx)+(dy*dy));
                    if(r2<=(neigh_dist*neigh_dist)){
                        nn[k]++;  
                        nn[j]++;
                    }
                }
            }
          
		    double e2, e4, deno = 0.e0, nume=0.e0, nume_sur=0.e0;
		    for(i=0;i<TP;i++){               
		        e2 = (delta_x[i]*delta_x[i])+(delta_y[i]*delta_y[i]);
		    	e4 = e2*e2;
		    	deno+=e4;
                	nume+=e2;
                	//if(nn[i]<20) nume_sur+=e2;
			if(((nn[i]*1.e0)/(acos(-1.e0)*neigh_dist*neigh_dist))<0.8e0) nume_sur+=e2;
		    }
		    PR[ens_id]=(1.e0/(double_TP*deno));
            sur_contr[ens_id]=nume_sur;

        for(i=0;i<TP-1;i++){
            for(j=i+1;j<TP;j++){
                dx = x[j]-x[i];
                dy = y[j]-y[i];
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
                r2 = (dx*dx)+(dy*dy);

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
                    phi_rrr = 4.e0*eps_ij*(-(12.e0*13.e0*14.e0*s12*ri14/r)+(6.e0*7.e0*8.e0*s6*ri8/r)+(24.e0*c_4*r/s4));
                    phi_rrrr = 4.e0*eps_ij*((12.e0*13.e0*14.e0*15.e0*s12*ri16)-(6.e0*7.e0*8.e0*9.e0*s6*ri10)+(24.e0*c_4/s4));

                    Tri_deri_1 = ((phi_rrr*ri2/r)-(3.e0*phi_rr*ri2*ri2)+(3.e0*phi_r*r*ri6));
                    Tri_deri_2 =  ((phi_rr*ri2)-(phi_r*ri2/r));

                    Fou_deri_1 = ((phi_rrrr*ri2*ri2)-(6.e0*phi_rrr*r*ri6)+(15.e0*phi_rr*ri6)-(15.e0*phi_r*r*ri8));
                    Fou_deri_2 = ((phi_rrr*r*ri2*ri2)-(3.e0*phi_rr*ri2*ri2)+(3.e0*phi_r*r*ri6));
                    Fou_deri_3 = ((phi_rr*ri2)-(phi_r*r*ri2*ri2));

                    for(cor1=0;cor1<2;cor1++){
                    	for(cor2=0;cor2<2;cor2++){
                        	if(cor1==0)	{
                                delta_r_ij_1=dx;
                                si_ij_1 = delta_x[j]-delta_x[i];
                            }
                            if(cor1==1){
                                delta_r_ij_1=dy;
                                si_ij_1 = delta_y[j]-delta_y[i];
                            }	

                            if(cor2==0){
                                delta_r_ij_2=dx;
                                si_ij_2 = delta_x[j]-delta_x[i];			
                            }	
                            if(cor2==1){
                                delta_r_ij_2=dy;
                                si_ij_2 = delta_y[j]-delta_y[i];
                            }
                            B2[ens_id]+=(((phi_rr*ri2)-(phi_r*ri2/r))*delta_r_ij_1*si_ij_1*delta_r_ij_2*si_ij_2);
                            if(cor1==cor2) B2[ens_id]+=((phi_r/r)*si_ij_1*si_ij_2);
                            
                    	}	
                	}
                    
                    for(cor1=0;cor1<2;cor1++){
                        for(cor2=0;cor2<2;cor2++){
                            for(cor3=0;cor3<2;cor3++){
                                if(cor1==0)	{
                                    delta_r_ij_1=dx;
                                    si_ij_1 = delta_x[j]-delta_x[i];
                                }
                                if(cor1==1){
                                    delta_r_ij_1=dy;
                                    si_ij_1 = delta_y[j]-delta_y[i];
                                }	

                                if(cor2==0){
                                    delta_r_ij_2=dx;
                                    si_ij_2 = delta_x[j]-delta_x[i];			
                                }	
                                if(cor2==1){
                                    delta_r_ij_2=dy;
                                    si_ij_2 = delta_y[j]-delta_y[i];
                                }
                                if(cor3==0){
                                    delta_r_ij_3=dx;
                                    si_ij_3 = delta_x[j]-delta_x[i];			
                                }	
                                if(cor3==1){
                                    delta_r_ij_3=dy;
                                    si_ij_3 = delta_y[j]-delta_y[i];
                                }
                                B3[ens_id]+=((Tri_deri_1)*delta_r_ij_1*si_ij_1*delta_r_ij_2*si_ij_2*delta_r_ij_3*si_ij_3);

                                if(cor2==cor3) B3[ens_id]+=(Tri_deri_2*delta_r_ij_1*si_ij_1*si_ij_2*si_ij_2);
                                if(cor1==cor3) B3[ens_id]+=(Tri_deri_2*delta_r_ij_2*si_ij_2*si_ij_1*si_ij_1);
                                if(cor1==cor2) B3[ens_id]+=(Tri_deri_2*delta_r_ij_3*si_ij_3*si_ij_1*si_ij_1);
                                

                            }
                        }
                    }
                    for(cor1=0;cor1<2;cor1++){
                        for(cor2=0;cor2<2;cor2++){
                            for(cor3=0;cor3<2;cor3++){
                                for(cor4=0;cor4<2;cor4++){
                                    if(cor1==0)	{
                                        delta_r_ij_1=dx;
                                        si_ij_1 = delta_x[j]-delta_x[i];
                                    }
                                    if(cor1==1){
                                        delta_r_ij_1=dy;
                                        si_ij_1 = delta_y[j]-delta_y[i];
                                    }	

                                    if(cor2==0){
                                        delta_r_ij_2=dx;
                                        si_ij_2 = delta_x[j]-delta_x[i];			
                                    }	
                                    if(cor2==1){
                                        delta_r_ij_2=dy;
                                        si_ij_2 = delta_y[j]-delta_y[i];
                                    }
                                    if(cor3==0){
                                        delta_r_ij_3=dx;
                                        si_ij_3 = delta_x[j]-delta_x[i];			
                                    }	
                                    if(cor3==1){
                                        delta_r_ij_3=dy;
                                        si_ij_3 = delta_y[j]-delta_y[i];
                                    }
                                    if(cor4==0){
                                        delta_r_ij_4=dx;
                                        si_ij_4 = delta_x[j]-delta_x[i];			
                                    }	
                                    if(cor4==1){
                                        delta_r_ij_4=dy;
                                        si_ij_4 = delta_y[j]-delta_y[i];
                                    }
                                    B4[ens_id]+=(Fou_deri_1*delta_r_ij_1*si_ij_1*delta_r_ij_2*si_ij_2*delta_r_ij_3*si_ij_3*delta_r_ij_4*si_ij_4);
                                    if (cor3==cor4) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_1*si_ij_1*delta_r_ij_2*si_ij_2*si_ij_3*si_ij_3);
                                    if (cor1==cor3) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_4*si_ij_4*delta_r_ij_2*si_ij_2*si_ij_3*si_ij_3);
                                    if (cor2==cor4) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_1*si_ij_1*delta_r_ij_3*si_ij_3*si_ij_2*si_ij_2);
                                    if (cor1==cor4) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_3*si_ij_3*delta_r_ij_2*si_ij_2*si_ij_1*si_ij_1);
                                    if (cor2==cor3) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_1*si_ij_1*delta_r_ij_4*si_ij_4*si_ij_2*si_ij_2);
                                    if (cor1==cor2) B4[ens_id]+=((Fou_deri_2)*delta_r_ij_3*si_ij_3*delta_r_ij_4*si_ij_4*si_ij_1*si_ij_1);

                                    if (cor1==cor3 && cor2==cor4) B4[ens_id]+=(Fou_deri_3*si_ij_1*si_ij_1*si_ij_2*si_ij_2);
                                    if (cor2==cor3 && cor1==cor4) B4[ens_id]+=(Fou_deri_3*si_ij_1*si_ij_1*si_ij_2*si_ij_2);
                                    if (cor1==cor2 && cor3==cor4) B4[ens_id]+=(Fou_deri_3*si_ij_1*si_ij_1*si_ij_3*si_ij_3);
                                }


                            }
                        }
                    }


			
                }
            }
        }
            f_up=0;
            if((3.e0*B2[ens_id]*B4[ens_id])<(B3[ens_id]*B3[ens_id])) f_up = 1;
            if(((fabs(B2[ens_id]-eigenvalue))<1e-5)){
                fprintf(NQLM_File,"%E\t%E\t%E\t%E\t%E\t%E\t%d\n",eigenvalue,PR[ens_id],B2[ens_id],B3[ens_id],B4[ens_id],sur_contr[ens_id],f_up);
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
            free(nn);
        }
    }
    
    free(B2);
    free(B3);
    free(B4);
    free(PR);
    free(sur_contr);
    fclose(NQLM_File);
 
}
