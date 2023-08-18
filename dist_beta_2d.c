#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>


int main(){
    int i, j, k;
    
        int TP = 400;
    unsigned int no_config = 200000,dummy;
    char low_PR_NQLM_FileName[128];
    FILE *low_PR_NQLM;
    sprintf(low_PR_NQLM_FileName,"NQLM_of_%06u_ens_id_%d_particles_cut_out_cg.txt",no_config,TP);
    low_PR_NQLM = fopen(low_PR_NQLM_FileName,"rb");
    double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
    int dummy7, dummy8, dummy9;
    unsigned int load_configs=0, configs;
    while ( fscanf(low_PR_NQLM,"%lf %lf %lf %lf %lf %lf %d %d %d",&dummy1,&dummy2,&dummy3,&dummy4,&dummy5,&dummy6,&dummy7,&dummy8,&dummy9) != EOF ){
        load_configs++;
    }
    rewind(low_PR_NQLM);
    
    double *eigenvalue, *participation, *B2, *B3, *B4, *sur_contr;
    eigenvalue = malloc(load_configs* sizeof(*eigenvalue));
    participation = malloc(load_configs* sizeof(*participation));
    B2 = malloc(load_configs* sizeof(*B2));
    B3 = malloc(load_configs* sizeof(*B3));
    B4 = malloc(load_configs* sizeof(*B4));
    sur_contr = malloc(load_configs* sizeof(*sur_contr));
    
    for(configs=0;configs<load_configs;configs++){
        fscanf(low_PR_NQLM,"%lf %lf %lf %lf %lf %lf %d %d %d ",&(eigenvalue[configs]),&(participation[configs]),&(B2[configs]),&(B3[configs]),&(B4[configs]),&(sur_contr[configs]),&dummy7,&dummy8,&dummy9);
    }
    
    //----------------------- beta = B3/sqrt(3B2B4)  -------------------------------------------------
    printf("\nbeta distribution\n");
    double low_beta, high_beta ,beta;

    high_beta = B3[0]/sqrt(3.e0*B2[0]*B4[0]);
    low_beta = B3[0]/sqrt(3.e0*B2[0]*B4[0]);

    for(configs=1;configs<load_configs;configs++){
        beta = B3[configs]/sqrt(3.e0*B2[configs]*B4[configs]);
        if(beta>high_beta){
            k_max = configs;       
            high_beta = beta;
        }
        if((beta)<low_beta){
            k_min=configs;       
            low_beta = beta;
        }
    }
    printf("%06u\t%E\t\t%06u\t%E\t%E\t%06u\t%E\n",load_configs,high_beta,k_max,B3[k_max]/sqrt(3.e0*B2[k_max]*B4[k_max]),low_beta,k_min,B3[k_min]/sqrt(3.e0*B2[k_min]*B4[k_min]));
    
    int no_of_binning_beta;
    no_of_binning_beta = 50;
    
    double *freq_beta;
    int *no_of_config_beta;
    freq_beta = malloc(no_of_binning_beta * sizeof(*freq_beta));
    no_of_config_beta = malloc(no_of_binning_beta * sizeof(*no_of_config_beta));
    
    double width_beta;
    width_beta= (high_beta-low_beta)/(no_of_binning_beta*1.e0);
    printf("%E\n",width_beta);
    
    for(i=0;i<no_of_binning_beta;i++){
        freq_beta[i]=0.e0;
        no_of_config_beta[i]=0;
    }
  
    for(configs=0;configs<load_configs;configs++){
	beta = B3[configs]/sqrt(3.e0*B2[configs]*B4[configs]);
	if(B4[configs]>1e-15){
        k= (int) ((beta-low_beta)/width_beta);
	freq_beta[k]+=1.e0;
        no_of_config_beta[k]++;
	}
	
    }
   
    char beta_distribution_FileName[138];
    FILE *beta_distribution;
    sprintf(beta_distribution_FileName,"%d_beta_distribution_linear_of_%d_particles_CG.txt",load_configs,TP);
    beta_distribution= fopen(beta_distribution_FileName,"w");
    for(i=0;i<no_of_binning_beta;i++){
        if(no_of_config_beta[i]!=0){
            freq_beta[i]=freq_beta[i]/(width_beta*load_configs*1.e0);
            fprintf(beta_distribution,"%E\t%E\n",(i+0.5e0)*width_beta+low_beta,freq_beta[i]);
        }    
    }
    
    double mean_beta =0.e0 , std_beta = 0.e0;
    for(configs=0;configs<load_configs;configs++){
    	mean_beta+=(B3[configs]/sqrt(3.e0*B2[configs]*B4[configs]));
    }
    mean_beta/=(load_configs*1.e0);
    for(configs=0;configs<load_configs;configs++){
    	std_beta+=(((B3[configs]/sqrt(3.e0*B2[configs]*B4[configs]))*(B3[configs]/sqrt(3.e0*B2[configs]*B4[configs])))/(load_configs*1.e0));
    }
    std_beta-=(mean_beta*mean_beta);
    printf("beta mean = %E\t std_devi = %E\n",mean_beta,std_beta);
    fclose(beta_distribution);
    free(freq_beta);
    free(no_of_config_beta);
    //----------------------------------------------------------------------------------------------------------------------
    
    

      free(eigenvalue);
    free(participation);
    free(B2);
    free(B3);
    free(sur_contr);
    free(B4);
}
