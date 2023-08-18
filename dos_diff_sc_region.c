#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>

#define TP (400)
#define double_TP (400.e0)
int main(){
    int i, j, k;
    unsigned int no_config = 606000;
    
    char parti_surf_FileName[128];
    FILE *parti_surf;
    sprintf(parti_surf_FileName,"parti_surf_of_%06u_ens_id_%d_particles_cut_out_dynamics_damping_new.txt",no_config,TP);
    parti_surf = fopen(parti_surf_FileName,"rb");
    double dummy1, dummy2, dummy3, dummy4;
    int dummy5;
    unsigned int avail_configs=0, configs;
    while ( fscanf(parti_surf,"%lf %lf %lf %lf %d",&dummy1,&dummy2,&dummy3,&dummy4,&dummy5) != EOF ){
        avail_configs++;
    }
    rewind(parti_surf);
    
    double *eigenvalue,*participation_tot, *B2_tot, *sur_contr_tot;
    eigenvalue = malloc(avail_configs* sizeof(*eigenvalue));
    participation_tot = malloc(avail_configs* sizeof(*participation_tot));
    B2_tot = malloc(avail_configs* sizeof(*B2_tot));
    sur_contr_tot = malloc(avail_configs* sizeof(*sur_contr_tot));
    
    for(configs=0;configs<avail_configs;configs++){
        fscanf(parti_surf,"%lf %lf %lf %lf %d",&(eigenvalue[configs]),&(participation_tot[configs]),&(B2_tot[configs]),&(sur_contr_tot[configs]),&dummy5);
    }
    fclose(parti_surf);


    
    //###########################################################################
    double low_sur_contr = 0.e0, high_sur_contr = 0.245e0;
    unsigned int load_configs=0;
    for(configs=0;configs<avail_configs;configs++){
        if(sur_contr_tot[configs]>low_sur_contr && sur_contr_tot[configs]<=high_sur_contr)    load_configs++;
    }
    printf("%06u\n",load_configs);
    double *omega, *participation;
    omega = malloc(load_configs* sizeof(*omega));
    participation = malloc(load_configs* sizeof(*participation));
    load_configs=0;
    for(configs=0;configs<avail_configs;configs++){
        if(sur_contr_tot[configs]>low_sur_contr && sur_contr_tot[configs]<=high_sur_contr){
            omega[load_configs]=sqrt(eigenvalue[configs]);
            participation[load_configs]=participation_tot[configs];
            load_configs++;
        }
    }
    
    //###############################################################################
    //////////////////////////////// distribution of vibrational frequencies ////////////////////////////////////////////////////
    char omega_distribution_FileName[128];
    FILE *omega_distribution;
    sprintf(omega_distribution_FileName,"vibrational_density_from_%.2lf_to_%.2lf_%06u_configs_%d_TP_40_bin_new.txt",low_sur_contr,high_sur_contr,load_configs,TP);
    omega_distribution= fopen(omega_distribution_FileName,"w");
    
   
    unsigned int k_max, k_min;
    double high_omega, low_omega, start_omega, end_omega, factor_omega;
    high_omega = omega[0];
    low_omega = omega[0];

    for(configs=1;configs<load_configs;configs++){
        if(omega[configs]>high_omega){
            k_max = configs;       
            high_omega = omega[configs];
        }
        if(omega[configs]<low_omega){
            k_min=configs;       
            low_omega = omega[configs];
        }
    }
    printf("%06u\t%E\t\t%06u\t%E\t%E\t%06u\t%E\n",load_configs,high_omega,k_max,omega[k_max],low_omega,k_min,omega[k_min]);
    int no_of_binning_omega;
    no_of_binning_omega = 40;
    start_omega = low_omega-0.01e0;
    end_omega = high_omega+0.01e0;
    factor_omega = pow(end_omega/start_omega,1.e0/(no_of_binning_omega*1.e0));
    double *freq_omega, *cdf;
    int *no_config_omega;
    freq_omega = malloc(no_of_binning_omega * sizeof(*freq_omega));
    cdf = malloc(no_of_binning_omega * sizeof(*cdf));
    no_config_omega = malloc(no_of_binning_omega * sizeof(*no_config_omega));
    
    double *width_omega, *mid_omega;
    width_omega = malloc(no_of_binning_omega * sizeof(*width_omega));
    mid_omega = malloc(no_of_binning_omega * sizeof(*mid_omega));
    
    
    for(i=0;i<no_of_binning_omega;i++){
        width_omega[i]=start_omega*(factor_omega-1.e0)*pow(factor_omega,i*1.e0);
        mid_omega[i]=start_omega*pow(factor_omega,i+0.5e0);
        freq_omega[i]=0.e0;
        cdf[i]=0.e0;
        no_config_omega[i]=0;
    }
    for(configs=0;configs<load_configs;configs++){
        k= (int) ( log(omega[configs]/start_omega)/log(factor_omega) );
        freq_omega[k]+=1.e0;
        no_config_omega[k]++;
    }

    
    for(i=0;i<no_of_binning_omega;i++){
        freq_omega[i]=freq_omega[i]/(width_omega[i]*load_configs*1.e0);
        for(j=0;j<=i;j++){
             cdf[i]+=(freq_omega[j]*width_omega[j]);
        }
        fprintf(omega_distribution,"%E\t%E\t%E\t%d\n",mid_omega[i],freq_omega[i],cdf[i],no_config_omega[i]);
    }
    
    fclose(omega_distribution);
    free(freq_omega);
    free(cdf);
    free(no_config_omega);
    
    free(mid_omega);
    free(width_omega);
    free(omega);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    free(eigenvalue);
    free(participation_tot);
    free(B2_tot);
    free(sur_contr_tot);
    
}
