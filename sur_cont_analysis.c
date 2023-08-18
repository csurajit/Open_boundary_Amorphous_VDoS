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
    unsigned int load_configs=0, configs;
    while ( fscanf(parti_surf,"%lf %lf %lf %lf %d",&dummy1,&dummy2,&dummy3,&dummy4,&dummy5) != EOF ){
        load_configs++;
    }
    rewind(parti_surf);
    
    double *omega, *participation, *B2, *pr_sur;
    omega = malloc(load_configs* sizeof(*omega));
    participation = malloc(load_configs* sizeof(*participation));
    B2 = malloc(load_configs* sizeof(*B2));
    pr_sur = malloc(load_configs* sizeof(*pr_sur));
    
    for(configs=0;configs<load_configs;configs++){
        fscanf(parti_surf,"%lf %lf %lf %lf %d",&(dummy1),&(participation[configs]),&(B2[configs]),&(pr_sur[configs]),&dummy5);
    }
    fclose(parti_surf);

    
    //------------------ surface_contribution ------------------------------
    printf("\nSURFACE CONTRIBUTION\n");
    double high_sur_contr, low_sur_contr;
    int k_max, k_min;
    high_sur_contr = pr_sur[0];
    low_sur_contr = pr_sur[0];

    for(configs=1;configs<load_configs;configs++){
        if(pr_sur[configs]>high_sur_contr){
            k_max = configs;       
            high_sur_contr = pr_sur[configs];
        }
        if(pr_sur[configs]<low_sur_contr){
            k_min=configs;       
            low_sur_contr = pr_sur[configs];
        }
    }   printf("%06u\t%E\t\t%06u\t%E\t%E\t%06u\t%E\n",load_configs,high_sur_contr,k_max,pr_sur[k_max],low_sur_contr,k_min,pr_sur[k_min]);
    int no_of_binning_sur_contr;
    no_of_binning_sur_contr = 60;
    double *freq_sur_contr;
    freq_sur_contr = malloc(no_of_binning_sur_contr * sizeof(*freq_sur_contr));
    double width_sur_contr;
    width_sur_contr= (high_sur_contr-low_sur_contr)/(no_of_binning_sur_contr*1.e0);
    printf("%E\n",width_sur_contr);
    for(i=0;i<no_of_binning_sur_contr;i++){
        freq_sur_contr[i]=0.e0;
    }
    for(configs=0;configs<load_configs;configs++){
        k= (int) ((pr_sur[configs]-low_sur_contr)/width_sur_contr);
        freq_sur_contr[k]+=1.e0;
    }
    char surface_contribution_distribution_FileName[128];
    FILE *surface_contribution_distribution;
    sprintf(surface_contribution_distribution_FileName,"%d_surface_contribution_distribution_linear_%d_TP.txt",configs,TP);
    surface_contribution_distribution= fopen(surface_contribution_distribution_FileName,"w");
    for(i=0;i<no_of_binning_sur_contr;i++){
        freq_sur_contr[i]=freq_sur_contr[i]/(width_sur_contr*load_configs*1.e0);
        fprintf(surface_contribution_distribution,"%E\t%E\n",(i+0.5e0)*width_sur_contr+low_sur_contr,freq_sur_contr[i]);
    }
    double mean_sur_contr =0.e0 , std_sur_contr = 0.e0;
    for(configs=0;configs<load_configs;configs++){
    	mean_sur_contr+=pr_sur[configs];
    }
    mean_sur_contr/=(load_configs*1.e0);
    for(configs=0;configs<load_configs;configs++){
    	std_sur_contr+=(((pr_sur[configs]*pr_sur[configs])/(load_configs*1.e0)));
    }
    std_sur_contr-=(mean_sur_contr*mean_sur_contr);
    printf("Surface contribution mean = %E\t std_devi = %E\n",mean_sur_contr,std_sur_contr);
    fclose(surface_contribution_distribution);
    free(freq_sur_contr);
    //--------------------------------------------------------------
    
    free(omega);
    free(participation);
    free(pr_sur);
    free(B2);
}
