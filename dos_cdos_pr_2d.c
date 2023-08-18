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
//int dimN=TP*dim;
int dimN=103;
int main(){
    
	unsigned int no_config, config, load_config, dummy,k_max, k_min;
	no_config =500000;
   
    int i,k,j,l;
	double *omega;   
    omega = malloc(no_config*dimN* sizeof(*omega));
    //********** Open Eigen value File ************
    load_config = 0; 
    config = 0;
    
    for(config=0;config<no_config;config++){
        char eigenvalue_FileName[128];
        FILE *eigenvalueFile;
        sprintf(eigenvalue_FileName,"eigenvalues_of_%d_particles_%06u_ens_id.txt",TP,config;
        //sprintf(eigenvalue_FileName,"../cg_eigenvalues_%d/eigenvalues_of_%06u_particles.txt",TP,config);
        eigenvalueFile = fopen(eigenvalue_FileName,"rb");

    	if(eigenvalueFile!=NULL){
		 	double *lambda, *sqrt_lambda, *pr;
	   		lambda = malloc(dimN * sizeof(*lambda));
	   		sqrt_lambda = malloc(dimN * sizeof(*sqrt_lambda));
	   		pr = malloc(dimN* sizeof(*pr));
			for(k=0;k<dimN;k++){
				fscanf(eigenvalueFile,"%lf %lf %lf",&(lambda[k]),&(sqrt_lambda[k]),&(pr[k]));	
			}
    		
	    	if(sqrt_lambda[3]>1e-3){
				for(i=3;i<dimN;i++){
		        		omega[load_config*dimN+i] = sqrt_lambda[i];
		     	}
            	load_config++;
	    	}
	    
	    	fclose(eigenvalueFile);
			free(lambda);
			free(sqrt_lambda);
			free(pr);
		}
    }
    printf("%06u\n",load_config);
    for(i=0;i<load_config;i++){
        for(j=0;j<3;j++){
            omega[i*dimN+j]=0.e0;
        }
    }
      
    double high_freq, low_freq;
    high_freq = omega[0];
    low_freq = omega[3];
    for(dummy=0;dummy<load_config;dummy++){
    for(l=3;l<dimN;l++){
        if(omega[dummy*dimN+l]>high_freq){
            k_max = dummy*dimN+l;       
            high_freq = omega[dummy*dimN+l];
        }
        if(omega[dummy*dimN+l]<low_freq){
            k_min=dummy*dimN+l;       
            low_freq = omega[dummy*dimN+l];
        }
        }
    }
    printf("%E\t\t%06u\t%E\t%E\t%06u\t%E\n",high_freq,k_max,omega[k_max],low_freq,k_min,omega[k_min]);
    
    double start, end, factor;
    int no_of_binning;
    start = low_freq-0.001e0;
    end = high_freq+0.001e0;
    no_of_binning = 100;
    factor = pow(end/start,1.e0/(no_of_binning*1.e0));
    printf("%E\n",factor);    
   //............... Data Files ....................
	char density_of_states_FileName[128];
	FILE *density_of_states;
	sprintf(density_of_states_FileName,"density_of_states_100_freq_%d_dim_%d_particles_%06u_configs_at_%d_bin.txt",dim,TP,load_config,no_of_binning);
	density_of_states = fopen(density_of_states_FileName,"w");
	//...............................................   
        
	double *dos,*cdf, *bin_width, *mid_omega;
	int *bin_config;

	dos = malloc(no_of_binning * sizeof(*dos));
	cdf = malloc(no_of_binning * sizeof(*cdf));
	bin_width = malloc(no_of_binning * sizeof(*bin_width));
	mid_omega = malloc(no_of_binning * sizeof(*mid_omega));
        
	bin_config = malloc(no_of_binning * sizeof(*bin_config));
	for(i=0;i<no_of_binning;i++){
		bin_width[i]=start*(factor-1.e0)*pow(factor,i*1.e0);
		mid_omega[i]=start*pow(factor,i+0.5e0);
		dos[i]=0.e0;
		cdf[i]=0.e0;
		bin_config[i]=0;
	}
    for(dummy=0;dummy<load_config;dummy++){
        for(i=3;i<dimN;i++){
            k= (int) ( log(omega[dummy*dimN+i]/start)/log(factor) );
            dos[k]+=1.e0;
            bin_config[k]++;
        }
        
    }
    for(i=0;i<no_of_binning;i++){
	for(dummy=0;dummy<load_config;dummy++){
		for(j=3;j<dimN;j++){
			if(omega[dummy*dimN+j]<mid_omega[i]) cdf[i]+=1.e0;
		}
	}
    }
    for(i=0;i<no_of_binning;i++){
        dos[i]=dos[i]/(bin_width[i]*(load_config*(dimN-3)));
	cdf[i]=cdf[i]/(load_config*(dimN-3)*1.e0);      
        fprintf(density_of_states,"%E\t%E\t%d\t%E\n",mid_omega[i],dos[i],bin_config[i],cdf[i]);
    }
    
   

    free(omega);
    fclose(density_of_states);
    
    free(cdf);
    free(dos);
    free(bin_width);
    free(mid_omega);
    free(bin_config);
}
 
