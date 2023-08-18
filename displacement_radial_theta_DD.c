#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>
clock_t start, end;
double  cpu_time_used;


#define TP (400)
#define double_TP (400.e0)

double r_cm,cmx,cmy;
double *x, *y, *x_0, *y_0;
int *tag;


int main(){
    start = clock();  // Initializing clock
    unsigned int ens_id,load_config=0,starting = 100000, max = 200000;
    int i, j, k, n;
    
    
    double pi = acos(-1.e0);
    int no_of_radial_bin = 100;
     double dr = 28.e0/(no_of_radial_bin*1.e0);
    double displacement_radial[no_of_radial_bin],displacement_theta[no_of_radial_bin],density[no_of_radial_bin],theta,mod_displacement_radial[no_of_radial_bin],mod_displacement_theta[no_of_radial_bin];

    double no_of_particle[no_of_radial_bin];

    for(j=0;j<no_of_radial_bin;j++) {
    	no_of_particle[j]=0.e0;
    	displacement_radial[j]=0.e0;
    	displacement_theta[j]=0.e0;
    	density[j]=0.e0;
    
	mod_displacement_radial[j]=0.e0;
        mod_displacement_theta[j]=0.e0;
	
    }
    
    for(ens_id=starting;ens_id<max;ens_id++){
    //printf("%06u\n",ens_id);
    //--------------- Input Data Files --------------
        char dataFileName[128];
        FILE *dataFile;
        sprintf(dataFileName,"min_configs_%d_TP_%06u_ens_id.txt",TP,ens_id);
        dataFile=fopen(dataFileName,"rb");
        
        char data1FileName[128];
        FILE *data1File;
        sprintf(data1FileName,"%d_particle_%06u_cut_liquid_configuration.txt",TP,ens_id);
        data1File=fopen(data1FileName,"rb");
        
        if(dataFile==NULL || data1File==NULL) {
                printf("ERROR No File at :'%s'\n",dataFileName);
        }
        else{
        x_0= malloc(TP * sizeof(*x_0));
        y_0= malloc(TP * sizeof(*y_0));
        x= malloc(TP * sizeof(*x));
        y= malloc(TP * sizeof(*y));
        tag = malloc (TP * sizeof(*tag));
 
        for(i=0;i<TP;i++){
            fscanf(dataFile,"%d %lf %lf",&(j),&(x[i]),&(y[i]));
        }
        fclose(dataFile);
        for(i=0;i<TP;i++){
            fscanf(data1File,"%d %lf %lf",&(tag[i]),&(x_0[i]),&(y_0[i]));
        }
        fclose(data1File);
        
    	cmx = 0.e0;
    	cmy = 0.e0;
    	for(i=0;i<TP;i++){
    		cmx+=x_0[i];
    		cmy+=y_0[i];
    	}
    	cmx/=double_TP;
    	cmy/=double_TP;
    	r_cm = sqrt((cmx*cmx)+(cmy*cmy));  
    	double delta_x, delta_y;
    
    	for(i=0;i<TP;i++){
    		delta_x = x_0[i]-cmx;
    		delta_y = y_0[i]-cmy;
    		theta = atan(delta_y/delta_x);
    	
    		if(delta_x<0.e0 && delta_y>0.e0) theta=theta+pi;
        	if(delta_x<0.e0 && delta_y<0.e0) theta=theta+pi;
        	if(delta_x>0.e0 && delta_y<0.e0) theta=theta+(2.e0*pi);
    	
    		n = (int) (sqrt(delta_x*delta_x + delta_y*delta_y)/dr);
    		no_of_particle[n]+=1.e0;
    		displacement_radial[n]+= (((x[i]-x_0[i])*cos(theta))+((y[i]-y_0[i])*sin(theta)));
    		displacement_theta[n]+= (((y[i]-y_0[i])*cos(theta))-((x[i]-x_0[i])*sin(theta)));
		
		mod_displacement_radial[n]+= fabs((((x[i]-x_0[i])*cos(theta))+((y[i]-y_0[i])*sin(theta))));
                mod_displacement_theta[n]+= fabs((((y[i]-y_0[i])*cos(theta))-((x[i]-x_0[i])*sin(theta))));

    	}
    
  
    	end = clock();   //  Ending clock
    	cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    
 
    	free(x);
    	free(y);
    	free(x_0);
    	free(y_0);
    	free(tag);
    	load_config++;
        }
    }
    printf("%d\n",load_config);
    char displacement_field_FileName[128];
    FILE  *displacement_field; 
    sprintf(displacement_field_FileName,"displacement_field_of_%d_particles_of_%06u_ens_id_%.2lf_bin_DD.txt",TP,load_config,dr);
    displacement_field=fopen(displacement_field_FileName,"w");
    for(j=0;j<no_of_radial_bin;j++) {
    	displacement_radial[j]/=(load_config*1.e0);
    	displacement_theta[j]/=(load_config*1.e0);
	mod_displacement_radial[j]/=(load_config*1.e0);
        mod_displacement_theta[j]/=(load_config*1.e0);
	no_of_particle[j]/=(load_config*1.e0);
    	density[j]=(no_of_particle[j])/(pi*dr*dr*(((j+1)*(j+1))-(j*j)));
    	fprintf(displacement_field,"%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",(j+0.5e0)*dr,displacement_radial[j],displacement_theta[j],mod_displacement_radial[j],mod_displacement_theta[j],no_of_particle[j],density[j]);
    }
    fclose(displacement_field);
}
