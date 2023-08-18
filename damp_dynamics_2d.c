#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>
/*#ifdef __INTEL_COMPILER
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
#endif*/
clock_t start, end;
double  cpu_time_used;

#define gam (0.1e0)
#define dim (2)
#define X_COMP (0)
#define Y_COMP (1)
static double double_TP;
int TP,t_step;
double te,r_cm,cmvx,cmvy,dt = 0.005e0;
double *x, *y, *vx, *vy, *fx, *fy;
int *tag;
//============ MODEL ===================
#define sig_AA (1.e0)  // LJ PARAMETER
#define eps_AA (1.e0)  // LJ PARAMETER
#define sig_AB (0.8e0*sig_AA)  // LJ PARAMETER
#define eps_AB (1.5e0*eps_AA)  // LJ PARAMETER
#define sig_BB (0.88e0*sig_AA)  // LJ PARAMETER
#define eps_BB (0.5e0*eps_AA)  // LJ PARAMETER
#define r_cut_AA (2.5e0*sig_AA)
#define r_cut_AB (2.5e0*sig_AB)
#define r_cut_BB (2.5e0*sig_BB)
#define c_0 (0.040490237952e0)
#define c_2 (-0.00970155098112e0)
#define c_4 (0.0006201261686784e0)
//...................

double pe,vir,ssum,ke,T;
//============ Neighbour List ==========
#define neigh_width  (0.3e0)
#define cfn (2.5e0+neigh_width) // cut off for finding neighbours
int *n_list[2500]; 		//id of neighbours
int *nn; 		//number of neighbours
// FOR HESSIAN
double *upTriH;
double *eigenValues, *eigenVectors, *pr;
int dimN  ,dimNsq;

#include "verlet_2d.c"
#include "force_nn_2d.c"
int main(int argc, char *argv[]){
	int i,j,k;

    //If an command line arg is not provided - error & exit
    if(argc!=2) {
    fprintf(stderr, "Ensemble ID is NOT provided - exiting.\n");
    exit(EXIT_FAILURE);
    }
    
    // Set the ensemble
    unsigned int core_number;
    sscanf(argv[1],"%u",&core_number);
    
    unsigned int ens_id;
    unsigned int starting, max;

    starting= core_number*3125;
    max = starting+3125;
for(ens_id=starting;ens_id<max;ens_id++){
    //printf("%06u\n",ens_id);
    //--------------- Input Data Files --------------
        char dataFileName[128];
        FILE *dataFile;
        sprintf(dataFileName,"../Data/%06u_configuration.txt",ens_id);
        dataFile=fopen(dataFileName,"rb");
        if(dataFile==NULL) {
                printf("ERROR No File at :'%s'\n",dataFileName);
        }
        else{
            int type;
            double rx, ry;
		double dist_x, dist_y, dist_r, dist_x_i, dist_x_j, dist_y_i, dist_y_j, dist_r_i, dist_r_j;
        int *Atag;
        double *Ax, *Ay;
        Ax= malloc(2500 * sizeof(*Ax));
        Ay= malloc(2500 * sizeof(*Ay));
        Atag= malloc(2500* sizeof(*Atag));
        for(i=0;i<2500;i++){
            fscanf(dataFile,"%d %lf %lf",&(Atag[i]),&(Ax[i]),&(Ay[i]));
        }
        fclose(dataFile);
        
        //---------------------- Cut out -------------------------------
        double rad = 25.e0;
        double cmx = 0.e0;
        double cmy = 0.e0;
        for(i=0;i<2500;i++){
            cmx+=Ax[i];
            cmy+=Ay[i];
        }
        
        cmx/=(2500*1.e0);
        cmy/=(2500*1.e0);
        //printf("%E\t%E\n",cmx,cmy);
        k=0;
        for(i=0;i<2500;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y);
            if(dist_r<=(rad*rad)) k++;
        }
        //printf("%d\n",k);
        if(k>=900){
        x= malloc(k * sizeof(*x));
        y= malloc(k * sizeof(*y));
        tag= malloc(k * sizeof(*tag));
        k=0;
        for(i=0;i<2500;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y);
            if(dist_r<=(rad*rad)){
                x[k]=Ax[i];
                y[k]=Ay[i];
                tag[k]=Atag[i];
                k++;
            }
        }
        //------------------------ Sorting ------------------------------
        for(j=1;j<k;j++){

            rx=x[j];
            ry=y[j];
            type=tag[j];

            i=j-1;

            dist_x_j = x[j]-cmx;
            dist_y_j = y[j]-cmy;
            dist_r_j = (dist_x_j*dist_x_j)+(dist_y_j*dist_y_j);
            
            dist_x_i = x[i]-cmx;
            dist_y_i = y[i]-cmy;
            dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i);

            while(i>=0 && dist_r_i>dist_r_j){
                x[i+1]=x[i];
                y[i+1]=y[i];
                tag[i+1]=tag[i];
                i--;
                dist_x_i = x[i]-cmx;
                dist_y_i = y[i]-cmy;
                dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i);
            }
            x[i+1]=rx;
            y[i+1]=ry;
            tag[i+1]=type;
        }
        //--------------------------------------------------------------
        
        //printf("%d\n",k);
        TP=900;
        double_TP = TP*1.e0;
        dimN = dim*TP;
        dimNsq = dimN*dimN;
/*	char liquidFileName[128];
	FILE *liquid;
	sprintf(liquidFileName,"../liquid_configs_900/%06u_liquid_configuration.txt",ens_id);
	liquid = fopen(liquidFileName,"w");
	for(i=0;i<TP;i++){
		fprintf(liquid,"%d\t%.16lf\t%.16lf\n",tag[i],x[i],y[i]);
	}	
	fclose(liquid);*/
    
    //-------------- Memory Allocation --------------
    fx= malloc(TP * sizeof(*fx));
    fy= malloc(TP * sizeof(*fy));
    vx= malloc(TP * sizeof(*vx));
    vy= malloc(TP * sizeof(*vy));
    
    //###### Neighbour list#######
    nn= malloc(TP * sizeof(*nn));
    for (i=0; i<2500; i++){
        n_list[i] = (int *)malloc(TP * sizeof(int));
    }
    
    //============================ For Hessian Analysis ====================
    
    upTriH = malloc(dimNsq * sizeof(*upTriH));
    eigenValues = malloc(dimN*sizeof(*eigenValues));
    eigenVectors = malloc(dimNsq*sizeof(*eigenVectors));
    pr = malloc(dimN*sizeof(*pr));
    
    for(i=0;i<TP;i++){
	vx[i]=0.e0;
	vy[i]=0.e0;
    }
  char configFileName[128];
        FILE *config;
        sprintf(configFileName,"min_configs_%d_TP_%06u_ens_id.txt",TP,ens_id);
        config=fopen(configFileName,"w");
    
    //#********************************************** fixing things *******************************************************
            cmx = 0.e0;
            cmy = 0.e0;
            for(i=0;i<TP;i++){
                cmx+=x[i];
                cmy+=y[i];
            }
            cmx/=double_TP;
            cmy/=double_TP;
            r_cm = sqrt((cmx*cmx)+(cmy*cmy));
            //#********************************************************************************************************************** 
            double max_D, delta_x,delta_y, delta_r;
            int neigh_counter;
            max_D=0.e0;
            neigh_counter=0;
            verlet_2d();  // To calculate force at t=0;
            force_nn_2d();  //  force at t=0
		double typicalGrad;
           /* for (typicalGrad = 0.0,i=0;i<TP;i++)
        		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i];
    		//typicalGrad = sqrt(typicalGrad/DOUBLE_N);
    		typicalGrad = sqrt(typicalGrad);
            printf("%E\t%E\n",pe,typicalGrad);*/

            ssum=0.e0;
            for(i=0;i<TP;i++){
                ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
            }
            ke=0.5e0*ssum;
            te=ke+pe;
            T = (0.5e0*ssum)/((TP-1)*1.e0);
            for(t_step=1;t_step<=100000;t_step++){
                //------------------- Update velocities from force at t-1 th time---------
                for(i=0;i<TP;i++){
                    vx[i]+= 0.5e0*(fx[i]-gam*vx[i])*dt;
                    vy[i]+= 0.5e0*(fy[i]-gam*vy[i])*dt;
                }
                //------------------------------------------------------------------------
                //------------------ Position Update -------------------------------------
                for(i=0;i<TP;i++){
                    x[i]+= vx[i]*dt;
                    y[i]+= vy[i]*dt;
                    delta_x = (vx[i]*dt);
                    delta_y = (vy[i]*dt);
                    delta_r = (delta_x*delta_x)+(delta_y*delta_y);
                    if(delta_r>max_D)   max_D= delta_r;
                }
                neigh_counter++;
                if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
                    verlet_2d();
                    max_D=0.e0;
                    neigh_counter=0;
                }
                force_nn_2d(); // Force at t th time
                //------------------- Update velocities from force at t th time---------
                ssum=0.e0;
                for(i=0;i<TP;i++){
                    vx[i]+= 0.5e0*(fx[i]-gam*vx[i])*dt;
                    vy[i]+= 0.5e0*(fy[i]-gam*vy[i])*dt;
                    ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
                }
                ke=0.5e0*ssum;
                te=ke+pe;
                T = (0.5e0*ssum)/((TP-1)*1.e0);
                if(t_step%100==0){
                    cmvx=0.e0;
                    cmvy=0.e0;
                    for(i=0;i<TP;i++){
                        cmvx+=vx[i];
                        cmvy+=vy[i];
                    }
                    //... making cm velocities zero  ....
                    for (i=0; i<TP; i++) {
                        vx[i]-= (cmvx/double_TP);
                        vy[i]-= (cmvy/double_TP);
                    }
                }
                //#*****************************************************
                //#*****************************************************
            }
           
            for (typicalGrad = 0.0,i=0;i<TP;i++)
        		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i];
    		//typicalGrad = sqrt(typicalGrad/DOUBLE_N);
    		typicalGrad = sqrt(typicalGrad);
            printf("%06u\t%E\t%E\n",ens_id,pe,typicalGrad);
            for(i=0;i<TP;i++){
                //if(nn[i]!=0){
                    fprintf(config,"%d\t%.16lf\t%.16lf\n",tag[i],x[i],y[i]);
                //}
            }
            fclose(config);
    
    
    //+++++++++++++++++++++++++++++++++++++++++ CALCULATE HESSIAN ELEMENTS ++++++++++++++++++++++++++++++++++++++++++
    int indx, indy;
    double xi, yi, dx, dy, r2, s2, s4, s6, s12, ri2, ri8, ri14, ri10, ri16, phi, si;
	int tag_i, tag_j;
    double hesXY;

    for(i=0;i<dimNsq;i++){
        upTriH[i]=0.e0;
    }
    //................................
        for(i=0;i<TP;i++){
            indx = dimN*((2*i)+X_COMP);
            indy = dimN*((2*i)+Y_COMP);
            xi = x[i];
            yi = y[i];
            tag_i=tag[i];
            
            for(j=i+1;j<TP;j++){
                tag_j = tag[j];
                dx = x[j]-xi;
                dy = y[j]-yi;

                r2=(dx*dx)+(dy*dy);
                if(tag_i==1 && tag_j==1 && r2<(r_cut_AA*r_cut_AA)){
                    s2= sig_AA*sig_AA;
                    ri2 = 1.e0/r2;
                    s4 = s2*s2;
                    s6 = s4*s2;
                    s12= s6*s6;
                    ri8 = ri2*ri2*ri2*ri2;
                    ri10 = ri8*ri2;
                    ri14 = ri10*ri2*ri2;
                    ri16 = ri8*ri8;
                    phi = 4.e0*eps_AA*((168.e0*(s12*ri16))-(48.e0*(s6*ri10))+(8.e0*(c_4/s4)));
                    si = 4.e0*eps_AA*((-12.e0*(s12*ri14))+(6.e0*(s6*ri8))+(2.e0*(c_2/s2))+(4.e0*c_4*(r2/s4)));

                    hesXY = - (phi*dx*dy);
                    //--------- XX, XY -------------------------------------
                    upTriH[indx+(2*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(2*j)+Y_COMP] += hesXY;
                    //----------- ends here --------------------------------
                    //--------- YX, YY -------------------------------------
                    upTriH[indy+(2*j)+X_COMP] += hesXY;
                    upTriH[indy+(2*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    //----------- ends here --------------------------------
                }
                if(tag_i==0 && tag_j==0 && r2<(r_cut_BB*r_cut_BB)){
                    s2= sig_BB*sig_BB;
                    ri2 = 1.e0/r2;
                    s4 = s2*s2;
                    s6 = s4*s2;
                    s12= s6*s6;
                    ri8 = ri2*ri2*ri2*ri2;
                    ri10 = ri8*ri2;
                    ri14 = ri10*ri2*ri2;
                    ri16 = ri8*ri8;
                    phi = 4.e0*eps_BB*((168.e0*(s12*ri16))-(48.e0*(s6*ri10))+(8.e0*(c_4/s4)));
                    si = 4.e0*eps_BB*((-12.e0*(s12*ri14))+(6.e0*(s6*ri8))+(2.e0*(c_2/s2))+(4.e0*c_4*(r2/s4)));

                    hesXY = - (phi*dx*dy);
                    //--------- XX, XY -------------------------------------
                    upTriH[indx+(2*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(2*j)+Y_COMP] += hesXY;
                    //----------- ends here --------------------------------
                    //--------- YX, YY -------------------------------------
                    upTriH[indy+(2*j)+X_COMP] += hesXY;
                    upTriH[indy+(2*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    //----------- ends here --------------------------------
                }
                if(tag_i!=tag_j && r2<(r_cut_AB*r_cut_AB)){
                    s2= sig_AB*sig_AB;
                    ri2 = 1.e0/r2;
                    s4 = s2*s2;
                    s6 = s4*s2;
                    s12= s6*s6;
                    ri8 = ri2*ri2*ri2*ri2;
                    ri10 = ri8*ri2;
                    ri14 = ri10*ri2*ri2;
                    ri16 = ri8*ri8;
                    phi = 4.e0*eps_AB*((168.e0*(s12*ri16))-(48.e0*(s6*ri10))+(8.e0*(c_4/s4)));
                    si = 4.e0*eps_AB*((-12.e0*(s12*ri14))+(6.e0*(s6*ri8))+(2.e0*(c_2/s2))+(4.e0*c_4*(r2/s4)));

                    hesXY = - (phi*dx*dy);
                    //--------- XX, XY -------------------------------------
                    upTriH[indx+(2*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(2*j)+Y_COMP] += hesXY;
                    //----------- ends here --------------------------------
                    //--------- YX, YY -------------------------------------
                    upTriH[indy+(2*j)+X_COMP] += hesXY;
                    upTriH[indy+(2*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    //----------- ends here --------------------------------
                }
            }

        }

        int indi,indi1,indi2N1;
        for (i=0; i<TP; i++){
            indi = (dimN+1)*(2*i);
            indi1 = indi + 1;
            indi2N1 = indi+dimN+1;
            for (j=1; j<(TP-i); j++){
                upTriH[indi] -= upTriH[indi+(2*j)];
                upTriH[indi1] -= upTriH[indi1+(2*j)];
                upTriH[indi2N1] -= upTriH[indi2N1+(2*j)];
            }
            for (j=0; j<i; j++){
                upTriH[indi] -= upTriH[indi-(dimN*2*(j+1))];
                upTriH[indi1] -= upTriH[indi1-(dimN*2*(j+1))];
                upTriH[indi2N1] -= upTriH[indi2N1-(dimN*2*(j+1))];
            }
        }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ Diagonalize Hessian Matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    double  value;
                
    int numVals = dimN;
    char jobz = 'V'; //for getting eigenvectors too, use V. If only eigenvalues required, use 'N'
    char range = 'A'; //for getting only some eigenvalues/vectors, use 'I'; for all, use 'A'
    char uplo = 'L'; //for considering the upper 'U' or lower 'L' triangular part.
    int dof;
    int returnStat;
    int *found;
    int fromEigenValue = 1;
    
    int uptoEigenValue = fromEigenValue + numVals - 1;
    double dummy1 = 0.0, dummy2 = 0.0;
    int *isuppz;
    int sizeOfMatrix= dimN;
    dof = sizeOfMatrix;
    found = malloc(sizeof(*found));
    if ( (jobz == 'V') || (range == 'A') ){
        isuppz = malloc(2*dof*sizeof(*isuppz));
    } else {
        isuppz = malloc(sizeof(*isuppz));
    }
    static const double TOL             = 1e-20;
    returnStat = LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplo, dof, upTriH, dof, dummy1, dummy2, fromEigenValue, uptoEigenValue, TOL, found, eigenValues, eigenVectors, dof, isuppz);
    int i,j;
    double e2, e4;
    double nume = 0.e0, deno = 0.e0;
    for(i=0;i<dimN;i++){
        for(j=0;j<TP;j++){
            e2 = (eigenVectors[(i*dimN)+dim*j]*eigenVectors[(i*dimN)+dim*j])+(eigenVectors[(i*dimN)+(dim*j)+1]*eigenVectors[(i*dimN)+(dim*j)+1]);
            e4 = e2*e2;
            nume+=e2;
            deno+=e4;
        }
        pr[i] = (1.e0)/(TP*deno);
    }
    free(found);
    free(isuppz);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    char eigenvalues_FileName[128];
    FILE  *eigenvalues; 
    sprintf(eigenvalues_FileName,"eigenvalues_of_%d_particles_%06u_ens_id.txt",TP,ens_id);
    eigenvalues=fopen(eigenvalues_FileName,"w");
    for(k=0;k<dimN;k++){
        fprintf(eigenvalues,"%E\t%E\t%E\n",eigenValues[k],sqrt(eigenValues[k]),pr[k]);
    }
    fclose(eigenvalues);

    char first_non_zero_eigenvector_FileName[128];
    FILE  *first_non_zero_eigenvector; 
    sprintf(first_non_zero_eigenvector_FileName,"first_non_zero_mode_of_%d_particles_%06u_ens_id.txt",TP,ens_id);
    first_non_zero_eigenvector=fopen(first_non_zero_eigenvector_FileName,"w");       
    for(i=0;i<dimN;i++){
	if(eigenValues[i]>1e-5){
		for(j=0;j<TP;j++){
    			fprintf(first_non_zero_eigenvector,"%.16lf\t%.16lf\t%.16lf\t%.16lf\n",x[j],y[j],eigenVectors[(i*dimN)+(dim*j)],eigenVectors[(i*dimN)+(dim*j)+1]);
       		}
		break;
    	}
    }
    fclose(first_non_zero_eigenvector);
    
    
            free(Ax);
            free(Ay);
            free(Atag);
    free(x);
    free(y);
    free(fx);
    free(fy);
    free(vx);
    free(vy);
    free(tag);
    free(nn);
    for(i=0;i<2500;i++) free(n_list[i]);
    free(eigenValues);
    free(eigenVectors);
    free(pr);
    free(upTriH);
    
    }
    }
    }
}
