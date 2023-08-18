#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<string.h>
#ifdef __INTEL_COMPILER
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
#endif
clock_t start, end;
double  cpu_time_used;

#define gam (0.1e0)
#define dim (3)
#define X_COMP (0)
#define Y_COMP (1)
#define Z_COMP (2)
static double double_TP;
int TP,t_step;
double te,r_cm,cmvx,cmvy,cmvz,dt = 0.005e0;
double *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz;
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
int *n_list[5832]; 		//id of neighbours
int *nn; 		//number of neighbours
// FOR HESSIAN
double *upTriH;
double *eigenValues, *eigenVectors, *pr;
int dimN  ,dimNsq;

#include "verlet_3d.c"
#include "force_nn_3d.c"
int main(int argc, char *argv[]){
	int i,j,k;
    start = clock();  // Initializing clock
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

    starting= core_number*50;
    max = starting+50;
for(ens_id=starting;ens_id<max;ens_id++){
    //printf("%06u\n",ens_id);
    //--------------- Input Data Files --------------
        char dataFileName[128];
        FILE *dataFile;
        sprintf(dataFileName,"%06u_configuration.txt",ens_id);
        dataFile=fopen(dataFileName,"rb");
        if(dataFile==NULL) {
                printf("ERROR No File at :'%s'\n",dataFileName);
        }
        else{
            int type;
            double rx, ry, rz, velo_x, velo_y, velo_z;
		double dist_x, dist_y, dist_z, dist_r, dist_x_i, dist_x_j, dist_y_i, dist_y_j, dist_z_i, dist_z_j, dist_r_i, dist_r_j;
        int *Atag;
        double *Ax, *Ay, *Az;
        Ax= malloc(5832 * sizeof(*Ax));
        Ay= malloc(5832 * sizeof(*Ay));
        Az= malloc(5832 * sizeof(*Az));
        Atag= malloc(5832* sizeof(*Atag));
            fscanf(dataFile,"%*[^\n]\n");
        for(i=0;i<5832;i++){
            fscanf(dataFile,"%d %lf %lf %lf\n",&(Atag[i]),&(Ax[i]),&(Ay[i]),&(Az[i]));
        }
        fclose(dataFile);
        //---------------------- Cut out -------------------------------
        double rad = 18.e0;
        double cmx = 0.e0;
        double cmy = 0.e0;
        double cmz = 0.e0;
        for(i=0;i<5832;i++){
            cmx+=Ax[i];
            cmy+=Ay[i];
            cmz+=Az[i];
        }
        
        cmx/=(5832*1.e0);
        cmy/=(5832*1.e0);
        cmz/=(5832*1.e0);
        //printf("%E\t%E\n",cmx,cmy);
        k=0;
        for(i=0;i<5832;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_z = Az[i]-cmz;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z);
            if(dist_r<=(rad*rad)) k++;
        }
        //printf("%d\n",k);
        if(k>=4096){
        x= malloc(k * sizeof(*x));
        y= malloc(k * sizeof(*y));
        z= malloc(k * sizeof(*z));
        vx= malloc(k * sizeof(*vx));
        vy= malloc(k * sizeof(*vy));
        vz= malloc(k * sizeof(*vz));
        tag= malloc(k * sizeof(*tag));
        k=0;
        for(i=0;i<5832;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_z = Az[i]-cmz;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z);;
            if(dist_r<=(rad*rad)){
                x[k]=Ax[i];
                y[k]=Ay[i];
                z[k]=Az[i];
                tag[k]=Atag[i];
                k++;
            }
        }
        //------------------------ Sorting ------------------------------
        for(j=1;j<k;j++){

            rx=x[j];
            ry=y[j];
            rz=z[j];
            type=tag[j];

            i=j-1;

            dist_x_j = x[j]-cmx;
            dist_y_j = y[j]-cmy;
            dist_z_j = z[j]-cmz;
            dist_r_j = (dist_x_j*dist_x_j)+(dist_y_j*dist_y_j)+(dist_z_j*dist_z_j);
            
            dist_x_i = x[i]-cmx;
            dist_y_i = y[i]-cmy;
            dist_z_i = z[i]-cmz;
            dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i)+(dist_z_i*dist_z_i);

            while(i>=0 && dist_r_i>dist_r_j){
                x[i+1]=x[i];
                y[i+1]=y[i];
                z[i+1]=z[i];
                tag[i+1]=tag[i];
                i--;
                dist_x_i = x[i]-cmx;
                dist_y_i = y[i]-cmy;
                dist_z_i = z[i]-cmz;
                dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i)+(dist_z_i*dist_z_i);
            }
            x[i+1]=rx;
            y[i+1]=ry;
            z[i+1]=rz;
            tag[i+1]=type;
        }
        //--------------------------------------------------------------
        
        //printf("%d\n",k);
        TP=4096;
        double_TP = TP*1.e0;
        dimN = dim*TP;
        dimNsq = dimN*dimN;
        for(i=0;i<TP;i++){
            vx[i]=0.e0;
            vy[i]=0.e0;
            vz[i]=0.e0;
        }

    
    //-------------- Memory Allocation --------------
    fx= malloc(TP * sizeof(*fx));
    fy= malloc(TP * sizeof(*fy));
    fz= malloc(TP * sizeof(*fz));
    
    //###### Neighbour list#######
    nn= malloc(TP * sizeof(*nn));
    for (i=0; i<5832; i++){
        n_list[i] = (int *)malloc(TP * sizeof(int));
    }
    //============================ For Hessian Analysis ====================
    upTriH = malloc(dimNsq * sizeof(*upTriH));
    eigenValues = malloc(dimN*sizeof(*eigenValues));
    eigenVectors = malloc(dimNsq*sizeof(*eigenVectors));
    pr = malloc(dimN*sizeof(*pr));
    
    
  char configFileName[128];
  FILE *config;
  sprintf(configFileName,"%d_particles_min_configuration_%06u_ens_id.txt",TP,ens_id);
  config=fopen(configFileName,"w");
            //#********************************************************************************************************************** 
            double max_D, delta_x,delta_y, delta_z, delta_r;
            int neigh_counter;
            max_D=0.e0;
            neigh_counter=0;
            verlet_3d();  // To calculate force at t=0;
            force_nn_3d();  //  force at t=0
		double typicalGrad;

            ssum=0.e0;
            for(i=0;i<TP;i++){
                ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i])+(vz[i]*vz[i]));
            }
            ke=0.5e0*ssum;
            te=ke+pe;
            T = (ssum)/((TP-1)*3.e0);
            for(t_step=1;t_step<=150000;t_step++){
                //------------------- Update velocities from force at t-1 th time---------
                for(i=0;i<TP;i++){
                    vx[i]+= 0.5e0*(fx[i]-gam*vx[i])*dt;
                    vy[i]+= 0.5e0*(fy[i]-gam*vy[i])*dt;
                    vz[i]+= 0.5e0*(fz[i]-gam*vz[i])*dt;
                }
                //------------------------------------------------------------------------
                //------------------ Position Update -------------------------------------
                for(i=0;i<TP;i++){
                    x[i]+= vx[i]*dt;
                    y[i]+= vy[i]*dt;
                    z[i]+= vz[i]*dt;
                    delta_x = (vx[i]*dt);
                    delta_y = (vy[i]*dt);
                    delta_z = (vz[i]*dt);
                    delta_r = (delta_x*delta_x)+(delta_y*delta_y)+(delta_z*delta_z);
                    if(delta_r>max_D)   max_D= delta_r;
                }
                neigh_counter++;
                if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
                    verlet_3d();
                    max_D=0.e0;
                    neigh_counter=0;
                }
                force_nn_3d(); // Force at t th time
                //------------------- Update velocities from force at t th time---------
                ssum=0.e0;
                for(i=0;i<TP;i++){
                    vx[i]+= 0.5e0*(fx[i]-gam*vx[i])*dt;
                    vy[i]+= 0.5e0*(fy[i]-gam*vy[i])*dt;
                    vz[i]+= 0.5e0*(fz[i]-gam*vz[i])*dt;
                    ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i])+(vz[i]*vz[i]));
                }
                ke=0.5e0*ssum;
                te=ke+pe;
                T = (ssum)/((TP-1)*3.e0);
                if(t_step%100==0){
                    cmvx=0.e0;
                    cmvy=0.e0;
                    cmvz=0.e0;
                    for(i=0;i<TP;i++){
                        cmvx+=vx[i];
                        cmvy+=vy[i];
                        cmvz+=vz[i];
                    }
                    //... making cm velocities zero  ....
                    for (i=0; i<TP; i++) {
                        vx[i]-= (cmvx/double_TP);
                        vy[i]-= (cmvy/double_TP);
                        vz[i]-= (cmvz/double_TP);
                    }
                }
                //#*****************************************************
                //#*****************************************************
            }
           
            for (typicalGrad = 0.0,i=0;i<TP;i++)
        		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i];
    		//typicalGrad = sqrt(typicalGrad/DOUBLE_N);
    		typicalGrad = sqrt(typicalGrad);
            printf("%06u\t%E\t%E\n",ens_id,pe,typicalGrad);
            for(i=0;i<TP;i++){
                //if(nn[i]!=0){
                    fprintf(config,"%d\t%.16lf\t%.16lf\t%.16lf\n",tag[i],x[i],y[i],z[i]);
                //}
            }
            fclose(config);
    
    
    //+++++++++++++++++++++++++++++++++++++++++ CALCULATE HESSIAN ELEMENTS ++++++++++++++++++++++++++++++++++++++++++
    int indx, indy, indz;
    double xi, yi, zi, dx, dy, dz, r2, s2, s4, s6, s12, ri2, ri8, ri14, ri10, ri16, phi, si;
	int tag_i, tag_j;
    double hesXY, hesXZ, hesYZ;

    for(i=0;i<dimNsq;i++){
        upTriH[i]=0.e0;
    }
    //................................
        for(i=0;i<TP;i++){
            indx = dimN*((3*i)+X_COMP);
            indy = dimN*((3*i)+Y_COMP);
            indz = dimN*((3*i)+Z_COMP);
            xi = x[i];
            yi = y[i];
            zi = z[i];
            tag_i=tag[i];
            
            for(j=i+1;j<TP;j++){
                tag_j = tag[j];
                dx = x[j]-xi;
                dy = y[j]-yi;
                dz = z[j]-zi;

                r2=(dx*dx)+(dy*dy)+(dz*dz);
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
                    hesXZ = -(phi*dx*dz);
                    hesYZ = -(phi*dy*dz);
                    //--------- XX, XY, XZ -------------------------------------
                    upTriH[indx+(3*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(3*j)+Y_COMP] += hesXY;
                    upTriH[indx+(3*j)+Z_COMP] += hesXZ;
                    //----------- ends here --------------------------------
                    //--------- YX, YY, YZ -------------------------------------
                    upTriH[indy+(3*j)+X_COMP] += hesXY;
                    upTriH[indy+(3*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    upTriH[indy+(3*j)+Z_COMP] += hesYZ;
                    //----------- ends here --------------------------------
                    //--------- ZX, ZY, ZZ -------------------------------------
                    upTriH[indz+(3*j)+X_COMP] += hesXZ;
                    upTriH[indz+(3*j)+Y_COMP] += hesYZ;
                    upTriH[indz+(3*j)+Z_COMP] += ((-phi*dz*dz)-si);
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
                    hesXZ = -(phi*dx*dz);
                    hesYZ = -(phi*dy*dz);
                    //--------- XX, XY, XZ -------------------------------------
                    upTriH[indx+(3*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(3*j)+Y_COMP] += hesXY;
                    upTriH[indx+(3*j)+Z_COMP] += hesXZ;
                    //----------- ends here --------------------------------
                    //--------- YX, YY, YZ -------------------------------------
                    upTriH[indy+(3*j)+X_COMP] += hesXY;
                    upTriH[indy+(3*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    upTriH[indy+(3*j)+Z_COMP] += hesYZ;
                    //----------- ends here --------------------------------
                    //--------- ZX, ZY, ZZ -------------------------------------
                    upTriH[indz+(3*j)+X_COMP] += hesXZ;
                    upTriH[indz+(3*j)+Y_COMP] += hesYZ;
                    upTriH[indz+(3*j)+Z_COMP] += ((-phi*dz*dz)-si);
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
                    hesXZ = -(phi*dx*dz);
                    hesYZ = -(phi*dy*dz);
                    //--------- XX, XY, XZ -------------------------------------
                    upTriH[indx+(3*j)+X_COMP] += ((-phi*dx*dx)-si);
                    upTriH[indx+(3*j)+Y_COMP] += hesXY;
                    upTriH[indx+(3*j)+Z_COMP] += hesXZ;
                    //----------- ends here --------------------------------
                    //--------- YX, YY, YZ -------------------------------------
                    upTriH[indy+(3*j)+X_COMP] += hesXY;
                    upTriH[indy+(3*j)+Y_COMP] += ((-phi*dy*dy)-si);
                    upTriH[indy+(3*j)+Z_COMP] += hesYZ;
                    //----------- ends here --------------------------------
                    //--------- ZX, ZY, ZZ -------------------------------------
                    upTriH[indz+(3*j)+X_COMP] += hesXZ;
                    upTriH[indz+(3*j)+Y_COMP] += hesYZ;
                    upTriH[indz+(3*j)+Z_COMP] += ((-phi*dz*dz)-si);
                    //----------- ends here --------------------------------
                }
            }

        }

            int indi,indi1,indi2N1,indi2,indi2N2,indi3N1;
            for (i=0; i<TP; i++){
                indi = (dimN+1)*(3*i);
                indi1 = indi + 1;
                indi2 = indi1 +1;
                indi2N1 = indi+dimN+1;
                indi2N2 = indi2N1+1;
                indi3N1 = indi+(2*(dimN+1));
                for (j=1; j<(TP-i); j++){
                    upTriH[indi] -= upTriH[indi+(3*j)];
                    upTriH[indi1] -= upTriH[indi1+(3*j)];
                    upTriH[indi2] -= upTriH[indi2+(3*j)];
                    upTriH[indi2N1] -= upTriH[indi2N1+(3*j)];
                    upTriH[indi2N2] -= upTriH[indi2N2+(3*j)];
                    upTriH[indi3N1] -= upTriH[indi3N1+(3*j)];
                }
                for (j=0; j<i; j++){
                    upTriH[indi] -= upTriH[indi-(dimN*3*(j+1))];
                    upTriH[indi1] -= upTriH[indi1-(dimN*3*(j+1))];
                    upTriH[indi2] -= upTriH[indi2-(dimN*3*(j+1))];
                    upTriH[indi2N1] -= upTriH[indi2N1-(dimN*3*(j+1))];
                    upTriH[indi2N2] -= upTriH[indi2N2-(dimN*3*(j+1))];
                    upTriH[indi3N1] -= upTriH[indi3N1-(dimN*3*(j+1))];
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
            e2 = (eigenVectors[(i*dimN)+dim*j]*eigenVectors[(i*dimN)+dim*j])+(eigenVectors[(i*dimN)+(dim*j)+1]*eigenVectors[(i*dimN)+(dim*j)+1])+(eigenVectors[(i*dimN)+(dim*j)+2]*eigenVectors[(i*dimN)+(dim*j)+2]);
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
            fprintf(first_non_zero_eigenvector,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",tag[j],x[j],y[j],z[j],eigenVectors[(i*dimN)+dim*j],eigenVectors[(i*dimN)+(dim*j)+1],eigenVectors[(i*dimN)+(dim*j)+2]);
       		}
		break;
    	}
    }
    fclose(first_non_zero_eigenvector);
    
        free(Ax);
        free(Ay);
        free(Az);
        free(Atag);
    free(x);
    free(y);
    free(z);
    free(fx);
    free(fy);
    free(fz);
    free(vx);
    free(vy);
    free(vz);
    free(tag);
    free(nn);
    for(i=0;i<5832;i++) free(n_list[i]);
    free(eigenValues);
    free(eigenVectors);
    free(pr);
    free(upTriH);
        }
    }
    }
}
