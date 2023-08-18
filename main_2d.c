#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>

clock_t start, end;
double  cpu_time_used;

#define TP (2500)
#define double_TP (2500.e0)
#define dt (5e-3)
#define rate (1e-4)
#define T_steps (200000)

//...bi-desperse 2D..
#define n_A (int) (double_TP*(65.e0/100.e0))
#define n_B (TP-n_A)
int *tag;
//...................
//...LJ parameters...
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
// //.................................
double *x, *y, *rx, *ry, *fy, *fx, *vx, *vy;
double rho, length, half_length;
double T_0, T, P_0, P;
double ke, pe, te, ssum, vir, cmvx, cmvy, temp, cmx, cmy;
//.................................

//.... Neighbour list .....
int *nn, *n_list[TP];
#define neigh_width  (0.3e0)
#define cfn (2.5e0+neigh_width) // cut off for finding neighbours

//........................
//********** cell list *******
double r_c;
int n_c,tn_c;
int *nc, *cells;
int *np_c[TP];
//****************************
//.. Thermostat & Barrostat (Berendsen)..
#define tao (5.e0)
#define beta (50.e0)
double lambda, mu;
//.......................................

#include "uniform_random.h"
#include "initialize_2d.c"
#include "cell_pbc_2d.c"
#include "verlet_pbc_2d.c"
#include "force_pbc_nn_2d.c"

int main(int argc, char *argv[]) {
        // If an command line arg is not provided - error & exit
    if(argc!=2) {
        fprintf(stderr, "Ensemble ID is NOT provided - exiting.\n");
        exit(EXIT_FAILURE);
    }
// Set the ensemble
        unsigned int core_number;
        sscanf(argv[1],"%u",&core_number);

    unsigned int ens_id;
        unsigned int starting, max;
        starting = core_number*120;
        max = starting +120;
for(ens_id=starting;ens_id<max;ens_id++){

    start = clock();  // Initializing clock
    //---------------------------------------------- Open files ------------------------------------------//

   char configuration_FileName[128];
   FILE *configuration;
    sprintf(configuration_FileName,"Data/%06u_configuration.txt",ens_id);
   configuration= fopen(configuration_FileName,"w");

    uint64_t seed = (uint64_t) ens_id;
   // printf("PRNG Seed -> %lu\n", seed);
    initRanq1(seed );
    
    int i, j, t, t_step, k;
    int tag_i, tag_j;
    double xi,yi;
    double dx,dy,r2;
    //*****************
    double delta_x, delta_y, delta_r;
    int neigh_counter=0;
    double max_D=0.e0;
    //----------- Arrays -------------------------------
    x= malloc(TP * sizeof(*x));
    y= malloc(TP * sizeof(*y));
    rx= malloc(TP * sizeof(*rx));
    ry= malloc(TP * sizeof(*ry));
    fx= malloc(TP * sizeof(*fx));
    fy= malloc(TP * sizeof(*fy));
    
    vx= malloc(TP * sizeof(*vx));
    vy= malloc(TP * sizeof(*vy));
    tag= malloc(TP * sizeof(*tag));
    //###### Neighbour list#######
    nn= malloc(TP * sizeof(*nn));
    for (i=0; i<TP; i++){
        n_list[i] = (int *)malloc(TP * sizeof(int));
    }
    //-------------------------------------------------
//  ============================================================== MAIN PROGRAM ========================================================================================   
    //##################################################### Equilibriate at temp 1 ####################################################################
    initialize_2d(); // starting square lattice
    //*******************************
   r_c=2.8e0 ; //Upper cut-off to for cell list will be adjusted latter
   double rc_rem;
   rc_rem = fmod(length, r_c);
   r_c = r_c + rc_rem/((int)(length/r_c));
   rc_rem = fmod(length, r_c);
   n_c = (int)(length/r_c);
   tn_c=n_c*n_c;
   printf("%lf\t%lf\t%d\t%d\n",length,r_c,n_c,tn_c);
   nc= malloc(tn_c * sizeof(*nc));
   cells = malloc(4 * sizeof(cells));
   for (i=0; i<TP; i++){
        np_c[i] = (int *)malloc(TP * sizeof(int));
   }
  //***************************************
    T_0=1.e0;
    //************ Velocity rescale to get desired temp (can equilibriate) ******
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    temp = sqrt(T_0/T);
    for(i=0;i<TP;i++){
        vx[i]*= temp;
        vy[i]*= temp;
    }
    //****************************************************************************
    cell_pbc_2d();
    force_pbc_nn_2d();
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    ke=0.5e0*ssum;
    te=ke+pe;
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    P = (T*rho)+(vir/(2.e0*length*length));
    t_step=0;
    
    printf("%d\t%lf\t%lf\t%lf\n",t_step,pe,ke,te);
    printf("%d\t%lf\n",t_step,T);
    printf("%d\t%lf\n",t_step,P);
    
    for(t_step=1;t_step<=10000;t_step++){
        //------------------- Update velocities from force at t-1 th time---------
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
        }
        //------------------------------------------------------------------------
         //------------------ Position Update -------------------------------------
        for(i=0;i<TP;i++){
            x[i]+= vx[i]*dt;
            y[i]+= vy[i]*dt;
            //========= PBC ============
            if(x[i]>=half_length) x[i]-= length;
            else if(x[i]<-half_length)   x[i]+= length;
            if(y[i]>=half_length) y[i]-= length;
            else if(y[i]<-half_length)   y[i]+= length;
            //==========================
        
        //------------------------------------------------------------------------
            delta_x = (vx[i]*dt);
            delta_y = (vy[i]*dt);
            delta_r = (delta_x*delta_x)+(delta_y*delta_y);
            if(delta_r>max_D)   max_D= delta_r;
        }
        neigh_counter++;
        if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
            cell_pbc_2d();
            max_D=0.e0;
            neigh_counter=0;
            //printf("%d\n",t_step);
        }
        force_pbc_nn_2d(); // Force at t th time
        //------------------- Update velocities from force at t th time---------
        ssum=0.e0;
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
            ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
        }
        //------------------------------------------------------------------------
        //************** Thermostat ***************
        T = (0.5e0*ssum)/((TP-1)*1.e0);
        lambda=sqrt(1+((dt/tao)*((T_0/T)-1.e0)));
        for(i=0;i<TP;i++){
            vx[i]*=lambda;
            vy[i]*=lambda;
        }   
        //***************************************** 
    }
    printf("Initialization Done\n");
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    ke=0.5e0*ssum;
    te=ke+pe;
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    P = (T*rho)+(vir/(2.e0*length*length));
    printf("pe is %lf\tke is %lf\tte is%lf\n",pe,ke,te);
    printf("temp is %lf\n",T);
    printf("press is %lf\n",P);
    cmvx=0.e0; cmvy = 0.e0; cmx=0.e0; cmy=0.e0;
    for(i=0;i<TP;i++){
        cmvx+=vx[i];
        cmvy+=vy[i];
        cmx+=x[i];
        cmy+=y[i];
    }
    cmvx=cmvx/double_TP;
    cmvy=cmvy/double_TP;
    cmx=cmx/double_TP;
    cmy=cmy/double_TP;
    printf("CM %E\t%E\t%E\t%E\n",cmvx,cmvy,cmx,cmy);
    //######################################################################################################################################################
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ T=0.55, equlibriate $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    T_0=0.55e0;
    //************ Velocity rescale to get desired temp (can equilibriate) ******
    max_D=0.e0;
    neigh_counter=0;
    cell_pbc_2d();  // To calculate force at t=0;
    force_pbc_nn_2d();  //  force at t=0
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    temp = sqrt(T_0/T);
    for(i=0;i<TP;i++){
        vx[i]*= temp;
        vy[i]*= temp;
    }
    //======================================================
    for(t=1;t<=T_steps;t++){
        //------------------- Update velocities from force at t-1 th time---------
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
        }
        //------------------------------------------------------------------------
        //------------------ Position Update -------------------------------------
        for(i=0;i<TP;i++){
            x[i]+= vx[i]*dt;
            y[i]+= vy[i]*dt;
            //========= PBC ============
            if(x[i]>=half_length) x[i]-= length;
            else if(x[i]<-half_length)   x[i]+= length;
            if(y[i]>=half_length) y[i]-= length;
            else if(y[i]<-half_length)   y[i]+= length;
            //==========================
            
            delta_x = (vx[i]*dt);
            delta_y = (vy[i]*dt);
            delta_r = (delta_x*delta_x)+(delta_y*delta_y);
            if(delta_r>max_D)   max_D= delta_r;
        }
        neigh_counter++;
        if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
            cell_pbc_2d();
            max_D=0.e0;
            neigh_counter=0;
        }
        force_pbc_nn_2d(); // Force at t th time
        //------------------- Update velocities from force at t th time---------
        ssum=0.e0;
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
            ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
        }
        //------------------------------------------------------------------------
        //######### Thermostat ###########
        T = (0.5e0*ssum)/((TP-1)*1.e0);
        lambda=sqrt(1+((dt/tao)*((T_0/T)-1.e0)));
        for(i=0;i<TP;i++){
            vx[i]*=lambda;
            vy[i]*=lambda;
        }   
        //#################################
        if(t%10000==0){
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
    }
    printf("equilibriation done\n");
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    ke=0.5e0*ssum;
    te=ke+pe;
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    P = (T*rho)+(vir/(2.e0*length*length));
    printf("pe is %lf\tke is %lf\tte is%lf\n",pe,ke,te);
    printf("temp is %lf\n",T);
    printf("press is %lf\n",P);
    cmvx=0.e0; cmvy = 0.e0; cmx=0.e0; cmy=0.e0;
    for(i=0;i<TP;i++){
        cmvx+=vx[i];
        cmvy+=vy[i];
        cmx+=x[i];
        cmy+=y[i];
    }
    cmvx=cmvx/double_TP;
    cmvy=cmvy/double_TP;
    cmx=cmx/double_TP;
    cmy=cmy/double_TP;
    printf("CM %E\t%E\t%E\t%E\n",cmvx,cmvy,cmx,cmy);
 
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Cool to T=0.2 at a small rate (can't equilibriate anymore) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    printf("Cooling Start\n");
    max_D=0.e0;
    neigh_counter=0;
    cell_pbc_2d();  // To calculate force at t=0;
    force_pbc_nn_2d();  //  force at t=0
    while(T_0>0.2e0){
        T_0 -= (dt*rate);
        //------------------- Update velocities from force at t-1 th time---------
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
        }
        //------------------------------------------------------------------------
        //------------------ Position Update -------------------------------------
       for(i=0;i<TP;i++){
            x[i]+= vx[i]*dt;
            y[i]+= vy[i]*dt;
            //========= PBC ============
            if(x[i]>=half_length) x[i]-= length;
            else if(x[i]<-half_length)   x[i]+= length;
            if(y[i]>=half_length) y[i]-= length;
            else if(y[i]<-half_length)   y[i]+= length;
            //==========================
            delta_x = (vx[i]*dt);
            delta_y = (vy[i]*dt);
            delta_r = (delta_x*delta_x)+(delta_y*delta_y);
            if(delta_r>max_D)   max_D= delta_r;
        }
        neigh_counter++;
        if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
            cell_pbc_2d();
            max_D=0.e0;
            neigh_counter=0;
        }
        force_pbc_nn_2d(); // Force at t th time
        //------------------- Update velocities from force at t th time---------
        ssum=0.e0;
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
            ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
        }
        //------------------------------------------------------------------------
       //************** Thermostat ***************
        T = (0.5e0*ssum)/((TP-1)*1.e0);
        lambda=sqrt(1+((dt/tao)*((T_0/T)-1.e0)));
        for(i=0;i<TP;i++){
            vx[i]*=lambda;
            vy[i]*=lambda;
        }   
        //***************************************** 
    }
    printf("Cooling Done\n");
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    ke=0.5e0*ssum;
    te=ke+pe;
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    P = (T*rho)+(vir/(2.e0*length*length));
    printf("pe is %lf\tke is %lf\tte is%lf\n",pe,ke,te);
    printf("temp is %lf\n",T);
    printf("press is %lf\n",P);
    cmvx=0.e0; cmvy = 0.e0; cmx=0.e0; cmy=0.e0;
    for(i=0;i<TP;i++){
        cmvx+=vx[i];
        cmvy+=vy[i];
        cmx+=x[i];
        cmy+=y[i];
    }
    cmvx=cmvx/double_TP;
    cmvy=cmvy/double_TP;
    cmx=cmx/double_TP;
    cmy=cmy/double_TP;
    printf("CM %E\t%E\t%E\t%E\n",cmvx,cmvy,cmx,cmy);
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //**************** Fix drift **************
    cmvx=0.e0;
    cmvy=0.e0;
    for(i=0;i<TP;i++){
        cmvx+=vx[i];
        cmvy+=vy[i];
    }
    for (i=0; i<TP; i++) {
        vx[i]-= (cmvx/double_TP);
        vy[i]-= (cmvy/double_TP);
    }
    //******************************************
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NPT at t=0.2 p=0.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    printf("NPT Starts\n");
    T_0=0.2e0;
    P_0=0.e0;
    t_step=0;
    verlet_pbc_2d();
    force_pbc_nn_2d();
    max_D=0.e0;
    neigh_counter=0;
    for(t_step=1;t_step<=10000;t_step++){
        //------------------- Update velocities from force at t-1 th time---------
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
        }
        //------------------------------------------------------------------------
        //------------------ Position Update -------------------------------------
        for(i=0;i<TP;i++){
            x[i]+= vx[i]*dt;
            y[i]+= vy[i]*dt;
            //========= PBC ============
            if(x[i]>=half_length) x[i]-= length;
            else if(x[i]<-half_length)   x[i]+= length;
            if(y[i]>=half_length) y[i]-= length;
            else if(y[i]<-half_length)   y[i]+= length;
            //==========================
            delta_x = (vx[i]*dt);
            delta_y = (vy[i]*dt);
            delta_r = (delta_x*delta_x)+(delta_y*delta_y);
            if(delta_r>max_D)   max_D= delta_r;
        }
        neigh_counter++;
        if(2.e0*(sqrt(max_D)*(double)neigh_counter)>neigh_width){
            verlet_pbc_2d();
            max_D=0.e0;
            neigh_counter=0;
            //printf("%d\n",t_step);
        }
        //------------------------------------------------------------------------
        force_pbc_nn_2d(); // Force at t th time
        //------------------- Update velocities from force at t th time---------
        ssum=0.e0;
        for(i=0;i<TP;i++){
            vx[i]+= 0.5e0*fx[i]*dt;
            vy[i]+= 0.5e0*fy[i]*dt;
            ssum+= ((vx[i]*vx[i])+(vy[i]*vy[i]));
        }
        //------------------------------------------------------------------------
        T = (0.5e0*ssum)/((TP-1)*1.e0);
        P = (T*rho)+((1.e0/(length*length))*(0.5e0*vir));
        //********* THRMOSTAT & BARROSTAT ********************
        lambda=sqrt(1+((dt/tao)*((T_0/T)-1.e0)));
        mu =sqrt(1-((dt/beta)*(P_0-P)));
        for(i=0;i<TP;i++){
            vx[i]*=lambda;
            vy[i]*=lambda;
        }
        for(i=0;i<TP;i++){
            x[i]*=mu;
            y[i]*=mu;
        }
        length*=mu;
        half_length=length/2.e0;
        rho=double_TP/(length*length);
        //****************************************************
        verlet_pbc_2d();
        max_D=0.e0;
        neigh_counter=0;
        force_pbc_nn_2d();
    }
    printf("NPT Done\n");
    for(i=0;i<TP;i++){
        fprintf(configuration,"%d\t%.16lf\t%.16lf\n",tag[i],x[i],y[i]);
        fprintf(velocity,"%.16lf\t%.16lf\n",vx[i],vy[i]);
    }
    ssum=0.e0;
    for(i=0;i<TP;i++){
        ssum+=((vx[i]*vx[i])+(vy[i]*vy[i]));
    }
    ke=0.5e0*ssum;
    te=ke+pe;
    T = (0.5e0*ssum)/((TP-1)*1.e0);
    P = (T*rho)+(vir/(2.e0*length*length));
    printf("pe is %lf\tke is %lf\tte is%lf\n",pe,ke,te);
    printf("temp is %lf\n",T);
    printf("press is %lf\n",P);
    cmvx=0.e0; cmvy = 0.e0; cmx=0.e0; cmy=0.e0;
    for(i=0;i<TP;i++){
        cmvx+=vx[i];
        cmvy+=vy[i];
        cmx+=x[i];
        cmy+=y[i];
    }
    cmvx=cmvx/double_TP;
    cmvy=cmvy/double_TP;
    cmx=cmx/double_TP;
    cmy=cmy/double_TP;
    printf("CM %E\t%E\t%E\t%E\n",cmvx,cmvy,cmx,cmy);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//    ========================================================================== END ===================================================================================
    end = clock();   //  Ending clock
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    printf("Time taken: %g sec\n",cpu_time_used);
    free(x);
    free(y);
    free(rx);
    free(ry);

    free(fx);
    free(fy);

    free(vx);
    free(vy);
    free(tag);
    free(nn);
    for(i=0;i<TP;i++){
        free(n_list[i]);
        free(np_c[i]);
    }
    free(nc);
    free(cells);
    
    fclose(configuration);
    fclose(velocity);
 
}
}
