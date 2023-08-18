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

enum {
    X_COMP              = 0,
    Y_COMP              = 1,
    DIM                 = 2,
};
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
//-------------------- MINIMIZER SHIT ------------------------------------------------
enum {
    MAX_ITERATIONS          = 500000,
    LINMIN_MAX_ITERATIONS   = 120,
    MAX_BASIS               = 1500
};
static const double TOL             = 1e-20;
// factors for growing and shrinking the interval- don't change
static const double LINMIN_G1       = 2.0;
static const double LINMIN_G2       = 1.25;
static const double LINMIN_G3       = 0.5;
static const double LAST_X_DEFAULT  = 0.0001;
static const double LANCZOS_TOL     = 1.0e-6;
static int restart;		// whether to restart optimizer - fresh cg directions 
static double lastx;		// keeps track of typical step length 
static double gtyp;		// stores the rms gradient for linmin 
double *p, *q, *g, *h, *mosesX, *mosesY, *scratch;
//-------------------------- GLOBAL VARIABLES -----------------------------------------
int N, dimN  ,dimNsq;
double DOUBLE_N;
static double pe;			    /*	potential energy		*/
//static double PE0;			    /*	initial potential energy		*/
static double virial;
static double pressure;
static double *rx;			/*	x component of position		*/
static double *ry;			/*	y component of position		*/
static double *fx;			/*	x component of force		*/
static double *fy;			/*	y component of force		*/
static double typicalGrad;
static int *type;           /*  particle tag                */
//--------------------- FOR NEBZ LIST -------------------------------------------------
#define neigh_width  (0.3e0)
#define cfn (2.5e0+neigh_width) // cut off for finding neighbours
int *nebz[2500]; 		//maximum we give here MAX_NEBZ nebz.
int *numOfNebz; 		//number of elements in nebz[N][MAX_NEBZ]
static double maxD,listCutOffSqrd; 	//for updating nebz list.
static int nebListCounter; 		//for counting how many steps have took place for nebz list updating.
//------------------------------------------------------------------------------------

double *upTriH;
double *eigenValues, *eigenVectors, *pr;

static unsigned numOfConfigs;

static void updateNebzLists(){
    int i,j;
    double dx,dy,r2;
    double rxi,ryi;

    nebListCounter = 0;
    maxD = 0.0;
    //set all neb lists to zero
    for (i=0; i<N; i++)
        numOfNebz[i] = 0;

    for (i=0;i<N-1;i++){
        rxi = rx[i];
        ryi = ry[i];
        for (j=i+1; j<N; j++){
            dx = rx[j] - rxi;
            dy = ry[j] - ryi;

            r2 = ( dx*dx + dy*dy );
            if (r2 < cfn*cfn){
                nebz[i][numOfNebz[i]] = j;
                nebz[j][numOfNebz[j]] = i;
                numOfNebz[i]++;
                numOfNebz[j]++;
            }
        }
    }
}

static void calculateForces(){
    int i,j,m,k,type_I, type_J;
    double dx,dy,rxi,ryi;
    double r2,s2,r2i,r6i;

    double grad;
    double sum1 = 0.0;
    double sum2 = 0.0;
    //first set 'em all to zero
    pe = 0.0;
    virial = 0.0;
    for (i=0;i<N;i++){
        fx[i] = 0.0;
        fy[i] = 0.0;
    }
    for (i=0;i<N-1;i++){
        rxi = rx[i];
        ryi = ry[i];

        type_I = type[i];
        m = numOfNebz[i];
        for (k=0;k<m;k++){
            j = nebz[i][k];
            type_J = type[j];
            if (j > i){
                dx = rx[j] - rxi;
                dy = ry[j] - ryi;

                r2 = ( dx*dx + dy*dy );
                if(type_I==1 && type_J==1 && r2<(r_cut_AA*r_cut_AA)){
                    s2= sig_AA*sig_AA;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AA*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad = 4.e0*eps_AA*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad*dx;
                    sum2 = grad*dy;

                    fx[i] -= sum1;
                    fy[i] -= sum2;


                    fx[j] += sum1;
                    fy[j] += sum2;
                }
                if(type_I==0 && type_J==0 && r2<(r_cut_BB*r_cut_BB)){
                    s2= sig_BB*sig_BB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_BB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad = 4.e0*eps_BB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad*dx;
                    sum2 = grad*dy;

                    fx[i] -= sum1;
                    fy[i] -= sum2;


                    fx[j] += sum1;
                    fy[j] += sum2;
                }
                if(type_I!=type_J && r2<(r_cut_AB*r_cut_AB)){
                    s2= sig_AB*sig_AB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad = 4.e0*eps_AB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad*dx;
                    sum2 = grad*dy;

                    fx[i] -= sum1;
                    fy[i] -= sum2;


                    fx[j] += sum1;
                    fy[j] += sum2;
                }
            }
        }
    }
    for (typicalGrad = 0.0,i=0;i<N;i++)
        typicalGrad += fx[i]*fx[i] + fy[i]*fy[i];
    //typicalGrad = sqrt(typicalGrad/DOUBLE_N);
    typicalGrad = sqrt(typicalGrad);
    return;
}


static void grad( double *point, double *res ){
    int i,k,j,type_I,type_J,m;
    double rxi, ryi, dx, dy;
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    double r2,s2,r2i,r6i;
    double grad1;
    pe = 0.0;

    //set to zero
    for ( i=0; i<dimN; i++ )
        res[i] = 0.0;
    // -------------------------------------  force calculation   -------------------------------/
    for (i=0;i<N-1;i++){
        rxi = point[ DIM*i + X_COMP ];
        ryi = point[ DIM*i + Y_COMP ];
        type_I = type[i];
        m = numOfNebz[i];
        for ( k=0; k<m; k++ ){
            j = nebz[i][k];
            type_J=type[j];
            if (j > i){
                dx = point[ DIM*j + X_COMP ] - rxi;
                dy = point[ DIM*j + Y_COMP ] - ryi;
                r2 = ( dx*dx + dy*dy );

                if(type_I==1 && type_J==1 && r2<(r_cut_AA*r_cut_AA)){
                    s2= sig_AA*sig_AA;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AA*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad1 = 4.e0*eps_AA*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad1*dx;
                    sum2 = grad1*dy;

                    res[(DIM*i) + X_COMP] += sum1; 	//the sign is always opposite of force
                    res[(DIM*i) + Y_COMP] += sum2;

                    res[(DIM*j) + X_COMP] -= sum1;
                    res[(DIM*j) + Y_COMP] -= sum2;
                }
                if(type_I==0 && type_J==0 && r2<(r_cut_BB*r_cut_BB)){
                    s2= sig_BB*sig_BB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_BB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad1 = 4.e0*eps_BB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad1*dx;
                    sum2 = grad1*dy;

                    res[(DIM*i) + X_COMP] += sum1;     //the sign is always opposite of force
                    res[(DIM*i) + Y_COMP] += sum2;

                    res[(DIM*j) + X_COMP] -= sum1;
                    res[(DIM*j) + Y_COMP] -= sum2;
                }
                if(type_I!=type_J && r2<(r_cut_AB*r_cut_AB)){
                    s2= sig_AB*sig_AB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    grad1 = 4.e0*eps_AB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    sum1 = grad1*dx;
                    sum2 = grad1*dy;

                    res[(DIM*i) + X_COMP] += sum1;     //the sign is always opposite of force
                    res[(DIM*i) + Y_COMP] += sum2;

                    res[(DIM*j) + X_COMP] -= sum1;
                    res[(DIM*j) + Y_COMP] -= sum2;
                }
            }
        }
    }
}

//----------------------------------------------/
static void fix(double *point){		  ////fix the boundary condition in real space coordinate..............
    int i, xIndex, yIndex;
    double tmpd, dx, dy, halfLength;
    maxD = 0.0;
    for (i=0;i<N;i++){
        xIndex = DIM*i + X_COMP;
        yIndex = xIndex + Y_COMP;

        dx = point[xIndex] - rx[i];
        dy = point[yIndex] - ry[i];
        // mess due to LEES EDWARDS periodic boundary conditions
        tmpd = dx*dx + dy*dy;
        if ( tmpd > maxD )
            maxD = tmpd;

    }
    if ( 2.0*sqrt(maxD) > neigh_width ){
        //printf("Cell list update\n");
        for (i=0;i<N;i++){
            rx[i] = point[DIM*i + X_COMP];
            ry[i] = point[DIM*i + Y_COMP];
        }
        updateNebzLists();
    }
    //that's it.
    return;
}
static void initializeMinimizer(int start){
    int i;
    if ( start==0 )
        lastx = LAST_X_DEFAULT;
    for (i=0;i<dimN;i++){
        if ( restart !=2 )
            g[i] = -q[i];
        q[i] = h[i] = g[i];
    }
    restart = 0;
    return;
}

static double prod(double x, double *gradScratch){		// used inside the linmin function change the position p[i] nth to p[i] n+1 th and calculate
    int i;
    double s;
    for (i=0;i<dimN;i++)
        scratch[i] = p[i] + x*q[i];
    fix(scratch);
    grad(scratch, gradScratch);
    for (s = 0.0, i=0; i<dimN; i++)
        s += gradScratch[i]*q[i];
    return s;
}

static double linmin(){
    double x,y,s,t,m,tmpd,step;
    double *gx,*gy,*gUnused;
    int it,i;

    gx = mosesX; gy = mosesY;

    x = lastx/gtyp;
    s = prod(x, gx);
    //printf("%lf\n",s);			//result goes into gx
    it = 0;
    if ( s<0.0 ){
        do{
            y = x*LINMIN_G1;
            t = prod(y, gy);
            if ( t >= 0.0 ) break ;
            x = y; s = t; gUnused = gx; gx = gy; gy = gUnused;
            it++;
        }while (it < LINMIN_MAX_ITERATIONS);
    }
    else if ( s>0 ){
        do{
            y = x*LINMIN_G3;
            t = prod(y, gy);
            if ( t <= 0.0 ) break ;
            x = y ; s = t ; gUnused = gx ; gx = gy ; gy = gUnused;
            it++;
        }while (it < LINMIN_MAX_ITERATIONS);
    }
    else{		// hole in one s = 0.0 
        t = 1.0; y = x;
    }

    if ( it >= LINMIN_MAX_ITERATIONS){
        printf("Warning! linmin overran, ");
        tmpd = prod(0.0, gy);
        if ( tmpd > 0.0 )
            restart = 1;
    }

    if ( s < 0.0 ) s = - s;
    if ( t < 0.0 ) t = - t;
    m = ( s + t );
    s /= m; t /= m;

    m =  s * y + t * x;
    // evaluate the step length, not that it necessarily means anything 
    for ( step=0.0, i=0; i<dimN; i++ ){
        tmpd = m*q[i];
        p[i] += tmpd; 			// this is the point where the parameter vector steps 
        step += fabs(tmpd);
        q[i] = s*gy[i] + t*gx[i];	// send back the estimated gradient in q (NB not like linmin) 
    }
    fix(p);
    lastx = m*LINMIN_G2*gtyp;
    return ( step / ((double)dimN) );
}

static int min(){
    int it,i;
    double gg, gam, dgg ,tmpd,maxGradSqr;

    lastx = LAST_X_DEFAULT;
    restart = 0;
    updateNebzLists();
    for (i=0;i<N;i++){
        p[DIM*i+X_COMP] = rx[i];
        p[DIM*i+Y_COMP] = ry[i];
    }

    grad(p, q);
    initializeMinimizer(1);
    for (it = 0; it < MAX_ITERATIONS; it++){
        maxGradSqr = 0.0; gg=0.0;
        for (i=0;i<dimN;i++){
            tmpd = g[i]*g[i];  
            gg += tmpd;
            if (tmpd > maxGradSqr)    //To check when to stop
                maxGradSqr = tmpd;  
        }

        gtyp = sqrt( gg / ((double)dimN) ); 

        if ( maxGradSqr < TOL ){   //To stop
            for (i=0;i<N;i++){
                rx[i] = p[DIM*i+X_COMP];
                ry[i] = p[DIM*i+Y_COMP];
            }
            updateNebzLists();
            calculateForces();
            return 1;
        }

        //step = linmin();
        linmin();
        if (restart)
            grad(p, q);
        if (restart){
            printf("\nRestarting optimizer1\n");
            initializeMinimizer(0);
        }
        else{
            for ( dgg=0.0, i=0; i<dimN; i++)
                dgg += ( q[i]+g[i] )*q[i];
            gam = dgg/gg;
            if (gam > 3.0){
                //printf("gamma was too big, restarting optimizer\n");
                restart = 1;
                initializeMinimizer(0);
            }
            else{//this is my addition, since <gam> can't be too large, it fucks everything up.
                for ( tmpd=0.0, i=0; i<dimN; i++){
                    g[i] = -q[i];                // g stores (-) the most recent gradient 
                    q[i] = h[i] = g[i] + gam*h[i];
                    // h stores q, the current line direction 
                    // check that the inner product of gradient and line search is < 0 
                    tmpd -= q[i]*g[i];
                }
                if ( tmpd > 0.0 ){
                    restart = 2; 	// signifies that g[j] = -q[j] is already done 
                    printf("Restarting optimizer (2)\n");
                    initializeMinimizer(0);
                }
            }
        }
    }
    return 0;
}

static void calculateHessian(){
    int indx, indy;
    int i, j, k, n, tag_i, tag_j;
    double xi, yi, dx, dy, r2, s2, s4, s6, s12, ri2, ri8, ri14, ri10, ri16, phi, si;
    double hesXY;

    for(i=0;i<dimNsq;i++){
        upTriH[i]=0.e0;
    }
    //................................
    for(i=0;i<N;i++){
        indx = dimN*((2*i)+X_COMP);
        indy = dimN*((2*i)+Y_COMP);
        xi = rx[i];
        yi = ry[i];
        tag_i=type[i];
        
        for(j=i+1;j<N;j++){
            tag_j = type[j];
            dx = rx[j]-xi;
            dy = ry[j]-yi;

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
    for (i=0; i<N; i++){
        indi = (dimN+1)*(2*i);
        indi1 = indi + 1;
        indi2N1 = indi+dimN+1;
        for (j=1; j<(N-i); j++){
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
}



static void calculateEigenValues(){
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

        for(j=0;j<N;j++){
            e2 = (eigenVectors[(i*dimN)+DIM*j]*eigenVectors[(i*dimN)+DIM*j])+(eigenVectors[(i*dimN)+(DIM*j)+1]*eigenVectors[(i*dimN)+(DIM*j)+1]);
            e4 = e2*e2;
            nume+=e2;
            deno+=e4;
        }
        pr[i] = (1.e0)/(N*deno);
    }
}


int main(){

    int i,j,k, tag;
    double x, y, dist_x, dist_y, dist_r, dist_x_i, dist_y_i, dist_r_i, dist_x_j, dist_y_j, dist_r_j;
    int dummy1;
    double dummy2, dummy3;

    start = clock();  // Initializing clock
    //If an command line arg is not provided - error & exit
    
	unsigned int starting, max, ens_id;
    starting= 100000;
    max = 200000;
    for(ens_id=starting;ens_id<max;ens_id++){
	printf("ens_id->%06u\n",ens_id);
        char dataFileName[128];
        FILE *dataFile;
        sprintf(dataFileName,"%06u_configuration.txt",ens_id);
        dataFile=fopen(dataFileName,"rb");
        if(dataFile==NULL) {
                printf("ERROR No File at :'%s'\n",dataFileName);
        }
        k = 0;
        while ( fscanf(dataFile,"%d %lf %lf",&dummy1,&dummy2,&dummy3) != EOF ){
            k++;
        }
        rewind(dataFile);

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
        double rad = 30.e0;
        double cmx = 0.e0;
        double cmy = 0.e0;
        for(i=0;i<2500;i++){
            cmx+=Ax[i];
            cmy+=Ay[i];
        }
        
        cmx/=(2500*1.e0);
        cmy/=(2500*1.e0);
        printf("%E\t%E\n",cmx,cmy);
        k=0;
        for(i=0;i<2500;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y);
            if(dist_r<=(rad*rad)) k++;        
        }
	printf("%d\n",k);
        if(k>=400){
        rx= malloc(k * sizeof(*rx));
        ry= malloc(k * sizeof(*ry));
        type= malloc(k * sizeof(*type));
        k=0;
        for(i=0;i<2500;i++){
            dist_x = Ax[i]-cmx;
            dist_y = Ay[i]-cmy;
            dist_r = (dist_x*dist_x)+(dist_y*dist_y);
            if(dist_r<=(rad*rad)){
                rx[k]=Ax[i];
                ry[k]=Ay[i];
                type[k]=Atag[i];
                k++;
            }        
        }
        //------------------------ Sorting ------------------------------
        for(j=1;j<k;j++){

            x=rx[j];
            y=ry[j];
            tag=type[j];

            i=j-1;

            dist_x_j = rx[j]-cmx;
            dist_y_j = ry[j]-cmy;
            dist_r_j = (dist_x_j*dist_x_j)+(dist_y_j*dist_y_j);
            
            dist_x_i = rx[i]-cmx;
            dist_y_i = ry[i]-cmy;
            dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i);

            while(i>=0 && dist_r_i>dist_r_j){
                rx[i+1]=rx[i];
                ry[i+1]=ry[i];
                type[i+1]=type[i];
                i--;
                dist_x_i = rx[i]-cmx;
                dist_y_i = ry[i]-cmy;
                dist_r_i = (dist_x_i*dist_x_i)+(dist_y_i*dist_y_i);
            }
            rx[i+1]=x;
            ry[i+1]=y;
            type[i+1]=tag;
        }
        //--------------------------------------------------------------
        
        printf("%d\n",k);
        N=400;
        DOUBLE_N = N*1.e0;
        dimN = DIM*N;
        dimNsq = dimN*dimN;
        
        numOfNebz= malloc(N* sizeof(*numOfNebz));
        for(k=0;k<2500;k++){
            nebz[k]= malloc(N* sizeof(*nebz));
        }    

        fx= malloc(N * sizeof(*fx));
        fy= malloc(N * sizeof(*fy));

        p = malloc(dimN * sizeof(*p));
        q = malloc(dimN * sizeof(*q));
        g = malloc(dimN * sizeof(*g));
        h = malloc(dimN * sizeof(*h));
        mosesX = malloc(dimN * sizeof(*mosesX));
        mosesY = malloc(dimN * sizeof(*mosesY));
        scratch = malloc(dimN * sizeof(*scratch));

        upTriH = malloc(dimNsq * sizeof(*upTriH));
        eigenValues = malloc(dimN*sizeof(*eigenValues));
        eigenVectors = malloc(dimNsq*sizeof(*eigenVectors));
        pr = malloc(dimN*sizeof(*pr));
        

        
        updateNebzLists();
        calculateForces();
        printf("%E\t%E\n",pe,typicalGrad);
        
        min();
        
        calculateForces();
        printf("%E\t%E\n",pe,typicalGrad);

        char cg_config_FileName[128];
        FILE  *cg_config; 
        sprintf(cg_config_FileName,"cg_config_of_%d_particles_%06u_ens_id.txt",ens_id);
        cg_config=fopen(cg_config_FileName,"w");
        for(k=0;k<N;k++){
            fprintf(cg_config,"%d\t%.16lf\t%.16lf\n",type[k],rx[k],ry[k]);
        }
        fclose(cg_config);

        calculateHessian();
        calculateEigenValues();

        char eigenvalues_FileName[128];
        FILE  *eigenvalues; 
        sprintf(eigenvalues_FileName,"eigenvalues_of_%d_particles_%06u_ens_id_CG.txt",TP,ens_id);
        eigenvalues=fopen(eigenvalues_FileName,"w");
        for(k=0;k<dimN;k++){
            fprintf(eigenvalues,"%E\t%E\t%E\n",eigenValues[k],sqrt(eigenValues[k]),pr[k]);
        }
        fclose(eigenvalues);
        
	char eigenvectors_FileName[128];
        FILE  *eigenvectors; 
        sprintf(eigenvectors_FileName,"first_non_zero_mode_of_%d_particles_%06u_ens_id_CG.txt",TP,ens_id);
        eigenvectors=fopen(eigenvectors_FileName,"w");
        
	i=3;

        for(j=0;j<N;j++){
            fprintf(eigenvectors,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",type[j],rx[j],ry[j],eigenVectors[(i*dimN)+DIM*j],eigenVectors[(i*dimN)+(DIM*j)+1]);
	 } 
        fclose(eigenvectors);

        
        free(rx);
        free(ry);
        free(type);
        free(numOfNebz);
        for(k=0;k<2500;k++){
            free(nebz[k]);
        }
        free(fx);
        free(fy);
        free(p);
        free(q);
        free(g);
        free(h);
        free(mosesX);
        free(mosesY);
        free(scratch);
        free(eigenValues);
        free(eigenVectors);
        free(pr);
	free(upTriH);
        }
        free(Ax);
        free(Ay);
        free(Atag);
        

    }
    end = clock();   //  Ending clock
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    printf("Time taken: %g sec\n",cpu_time_used);
}
