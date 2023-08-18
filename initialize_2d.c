void initialize_2d(){
    int np, i, j, n;
    double lc;
    rho = 1.2e0;
    np = 1;
    while(np*np<TP) np++;
    printf("np is %d\n",np);
    length = pow((double_TP/rho),0.5e0);
    half_length = length/2.e0;
    lc = length/(np*1.e0);
    printf("Initial length is %lf\n",length);
    printf("Initial half_length is %lf\n",half_length);
    printf("lc is %lf\n",lc);

    for(j=0;j<np;j++){
        for(i=0;i<np;i++){
            n=i+(j*np);
            x[n]=i*lc-(half_length);
            y[n]=j*lc-(half_length);
        }
    }
    for(n=0;n<TP;n++){
        if(n<n_A){
            tag[n]=1;
        }
        else{
            tag[n]=0;
        }
    }
    for(i=0;i<100000;i++){
        int p=(int)(dblRanq1()*TP);
        int s=(int)(dblRanq1()*TP);
        n=tag[p];
        tag[p]=tag[s];
        tag[s]=n;
    }
    int A=0;
    int B=0;
    for(i=0;i<TP;i++){
        if(tag[i]==1){
            A+=1;
        }
        if(tag[i]==0){
            B+=1;
        }
    }
    
    printf("No of A particle %d\tNo of B particle is%d\t total is%d\n",A,B,A+B);
    //-------------------------------------------------------

    //------------------------------------------------------
    for(n=0;n<TP;n++){
        vx[n]=dblRanq1()-0.50;
        vy[n]=dblRanq1()-0.50;
    }
    cmvx=0.e0;
    cmvy=0.e0;
    for(n=0;n<TP;n++){
        cmvx+=vx[n];
        cmvy+=vy[n];
    }
    //... making cm velocities zero  ....
    for (n=0; n<TP; n++) {
        vx[n]-= (cmvx/double_TP);
        vy[n]-= (cmvy/double_TP);
    }
    //...... check ......
    cmvx=0.e0;
    cmvy=0.e0;
    ssum=0.e0;
    for(n=0;n<TP;n++){
        cmvx+=vx[n];
        cmvy+=vy[n];
        ssum +=((vx[n]*vx[n])+(vy[n]*vy[n]));
    }
    T = (0.5e0*ssum)/(TP-1);
    printf("Begin with %E\t%E\t%E\n",cmvx,cmvy,T);
}
