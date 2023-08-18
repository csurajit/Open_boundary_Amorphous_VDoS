void verlet_pbc_2d(){
    int i, j;
    double xi, yi, dx, dy, dz, dr;
    for(i=0;i<TP;i++){
        nn[i]=0;
    }

    for(i=0;i<TP-1;i++){
        xi = x[i];
        yi = y[i];
        for(j=i+1;j<TP;j++){
            dx = x[j]-xi;
            dy = y[j]-yi;
            if(dx>=half_length)  dx=dx-length;
            else if(dx< -half_length)  dx=dx+length;
            if(dy>=half_length)  dy=dy-length;
            else if(dy< -half_length)  dy=dy+length;
            dr = (dx*dx)+(dy*dy);
            if(dr<cfn*cfn){
                n_list[i][nn[i]] = j;
                n_list[j][nn[j]] = i;
                nn[i]++;
                nn[j]++;
            }
        }
    }

}
