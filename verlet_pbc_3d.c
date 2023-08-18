void verlet_pbc_3d(){
    int i, j;
    double xi, yi, zi, dx, dy, dz, dr;
    for(i=0;i<TP;i++){
        nn[i]=0;
    }

    for(i=0;i<TP-1;i++){
        xi = x[i];
        yi = y[i];
        zi = z[i];
        for(j=i+1;j<TP;j++){
            dx = x[j]-xi;
            dy = y[j]-yi; 
            dz = z[j]-zi;
            if(dx>=half_length)  dx=dx-length;
            else if(dx< -half_length)  dx=dx+length;
            if(dy>=half_length)  dy=dy-length;
            else if(dy< -half_length)  dy=dy+length;
            if(dz>=half_length)  dz=dz-length;
            else if(dz< -half_length)  dz=dz+length;
            dr = (dx*dx)+(dy*dy)+(dz*dz);
            if(dr<cfn*cfn){
                n_list[i][nn[i]] = j;
                n_list[j][nn[j]] = i;
                nn[i]++;
                nn[j]++;
            }
        }
    }

}
