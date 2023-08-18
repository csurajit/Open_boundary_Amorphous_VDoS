void force_pbc_nn_3d(){
    int tag_i, tag_j;
    double xi,yi,zi,r2i,r6i,r;
    double dx,dy,dz,r2,s2,fr,store_x,store_y,store_z;
    int i,k,j,n;

    for(i=0;i<TP;i++){
        fx[i]=0.e0;
        fy[i]=0.e0;
        fz[i]=0.e0;
    }
    pe = 0.e0;
    vir = 0.e0;

    double Tsig_xx=0.e0;
    double Tsig_xy=0.e0;
    double Tsig_yy=0.e0;
    double Tsig_zz=0.e0;
    double Tsig_xz=0.e0;
    double Tsig_yz=0.e0;
    for(i=0;i<TP-1;i++){
        xi=x[i];
        yi=y[i];
        zi=z[i];
        tag_i= tag[i];
        n=nn[i];
        for(k=0;k<n;k++){
            j=n_list[i][k];
            tag_j = tag[j];
            if(j>i){
                dx=x[j]-xi;
                dy=y[j]-yi;
                dz=z[j]-zi;
                if(dx>=half_length)  dx=dx-length;
                else if(dx< -half_length)  dx=dx+length;
                if(dy>=half_length)  dy=dy-length;
                else if(dy< -half_length)  dy=dy+length;
                if(dz>=half_length)  dz=dz-length;
                else if(dz< -half_length)  dz=dz+length;

                r2=(dx*dx)+(dy*dy)+(dz*dz);
                if(tag_i==1 && tag_j==1 && r2<(r_cut_AA*r_cut_AA)){
                    s2= sig_AA*sig_AA;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AA*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    fr = 4.e0*eps_AA*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    
                    store_x=fr*dx;
                    fx[j] += store_x;  //
                    fx[i] -= store_x;  //
                    store_y=fr*dy;
                    fy[j] += store_y;  //Final force calculation componantwise
                    fy[i] -= store_y;  //
                    store_z=fr*dz;
                    fz[j] += store_z;  //
                    fz[i] -= store_z;  //

                    Tsig_xx+=store_x*dx;
                    Tsig_yy+=store_y*dy;
                    Tsig_xy+=store_x*dy;
                    Tsig_zz+=store_z*dz;
                    Tsig_yz+=store_y*dz;
                    Tsig_xz+=store_x*dz;

                    
                }
                if (tag_i==0 && tag_j==0 && r2<(r_cut_BB*r_cut_BB)){
                    s2= sig_BB*sig_BB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_BB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    fr = 4.e0*eps_BB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    
                    store_x=fr*dx;
                    fx[j] += store_x;  //
                    fx[i] -= store_x;  //
                    store_y=fr*dy;
                    fy[j] += store_y;  //Final force calculation componantwise
                    fy[i] -= store_y;  //
                    store_z=fr*dz;
                    fz[j] += store_z;  //
                    fz[i] -= store_z;  //

                    Tsig_xx+=store_x*dx;
                    Tsig_yy+=store_y*dy;
                    Tsig_xy+=store_x*dy;
                    Tsig_zz+=store_z*dz;
                    Tsig_yz+=store_y*dz;
                    Tsig_xz+=store_x*dz;

                    
                }
                if (tag_i!=tag_j && r2<(r_cut_AB*r_cut_AB)){
                    s2= sig_AB*sig_AB;
                    r2i = 1.e0/r2;
                    r6i = r2i*r2i*r2i*s2*s2*s2;
                    pe+=4.e0*eps_AB*((r6i*(r6i-1.e0))+c_0+(c_2*r2/s2)+(c_4*r2*r2/(s2*s2)));
                    fr = 4.e0*eps_AB*(12.e0*r2i*r6i*(r6i-0.5e0)-(2.e0*c_2/s2)-4.e0*c_4*r2/(s2*s2));
                    
                    store_x=fr*dx;
                    fx[j] += store_x;  //
                    fx[i] -= store_x;  //
                    store_y=fr*dy;
                    fy[j] += store_y;  //Final force calculation componantwise
                    fy[i] -= store_y;  //
                    store_z=fr*dz;
                    fz[j] += store_z;  //
                    fz[i] -= store_z;  //

                    Tsig_xx+=store_x*dx;
                    Tsig_yy+=store_y*dy;
                    Tsig_xy+=store_x*dy;
                    Tsig_zz+=store_z*dz;
                    Tsig_yz+=store_y*dz;
                    Tsig_xz+=store_x*dz;

                    
                }
            }    
        }
    }
    vir = Tsig_xx+Tsig_yy+Tsig_zz;
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
}
