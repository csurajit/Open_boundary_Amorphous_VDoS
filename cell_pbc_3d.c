void cell_pbc_3d(){
    
    int i,j,k,n,m,ix,iy,iz;
    double dx,dy,dz,dr;
    
    int currentCell,targetCell, numberInCurrent,numberInTarget;
    int c,x1, x2, y1, y2,z1,z2;
    //.................. Neighbours outputs are Initialising to 0 ...................
    for(i=0;i<TP;i++){
        nn[i]=0; 
    } 
    //...............................................................................
    //........ initialising number of particle in each cells is 0..............
    for(m=0;m<tn_c;m++){
        nc[m]=0;           
    }
    //........................................................................
    //################################ CELL LIST + VERLET LIST ######################################################

    //....................... Identifying number of particles for each cell....................
    for(n=0;n<TP;n++){        
        m=(int)((x[n]+half_length)/r_c)+(n_c*(int)((y[n]+half_length)/r_c))+(n_c*n_c*(int)((z[n]+half_length)/r_c));
        if(m<0 || m>tn_c){
            printf("wrong\n%d\n",m);
            printf("%lf\t%lf\n",x[n],y[n]);
            exit(EXIT_FAILURE);
        }
        np_c[m][nc[m]]=n;    //np_c[8][1]=2,, 1st particle in the 8th cell is 2
        nc[m]=nc[m]+1;         //no of particle in mth cell
        
    }
    //printf("%lf\t%lf\t%d\t%d\n",L,r_c,n_c,m);
    //...................................................................................................
    for(iz=0;iz<n_c;iz++){
        for(iy=0;iy<n_c;iy++){
            for(ix=0;ix<n_c;ix++){
                currentCell=ix+(iy*n_c)+(iz*n_c*n_c);
                numberInCurrent=nc[currentCell];

                //................In Current cell....................
                for(k=0;k<numberInCurrent-1;k++){
                    for(n=k+1;n<numberInCurrent;n++){
                        i=np_c[currentCell][k];
                        j=np_c[currentCell][n];
                        dx=x[i]-x[j];
                        dy=y[i]-y[j];
                        dz=z[i]-z[j];
                        dr=((dx*dx)+(dy*dy)+(dz*dz));
                        if(dr<(cfn*cfn)){      
                            n_list[i][nn[i]]=j;            
                            n_list[j][nn[j]]=i;            
                            nn[i]+=1; 
                            nn[j]+=1;
                        }
                    }
                }
                //...............................................................
                //###### To find the Neighbouring 13 cells For PBC ###############
                x1=(ix+1);
                x2=(ix-1);
                y1=(iy+1);
                y2=(iy-1);
                z1=(iz+1);
                z2=(iz-1);
                
                //... For PBC.....
                if(ix==(n_c-1)){
                    x1=0;
                }
                if(ix==0){
                    x2=n_c-1;
                }
                if(iy==(n_c-1)){
                    y1=0;
                }
                if(iy==0){
                    y2=n_c-1;
                }
                if(iz==(n_c-1)){
                    z1=0;
                }
                if(iz==0){
                    z2=n_c-1;
                }
                
                cells[0]=ix+(iy*n_c)+(z1*n_c*n_c);
                cells[1]=x1+(iy*n_c)+(z1*n_c*n_c);
                cells[2]=x2+(iy*n_c)+(z1*n_c*n_c);
                cells[3]=ix+(y1*n_c)+(z1*n_c*n_c);
                cells[4]=ix+(y2*n_c)+(z1*n_c*n_c);
                cells[5]=x1+(y1*n_c)+(z1*n_c*n_c);
                cells[6]=x1+(y2*n_c)+(z1*n_c*n_c);
                cells[7]=x2+(y1*n_c)+(z1*n_c*n_c);
                cells[8]=x2+(y2*n_c)+(z1*n_c*n_c);
                
                cells[9]=ix+(y1*n_c)+(iz*n_c*n_c);
                
                cells[10]=x2+(iy*n_c)+(iz*n_c*n_c);
                cells[11]=x2+(y1*n_c)+(iz*n_c*n_c);
                cells[12]=x2+(y2*n_c)+(iz*n_c*n_c);
                //###################################################
                for(c=0;c<13;c++){
                    targetCell=cells[c];

                    numberInTarget=nc[targetCell];

                    for(k=0;k<numberInCurrent;k++){
                        for(n=0;n<numberInTarget;n++){
                            i=np_c[currentCell][k];
                            j=np_c[targetCell][n];

                            dx=x[i]-x[j];
                            dy=y[i]-y[j];
                            dz=z[i]-z[j];
                            if(dx>=half_length)
                                dx=dx-length;
                            if(dx<-half_length)
                                dx=dx+length;
                            
                            if(dy>=half_length)
                                dy=dy-length;
                            if(dy<-half_length)
                                dy=dy+length;
                            
                            if(dz>=half_length)
                                dz=dz-length;
                            if(dy<-half_length)
                                dy=dy+length;
                            
                            dr=((dx*dx)+(dy*dy)+(dz*dz));
                            
                            if(dr<(cfn*cfn)){      
                                n_list[i][nn[i]]=j;            
                                n_list[j][nn[j]]=i;           
                                nn[i]+=1; 
                                nn[j]+=1;
                            }
                        }
                    }
                }
                //................................................................

            }
        }
    }
    
    //###################################################################################################
    
    for(i=0;i<TP;i++)
    {
        for(j=0;j<nn[i]-1;j++)
        {
            for(k=j+1;k<nn[i];k++)
            {
                if(n_list[i][j]>n_list[i][k])
                {
                    n=n_list[i][j];
                    n_list[i][j]=n_list[i][k];
                    n_list[i][k]=n;
                }
            }
        }
    }
    
}
