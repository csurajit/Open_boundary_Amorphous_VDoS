void cell_pbc_2d(){
    
    int i,j,k,n,m,ix,iy;
    double dx,dy,dz,dr;
    
    int currentCell,targetCell, numberInCurrent,numberInTarget;
    int c,x1, x2, y1, y2;
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
        m=(int)((x[n]+half_length)/r_c)+(n_c*(int)((y[n]+half_length)/r_c));
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
        for(iy=0;iy<n_c;iy++){
            for(ix=0;ix<n_c;ix++){
                currentCell=ix+(iy*n_c);
                numberInCurrent=nc[currentCell];

                //................In Current cell....................
                for(k=0;k<numberInCurrent-1;k++){
                    for(n=k+1;n<numberInCurrent;n++){
                        i=np_c[currentCell][k];
                        j=np_c[currentCell][n];
                        dx=x[i]-x[j];
                        dy=y[i]-y[j];
                        
                        dr=((dx*dx)+(dy*dy));
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
                cells[0]=x1+(iy*n_c);
                cells[1]=x1+(y1*n_c);
                cells[2]=ix+(y1*n_c);
                cells[3]=x2+(y1*n_c);
                
                //###################################################
                for(c=0;c<4;c++){
                    targetCell=cells[c];

                    numberInTarget=nc[targetCell];

                    for(k=0;k<numberInCurrent;k++){
                        for(n=0;n<numberInTarget;n++){
                            i=np_c[currentCell][k];
                            j=np_c[targetCell][n];

                            dx=x[i]-x[j];
                            dy=y[i]-y[j];
                            if(dx>=half_length)
                                dx-=length;
                            if(dx<-half_length)
                                dx+=length;
                            
                            if(dy>=half_length)
                                dy-=length;
                            if(dy<-half_length)
                                dy+=length;
                            
                            dr=((dx*dx)+(dy*dy));
                            
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
