
function [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(sizeG,sizeI)

  
if(length(sizeG)==2)
    % Determine grid spacing
    dx=sizeI(1)/sizeG(1);
    dy=sizeI(2)/sizeG(2);

    % Calculate te grid coordinates (make the grid)
    [X_1,Y_1]=ndgrid(1:dx:sizeI(1),1:dy:sizeI(2));
    X_1=fix(X_1);
    Y_1=fix(Y_1);
    [X_2,Y_2]=ndgrid(dx:dx:sizeI(1),dy:dy:sizeI(2));
    X_2=ceil(X_2);
    Y_2=ceil(Y_2);

     O_trans_X(:,:,1)=X_1;
     O_trans_X(:,:,2)=X_2;
     O_trans_Y(:,:,1)=Y_1;
     O_trans_Y(:,:,2)=Y_2;
else
    % Determine grid spacing
%     dx=ceil(sizeI(1)/sizeG(1));
%     dy=ceil(sizeI(2)/sizeG(2));
%     dz=ceil(sizeI(3)/sizeG(3));
    dx=sizeG(1);
    dy=sizeG(2);
    dz=sizeG(3);
    nx=fix(sizeI(1)/dx);
    ny=fix(sizeI(2)/dy);
    nz=fix(sizeI(3)/dz);
    resid_x=mod(sizeI(1),dx);
    resid_y=mod(sizeI(2),dy);
    resid_z=mod(sizeI(3),dz);
    
    if resid_x==0
        x_1=(1:dx:(sizeI(1)-dx+1));
        x_2=(dx:dx:sizeI(1));
    else
        star=fix(resid_x/2)+1;
        x_1=(star:dx:(sizeI(1)-dx+1));
        x_2=((star+dx-1):dx:sizeI(1));
    end
    
    
    if resid_y==0
        y_1=(1:dy:(sizeI(2)-dy+1));
        y_2=(dy:dy:sizeI(2));
    else
        star=fix(resid_y/2)+1;
        y_1=(star:dy:(sizeI(2)-dy+1));
        y_2=((star+dy-1):dy:sizeI(2));
    end
    
        
    if resid_z==0
        z_1=(1:dz:(sizeI(3)-dz+1));
        z_2=(dz:dz:sizeI(3));
    else
        star=fix(resid_z/2)+1;
        z_1=(star:dz:(sizeI(3)-dz+1));
        z_2=((star+dz-1):dz:sizeI(3));
    end
             
    % Calculate te grid coordinates (make the grid)
    [X_1,Y_1,Z_1]=ndgrid(x_1,y_1,z_1);
    X_1=fix(X_1);
    Y_1=fix(Y_1);
    Z_1=fix(Z_1);

    [X_2,Y_2,Z_2]=ndgrid(x_2,y_2,z_2);
    X_2=ceil(X_2);
    Y_2=ceil(Y_2);
    Z_2=fix(Z_2);

    O_trans_X(:,:,:,1)=X_1;
    O_trans_X(:,:,:,2)=X_2;
    O_trans_Y(:,:,:,1)=Y_1;
    O_trans_Y(:,:,:,2)=Y_2;
    O_trans_Z(:,:,:,1)=Z_1;
    O_trans_Z(:,:,:,2)=Z_2;
end
end