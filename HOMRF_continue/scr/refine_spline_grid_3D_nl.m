function [Knots_n] = refine_spline_grid_3D_nl(Knots, O_trans_X_or,O_trans_Y_or,O_trans_Z_or, O_trans_X,O_trans_Y,O_trans_Z)
    %[Knots_n] = refine_linear_grid_3d(Knots, old_spacing, ds, volsz,
    %  {upsampling_type}
    ksz_new=size(O_trans_X);
    Kn=zeros(ksz_new);
    x_or=O_trans_X_or(:,1,1,2)-O_trans_X_or(:,1,1,1);
    x_or=O_trans_X_or(:,1,1,1)+x_or./2;
    y_or=O_trans_Y_or(1,:,1,2)-O_trans_Y_or(1,:,1,1);
    y_or=O_trans_Y_or(1,:,1,1)+y_or./2;
    z_or=O_trans_Z_or(1,1,:,2)-O_trans_Z_or(1,1,:,1);
    z_or=O_trans_Z_or(1,1,:,1)+z_or./2;
    [x2,y2,z2]=meshgrid(y_or,x_or,z_or);


    x_next=O_trans_X(:,1,1,2)-O_trans_X(:,1,1,1);
    x_next=O_trans_X(:,1,1,1)+x_next./2;
    y_next=O_trans_Y(1,:,1,2)-O_trans_Y(1,:,1,1);
    y_next=O_trans_Y(1,:,1,1)+y_next./2;
    z_next=O_trans_Z(1,1,:,2)-O_trans_Z(1,1,:,1);
    z_next=O_trans_Z(1,1,:,1)+z_next./2;
    [x1,y1,z1]=meshgrid(y_next,x_next,z_next);


    Kn(:,:,:,1)=interp3(x2,y2,z2,Knots(:,:,:,1),x1,y1,z1,'linear',0);
    Kn(:,:,:,2)=interp3(x2,y2,z2,Knots(:,:,:,2),x1,y1,z1,'linear',0);
    Kn(:,:,:,3)=interp3(x2,y2,z2,Knots(:,:,:,3),x1,y1,z1,'linear',0);
    Knots_n=Kn;





end
