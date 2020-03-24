function Knots_n_rsz_smooth=filter_3D(Knots_n_rsz,grid_space)

X_d=Knots_n_rsz(:,:,:,1);
Y_d=Knots_n_rsz(:,:,:,2);
Z_d=Knots_n_rsz(:,:,:,3);
if nargin==1
    Hsize=[3,3,3];
    Hsmooth=fspecial3_v1('gaussian',Hsize);  
else
    Hsize=fix(grid_space/2)*2+1;
    Hsmooth=fspecial3_v1('gaussian',Hsize);
end
X_d_smooth=imfilter(X_d,Hsmooth);
Y_d_smooth=imfilter(Y_d,Hsmooth);
Z_d_smooth=imfilter(Z_d,Hsmooth);
Knots_n_rsz_smooth(:,:,:,1)=X_d_smooth;
Knots_n_rsz_smooth(:,:,:,2)=Y_d_smooth;
Knots_n_rsz_smooth(:,:,:,3)=Z_d_smooth;