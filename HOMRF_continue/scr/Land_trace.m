function landregis=Land_trace(Tmcmc_rsz,Movxyz,init_size)
Tmcmc_rsz=round(Tmcmc_rsz);
DIR_coordinate_x=Tmcmc_rsz(:,:,:,1);
DIR_coordinate_y=Tmcmc_rsz(:,:,:,2);
DIR_coordinate_z=Tmcmc_rsz(:,:,:,3);

landregis=zeros(300,3);
numlandmark=zeros(300,1);

for i=1:300
    numlandmark(i)=(Movxyz(i,3)-1)*init_size(1)*init_size(2)+(Movxyz(i,2)-1)*init_size(1)+Movxyz(i,1);
end

Movindex=numlandmark; 

for j=1:300   
    landregis(j,1)=Movxyz(j,1)+DIR_coordinate_x(Movindex(j));
    landregis(j,2)=Movxyz(j,2)+DIR_coordinate_y(Movindex(j));
    landregis(j,3)=Movxyz(j,3)+DIR_coordinate_z(Movindex(j));
end

end