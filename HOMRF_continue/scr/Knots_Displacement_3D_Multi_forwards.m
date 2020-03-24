

function Knots=Knots_Displacement_3D_Multi_forwards(Ln,sx,sy,sz,griddim,level,quant)
sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);

n_each_layer=sx*sy;
hz=floor(sz/2);
hx=floor(sx/2);
hy=floor(sy/2);



r=griddim(1);
c=griddim(2);
l=griddim(3);


Knots=zeros(r,c,l,3);


regisMovnum=reshape(Ln,[r,c,l]);


Tx=zeros(r,c,l);Ty=zeros(r,c,l);Tz=zeros(r,c,l);
for j=1:c
    for i=1:r
        for k=1:l
            orpz=fix((regisMovnum(i,j,k)-1)/n_each_layer)-hz;
            intermediate=mod(regisMovnum(i,j,k),n_each_layer);
            if intermediate==0
              intermediate=n_each_layer;
            end
            orpy=floor((intermediate-0.5)/sx)-hy;
            orpx=mod((intermediate-1),sx)-hx;
            lx=orpx*sample_space_x;
            ly=orpy*sample_space_y;
            lz=orpz*sample_space_z;
            
            Tz(i,j,k)=lz;
            Ty(i,j,k)=ly;
            Tx(i,j,k)=lx;
        end
    end
end
Knots(:,:,:,1)=Tx;
Knots(:,:,:,2)=Ty;
Knots(:,:,:,3)=Tz;

end