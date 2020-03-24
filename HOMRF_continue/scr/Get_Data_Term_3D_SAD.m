%use the sad metric
function unary=Get_Data_Term_3D_SAD(Mov,Fix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant)

[r,c,s]=size(Mov);
sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);
min_max=[min(previous,[],1)',max(previous,[],1)'];
min_max=abs(round(min_max));
shx_or=max(min_max(1,:));
shy_or=max(min_max(2,:));
shz_or=max(min_max(3,:));
sx_or=max(min_max(1,:))*2+1;
sy_or=max(min_max(2,:))*2+1;
sz_or=max(min_max(3,:))*2+1;
add_hx=shx_or+labels.hx*sample_space_x;
add_hy=shy_or+labels.hy*sample_space_y;
add_hz=shz_or+labels.hz*sample_space_z;
size_orig=[r+(add_hx)*2+1,c+(add_hy)*2+1,s+(add_hz)*2+1];
crop_v=[1+add_hx,r+add_hx;
        1+add_hy,c+add_hy;
        1+add_hz,s+add_hz;];
    weight=100;
    volmovl= uncrop_data_xp_3d(Mov, crop_v, size_orig,weight);
    volfixl= uncrop_data_xp_3d(Fix, crop_v, size_orig,weight);
    
labels_or.sx=sx_or;
labels_or.sy=sy_or;
labels_or.sz=sz_or;
labels_or.n_each_layer=sx_or*sy_or;
labels_or.nlabels=sx_or*sy_or*sz_or;
labels_or.hx=fix(sx_or/2);
labels_or.hy=fix(sy_or/2);
labels_or.hz=fix(sz_or/2);
labels_or_index=trans_labels(previous,sx_or,sy_or,sz_or);

dim_r=griddim(1);
dim_c=griddim(2);
dim_s=griddim(3);

[x_index,y_index,z_index]=Label_Coordinate(labels);
[x_index_or,y_index_or,z_index_or]=Label_Coordinate(labels_or);

unary=cell(dim_r*dim_c*dim_s,1);
parfor num=1:dim_r*dim_c*dim_s
    unary_cur=zeros(1,labels.nlabels);
    previous_label=labels_or_index(num);
    [i,j,k]=ind2sub([dim_r,dim_c,dim_s],num);
    i_1=O_trans_X(i,j,k,1)+add_hx;
    i_2=O_trans_X(i,j,k,2)+add_hx;
    j_1=O_trans_Y(i,j,k,1)+add_hy;
    j_2=O_trans_Y(i,j,k,2)+add_hy;
    k_1=O_trans_Z(i,j,k,1)+add_hz;
    k_2=O_trans_Z(i,j,k,2)+add_hz;
    x_width=i_2-i_1+1;
    y_length=j_2-j_1+1;
    z_hight=k_2-k_1+1;
    v_zone=x_width*y_length*z_hight;
    for nlb=1:labels.nlabels
        x_mov=x_index(nlb)*sample_space_x+x_index_or(previous_label);
        y_mov=y_index(nlb)*sample_space_y+y_index_or(previous_label);
        z_mov=z_index(nlb)*sample_space_z+z_index_or(previous_label);
        M=volmovl(i_1+x_mov:i_2+x_mov,j_1+y_mov:j_2+y_mov,k_1+z_mov:k_2+z_mov);
        F=volfixl(i_1:i_2,j_1:j_2,k_1:k_2);
        sad=abs(M-F);
        lcc=sum(sad(:))/v_zone;
        unary_cur(nlb)=lcc;
            
    end
    unary{num}=unary_cur;
end
   unary=cell2mat(unary);
end


