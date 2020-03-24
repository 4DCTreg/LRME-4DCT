
function MIND_cache=Get_Data_Term_3D_MIND_cache(Fix_mind,O_trans_all,nlevel,grid_space)
MIND_cache=cell(nlevel);

for idx=nlevel:-1:2
    cur_grid_space=grid_space(idx,:);
    O_trans=O_trans_all{idx};
    O_trans_X=O_trans(:,:,:,:,1);
    O_trans_Y=O_trans(:,:,:,:,2);
    O_trans_Z=O_trans(:,:,:,:,3);
    griddim=[size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
    dim_r=griddim(1);dim_c=griddim(2);dim_s=griddim(3);
    F=zeros(dim_r*dim_c*dim_s,prod(cur_grid_space),R);
    for R_index=1:R     
        volfixl= Fix_mind(:,:,:,R_index);
        for num=1:dim_r*dim_c*dim_s  
            [i,j,k]=ind2sub([dim_r,dim_c,dim_s],num);
            i_1=O_trans_X(i,j,k,1);
            i_2=O_trans_X(i,j,k,2);
            j_1=O_trans_Y(i,j,k,1);
            j_2=O_trans_Y(i,j,k,2);
            k_1=O_trans_Z(i,j,k,1);
            k_2=O_trans_Z(i,j,k,2);
            F_cur=volfixl(i_1:i_2,j_1:j_2,k_1:k_2);
            F(num,:,R_index)=F_cur(:)';
        end
    end
    MIND_cache{idx}=F;
end


