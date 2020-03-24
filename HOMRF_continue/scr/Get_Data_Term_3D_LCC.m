% use the lcc metric
function unary=Get_Data_Term_3D_LCC(Mov,Fix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,cache,metric,pix_resolution,fixed_mask,quant)

sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);
min_max=[min(previous,[],1)',max(previous,[],1)'];
min_max=abs(round(min_max));
sx_or=max(min_max(1,:))*2+1;
sy_or=max(min_max(2,:))*2+1;
sz_or=max(min_max(3,:))*2+1;
labels_or.sx=sx_or;
labels_or.sy=sy_or;
labels_or.sz=sz_or;
labels_or.n_each_layer=sx_or*sy_or;
labels_or.nlabels=sx_or*sy_or*sz_or;
labels_or.hx=fix(sx_or/2);
labels_or.hy=fix(sy_or/2);
labels_or.hz=fix(sz_or/2);
labels_or_index=trans_labels(previous,sx_or,sy_or,sz_or);
[r,c,s]=size(Mov);
dim_r=griddim(1);
dim_c=griddim(2);
dim_s=griddim(3);
unary=zeros(dim_r*dim_c*dim_s,labels.sx*labels.sy*labels.sz);
[x_index,y_index,z_index]=Label_Coordinate(labels);
[x_index_or,y_index_or,z_index_or]=Label_Coordinate(labels_or);

for it=1:sx_or*sy_or*sz_or
    cur_index=find(labels_or_index==it);
    if isempty(cur_index)
        continue;
    end
for nlb=1:labels.nlabels
    volmov=ones(r,c,s);%%%set
    x_mov=-x_index(nlb)*sample_space_x-x_index_or(it);
    y_mov=-y_index(nlb)*sample_space_y-y_index_or(it);
    z_mov=-z_index(nlb)*sample_space_z-z_index_or(it);
    if 1+x_mov>=1
        x_init=1+x_mov;
        x_end=r;
        x_init_bf=1;
        x_end_bf=x_end-x_init+1;
    else
        x_init=1;
        x_end=r+x_mov;
        x_end_bf=r;
        x_init_bf=r-(x_end-x_init);
    end
    if 1+y_mov>=1
        y_init=1+y_mov;
        y_end=c;
        y_init_bf=1;
        y_end_bf=y_end-y_init+1;
    else
        y_init=1;
        y_end=c+y_mov;
        y_end_bf=c;
        y_init_bf=c-(y_end-y_init);
    end
    if 1+z_mov>=1
        z_init=1+z_mov;
        z_end=s;
        z_init_bf=1;
        z_end_bf=z_end-z_init+1;
    else
        z_init=1;
        z_end=s+z_mov;
        z_end_bf=s;
        z_init_bf=s-(z_end-z_init);
    end
    volmov(x_init:x_end,y_init:y_end,z_init:z_end)=Mov(x_init_bf:x_end_bf,y_init_bf:y_end_bf,z_init_bf:z_end_bf);
    volfix=Fix;
    em_opts = [];
    if ~isempty(cache)
        em_opts.cache = cache{1};
    end
    LCC = eval_metric_xp(metric, squeeze(volfix(:, :, :, 1)), squeeze(volmov(:,:,:, 1, 1)), pix_resolution, em_opts, fixed_mask);
    for n=1:length(cur_index)
        [i,j,k]=ind2sub([dim_r,dim_c,dim_s],cur_index(n));
                i_1=O_trans_X(i,j,k,1);
                i_2=O_trans_X(i,j,k,2);
                j_1=O_trans_Y(i,j,k,1);
                j_2=O_trans_Y(i,j,k,2);
                k_1=O_trans_Z(i,j,k,1);
                k_2=O_trans_Z(i,j,k,2);
                x_width=i_2-i_1+1;
                y_length=j_2-j_1+1;
                z_hight=k_2-k_1+1;
                v_zone=x_width*y_length*z_hight;
                sad_cg=LCC(i_1:i_2,j_1:j_2,k_1:k_2);
                sad_cg=sum(sad_cg(:))/v_zone;
               unary(cur_index(n),nlb)=sad_cg.*50;

    end
end
end
end


