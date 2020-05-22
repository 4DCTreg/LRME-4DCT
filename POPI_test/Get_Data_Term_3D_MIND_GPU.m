
function unary_tot=Get_Data_Term_3D_MIND_GPU(Mov_mind,Fix_mind,mov_maks,fix_mask,griddim,O_trans_X,...
    O_trans_Y,O_trans_Z,labels,level,previous,quant,cur_grid_space)
weight = 100;
[r,c,s,R]=size(Mov_mind);
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
labels_or.sx=sx_or;
labels_or.sy=sy_or;
labels_or.sz=sz_or;
labels_or.n_each_layer=sx_or*sy_or;
labels_or.nlabels=sx_or*sy_or*sz_or;
labels_or.hx=fix(sx_or/2);
labels_or.hy=fix(sy_or/2);
labels_or.hz=fix(sz_or/2);
labels_or_index=trans_labels(previous,sx_or,sy_or,sz_or);
[x_index_or,y_index_or,z_index_or]=Label_Coordinate(labels_or);


control_size = size(O_trans_X);
grid_space = cur_grid_space;
grid_size = grid_space(1)*grid_space(2)*grid_space(3);
N_control_size = control_size(1)*control_size(2)*control_size(3);
% filt = ones(N_control_size,R,grid_space(1),grid_space(2),grid_space(3),'single');
x_st = O_trans_X(1,1,1,1);
x_end = O_trans_X(end,1,1,2);
y_st = O_trans_Y(1,1,1,1);
y_end = O_trans_Y(1,end,1,2);
z_st = O_trans_Z(1,1,1,1);
z_end = O_trans_Z(1,1,end,2);

Fix_sum_all = zeros(N_control_size,R);
Mov_sum_all = zeros(N_control_size,R);

division_x = grid_space(1)*ones(1,griddim(1));
division_y = grid_space(2)*ones(1,griddim(2));
division_z = grid_space(3)*ones(1,griddim(3));
Fix_crop = Fix_mind(x_st:x_end,y_st:y_end,z_st:z_end,:);
Mov_crop = Mov_mind(x_st:x_end,y_st:y_end,z_st:z_end,:);

for channel = 1:R
    division_fix = mat2cell(Fix_crop(:,:,:,channel).^2,division_x,division_y,division_z);
    
    sum_fix = cellfun(@(x)sum(x(:)),division_fix);
    Fix_sum_all(:,channel) = reshape(sum_fix,N_control_size,1);

end
Fix_sum_all = mean(Fix_sum_all,2);
Mov_sum_all = mean(Mov_sum_all,2);
Fix_sum_all = reshape(Fix_sum_all,1,1,1,N_control_size);
Mov_sum_all = reshape(Mov_sum_all,1,1,1,N_control_size);
Fix_sum = repmat(Fix_sum_all,[labels.sx,labels.sy,labels.sz,1])./grid_size;

refBlock = zeros(N_control_size,R,grid_space(1),grid_space(2),grid_space(3));
tarSearch = zeros(N_control_size,R,grid_space(1)+2*labels.hx,grid_space(2)+2*labels.hy,grid_space(3)+2*labels.hz);
toc

tic
for R_index = 1:R
    volmovl= single(uncrop_data_xp_3d(Mov_mind(:,:,:,R_index), crop_v, size_orig,weight));
    volfixl= single(uncrop_data_xp_3d(Fix_mind(:,:,:,R_index), crop_v, size_orig,weight));
for i = 1:griddim(1)
    for j = 1:griddim(2)
        for k = 1:griddim(3)
            num = sub2ind(griddim,i,j,k);
            previous_label=labels_or_index(num);
            x_start = O_trans_X(i,j,k,1);
            x_end = O_trans_X(i,j,k,2);
            y_start = O_trans_Y(i,j,k,1);
            y_end = O_trans_Y(i,j,k,2);
            z_start = O_trans_Z(i,j,k,1);
            z_end = O_trans_Z(i,j,k,2);
            
            refBlock(num,R_index,:,:,:) = volfixl(x_start+add_hx:x_end+add_hx,y_start+add_hy:y_end+add_hy,...
                z_start+add_hz:z_end+add_hz);
            tarSearch(num,R_index,:,:,:) = volmovl(x_start+x_index_or(previous_label)+add_hx-labels.hx:x_end+x_index_or(previous_label)+add_hx+labels.hx,...
                y_start+y_index_or(previous_label)+add_hy-labels.hy:y_end+y_index_or(previous_label)+add_hy+labels.hy,...
                z_start+z_index_or(previous_label)+add_hz-labels.hz:z_end+z_index_or(previous_label)+add_hz+labels.hz);
        end
    end
end
end
toc
m_filter = single(ones(N_control_size,grid_space(1),grid_space(2),grid_space(3)));
tarSearch = single(tarSearch);
refBlock = single(refBlock);
[gy_filter,tar_sum] = conv3d_gpu(tarSearch,refBlock,m_filter,N_control_size);
gy_filter = squeeze(mean(gy_filter,1));
Mov_sum = squeeze(mean(tar_sum,1));
Mov_sum = permute(Mov_sum,[2,3,4,1]);

gy_filter = permute(gy_filter,[2,3,4,1]);

distSSC = Mov_sum./(grid_size)-2.*gy_filter./(grid_size)+Fix_sum;

unary_tot = reshape(abs(distSSC),labels.sx*labels.sy*labels.sz,N_control_size);
unary_tot = sqrt(unary_tot');
end


