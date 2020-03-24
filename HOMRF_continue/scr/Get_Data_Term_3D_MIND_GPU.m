
function unary_tot=Get_Data_Term_3D_MIND_GPU(Mov_mind,Fix_mind,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant)


grid_space = [8 8 8];
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
dim_r=griddim(1);
dim_c=griddim(2);
dim_s=griddim(3);
[x_index,y_index,z_index]=Label_Coordinate(labels);
[x_index_or,y_index_or,z_index_or]=Label_Coordinate(labels_or);

weight=100;
unary_tot=zeros(dim_r*dim_c*dim_s,labels.nlabels);
O_trans_X_1=O_trans_X(:,:,:,1);
O_trans_X_2=O_trans_X(:,:,:,2);
O_trans_Y_1=O_trans_Y(:,:,:,1);
O_trans_Y_2=O_trans_Y(:,:,:,2);
O_trans_Z_1=O_trans_Z(:,:,:,1);
O_trans_Z_2=O_trans_Z(:,:,:,2);
nlabels=labels.nlabels;
O_trans_X_1_ex=kron(O_trans_X_1(:),ones(nlabels,1));
O_trans_X_2_ex=kron(O_trans_X_2(:),ones(nlabels,1));
O_trans_Y_1_ex=kron(O_trans_Y_1(:),ones(nlabels,1));
O_trans_Y_2_ex=kron(O_trans_Y_2(:),ones(nlabels,1));
O_trans_Z_1_ex=kron(O_trans_Z_1(:),ones(nlabels,1));
O_trans_Z_2_ex=kron(O_trans_Z_2(:),ones(nlabels,1));
x_index_ex=repmat(x_index,[dim_r*dim_c*dim_s,1]);
y_index_ex=repmat(y_index,[dim_r*dim_c*dim_s,1]);
z_index_ex=repmat(z_index,[dim_r*dim_c*dim_s,1]);
labels_or_index_ex=kron(labels_or_index,ones(nlabels,1));


blockSize_x = grid_space(1);
blockSize_y = grid_space(2);
blockSize_z = grid_space(3);
Npoints = prod(griddim);
mindDim = R;

refBlock = zeros(blockSize_x,blockSize_y,blockSize_z,mindDim,Npoints,'single');

searchSize_x = add_hx*2+1;
searchSize_y = add_hy*2+1;
searchSize_z = add_hz*2+1;
tarSearch = zeros(searchSize_x+blockSize_x-1,searchSize_y+blockSize_y-1,searchSize_z+blockSize_z-1,mindDim,Npoints,'single'); 
div = (blockSize_x.*blockSize_y*blockSize_z.*mindDim);
for R_index = 1:R
    volmovl= single(uncrop_data_xp_3d(Mov_mind(:,:,:,R_index), crop_v, size_orig,weight));
    volfixl= single(uncrop_data_xp_3d(Fix_mind(:,:,:,R_index), crop_v, size_orig,weight));
for i = 1:griddim(1)
    for j = 1:griddim(2)
        for k = 1:griddim(3)
            num = sub2ind(griddim,i,j,k);
            x_start = O_trans_X(i,j,k,1);
            x_end = O_trans_X(i,j,k,2);
            y_start = O_trans_Y(i,j,k,1);
            y_end = O_trans_Y(i,j,k,2);
            z_start = O_trans_Z(i,j,k,1);
            z_end = O_trans_Z(i,j,k,2);
            
            refBlock(:,:,:,R_index,num) = volfixl(x_start+add_hx:x_end+add_hx,y_start+add_hy:y_end+add_hy,...
                z_start+add_hz:z_end+add_hz);
            tarSearch(:,:,:,R_index,num) = volmovl(x_start+add_hx-add_hx:x_end+add_hx+add_hx,y_start+add_hy-add_hy:y_end+add_hy+add_hy,...
                z_start+add_hz-add_hz:z_end+add_hz+add_hz);
        end
    end
end
end
filt = ones(blockSize_x*blockSize_y,blockSize_z,mindDim,Npoints,'single');
gfilt = gpuArray(filt);
gRefBlock = gpuArray(reshape(refBlock,blockSize_x*blockSize_y,blockSize_z,[]));
% filt = ones(blockSize_x,blockSize_y,blockSize_z);
% gfilt = gpuArray(filt);
% tic
% for n = 1:Npoints
% gRefsum2 = imfilter(gRefBlock(:,:,:,n).^2,gfilt);
% end
% toc
% t = toc;
gRefSum2 = vl_nnconv(gRefBlock.^2,gfilt,[]);
grefSum2 = repmat(gRefSum2,[searchSize_x,searchSize_y,searchSize_z,1]);
grefSum2 = reshape(grefSum2,searchSize_x,searchSize_y,searchSize_z,Npoints);
gRefBlockr = double(gather(reshape(gRefBlock,blockSize_x,blockSize_y,blockSize_z,mindDim,Npoints)));
gTarSearch = gpuArray(reshape(tarSearch,(blockSize_x+searchSize_x-1)*(blockSize_y+searchSize_y-1),blockSize_z+searchSize_z-1,[]));
gTarSearchr = double(reshape(tarSearch,(blockSize_x+searchSize_x-1),(blockSize_y+searchSize_y-1),blockSize_z+searchSize_z-1,mindDim,Npoints));
y_filt = zeros(searchSize_x,searchSize_y,searchSize_z,Npoints);


    parfor n = 1:Npoints
        cur2 =zeros(searchSize_x,searchSize_y,searchSize_z);
        for nR = 1:R
            
            cur = imfilter(gTarSearchr(:,:,:,nR,n),gRefBlockr(:,:,:,nR,n));
            cur1 = cur(4:14,4:14,4:14);
            cur2 = cur1+cur2;

        end
        y_filt(:,:,:,n) = cur2/mindDim;
    end


filt = ones(blockSize_x,blockSize_y,blockSize_z);
gtargetSum2 = zeros(searchSize_x,searchSize_y,searchSize_z,Npoints);

    parfor n = 1:Npoints
        cur2 =zeros(searchSize_x,searchSize_y,searchSize_z);
        for nR = 1:R
            
            cur = imfilter(gTarSearchr(:,:,:,nR,n),filt);
            cur1 = cur(4:14,4:14,4:14);
            cur2 = cur1+cur2;

        end
        gtargetSum2(:,:,:,n) = cur2/mindDim;
    end


% tic
%     parfor n = 1:Npoints*R
% 
%             
%             cur = imfilter(gTarSearchr(:,:,:,n),gRefBlockr(:,:,:,n));
%             cur1 = cur(4:14,4:14,4:14);
%             y_filt(:,:,:,n)=cur1;
%      
% 
%      end
% 
% toc
% t = toc;

% gy_filt = vl_nnconv(gTarSearch,reshape(gRefBlock,blockSize_x*blockSize_y,blockSize_z,mindDim,Npoints),[]);
% y_filt =reshape(gy_filt,searchSize_x,searchSize_y,searchSize_z,Npoints);

% gtargetSum2 = vl_nnconv(gTarSearch.^2,gfilt,[]);
% gtargetSum2 = reshape(gtargetSum2,searchSize_x,searchSize_y,searchSize_z,Npoints);
gdistSSC = gtargetSum2./div-2.*y_filt./div+grefSum2./div;
distSSC = gather(gdistSSC);

unary_tot = reshape(distSSC,searchSize_x*searchSize_y*searchSize_z,Npoints);
unary_tot = unary_tot';
end


