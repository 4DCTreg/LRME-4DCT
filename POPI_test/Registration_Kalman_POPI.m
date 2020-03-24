% This is the tese of the POPI dataset
clear;

addpath(genpath('../DIRpackage/'));
addpath(genpath('../HOMRF_continue/MIND-SSC/'));
addpath(genpath('../HOMRF_continue/scr/'));

isoPTVpath =  '../Observation Data/POPI/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,10);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions; 
volmov = zeros(all_size);

obmethod = 'isoPTV';
init_size = info.Dimensions;
ir = init_size(1);
ic = init_size(2);
is = init_size(3);
[Xgrid,Ygrid,Zgrid] = ndgrid(1:ir,1:ic,1:is);
Xgrid = Xgrid.*spc_or(1);
Ygrid = Ygrid.*spc_or(2);
Zgrid = Zgrid.*spc_or(3);

T_ob = [];
for phase = 1:10
    if phase == 10
        LMname = [LMpath,'case00.txt'];
    else
        LMname = [LMpath,'case',num2str(phase),'0.txt'];
    end
    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
    if strcmp(obmethod,'isoPTV')
        if phase == 10
            fieldFilename = [isoPTVpath,'case0_0.mat'];
        else
            fieldFilename = [isoPTVpath,'case',num2str(phase),'_0.mat'];
        end
        load(fieldFilename);
        T_ob(:,:,:,:,phase) = Tptv;

    end 
end


min1 = min(min(squeeze(pts_all_phase(:,1,:))));
max1 = max(max(squeeze(pts_all_phase(:,1,:))));
min2 = min(min(squeeze(pts_all_phase(:,2,:))));
max2 = max(max(squeeze(pts_all_phase(:,2,:))));
min3 = min(min(squeeze(pts_all_phase(:,3,:))));
max3 = max(max(squeeze(pts_all_phase(:,3,:))));
d = [10, 10, 5];
min_max = [min1,max1;min2,max2;min3,max3];
crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];
           
nrsz = size(T_ob);
nrk = nrsz(1);nck = nrsz(2);nsk  = nrsz(3);
nx = [1:nrk];
ny = [1:nck]; 
nz = [1:nsk];
    [NX, NY, NZ] = ndgrid(nx,ny,nz);
    NX = NX.*spc_or(1);
    NY = NY.*spc_or(2);
    NZ = NZ.*spc_or(3);
    N = 10;

    Z = zeros(1,N);
    Px = zeros(nrk,nck,nsk,N,'single')+1; 
    Py = zeros(nrk,nck,nsk,N,'single')+1;
    Pz = zeros(nrk,nck,nsk,N,'single')+1;
    Q = 0.01;
    R = 0.01;

    H=ones(nrk,nck,nsk,'single');
    I=ones(nrk,nck,nsk,'single'); 
    
    xd_kf = zeros(nrk,nck,nsk,N,'single');
    yd_kf = zeros(nrk,nck,nsk,N,'single');
    zd_kf = zeros(nrk,nck,nsk,N,'single');
    
    xd_kf(:,:,:,1) = NX;
    yd_kf(:,:,:,1) = NY;
    zd_kf(:,:,:,1) = NZ;
    
    fix_all = zeros(nrk,nck,nsk,3,'single');
    fix_all(:,:,:,1) = xd_kf(:,:,:,1);
    fix_all(:,:,:,2) = yd_kf(:,:,:,1);
    fix_all(:,:,:,3) = zd_kf(:,:,:,1);
    filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
    TRE_homrf = zeros(10,1);
    TRE_kalman = zeros(10,1);
    TRE_isoPTV = zeros(10,1);
    STD_homrf = zeros(10,1);
    STD_kalman = zeros(10,1);
    STD_isoPTV = zeros(10,1);
    maskpath = 'H:/POPI_model/example/4DMask_metalmage/';
    LMpath = 'H:/POPI_model/example/4DLandmarks/';
    for phase = 1:10
        
        if phase == 10
            Movmaskname = [maskpath,'00-air-body-lungs.raw'];
            filename = [filepath,'00_P.raw'];
            LMname = [LMpath,'case00.txt'];
        else
            Movmaskname = [maskpath,num2str(phase),'0-air-body-lungs.raw'];
            filename = [filepath,num2str(phase),'0_P.raw'];
            LMname = [LMpath,'case',num2str(phase),'0.txt'];
        end
        Fixmaskname = [maskpath,'10-air-body-lungs.raw'];
        mov_mask = readrawPOPImask(Movmaskname,init_size);
        fix_mask = readrawPOPImask(Fixmaskname,init_size);
        mov_mask(mov_mask~=2) = 0;
        fix_mask(fix_mask~=2) = 0;
        mov_mask(mov_mask==2) = 1;
        fix_mask(fix_mask==2) = 1;
        spc = spc_or;
        init_size = info.Dimensions;
        
        volmov = readrawPOPImeta(filename,init_size);
        volfix = readrawPOPImeta('H:/POPI_model/example/4DCT_MetaImage/10_P.raw',init_size);
        volmov = volmov+1024;
        volfix = volfix+1024;
        volmov = img_thr(volmov, min(volmov(:)), max(volmov(:)), 1);
        volfix = img_thr(volfix, min(volfix(:)), max(volfix(:)), 1);

        
        
        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);
        pts_fix = round(pts_fix_or./koef);
        pts_mov = round(pts_mov_or./koef);
%         crop_v = [1,init_size(1);1, init_size(2);1,init_size(3)];

        volmov = crop_data(volmov, crop_v);
        volfix = crop_data(volfix, crop_v);
        mov_mask = crop_data(mov_mask,crop_v);
        fix_mask = crop_data(fix_mask,crop_v);

        useTop = 1;
        resize =1;
        if resize==1
            bszv = size(volmov);
            spc_tmp = [1, 1, 1];
            volfix = volresize(volfix, round(bszv .* spc .* spc_tmp), 1);
            volmov = volresize(volmov, round(bszv .* spc .* spc_tmp), 1);
            mov_mask = imresize3(mov_mask, round(bszv .* spc .* spc_tmp), 'nearest');
            fix_mask = imresize3(fix_mask, round(bszv .* spc .* spc_tmp), 'nearest');
            spc = [1,1,1] ./ spc_tmp;
        end
        parameter=homrf_get_POPI_parameter(useTop,spc);
        cur_grid_space = parameter.grid_space(1,:);
        
        [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
        griddim = [size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
        
        if phase == 1
            Knots_n = zeros(griddim(1),griddim(2),griddim(3),3);
        else
            Knots_n = Knots_n_kalman_sparse_rsz_K;
        end
 
        [Knots_n_tot,Knots_n]=register_homrf_kalman_POPI_mask(parameter,volmov,volfix,mov_mask,fix_mask,...
            init_size,pts_fix_or,pts_mov_or,Knots_n,bszv); 
    
        Knots_n_pre = Knots_n;

        if parameter.resize ==1
            Tptv_rsz = cat(4, volresize(Knots_n_tot(:,:,:,1), bszv), ...
       volresize(Knots_n_tot(:,:,:,2), bszv), volresize(Knots_n_tot(:,:,:,3), bszv));  
            voldef = volresize(volmov,bszv);
        end
        [~, Knots_n_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);

        HOMRFgriddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,1));
        HOMRFgriddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,2));
        HOMRFgriddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,3));
        displhomrf(:,1) = HOMRFgriddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,2) = HOMRFgriddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,3) = HOMRFgriddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMsHOMRF = displhomrf+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMsHOMRF-pts_mov_or)).^2,2)));
        TRE_homrf(phase) = tre;
        STD_homrf(phase) = std(sqrt( sum(( (warpedLMsHOMRF-pts_mov_or)).^2, 2) ));

        if phase > 1
            es_all = Knots_n_tot;

            es_mov_all = es_all+fix_all;
            move_all = T_ob(:,:,:,:,phase);
            move_all = move_all+fix_all;
            Fx_cur = es_mov_all(:,:,:,1)./xd_kf(:,:,:,phase-1);
            Fy_cur = es_mov_all(:,:,:,2)./yd_kf(:,:,:,phase-1);
            Fz_cur = es_mov_all(:,:,:,3)./zd_kf(:,:,:,phase-1);
            xd_pre = Fx_cur.*xd_kf(:,:,:,phase-1)+Q;
            yd_pre = Fy_cur.*yd_kf(:,:,:,phase-1)+Q;
            zd_pre = Fz_cur.*zd_kf(:,:,:,phase-1)+Q;
            
            p_pre_x = Fx_cur.*Px(:,:,:,phase-1).*Fx_cur+Q;
            p_pre_y = Fy_cur.*Py(:,:,:,phase-1).*Fy_cur+Q;
            p_pre_z = Fz_cur.*Pz(:,:,:,phase-1).*Fz_cur+Q;
            
            Kg_x = p_pre_x.*H.*1./(H.*p_pre_x.*H+R);
            Kg_y = p_pre_y.*H.*1./(H.*p_pre_y.*H+R);
            Kg_z = p_pre_z.*H.*1./(H.*p_pre_z.*H+R);
            
            ex = move_all(:,:,:,1)-H.*xd_pre;
            ey = move_all(:,:,:,2)-H.*yd_pre;
            ez = move_all(:,:,:,3)-H.*zd_pre;
            
            xd_kf(:,:,:,phase) = xd_pre+Kg_x.*ex;
            yd_kf(:,:,:,phase) = yd_pre+Kg_y.*ey;
            zd_kf(:,:,:,phase) = zd_pre+Kg_z.*ez;
            
            Px(:,:,:,phase) = (I-Kg_x.*H).*p_pre_x;
            Py(:,:,:,phase) = (I-Kg_y.*H).*p_pre_y;
            Pz(:,:,:,phase) = (I-Kg_z.*H).*p_pre_z;
        end
        
        Knots_n_kalman = [];      
        Knots_n_kalman(:,:,:,1) = xd_kf(:,:,:,phase)-fix_all(:,:,:,1);
        Knots_n_kalman(:,:,:,2) = yd_kf(:,:,:,phase)-fix_all(:,:,:,2);
        Knots_n_kalman(:,:,:,3) = zd_kf(:,:,:,phase)-fix_all(:,:,:,3);
            
        if parameter.resize==1
            Knots_n_kalman_rsz = cat(4, volresize(Knots_n_kalman(:,:,:,1), bszv), volresize(Knots_n_kalman(:,:,:,2), bszv), volresize(Knots_n_kalman(:,:,:,3), bszv)); 
            Tptv_rsz = cat(4, volresize(T_ob(:,:,:,1,phase), bszv), volresize(T_ob(:,:,:,2,phase), bszv), volresize(T_ob(:,:,:,3,phase), bszv)); 
        end
        
        [~,Knots_n_kalman_rsz_all] = uncrop_data(voldef,Knots_n_kalman_rsz,crop_v,init_size);
        [~,Tptv_rsz] = uncrop_data(voldef,Tptv_rsz,crop_v,init_size);
        griddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Tptv_rsz(:,:,:,1));
        griddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Tptv_rsz(:,:,:,2));
        griddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Tptv_rsz(:,:,:,3));
        displ(:,1) = griddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,2) = griddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,3) = griddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMs = displ+pts_fix_or;
        tre = mean(sqrt(sum(((warpedLMs-pts_mov_or)).^2,2)));
        TRE_isoPTV(phase) = tre;
        STD_isoPTV(phase) = std(sqrt( sum(( (warpedLMs-pts_mov_or)).^2, 2) ));

        Knots_n_kalman_sparse_rsz_K = dense_to_Knots(Knots_n_kalman,O_trans_X,O_trans_Y,O_trans_Z);

        griddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,1));
        griddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,2));
        griddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,3));
        displ(:,1) = griddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,2) = griddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,3) = griddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMs = displ+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMs-pts_mov_or)).^2,2)));
        TRE_kalman(phase) = tre;
        STD_kalman(phase) = std(sqrt( sum(( (warpedLMs-pts_mov_or)).^2, 2) ));
      
    end


