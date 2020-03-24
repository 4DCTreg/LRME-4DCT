% This is the same version as 'registration_kalman_v1'
clear;

addpath(genpath('../pTVreg-master/mutils/My/'));
addpath(genpath('../pTVreg-master/ptv'));
addpath(genpath('../HOMRF_groupwise/MIND-SSC/'));
addpath(genpath('../HOMRF_groupwise/scr/'));


T_ob = [];
isoPTVpath = 'H:/groupwise_registration/POPI_iso_continus/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,9);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions;
volmov = zeros(all_size);
dvfpath = 'H:/POPI_model/example/4DVF/Parametric/';
obmethod = 'isoPTV';
init_size = info.Dimensions;
        ir = init_size(1);
        ic = init_size(2);
        is = init_size(3);
        [Xgrid,Ygrid,Zgrid] = ndgrid(1:ir,1:ic,1:is);
        Xgrid = Xgrid.*spc_or(1);
        Ygrid = Ygrid.*spc_or(2);
        Zgrid = Zgrid.*spc_or(3);


for phase = 1:9
    LMname = [LMpath,'case',num2str(phase),'0.txt'];

    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
    if strcmp(obmethod,'isoPTV')
        
        fieldFilename = [isoPTVpath,'case',num2str(phase),'_0.mat'];
        load(fieldFilename);
        T_ob(:,:,:,:,phase) = Tptv;
    else
        if strcmp(obmethod,'parametric')
            VFname = [dvfpath,'VF1',num2str(phase),'.vf'];
            T_ob_cur = readvf(VFname);
            T_ob(:,:,:,:,phase) = T_ob_cur.*2;
        end
    end
        
    
end
if strcmp(obmethod,'parametric')
    nrsz = size(T_ob);
    nrk = nrsz(1);nck = nrsz(2);nsk  = nrsz(3);
    T_ob(:,:,:,:,1) = zeros(nrk,nck,nsk,3,'single');
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
           



    lu_orig=[0,0,0];

    N=9;



    filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
    TRE_homrf = zeros(10,1);
    TRE_kalman = zeros(10,1);
    TRE_isoPTV = zeros(10,1);

    for phase = 1:9
        spc = spc_or;
        init_size = info.Dimensions;
        filename = [filepath,num2str(phase),'0_P.raw'];
        volmov = readrawPOPImeta(filename,init_size);
        volfix = readrawPOPImeta('H:/POPI_model/example/4DCT_MetaImage/10_P.raw',init_size);
        volmov = volmov+1024;
        volfix = volfix+1024;
        volmov = img_thr(volmov, min(volmov(:)), max(volmov(:)), 1);
        volfix = img_thr(volfix, min(volfix(:)), max(volfix(:)), 1);

        LMpath = 'H:/POPI_model/example/4DLandmarks/';
        LMname = [LMpath,'case',num2str(phase),'0.txt'];
        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);
        pts_fix = round(pts_fix_or./koef);
        pts_mov = round(pts_mov_or./koef);
        crop_v = [1,init_size(1);1, init_size(2);1,init_size(3)];

        volmov = crop_data(volmov, crop_v);
        volfix = crop_data(volfix, crop_v);

        useTop = 0;
        resize =1;
        if resize==1
            bszv = size(volmov);
            spc_tmp = [0.5, 0.5, 0.5];
            volfix = volresize(volfix, round(bszv .* spc .* spc_tmp), 1);
            volmov = volresize(volmov, round(bszv .* spc .* spc_tmp), 1);
            spc = [1,1,1] ./ spc_tmp;
        end
        if resize==1
          
            Tptv_rsz = cat(4, volresize(T_ob(:,:,:,1,phase), bszv), volresize(T_ob(:,:,:,2,phase), bszv), volresize(T_ob(:,:,:,3,phase), bszv)); 
            voldef = volresize(volmov,bszv);
        end
        
       
%         [~,Knots_n_kalman_rsz_all] = uncrop_data(voldef,Knots_n_kalman_rsz,crop_v,init_size);
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

        
    end


