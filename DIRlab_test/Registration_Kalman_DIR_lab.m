
clear;
addpath(genpath('../DIRpackage/'));
addpath(genpath('../HOMRF_continue/MIND-SSC/'));
addpath(genpath('../HOMRF_continue/scr/'));

basepth = '../data_prj/dir_dataset/DIR_files/'; 
HOMRF_all = zeros(10,6);
Kalman_all = zeros(10,6);
isoPTV_all = zeros(10,6);
HOMRF_all_std = zeros(10,6);
Kalman_all_std = zeros(10,6);
isoPTV_all_std = zeros(10,6);
for idx = 1:10
%% load the observation value

    filePathiso_rsz = '../Observation Data/DIR-lab/';
    T_ob = [];   
    for phase = 0:5
        fieldFilename = [filePathiso_rsz,'case',num2str(idx),'/','T',num2str(phase),'_0','.mat'];
        load(fieldFilename);
        T_ob(:,:,:,:,phase+1) = Tptv;
    end

 %% set the kalman parameter
    
    nrsz = size(T_ob);
    nr = nrsz(1);nc = nrsz(2);ns  = nrsz(3);
    N=6;
    nx = [1:nr];
    ny = [1:nc]; 
    nz = [1:ns];
    [Xgrid, Ygrid, Zgrid] = ndgrid(nx,ny,nz);
    Z = zeros(1,N);
    Px = zeros(nr,nc,ns,N,'single')+1; 
    Py = zeros(nr,nc,ns,N,'single')+1;
    Pz = zeros(nr,nc,ns,N,'single')+1;
    Q = 0.25;
    R = 0.01;  
    H=ones(nr,nc,ns,'single');
    I=ones(nr,nc,ns,'single'); 
    xd_kf = zeros(nr,nc,ns,N,'single');
    yd_kf = zeros(nr,nc,ns,N,'single');
    zd_kf = zeros(nr,nc,ns,N,'single');
    xd_kf(:,:,:,1) = Xgrid;
    yd_kf(:,:,:,1) = Ygrid;
    zd_kf(:,:,:,1) = Zgrid;
    fix_all = zeros(nr,nc,ns,3,'single');
    fix_all(:,:,:,1) = xd_kf(:,:,:,1);
    fix_all(:,:,:,2) = yd_kf(:,:,:,1);
    fix_all(:,:,:,3) = zd_kf(:,:,:,1);
  
 %% HOMRF continus registration and kalman filter
    for phase = 0:5
        pts_struct = DIR_get_all_points_for_the_case(idx, basepth);
        [volmov, spc] = read_DIR_volume_4dCT(idx, phase, basepth);
        volmov = double(volmov);
%         seg_mask_mov = parenchyma_extration(volmov);
        pts_mov = pts_struct.smp{phase+1}.pts;
        pts_mov_max = pts_struct.smp{6}.pts;
        pts_fix = pts_struct.smp{1}.pts;
        [volfix, spc] = read_DIR_volume_4dCT(idx, 0, basepth);
        volfix = double(volfix); 
%         seg_mask_fix = parenchyma_extration(volfix);
        init_size=size(volmov);
        volmov=img_thr(volmov,80,900,1);
        volfix=img_thr(volfix,80,900,1);
        min_max1 = [ min(pts_mov_max, [], 1)', max(pts_mov_max, [], 1)'];
        min_max2 = [ min(pts_fix, [], 1)', max(pts_fix, [], 1)'];
        min_max = [ min(min_max1(:, 1), min_max2(:, 1)), max(min_max1(:, 2), min_max2(:, 2))];
        d = [10, 10, 5];
        crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];

        volmov = crop_data(volmov, crop_v);
        volfix = crop_data(volfix, crop_v);
        useTop = 0;

        spc_orig = spc;
        resize = 1;
        if resize==1
            bszv = size(volmov);
            spc_tmp = [1, 1, 1];
            volfix = volresize(volfix, round(bszv .* spc .* spc_tmp), 1);
            volmov = volresize(volmov, round(bszv .* spc .* spc_tmp), 1);
            spc = [1,1,1] ./ spc_tmp;
        end
        parameter=homrf_get_parameter(idx,useTop,spc); 
        parameter.nlevel = 1;
        cur_grid_space = parameter.grid_space(parameter.nlevel,:);
        
        [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
        
        griddim = [size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
        
        if phase == 0
            Knots_n = zeros(griddim(1),griddim(2),griddim(3),3);
        else
            Knots_n = Knots_n_kalman_sparse_rsz_K;
        end
 
        [Knots_n_tot,Knots_n]=register_homrf_kalman(parameter,volmov,volfix,crop_v,init_size,pts_fix,pts_mov,Knots_n,bszv); 
            
        if parameter.resize ==1
            Tptv_rsz = cat(4, volresize(Knots_n_tot(:,:,:,1), bszv), ...
       volresize(Knots_n_tot(:,:,:,2), bszv), volresize(Knots_n_tot(:,:,:,3), bszv));  
            voldef = volresize(volmov,bszv);
        end
        [~, Knots_n_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);
   
        [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints_homrf(pts_mov, pts_fix, Knots_n_rsz, spc_orig,parameter.resize, []);

        HOMRF_all(idx,phase+1) = TRE_phys;
        HOMRF_all_std(idx,phase+1) = TREstd_phys;
        
        if phase > 0
            T_iso = cat(4, volresize(T_ob(:,:,:,1,phase+1), bszv), ...
       volresize(T_ob(:,:,:,2,phase+1), bszv), volresize(T_ob(:,:,:,3,phase+1), bszv)); 
            [~, T_iso_rsz] = uncrop_data(voldef, T_iso, crop_v, init_size);
            [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_mov, pts_fix, T_iso_rsz, spc_orig, []);
            isoPTV_all(idx,phase+1) = TRE_phys;
            isoPTV_all_std(idx,phase+1) = TREstd_phys;
        end
    
 %% kalman filter
        if phase > 0
            es_all = Knots_n_tot;
            es_mov_all = es_all+fix_all;
            move_all = T_ob(:,:,:,:,phase+1);
            move_all = move_all+fix_all;
            Hx = ones(nr,nc,ns,'single');
            Hy = ones(nr,nc,ns,'single');
            Hz = ones(nr,nc,ns,'single');
            
            Fx_cur = es_mov_all(:,:,:,1)./xd_kf(:,:,:,phase);
            Fy_cur = es_mov_all(:,:,:,2)./yd_kf(:,:,:,phase);
            Fz_cur = es_mov_all(:,:,:,3)./zd_kf(:,:,:,phase);
            xd_pre = Fx_cur.*xd_kf(:,:,:,phase)+Q;
            yd_pre = Fy_cur.*yd_kf(:,:,:,phase)+Q;
            zd_pre = Fz_cur.*zd_kf(:,:,:,phase)+Q;
            
            p_pre_x = Fx_cur.*Px(:,:,:,phase).*Fx_cur+Q;
            p_pre_y = Fy_cur.*Py(:,:,:,phase).*Fy_cur+Q;
            p_pre_z = Fz_cur.*Pz(:,:,:,phase).*Fz_cur+Q;
            
            Kg_x = p_pre_x.*Hx.*1./(Hx.*p_pre_x.*Hx+R);
            Kg_y = p_pre_y.*Hy.*1./(Hy.*p_pre_y.*Hy+R);
            Kg_z = p_pre_z.*Hz.*1./(Hz.*p_pre_z.*Hz+R);
            
            ex = move_all(:,:,:,1)-Hx.*xd_pre;
            ey = move_all(:,:,:,2)-Hy.*yd_pre;
            ez = move_all(:,:,:,3)-Hz.*zd_pre;
            
            xd_kf(:,:,:,phase+1) = xd_pre+Kg_x.*ex;
            yd_kf(:,:,:,phase+1) = yd_pre+Kg_y.*ey;
            zd_kf(:,:,:,phase+1) = zd_pre+Kg_z.*ez;
            
            xd_kf_cur = xd_kf(:,:,:,phase+1);
            yd_kf_cur = yd_kf(:,:,:,phase+1);
            zd_kf_cur = zd_kf(:,:,:,phase+1);
            
            tx = ismissing(xd_kf_cur);
            fx = find(tx==1);
            if ~isempty(fx)
                xd_kf_cur(fx) = xd_kf_cur(fx+1);
            end
            ty = ismissing(yd_kf_cur);
            fy = find(ty==1);
            if ~isempty(fy)
                yd_kf_cur(fy) = yd_kf_cur(fy+1);
            end
            tz = ismissing(zd_kf_cur);
            fz = find(tz==1);
            if ~isempty(fz)
                zd_kf_cur(fz) = zd_kf_cur(fz+1);
            end
            
            xd_kf(:,:,:,phase+1) = xd_kf_cur;
            yd_kf(:,:,:,phase+1) = yd_kf_cur;
            zd_kf(:,:,:,phase+1) = zd_kf_cur;
            
            Px(:,:,:,phase+1) = (I-Kg_x.*Hx).*p_pre_x;
            Py(:,:,:,phase+1) = (I-Kg_y.*Hy).*p_pre_y;
            Pz(:,:,:,phase+1) = (I-Kg_z.*Hz).*p_pre_z;
        end
        
        Knots_n_kalman = [];  
        Knots_n_kalman(:,:,:,1) = xd_kf(:,:,:,phase+1)-fix_all(:,:,:,1);
        Knots_n_kalman(:,:,:,2) = yd_kf(:,:,:,phase+1)-fix_all(:,:,:,2);
        Knots_n_kalman(:,:,:,3) = zd_kf(:,:,:,phase+1)-fix_all(:,:,:,3);
            
        if parameter.resize==1
            Knots_n_kalman_rsz = cat(4, volresize(Knots_n_kalman(:,:,:,1), bszv), volresize(Knots_n_kalman(:,:,:,2), bszv), volresize(Knots_n_kalman(:,:,:,3), bszv));   
        end
        
        [~,Knots_n_kalman_rsz_all] = uncrop_data(voldef,Knots_n_kalman_rsz,crop_v,init_size);

        Knots_n_kalman_sparse_rsz_K = dense_to_Knots(Knots_n_kalman,O_trans_X,O_trans_Y,O_trans_Z);
        
        [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints_homrf(pts_mov, pts_fix, Knots_n_kalman_rsz_all, spc_orig, parameter.resize, []);
        Kalman_all(idx,phase+1) = TRE_phys;
        Kalman_all_std(idx,phase+1) = TREstd_phys;
   
    end

end
