
function unary_tot=Get_Data_Term_3D(volmov,volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant,metric,cur_grid_space)
% sigmas=[0.285,0.285,0.285];
% volmov = imgaussfilt3(volmov, sigmas);
% volfix = imgaussfilt3(volfix, sigmas);
if strcmp(metric, 'MIND')
    mind_volmov=MIND_descriptor(volmov);
    mind_volfix=MIND_descriptor(volfix);
    mind_volmov=mind_volmov.*100;
    mind_volfix=mind_volfix.*100;
    unary_tot=Get_Data_Term_3D_MIND(mind_volmov,mind_volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant);
%     unary_tot=Get_Data_Term_3D_MIND_parallel(mind_volmov,mind_volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant);
%     unary_tot=Get_Data_Term_3D_MIND_GPU(mind_volmov,mind_volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant);
%     unary_tot=Get_Data_Term_3D_MIND_filter(mind_volmov,mind_volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant,cur_grid_space);
else
    if strcmp(metric, 'SAD')
        unary_tot=Get_Data_Term_3D_SAD(volmov,volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,quant);
    else
        if strcmp(metric, 'LCC')
            
            Nd=3;
            Nch=1;
            fixed_mask=[];
            internal_dtype = 'CPU_double';
            loc_cc_approximate=0;
            scale_metric_param=1;
            cur_pix_resolution=[1,1,1];
            pix_resolution=[1,1,1];
            metric_lcc='loc_cc_fftn_gpu';
            metric_param = [1,1,1] * 2.1;
            border_mask=5;
            imsz=size(volmov);
            if ~isempty(border_mask) && border_mask > 0
                bmask = ones(imsz);
                bmask(1:border_mask, :, :) = 0; bmask(end-border_mask+1:end, :, :) = 0;
                bmask(:, 1:border_mask, :) = 0; bmask(:, end-border_mask+1:end, :) = 0; 
                if Nd == 3
                    bmask(:, :, 1:border_mask) = 0; bmask(:, :, end-border_mask+1:end) = 0;
                end
                if isempty(fixed_mask)
                    fixed_mask = bmask;
                else
                    fixed_mask = fixed_mask .* bmask;
                end
            end
            metric_param_pix = [];
            if strcmp(metric_lcc, 'loc_cc_fftn_gpu') || strcmp(metric_lcc, 'loc_cc_fftn_gpu_single') || strcmp(metric_lcc, 'loc_cc_fftn')|| strcmp(metric_lcc, 'loc_cc_fftn_single')
               if scale_metric_param
                   metric_param_pix = metric_param ./ cur_pix_resolution;
               else
                   metric_param_pix = metric_param ./ pix_resolution;
               end
            end

            cache = []; 
            if strcmp(metric_lcc, 'loc_cc_fftn_gpu') || strcmp(metric_lcc, 'loc_cc_fftn_gpu_single') || strcmp(metric_lcc, 'loc_cc_fftn') || strcmp(metric_lcc, 'loc_cc_fftn_single')
               cache = cell(Nch, 1);
                for ic = 1 : Nch
                    cache{ic} = create_loc_cc_fftn_cache(metric_lcc, metric_param_pix, volmov(:,:,:, ic, 1), volfix(:,:,:, ic, 1), internal_dtype, fixed_mask, 1e-4, loc_cc_approximate);
                end
            end
            unary_tot=Get_Data_Term_3D_LCC(volmov,volfix,griddim,O_trans_X,O_trans_Y,O_trans_Z,labels,level,previous,cache,metric_lcc,pix_resolution,fixed_mask,quant);
        end
    end
end
              
end


