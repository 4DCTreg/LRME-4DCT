function [gy_filter,tar_sum] = conv3d_gpu(img,filter,m_filter,groups)



save inter.mat img filter m_filter groups -v7

path = cd;
opean_path = ['CD',' ',path];
activate_env = ['conda activate Conv3D_GPU'];
file = ['python conv3d_gpu_cuda.py'];
tic
commd = [opean_path '&&' activate_env '&&' file];

state = system(commd);
toc
if state ~= 0
    error(['cmd error.']); 
end
load('results.mat');
conv = squeeze(array);
gy_filter = conv(7:12,:,:,:,:);
tar_sum = conv(1:6,:,:,:,:);

end

