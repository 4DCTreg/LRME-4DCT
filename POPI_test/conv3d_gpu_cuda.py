import torch
import torch.nn.functional as F
import numpy as np
import scipy.io as scio
path = 'inter.mat'
data = scio.loadmat(path)

target_m = data['img']
filter_m = data['filter']
group_m = data['groups']
m_filter_m = data['m_filter']
m_target_m = target_m**2

t_shape = target_m.shape
f_shape = filter_m.shape
m_f_shape = m_filter_m.shape
filter_m_1 = filter_m[:, 0, :, :, :]
filter_m_2 = filter_m[:, 1, :, :, :]
filter_m_3 = filter_m[:, 2, :, :, :]
filter_m_4 = filter_m[:, 3, :, :, :]
filter_m_5 = filter_m[:, 4, :, :, :]
filter_m_6 = filter_m[:, 5, :, :, :]

target_m_1 = target_m[:, 0, :, :, :]
target_m_2 = target_m[:, 1, :, :, :]
target_m_3 = target_m[:, 2, :, :, :]
target_m_4 = target_m[:, 3, :, :, :]
target_m_5 = target_m[:, 4, :, :, :]
target_m_6 = target_m[:, 5, :, :, :]

m_target_m_1 = m_target_m[:, 0, :, :, :]
m_target_m_2 = m_target_m[:, 1, :, :, :]
m_target_m_3 = m_target_m[:, 2, :, :, :]
m_target_m_4 = m_target_m[:, 3, :, :, :]
m_target_m_5 = m_target_m[:, 4, :, :, :]
m_target_m_6 = m_target_m[:, 5, :, :, :]



ngroup = int(group_m)
m_filter_tensor = torch.from_numpy(m_filter_m).view(m_f_shape[0],1,m_f_shape[1],m_f_shape[2],m_f_shape[3]).cuda()

target_tensor_1 = torch.from_numpy(target_m_1).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()

filter_tensor_1 = torch.from_numpy(filter_m_1).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_1 = F.conv3d(target_tensor_1, filter_tensor_1, groups=ngroup)
m_output3d_1 = F.conv3d(target_tensor_1**2, m_filter_tensor, groups=ngroup)

target_tensor_2 = torch.from_numpy(target_m_2).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()
filter_tensor_2 = torch.from_numpy(filter_m_2).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_2 = F.conv3d(target_tensor_2, filter_tensor_2, groups=ngroup)
m_output3d_2 = F.conv3d(target_tensor_2**2, m_filter_tensor, groups=ngroup)

target_tensor_3 = torch.from_numpy(target_m_3).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()
filter_tensor_3 = torch.from_numpy(filter_m_3).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_3 = F.conv3d(target_tensor_3, filter_tensor_3, groups=ngroup)
m_output3d_3 = F.conv3d(target_tensor_3**2, m_filter_tensor, groups=ngroup)

target_tensor_4 = torch.from_numpy(target_m_4).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()
filter_tensor_4 = torch.from_numpy(filter_m_4).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_4 = F.conv3d(target_tensor_4, filter_tensor_4, groups=ngroup)
m_output3d_4 = F.conv3d(target_tensor_4**2, m_filter_tensor, groups=ngroup)

target_tensor_5 = torch.from_numpy(target_m_5).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()
filter_tensor_5 = torch.from_numpy(filter_m_5).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_5 = F.conv3d(target_tensor_5, filter_tensor_5, groups=ngroup)
m_output3d_5 = F.conv3d(target_tensor_5**2, m_filter_tensor, groups=ngroup)

target_tensor_6 = torch.from_numpy(target_m_6).view(1,t_shape[0],t_shape[2],t_shape[3],t_shape[4]).cuda()
filter_tensor_6 = torch.from_numpy(filter_m_6).view(f_shape[0],1,f_shape[2],f_shape[3],f_shape[4]).cuda()
output3d_6 = F.conv3d(target_tensor_6, filter_tensor_6, groups=ngroup)
m_output3d_6 = F.conv3d(target_tensor_6**2, m_filter_tensor, groups=ngroup)

result1 = output3d_1.cpu().numpy()
result2 = output3d_2.cpu().numpy()
result3 = output3d_3.cpu().numpy()
result4 = output3d_4.cpu().numpy()
result5 = output3d_5.cpu().numpy()
result6 = output3d_6.cpu().numpy()

m_result1 = m_output3d_1.cpu().numpy()
m_result2 = m_output3d_2.cpu().numpy()
m_result3 = m_output3d_3.cpu().numpy()
m_result4 = m_output3d_4.cpu().numpy()
m_result5 = m_output3d_5.cpu().numpy()
m_result6 = m_output3d_6.cpu().numpy()



result_all = np.stack((m_result1, m_result2, m_result3, m_result4, m_result5, m_result6, result1, result2, result3, result4, result5, result6), axis=0)
save_fname = 'results.mat'

scio.savemat(save_fname, {'array': result_all})




