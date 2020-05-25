import torch
import torch.nn.functional as F
import numpy as np

def conv3D(img_1, filter_1, groups,channel,Nc, Nr,Ns, Nfc, Nfr, Nfs):
    tar_size = int(groups*Nc*Nr*Ns)
    f_size = int(groups*Nfc*Nfr*Nfs)
    img_1 = np.array(img_1)
    filter_1 = np.array(filter_1)
    groups = int(groups)
    Nc = int(Nc)
    Nr = int(Nr)
    Ns = int(Ns)
    Nfc = int(Nfc)
    Nfr = int(Nfr)
    Nfs = int(Nfs)
    filter_m_1 = filter_1[0:f_size]
    filter_m_1 = filter_m_1.reshape((groups, Nfc, Nfr, Nfs), order='F')
    filter_m_2 = filter_1[f_size:f_size*2]
    filter_m_2 = filter_m_2.reshape((groups, Nfc, Nfr, Nfs), order='F')
    filter_m_3 = filter_1[f_size*2:f_size*3]
    filter_m_3 = filter_m_3.reshape((groups, Nfc, Nfr, Nfs), order='F')
    filter_m_4 = filter_1[f_size*3:f_size*4]
    filter_m_4 = filter_m_4.reshape((groups, Nfc, Nfr, Nfs), order='F')
    filter_m_5 = filter_1[f_size*4:f_size*5]
    filter_m_5 = filter_m_5.reshape((groups, Nfc, Nfr, Nfs), order='F')
    filter_m_6 = filter_1[f_size*5:f_size*6]
    filter_m_6 = filter_m_6.reshape((groups, Nfc, Nfr, Nfs), order='F')

    target_m_1 = img_1[0:tar_size]
    target_m_1 = target_m_1.reshape((groups, Nr, Nc, Ns), order='F')
    target_m_2 = img_1[tar_size:tar_size*2]
    target_m_2 = target_m_2.reshape((groups, Nr, Nc, Ns), order='F')
    target_m_3 = img_1[tar_size*2:tar_size*3]
    target_m_3 = target_m_3.reshape((groups, Nr, Nc, Ns), order='F')
    target_m_4 = img_1[tar_size*3:tar_size*4]
    target_m_4 = target_m_4.reshape((groups, Nr, Nc, Ns), order='F')
    target_m_5 = img_1[tar_size*4:tar_size*5]
    target_m_5 = target_m_5.reshape((groups, Nr, Nc, Ns), order='F')
    target_m_6 = img_1[tar_size*5:tar_size*6]
    target_m_6 = target_m_6.reshape((groups, Nr, Nc, Ns), order='F')

    tar_filter = np.full((groups, 1, Nfr, Nfc, Nfs), 1)


    m_filter_tensor = torch.from_numpy(tar_filter).view(groups, 1, Nr, Nc, Ns).cuda()

    target_tensor_1 = torch.from_numpy(target_m_1).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_1 = torch.from_numpy(filter_m_1).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_1 = F.conv3d(target_tensor_1, filter_tensor_1, groups=groups)
    m_output3d_1 = F.conv3d(target_tensor_1 ** 2, m_filter_tensor, groups=groups)

    target_tensor_2 = torch.from_numpy(target_m_2).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_2 = torch.from_numpy(filter_m_2).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_2 = F.conv3d(target_tensor_2, filter_tensor_2, groups=groups)
    m_output3d_2 = F.conv3d(target_tensor_2 ** 2, m_filter_tensor, groups=groups)

    target_tensor_3 = torch.from_numpy(target_m_3).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_3 = torch.from_numpy(filter_m_3).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_3 = F.conv3d(target_tensor_3, filter_tensor_3, groups=groups)
    m_output3d_3 = F.conv3d(target_tensor_3 ** 2, m_filter_tensor, groups=groups)

    target_tensor_4 = torch.from_numpy(target_m_4).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_4 = torch.from_numpy(filter_m_4).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_4 = F.conv3d(target_tensor_4, filter_tensor_4, groups=groups)
    m_output3d_4 = F.conv3d(target_tensor_4 ** 2, m_filter_tensor, groups=groups)

    target_tensor_5 = torch.from_numpy(target_m_5).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_5 = torch.from_numpy(filter_m_5).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_5 = F.conv3d(target_tensor_5, filter_tensor_5, groups=groups)
    m_output3d_5 = F.conv3d(target_tensor_5 ** 2, m_filter_tensor, groups=groups)

    target_tensor_6 = torch.from_numpy(target_m_6).view(1, groups, Nr, Nc, Ns).cuda()
    filter_tensor_6 = torch.from_numpy(filter_m_6).view(groups, 1, Nfr, Nfc, Nfs).cuda()
    output3d_6 = F.conv3d(target_tensor_6, filter_tensor_6, groups=groups)
    m_output3d_6 = F.conv3d(target_tensor_6 ** 2, m_filter_tensor, groups=groups)

    return output3d_1, output3d_2, output3d_3, output3d_4, output3d_5, output3d_6,\
           m_output3d_1, m_output3d_2, m_output3d_3, m_output3d_4, m_output3d_5, m_output3d_6
