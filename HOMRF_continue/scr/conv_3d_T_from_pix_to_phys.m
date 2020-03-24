function Tphys = conv_3d_T_from_pix_to_phys(Tpix, spc)
    Tphys = Tpix;
    for j = 1 : size(Tpix, 5)
    for i = 1 : size(Tpix, 4)
        Tphys(:,:,:, i, j) = Tphys(:,:,:,i, j) * spc(i);
    end
    end
end