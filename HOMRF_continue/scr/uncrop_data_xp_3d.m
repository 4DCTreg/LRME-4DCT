function [vols] = uncrop_data_xp_3d(vols_in, crp, size_orig,weight)


        vols = ones(size_orig)*weight;

        vols(crp(1,1):crp(1,2), crp(2,1):crp(2,2), crp(3,1):crp(3,2)) = vols_in;
   
end