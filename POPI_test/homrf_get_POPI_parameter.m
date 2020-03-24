function parameter=homrf_get_POPI_parameter(useTop,spc)




        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[ 7 7 7 ];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 9;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=1;
      
        parameter.top_co= 1;
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;

    

end
