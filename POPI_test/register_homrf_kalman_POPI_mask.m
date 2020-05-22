function [Knots_n_tot,Knots_n]=register_homrf_kalman_POPI_mask(parameter,...
    volmov,volfix,mov_mask,fix_mask,init_size,pts_fix,pts_mov,Knots_n_pre,bszv)

nlevel=parameter.nlevel;
k_down=parameter.k_down;
x_tot=parameter.x_tot;
y_tot=parameter.y_tot;
z_tot=parameter.z_tot;
useTop=parameter.useTop;
grid_space=parameter.grid_space;
metric=parameter.metric;
labels_search=parameter.labels;
presmooth=parameter.presmooth;
quant=parameter.quant;
smooth_co_pre=parameter.smooth_co;
endl=parameter.end;
Tcoe_n=parameter.Tcoe_n;
top_co=parameter.top_co;
spc=parameter.spc;
spc_orig = spc;
spc = [1,1,1];
resize=parameter.resize;
dist_co=parameter.dist_co;
max_gs=grid_space(nlevel,:);
volsz = size(volmov);

if presmooth==1
    dist=Get_Smooth_tot(x_tot,y_tot,z_tot,spc);     
else
    dist=[];
end

grid_sz=zeros(nlevel,3);

for level=nlevel:-1:endl
    
    labels=DIR_Process_labels(labels_search,level);
        
    % Current Grid Size
    cur_grid_space=grid_space(level,:);
    
    smooth_co = prod(max_gs)/prod(cur_grid_space)*smooth_co_pre;
%     smooth_co = smooth_co_pre;
    

    % Make the registration control grid 
%     [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_homrf(cur_grid_space,size(volmov));
    [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
    % Control grid dimensions (in the image)
    griddim=[size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
    grid_sz(level,:)=griddim;
    
    if level==nlevel
        Knots_n = get_Knots_pre(Knots_n_pre,O_trans_X,O_trans_Y,O_trans_Z,griddim);
        previous = get_previous_displacement(Knots_n,griddim);
    else 
        Knots_n = refine_spline_grid_3D_nl(Knots_n, O_trans_X_or,...
            O_trans_Y_or,O_trans_Z_or, O_trans_X,O_trans_Y,O_trans_Z);
        previous=get_previous_displacement(Knots_n,griddim);
            
    end
    previous=round(previous);
    labelstot=get_labels(previous,x_tot,y_tot,z_tot);
    O_trans_X_or=O_trans_X;O_trans_Y_or=O_trans_Y;O_trans_Z_or=O_trans_Z;
        
    % Get the neighbor
     neighbor=Get_Neighbor_18(griddim(1),griddim(2),griddim(3));
%      neighbor=Get_Neighbor_6(griddim(1),griddim(2),griddim(3));
    % Get the unary term

    unary=Get_Data_Term_3D_mask(volmov,volfix,mov_mask,fix_mask,griddim,O_trans_X,O_trans_Y,...
        O_trans_Z,labels,level,previous,quant,metric,cur_grid_space);

    [Ln,labelstot_Jaco] = UGM_Decode_ICM_3D(labels,dist,dist_co,neighbor,...
        griddim,unary,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,...
        grid_space,spc,quant,smooth_co,Tcoe_n,top_co);

    Knots_min = Knots_Displacement_3D_Multi_forwards(Ln,labels.sx,labels.sy,...
        labels.sz,griddim,level,quant);        
    Knots_n = Knots_n+Knots_min;

    Knots_n_tot = refine_spline_grid_3d(Knots_n,[1,1,1], size(volmov), O_trans_X,O_trans_Y,O_trans_Z);


end

end