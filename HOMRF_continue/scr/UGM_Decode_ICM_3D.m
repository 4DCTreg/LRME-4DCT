function  [Ln,labelstot] = UGM_Decode_ICM_3D(labels,dist,dist_co,neighbor,griddim,nodePot,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,grid_space,spc,quant,smooth_co_tot,Tcoe_n_tot,top_co_tot)
num_neighbor=18;
[junk Ln] = min(nodePot,[],2);
nneighbor=size(neighbor,2);
sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);
% top_co=top_co_tot(level);
Tcoe_n=Tcoe_n_tot;
r=griddim(1);
c=griddim(2);
s=griddim(3);
x_h=fix(x_tot/2)+1;
y_h=fix(y_tot/2)+1;
z_h=fix(z_tot/2)+1;


% itsmooth=size(neighbor,2);

% top_co=spc(1)*spc(2)*spc(3)/(grid_space(level,1))^3;
top_co=top_co_tot(level);
% Tcoe_n=Tcoe_n_tot(level);

cur_grid_space_x=grid_space(level,1);
cur_grid_space_y=grid_space(level,2);
cur_grid_space_z=grid_space(level,3);
cur_pix_spc=[cur_grid_space_x,cur_grid_space_y,cur_grid_space_z];
% spc_tot=ones(itsmooth,3).*spc;

smooth_co=smooth_co_tot;
[x_index_tot,y_index_tot,z_index_tot]=Label_Coordinate_tot(x_tot,y_tot,z_tot);

z_index_tot=-z_index_tot;

x_index_tot_spc=x_index_tot.*spc(1);
y_index_tot_spc=y_index_tot.*spc(2);
z_index_tot_spc=z_index_tot.*spc(3);


labelindex=1:x_tot*y_tot*z_tot;
labelindex=reshape(labelindex,[x_tot,y_tot,z_tot]);



% This label index the data term
labelnum=zeros(r*c*s+1,1)+floor(labels.nlabels/2)+1;
labelstot_1=zeros(r*c*s+1,1)+floor(x_tot*y_tot*z_tot/2)+1;
labelstot_1(1:r*c*s)=labelstot;
labelstot=labelstot_1;
  
maxIter=50;

[x_index,y_index,z_index]=Label_Coordinate(labels);
x_index=x_index*sample_space_x;
y_index=y_index*sample_space_y;
z_index=z_index*sample_space_z;



[nNodes,maxStates] = size(nodePot);

nStates = maxStates;
y=labelnum;
done = 0;
iter=1;
while ~done&&(iter<maxIter)
    done = 1;
	y2 = y;
    for n = 1:nNodes
        
        x_mov_or_add=x_index+previous(n,1)+x_h;
        y_mov_or_add=y_index+previous(n,2)+y_h;
        z_mov_or_add=z_index+previous(n,3)+z_h;
        totlabel=zeros(1,nStates);
        for idx=1:nStates 
            totlabel(idx)=labelindex(x_mov_or_add(idx),y_mov_or_add(idx),z_mov_or_add(idx));
        end
    
        % Compute Node Potential
        pot = nodePot(n,1:nStates);

        % Find Neighbors
        edges=neighbor(n,:);
        if num_neighbor==6
            up_point=edges(1);
            down_point=edges(6);
            left_point=edges(2);
            right_point=edges(5);
            anter_point=edges(3);
            poster_point=edges(4);
        else
            if num_neighbor==18
                up_point=edges(3);
                down_point=edges(16);
                left_point=edges(7);
                right_point=edges(12);
                anter_point=edges(9);
                poster_point=edges(10);
            end
        end
                
        

            smooth_cur=zeros(nStates,3);
            x_index_nstate=x_index_tot_spc(totlabel);
            y_index_nstate=y_index_tot_spc(totlabel);
            z_index_nstate=z_index_tot_spc(totlabel);
            smooth_cur(:,1)=x_index_nstate;
            smooth_cur(:,2)=y_index_nstate;
            smooth_cur(:,3)=z_index_nstate;
            for e=1:nneighbor
                if isempty(dist)
                ecurr=labelstot(edges(e));
                ecurr_x=x_index_tot_spc(ecurr);
                ecurr_y=y_index_tot_spc(ecurr);
                ecurr_z=z_index_tot_spc(ecurr);
                smooth_e=repmat([ecurr_x,ecurr_y,ecurr_z],nStates,1);
                
%                 ep=sum(abs(smooth_cur-smooth_e),2);
                ep=sqrt(sum((smooth_cur-smooth_e).^2,2));
                ep(ep>=10)=10;
                ep=ep.*smooth_co;
                ep=ep';
                pot=pot+ep;
                else
                    ecurr=labelstot(edges(e));
                    ep=dist(ecurr,totlabel);
                    ep=ep.*smooth_co;
                    pot = pot +ep;
                end
             end


        % Multiply Edge Potentials

        etp=zeros(1,nStates);

        if useTop==1
            for tpidx=1:nStates
                ortotlabel=totlabel(tpidx);
                tporpx=x_index_tot_spc(ortotlabel);
                tporpy=y_index_tot_spc(ortotlabel);
                tporpz=z_index_tot_spc(ortotlabel);
    
                or_p_up_x=x_index_tot_spc(labelstot(up_point));
                or_p_up_y=y_index_tot_spc(labelstot(up_point));
                or_p_up_z=z_index_tot_spc(labelstot(up_point));
                or_p_up_z=(or_p_up_z+cur_grid_space_z*spc(3));
             

                or_p_down_x=x_index_tot_spc(labelstot(down_point));
                or_p_down_y=y_index_tot_spc(labelstot(down_point));
                or_p_down_z=z_index_tot_spc(labelstot(down_point));
                or_p_down_z=(or_p_down_z-cur_grid_space_z*spc(3));
        

                or_p_left_x=x_index_tot_spc(labelstot(left_point));
                or_p_left_y=y_index_tot_spc(labelstot(left_point));
                or_p_left_z=z_index_tot_spc(labelstot(left_point));
                or_p_left_y=(or_p_left_y-cur_grid_space_y*spc(2));
        
                or_p_right_x=x_index_tot_spc(labelstot(right_point));
                or_p_right_y=y_index_tot_spc(labelstot(right_point));
                or_p_right_z=z_index_tot_spc(labelstot(right_point));
                or_p_right_y=or_p_right_y+cur_grid_space_y*spc(2);
        

                or_p_anter_x=x_index_tot_spc(labelstot(anter_point));
                or_p_anter_y=y_index_tot_spc(labelstot(anter_point));
                or_p_anter_z=z_index_tot_spc(labelstot(anter_point));
                or_p_anter_x=or_p_anter_x-cur_grid_space_x*spc(1);
        
         

                or_p_poster_x=x_index_tot_spc(labelstot(poster_point));
                or_p_poster_y=y_index_tot_spc(labelstot(poster_point));
                or_p_poster_z=z_index_tot_spc(labelstot(poster_point)); 
                or_p_poster_x=or_p_poster_x+cur_grid_space_x*spc(1);
                
                [tpenergy] = homrf_jac(tporpx,tporpy,tporpz,or_p_up_x,or_p_up_y,or_p_up_z,...
    or_p_down_x,or_p_down_y,or_p_down_z,or_p_left_x,or_p_left_y,or_p_left_z,...
    or_p_right_x,or_p_right_y,or_p_right_z,or_p_anter_x,or_p_anter_y,or_p_anter_z,...
    or_p_poster_x,or_p_poster_y,or_p_poster_z,Tcoe_n,top_co);
                
                 
                etp(tpidx)=tpenergy;
        

            end
        end
        

        pot=pot+etp;
        % Assign to Maximum State
        [junk newY] = min(pot);
        if newY ~= y(n)
            y(n) = newY;
            labelstot(n)=totlabel(newY);
%             done = 0;
        end
    end
   iter=iter+1;
   changes=sum(y2~=y);
   if changes~=0
       done=0;
   else
       done=1;
   end
%    fprintf('logPot = %f, changes = %d\n',UGM_ConfigurationPotential_correct_TP(y,nodePot,edgePot,edgeEnds,TP_energy,C_tp),sum(y2~=y));
   
   fprintf('changes = %d, iter = %d\n',sum(y2~=y), iter);
%             fprintf('changes = %d\n',sum(y2~=y));
 end


Ln=y(1:nNodes);
end
