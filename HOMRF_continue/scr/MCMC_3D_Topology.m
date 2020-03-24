
function [Ln,labelstot]=MCMC_3D_Topology(labels,dist,dist_co,neighbor,griddim,unary,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,grid_space,spc,quant,smooth_co_tot,Tcoe_n_tot,top_co_tot)

sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);
r=griddim(1);
c=griddim(2);
s=griddim(3);
x_h=fix(x_tot/2)+1;
y_h=fix(y_tot/2)+1;
z_h=fix(z_tot/2)+1;

itsmooth=size(neighbor,2);

% top_co=spc(1)*spc(2)*spc(3)/(grid_space(level,1))^3;
top_co=top_co_tot(level);
Tcoe_n=Tcoe_n_tot(level);

cur_grid_space_x=grid_space(level,1);
cur_grid_space_y=grid_space(level,2);
cur_grid_space_z=grid_space(level,3);

spc_tot=ones(itsmooth,3).*spc;

smooth_co=smooth_co_tot(level);
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


% Main function

T=0.0001;
   
totenergyend=zeros(5000,1);

max_iter=2000;

niter=fix(max_iter*k_down^(-level+7));
randnewlabel=unidrnd(labels.nlabels,[r*c*s*niter,1]);


[x_index,y_index,z_index]=Label_Coordinate(labels);
x_index=x_index*sample_space_x;
y_index=y_index*sample_space_y;
z_index=z_index*sample_space_z;

for iter=1:niter                  %iteration
      totenergy=0;  
      
for num=1:r*c*s


      
    originalnum=labelnum(num);
    x_mov_or=x_index(originalnum);
    y_mov_or=y_index(originalnum);
    z_mov_or=z_index(originalnum);
    x_mov_or_add=x_mov_or+previous(num,1)+x_h;
    y_mov_or_add=y_mov_or+previous(num,2)+y_h;
    z_mov_or_add=z_mov_or+previous(num,3)+z_h;

    ortotlabel=labelindex(x_mov_or_add,y_mov_or_add,z_mov_or_add);
    

    % Get the node label random
    newnum=randnewlabel((iter-1)*r*c*s+num);
    x_mov_newnum=x_index(newnum);
    y_mov_newnum=y_index(newnum);
    z_mov_newnum=z_index(newnum);
    x_mov_new_add=x_mov_newnum+previous(num,1)+x_h;
    y_mov_new_add=y_mov_newnum+previous(num,2)+y_h;
    z_mov_new_add=z_mov_newnum+previous(num,3)+z_h;

    
    newtotlabel=labelindex(x_mov_new_add,y_mov_new_add,z_mov_new_add);

          
    label_neighbor=neighbor(num,:);

    label_neighbor_or=neighbor(num,:);
    
    if isempty(dist)
        
        x1_or=x_mov_or+previous(num,1);
        y1_or=y_mov_or+previous(num,2);
        z1_or=z_mov_or+previous(num,3);
        smooth_or=repmat([x1_or,y1_or,z1_or],itsmooth,1).*spc_tot;
        smooth_neighbor=zeros(itsmooth,3);
        
        x2_new=x_mov_newnum+previous(num,1);
        y2_new=y_mov_newnum+previous(num,2);
        z2_new=z_mov_newnum+previous(num,3);
        smooth_new=repmat([x2_new,y2_new,z2_new],itsmooth,1).*spc_tot;
        
        smooth_neighbor(:,1)=x_index_tot(labelstot(label_neighbor));
        smooth_neighbor(:,2)=y_index_tot(labelstot(label_neighbor));
        smooth_neighbor(:,3)=-z_index_tot(labelstot(label_neighbor));
        smooth_neighbor=smooth_neighbor.*spc_tot;
        orb_tot=sqrt(sum((smooth_or-smooth_neighbor).^2,2));
        orb_tot(orb_tot>=80)=80;
        newb_tot=sqrt(sum((smooth_new-smooth_neighbor).^2,2));
        newb_tot(newb_tot>=80)=80;
        orbenergy=sum(orb_tot)/(dist_co)*smooth_co;
        newbenergy=sum(newb_tot)/(dist_co)*smooth_co;

    else
        orbenergy=sum(dist(ortotlabel,labelstot(label_neighbor)))*smooth_co;
        newbenergy=sum(dist(newtotlabel,labelstot(label_neighbor)))*smooth_co;
    end
    
    
         if useTop==1&&(iter>=niter/2)

             
             tporpx=x_index_tot_spc(ortotlabel);
             tporpy=y_index_tot_spc(ortotlabel);
             tporpz=z_index_tot_spc(ortotlabel);
    
             tpnewpx=x_index_tot_spc(newtotlabel);
             tpnewpy=y_index_tot_spc(newtotlabel);
             tpnewpz=z_index_tot_spc(newtotlabel);


             or_p_up_x=x_index_tot_spc(labelstot(label_neighbor_or(3)));
             or_p_up_y=y_index_tot_spc(labelstot(label_neighbor_or(3)));
             or_p_up_z=z_index_tot_spc(labelstot(label_neighbor_or(3)));

             or_p_up_z=(or_p_up_z+cur_grid_space_z*spc(3));
             

            or_p_down_x=x_index_tot_spc(labelstot(label_neighbor_or(16)));
            or_p_down_y=y_index_tot_spc(labelstot(label_neighbor_or(16)));
            or_p_down_z=z_index_tot_spc(labelstot(label_neighbor_or(16)));

            or_p_down_z=(or_p_down_z-cur_grid_space_z*spc(3));
        

            or_p_left_x=x_index_tot_spc(labelstot(label_neighbor_or(7)));
            or_p_left_y=y_index_tot_spc(labelstot(label_neighbor_or(7)));
            or_p_left_z=z_index_tot_spc(labelstot(label_neighbor_or(7)));

            or_p_left_y=(or_p_left_y-cur_grid_space_y*spc(2));
        

            or_p_right_x=x_index_tot_spc(labelstot(label_neighbor_or(12)));
            or_p_right_y=y_index_tot_spc(labelstot(label_neighbor_or(12)));
            or_p_right_z=z_index_tot_spc(labelstot(label_neighbor_or(12)));
         
%          [or_p_right_x,or_p_right_y,or_p_right_z]=Get_Deformation(trans_label(label_neighbor_or(16)),labels,space_coefficient);
            or_p_right_y=or_p_right_y+cur_grid_space_y*spc(2);
        
%             or_p_anter_x=x_index_tot(label_trans(label_neighbor_or(13)));
%             or_p_anter_y=y_index_tot(label_trans(label_neighbor_or(13)));
%             or_p_anter_z=z_index_tot(label_trans(label_neighbor_or(13)));
            or_p_anter_x=x_index_tot_spc(labelstot(label_neighbor_or(9)));
            or_p_anter_y=y_index_tot_spc(labelstot(label_neighbor_or(9)));
            or_p_anter_z=z_index_tot_spc(labelstot(label_neighbor_or(9)));
%          [or_p_anter_x,or_p_anter_y,or_p_anter_z]=Get_Deformation(trans_label(label_neighbor_or(13)),labels,space_coefficient);
            or_p_anter_x=or_p_anter_x-cur_grid_space_x*spc(1);
        
         
%             or_p_poster_x=x_index_tot(label_trans(label_neighbor_or(14)));
%             or_p_poster_y=y_index_tot(label_trans(label_neighbor_or(14)));
%             or_p_poster_z=z_index_tot(label_trans(label_neighbor_or(14)));
            or_p_poster_x=x_index_tot_spc(labelstot(label_neighbor_or(10)));
            or_p_poster_y=y_index_tot_spc(labelstot(label_neighbor_or(10)));
            or_p_poster_z=z_index_tot_spc(labelstot(label_neighbor_or(10)));
        
%          [or_p_poster_x,or_p_poster_y,or_p_poster_z]=Get_Deformation(trans_label(label_neighbor_or(14)),labels,space_coefficient);
            or_p_poster_x=or_p_poster_x+cur_grid_space_x*spc(1);
         
         
         %Jaco1 point 1 4 5
        %Jaco_1=(p_right_x-orpx)*(orpy-p_anter_y)*(p_up_z-orpz)+(orpx-p_anter_x)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(orpz-p_anter_z)*(p_up_x-orpx)-(p_right_z-orpz)*(orpy-p_anter_y)*(p_up_x-orpx)-(p_right_y-orpy)*(orpx-p_anter_x)*(p_up_z-orpz)-(orpz-p_anter_z)*(p_up_y-orpy)*(p_right_x-orpx);
        Jaco_1=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
        Jaco_1=Jaco_1*top_co;
        if Jaco_1>=0
            orthird_1=0;
        else

            orthird_1=log(-Jaco_1+1)*Tcoe_n;
        end
        
        
        %Jaco2 point 1 3 5
        %Jaco_2=(orpx-p_left_x)*(orpy-p_anter_y)*(p_up_z-orpz)+(orpx-p_anter_x)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(orpz-p_anter_z)*(p_up_x-orpx)-(orpz-p_left_z)*(orpy-p_anter_y)*(p_up_x-orpx)-(orpy-p_left_y)*(orpx-p_anter_x)*(p_up_z-orpz)-(orpz-p_anter_z)*(p_up_y-orpy)*(orpx-p_left_x);
        Jaco_2=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
        Jaco_2=Jaco_2*top_co;
        if Jaco_2>=0
            orthird_2=0;
        else

            orthird_2=log(-Jaco_2+1)*Tcoe_n;
        end
        
        
        %Jaco3 point 2 4 5
        %Jaco_3=(p_right_x-orpx)*(orpy-p_anter_y)*(orpz-p_down_z)+(orpx-p_anter_x)*(orpy-p_down_y)*(p_right_z-orpz)+(p_right_y-orpy)*(orpz-p_anter_z)*(orpx-p_down_x)-(p_right_z-orpz)*(orpy-p_anter_y)*(orpx-p_down_x)-(p_right_y-orpy)*(orpx-p_anter_x)*(orpz-p_down_z)-(orpz-p_anter_z)*(orpy-p_down_y)*(p_right_x-orpx);
        Jaco_3=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
        Jaco_3=Jaco_3*top_co;
        
        if Jaco_3>=0
            orthird_3=0;
        else

            orthird_3=log(-Jaco_3+1)*Tcoe_n;
        end
        
        
        %Jaco4 point 2 3 5
        %Jaco_4=(orpx-p_left_x)*(orpy-p_anter_y)*(orpz-p_down_z)+(orpx-p_anter_x)*(orpy-p_down_y)*(orpz-p_left_z)+(orpy-p_left_y)*(orpz-p_anter_z)*(orpx-p_down_x)-(orpz-p_left_z)*(orpy-p_anter_y)*(orpx-p_down_x)-(orpy-p_left_y)*(orpx-p_anter_x)*(orpz-p_down_z)-(orpz-p_anter_z)*(orpy-p_down_y)*(orpx-p_left_x);
        Jaco_4=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
       Jaco_4=Jaco_4*top_co; 
        
        if Jaco_4>=0
            orthird_4=0;
        else

            orthird_4=log(-Jaco_4+1)*Tcoe_n;
        end
        
        
        %Jaco5 point 1 4 6 fff
        %Jaco_5=(p_right_x-orpx)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(p_up_x-orpx)-(p_right_z-orpz)*(p_poster_y-orpy)*(p_up_x-orpx)-(p_right_y-orpy)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(p_right_x-orpx);
        Jaco_5=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
        Jaco_5=Jaco_5*top_co;
        if Jaco_5>=0
            orthird_5=0;
        else
            orthird_5=log(-Jaco_5+1)*Tcoe_n;
        end
        
        
        
        %Jaco6 point 2 3 6
        %Jaco_6=(orpx-p_left_x)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(orpx-p_down_x)-(orpz-p_left_z)*(p_poster_y-orpy)*(orpx-p_down_x)-(orpy-p_left_y)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(orpx-p_left_x);
        Jaco_6=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
        Jaco_6=Jaco_6*top_co;
        if Jaco_6>=0
            orthird_6=0;
        else

            orthird_6=log(-Jaco_6+1)*Tcoe_n;
        end
        
        
        
        %Jaco7 point 2 4 6
        %Jaco_7=(p_right_x-orpx)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(orpx-p_down_x)-(p_right_z-orpz)*(p_poster_y-orpy)*(orpx-p_down_x)-(p_right_y-orpy)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(-p_right_x-orpx);

        Jaco_7=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
        Jaco_7=Jaco_7*top_co;
        
        
        if Jaco_7>=0
            orthird_7=0;
        else
            orthird_7=log(-Jaco_7+1)*Tcoe_n;
        end
        
        
        %Jaco8 point 1 3 6
        %Jaco_8=(orpx-p_left_x)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(p_up_x-orpx)-(orpz-p_left_z)*(p_poster_y-orpy)*(p_up_x-orpx)-(orpy-p_left_y)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(orpx-p_left_x);
        Jaco_8=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
        Jaco_8=Jaco_8*top_co;
        
        
        if Jaco_8>=0
            orthird_8=0;
        else
  
            orthird_8=log(-Jaco_8+1)*Tcoe_n;
        end
        
        orthirdtot=orthird_1+orthird_2+orthird_3+orthird_4+orthird_5+orthird_6+orthird_7+orthird_8;
        
%         new_p_up_x=x_index_tot(label_trans(label_neighbor_new(5)));
%         new_p_up_y=y_index_tot(label_trans(label_neighbor_new(5)));
%         new_p_up_z=z_index_tot(label_trans(label_neighbor_new(5)));
%         new_p_up_x=x_index_tot(labelstot(label_neighbor_new(3)));
%         new_p_up_y=y_index_tot(labelstot(label_neighbor_new(3)));
%         new_p_up_z=z_index_tot(labelstot(label_neighbor_new(3)));
        new_p_up_x=or_p_up_x;
        new_p_up_y=or_p_up_y;
        new_p_up_z=or_p_up_z;
        

%          new_p_up_z=new_p_up_z+cur_grid_space;
        
%          new_p_down_x=x_index_tot(label_trans(label_neighbor_new(22)));
%          new_p_down_y=y_index_tot(label_trans(label_neighbor_new(22)));
%          new_p_down_z=z_index_tot(label_trans(label_neighbor_new(22)));
%          new_p_down_x=x_index_tot(labelstot(label_neighbor_new(16)));
%          new_p_down_y=y_index_tot(labelstot(label_neighbor_new(16)));
%          new_p_down_z=z_index_tot(labelstot(label_neighbor_new(16)));
% %          [new_p_down_x,new_p_down_y,new_p_down_z]=Get_Deformation(trans_label(label_neighbor_new(22)),labels,space_coefficient);
%          new_p_down_z=new_p_down_z-cur_grid_space;
         new_p_down_x=or_p_down_x;
         new_p_down_y=or_p_down_y;
         new_p_down_z=or_p_down_z;
         
%          new_p_left_x=x_index_tot(label_trans(label_neighbor_new(11)));
%          new_p_left_y=y_index_tot(label_trans(label_neighbor_new(11)));
%          new_p_left_z=z_index_tot(label_trans(label_neighbor_new(11)));
%          new_p_left_x=x_index_tot(labelstot(label_neighbor_new(7)));
%          new_p_left_y=y_index_tot(labelstot(label_neighbor_new(7)));
%          new_p_left_z=z_index_tot(labelstot(label_neighbor_new(7)));
% %          [new_p_left_x,new_p_left_y,new_p_left_z]=Get_Deformation(trans_label(label_neighbor_new(11)),labels,space_coefficient);
%          new_p_left_y=new_p_left_y-cur_grid_space;
         new_p_left_x=or_p_left_x;
         new_p_left_y=or_p_left_y;
         new_p_left_z=or_p_left_z;
         
%          new_p_right_x=x_index_tot(label_trans(label_neighbor_new(16)));
%          new_p_right_y=y_index_tot(label_trans(label_neighbor_new(16)));
%          new_p_right_z=z_index_tot(label_trans(label_neighbor_new(16)));
%          new_p_right_x=x_index_tot(labelstot(label_neighbor_new(12)));
%          new_p_right_y=y_index_tot(labelstot(label_neighbor_new(12)));
%          new_p_right_z=z_index_tot(labelstot(label_neighbor_new(12)));
% %          [new_p_right_x,new_p_right_y,new_p_right_z]=Get_Deformation(trans_label(label_neighbor_new(16)),labels,space_coefficient);
%          new_p_right_y=new_p_right_y+cur_grid_space;
         new_p_right_x=or_p_right_x;
         new_p_right_y=or_p_right_y;
         new_p_right_z=or_p_right_z;
%          new_p_anter_x=x_index_tot(label_trans(label_neighbor_new(13)));
%          new_p_anter_y=y_index_tot(label_trans(label_neighbor_new(13)));
%          new_p_anter_z=z_index_tot(label_trans(label_neighbor_new(13)));
%          new_p_anter_x=x_index_tot(labelstot(label_neighbor_new(9)));
%          new_p_anter_y=y_index_tot(labelstot(label_neighbor_new(9)));
%          new_p_anter_z=z_index_tot(labelstot(label_neighbor_new(9)));
% %          [new_p_anter_x,new_p_anter_y,new_p_anter_z]=Get_Deformation(trans_label(label_neighbor_new(13)),labels,space_coefficient);
%          new_p_anter_x=new_p_anter_x-cur_grid_space;
         new_p_anter_x=or_p_anter_x;           
         new_p_anter_y=or_p_anter_y;
         new_p_anter_z=or_p_anter_z;
         
%          new_p_poster_x=x_index_tot(label_trans(label_neighbor_new(14)));
%          new_p_poster_y=y_index_tot(label_trans(label_neighbor_new(14)));
%          new_p_poster_z=z_index_tot(label_trans(label_neighbor_new(14)));
%          new_p_poster_x=x_index_tot(labelstot(label_neighbor_new(10)));
%          new_p_poster_y=y_index_tot(labelstot(label_neighbor_new(10)));
%          new_p_poster_z=z_index_tot(labelstot(label_neighbor_new(10)));
% %          [new_p_poster_x,new_p_poster_y,new_p_poster_z]=Get_Deformation(trans_label(label_neighbor_new(14)),labels,space_coefficient);
%          new_p_poster_x=new_p_poster_x+cur_grid_space;
         new_p_poster_x=or_p_poster_x;
         new_p_poster_y=or_p_poster_y;
         new_p_poster_z=or_p_poster_z;
         %Jaco1 point 1 4 5
        %Jaco_1=(p_right_x-newpx)*(newpy-p_anter_y)*(p_up_z-newpz)+(newpx-p_anter_x)*(p_up_y-newpy)*(p_right_z-newpz)+(p_right_y-newpy)*(newpz-p_anter_z)*(p_up_x-newpx)-(p_right_z-newpz)*(newpy-p_anter_y)*(p_up_x-newpx)-(p_right_y-newpy)*(newpx-p_anter_x)*(p_up_z-newpz)-(newpz-p_anter_z)*(p_up_y-newpy)*(p_right_x-newpx);
        Jaco_1=(tpnewpx-new_p_anter_x)*(new_p_right_y-tpnewpy)*(new_p_up_z-tpnewpz)+(tpnewpy-new_p_anter_y)*(new_p_right_z-tpnewpz)*(new_p_up_x-tpnewpx)+(new_p_right_x-tpnewpx)*(new_p_up_y-tpnewpy)*(tpnewpz-new_p_anter_z)-(tpnewpz-new_p_anter_z)*(new_p_right_y-tpnewpy)*(new_p_up_x-tpnewpx)-(tpnewpy-new_p_anter_y)*(new_p_right_x-tpnewpx)*(new_p_up_z-tpnewpz)-(new_p_right_z-tpnewpz)*(new_p_up_y-tpnewpy)*(tporpx-new_p_anter_x);
        Jaco_1=Jaco_1*top_co;
        if Jaco_1>=0
            newthird_1=0;
        else
     
            newthird_1=log(-Jaco_1+1)*Tcoe_n;
        end
        
        
        %Jaco2 point 1 3 5
        %Jaco_2=(newpx-p_left_x)*(newpy-p_anter_y)*(p_up_z-newpz)+(newpx-p_anter_x)*(p_up_y-newpy)*(newpz-p_left_z)+(newpy-p_left_y)*(newpz-p_anter_z)*(p_up_x-newpx)-(newpz-p_left_z)*(newpy-p_anter_y)*(p_up_x-newpx)-(newpy-p_left_y)*(newpx-p_anter_x)*(p_up_z-newpz)-(newpz-p_anter_z)*(p_up_y-newpy)*(newpx-p_left_x);
        Jaco_2=(tpnewpx-new_p_anter_x)*(tpnewpy-new_p_left_y)*(new_p_up_z-tpnewpz)+(tpnewpy-new_p_anter_y)*(tpnewpz-new_p_left_z)*(new_p_up_x-tpnewpx)+(tpnewpx-new_p_left_x)*(new_p_up_y-tpnewpy)*(tpnewpz-new_p_anter_z)-(tpnewpz-new_p_anter_z)*(tpnewpy-new_p_left_y)*(new_p_up_x-tpnewpx)-(tpnewpy-new_p_anter_y)*(tpnewpx-new_p_left_x)*(new_p_up_z-tpnewpz)-(tpnewpz-new_p_left_z)*(new_p_up_y-tpnewpy)*(tpnewpx-new_p_anter_x);
        Jaco_2=Jaco_2*top_co;
        if Jaco_2>=0
            newthird_2=0;
        else
         
            newthird_2=log(-Jaco_2+1)*Tcoe_n;
        end
        
        
        %Jaco3 point 2 4 5
        %Jaco_3=(p_right_x-newpx)*(newpy-p_anter_y)*(newpz-p_down_z)+(newpx-p_anter_x)*(newpy-p_down_y)*(p_right_z-newpz)+(p_right_y-newpy)*(newpz-p_anter_z)*(newpx-p_down_x)-(p_right_z-newpz)*(newpy-p_anter_y)*(newpx-p_down_x)-(p_right_y-newpy)*(newpx-p_anter_x)*(newpz-p_down_z)-(newpz-p_anter_z)*(newpy-p_down_y)*(p_right_x-newpx);
        Jaco_3=(tpnewpx-new_p_anter_x)*(new_p_right_y-tpnewpy)*(tpnewpz-new_p_down_z)+(tpnewpy-new_p_anter_y)*(new_p_right_z-tpnewpz)*(tpnewpx-new_p_down_x)+(new_p_right_x-tpnewpx)*(tpnewpy-new_p_down_y)*(tpnewpz-new_p_anter_z)-(tpnewpz-new_p_anter_z)*(new_p_right_y-tpnewpy)*(tpnewpx-new_p_down_x)-(tpnewpy-new_p_anter_y)*(new_p_right_x-tpnewpx)*(tpnewpz-new_p_down_z)-(new_p_right_z-tpnewpz)*(tpnewpy-new_p_down_y)*(tpnewpx-new_p_anter_x);
        Jaco_3=Jaco_3*top_co;
        
        if Jaco_3>=0
            newthird_3=0;
        else
           
            newthird_3=log(-Jaco_3+1)*Tcoe_n;
        end
        
        
        %Jaco4 point 2 3 5
        %Jaco_4=(newpx-p_left_x)*(newpy-p_anter_y)*(newpz-p_down_z)+(newpx-p_anter_x)*(newpy-p_down_y)*(newpz-p_left_z)+(newpy-p_left_y)*(newpz-p_anter_z)*(newpx-p_down_x)-(newpz-p_left_z)*(newpy-p_anter_y)*(newpx-p_down_x)-(newpy-p_left_y)*(newpx-p_anter_x)*(newpz-p_down_z)-(newpz-p_anter_z)*(newpy-p_down_y)*(newpx-p_left_x);
        Jaco_4=(tpnewpx-new_p_anter_x)*(tpnewpy-new_p_left_y)*(tpnewpz-new_p_down_z)+(tpnewpy-new_p_anter_y)*(tpnewpz-new_p_left_z)*(tpnewpx-new_p_down_x)+(tpnewpx-new_p_left_x)*(tpnewpy-new_p_down_y)*(tpnewpz-new_p_anter_z)-(tpnewpz-new_p_anter_z)*(tpnewpy-new_p_left_y)*(tpnewpx-new_p_down_x)-(tpnewpy-new_p_anter_y)*(tpnewpx-new_p_left_x)*(tpnewpz-new_p_down_z)-(tpnewpz-new_p_left_z)*(tpnewpy-new_p_down_y)*(tpnewpx-new_p_anter_x);
        Jaco_4=Jaco_4*top_co;
        
        if Jaco_4>=0
            newthird_4=0;
        else
         
            newthird_4=log(-Jaco_4+1)*Tcoe_n;
        end
        
        
        %Jaco5 point 1 4 6
        %Jaco_5=(p_right_x-newpx)*(p_poster_y-newpy)*(p_up_z-newpz)+(p_poster_x-newpx)*(p_up_y-newpy)*(p_right_z-newpz)+(p_right_y-newpy)*(p_poster_z-newpz)*(p_up_x-newpx)-(p_right_z-newpz)*(p_poster_y-newpy)*(p_up_x-newpx)-(p_right_y-newpy)*(p_poster_x-newpx)*(p_up_z-newpz)-(p_poster_z-newpz)*(p_up_y-newpy)*(p_right_x-newpx);
        Jaco_5=(new_p_poster_x-tpnewpx)*(new_p_right_y-tpnewpy)*(new_p_up_z-tpnewpz)+(new_p_poster_y-tpnewpy)*(new_p_right_z-tpnewpz)*(new_p_up_x-tpnewpx)+(new_p_right_x-tpnewpx)*(new_p_up_y-tpnewpy)*(new_p_poster_z-tpnewpz)-(new_p_poster_z-tpnewpz)*(new_p_right_y-tpnewpy)*(new_p_up_x-tpnewpx)-(new_p_poster_y-tpnewpy)*(new_p_right_x-tpnewpx)*(new_p_up_z-tpnewpz)-(new_p_right_z-tpnewpz)*(new_p_up_y-tpnewpy)*(new_p_poster_x-tpnewpx);
        Jaco_5=Jaco_5*top_co;
        if Jaco_5>=0
            newthird_5=0;
        else
        
            newthird_5=log(-Jaco_5+1)*Tcoe_n;
        end
        
        
        
        %Jaco6 point 2 3 6
        %Jaco_6=(newpx-p_left_x)*(p_poster_y-newpy)*(newpz-p_down_z)+(p_poster_x-newpx)*(newpy-p_down_y)*(newpz-p_left_z)+(newpy-p_left_y)*(p_poster_z-newpz)*(newpx-p_down_x)-(newpz-p_left_z)*(p_poster_y-newpy)*(newpx-p_down_x)-(newpy-p_left_y)*(p_poster_x-newpx)*(newpz-p_down_z)-(p_poster_z-newpz)*(newpy-p_down_y)*(newpx-p_left_x);
        Jaco_6=(new_p_poster_x-tpnewpx)*(tpnewpy-new_p_left_y)*(tpnewpz-new_p_down_z)+(new_p_poster_y-tpnewpy)*(tpnewpz-new_p_left_z)*(tpnewpx-new_p_down_x)+(tpnewpx-new_p_left_x)*(tpnewpy-new_p_down_y)*(new_p_poster_z-tpnewpz)-(new_p_poster_z-tpnewpz)*(tpnewpy-new_p_left_y)*(tpnewpx-new_p_down_x)-(new_p_poster_y-tpnewpy)*(tpnewpx-new_p_left_x)*(tpnewpz-new_p_down_z)-(tpnewpz-new_p_left_z)*(tpnewpy-new_p_down_y)*(new_p_poster_x-tpnewpx);
        Jaco_6=Jaco_6*top_co;
        if Jaco_6>=0
            newthird_6=0;
        else
      
            newthird_6=log(-Jaco_6+1)*Tcoe_n;
        end
        
        
        
        %Jaco7 point 2 4 6
        %Jaco_7=(p_right_x-newpx)*(p_poster_y-newpy)*(newpz-p_down_z)+(p_poster_x-newpx)*(newpy-p_down_y)*(p_right_z-newpz)+(p_right_y-newpy)*(p_poster_z-newpz)*(newpx-p_down_x)-(p_right_z-newpz)*(p_poster_y-newpy)*(newpx-p_down_x)-(p_right_y-newpy)*(p_poster_x-newpx)*(newpz-p_down_z)-(p_poster_z-newpz)*(newpy-p_down_y)*(p_right_x-newpx);
        Jaco_7=(new_p_poster_x-tpnewpx)*(new_p_right_y-tpnewpy)*(tpnewpz-new_p_down_z)+(new_p_poster_y-tpnewpy)*(new_p_right_z-tpnewpz)*(tpnewpx-new_p_down_x)+(new_p_right_x-tpnewpx)*(tpnewpy-new_p_down_y)*(new_p_poster_z-tpnewpz)-(new_p_poster_z-tpnewpz)*(new_p_right_y-tpnewpy)*(tpnewpx-new_p_down_x)-(new_p_poster_y-tpnewpy)*(new_p_right_x-tpnewpx)*(tpnewpz-new_p_down_z)-(new_p_right_z-tpnewpz)*(tpnewpy-new_p_down_y)*(new_p_poster_x-tpnewpx);
        Jaco_7=Jaco_7*top_co;
        if Jaco_7>=0
            newthird_7=0;
        else
        
            newthird_7=log(-Jaco_7+1)*Tcoe_n;
        end
        
        
        %Jaco8 point 1 3 6
        %Jaco_8=(newpx-p_left_x)*(p_poster_y-newpy)*(p_up_z-newpz)+(p_poster_x-newpx)*(p_up_y-newpy)*(newpz-p_left_z)+(newpy-p_left_y)*(p_poster_z-newpz)*(p_up_x-newpx)-(newpz-p_left_z)*(p_poster_y-newpy)*(p_up_x-newpx)-(newpy-p_left_y)*(p_poster_x-newpx)*(p_up_z-newpz)-(p_poster_z-newpz)*(p_up_y-newpy)*(newpx-p_left_x);
        Jaco_8=(new_p_poster_x-tpnewpx)*(tpnewpy-new_p_left_y)*(new_p_up_z-tpnewpz)+(new_p_poster_y-tpnewpy)*(tpnewpz-new_p_left_z)*(new_p_up_x-tpnewpx)+(tpnewpx-new_p_left_x)*(new_p_up_y-tpnewpy)*(new_p_poster_z-tpnewpz)-(new_p_poster_z-tpnewpz)*(tpnewpy-new_p_left_y)*(new_p_up_x-tpnewpx)-(new_p_poster_y-tpnewpy)*(tpnewpx-new_p_left_x)*(new_p_up_z-tpnewpz)-(tpnewpz-new_p_left_z)*(new_p_up_y-tpnewpy)*(new_p_poster_x-tpnewpx);
        Jaco_8=Jaco_8*top_co;
        
          if   Jaco_8>=0
               newthird_8=0;
          else
        
               newthird_8=log(-Jaco_8+1)*Tcoe_n;
          end
        
          newthirdtot=newthird_1+newthird_2+newthird_3+newthird_4+newthird_5+newthird_6+newthird_7+newthird_8;
         else
             orthirdtot=0;
             newthirdtot=0;
         end
         

       unaryweight=unary(num,:);

       orunaryweight=unaryweight(originalnum);

       newunaryweight=unaryweight(newnum);

               
         %Third
         orenergy=orunaryweight+orbenergy+orthirdtot;
  

         newenergy=newunaryweight+newbenergy+newthirdtot;


      diff=orenergy-newenergy;
            
        ratio=exp(diff/T);
         
          U=rand;
         if  (ratio > 1)
            labelnum(num)=newnum;
            labelstot(num)=newtotlabel;
            totenergy=totenergy+newenergy;

         else
             if ratio==1
                 labelnum(num)=labelnum(num);
                 labelstot(num)=ortotlabel;
                 totenergy=totenergy+orenergy;

             else
                 if ratio>U
                     labelnum(num)=newnum;
                     labelstot(num)=newtotlabel;
                     totenergy=totenergy+orenergy;

                 else
                     labelnum(num)=labelnum(num);
                     labelstot(num)=ortotlabel;
                     totenergy=totenergy+orenergy;

                 end
             end
         end
    
end
totenergyend(iter)=totenergy;

end
Ln=labelnum(1:r*c*s,1);

