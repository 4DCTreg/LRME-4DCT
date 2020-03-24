function [count_Jaco,count1]=calculate_Jaco(neighbor,griddim,level,x_tot,y_tot,z_tot,labelstot,grid_space,spc)

r=griddim(1);
c=griddim(2);
s=griddim(3);
count_Jaco=zeros(r*c*s,8);
top_co=spc(1)*spc(2)*spc(3)/(grid_space(level,1))^3;

cur_grid_space_x=grid_space(level,1);
cur_grid_space_y=grid_space(level,2);
cur_grid_space_z=grid_space(level,3);
[x_index_tot,y_index_tot,z_index_tot]=Label_Coordinate_tot(x_tot,y_tot,z_tot);
z_index_tot=-z_index_tot;
x_index_tot_spc=x_index_tot.*spc(1);
y_index_tot_spc=y_index_tot.*spc(2);
z_index_tot_spc=z_index_tot.*spc(3);

for index=1:r*c*s
    ortotlabel=labelstot(index);
    label_neighbor_or=neighbor(index,:);
    tporpx=x_index_tot_spc(ortotlabel);
    tporpy=y_index_tot_spc(ortotlabel);
    tporpz=z_index_tot_spc(ortotlabel);
    
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
         
            or_p_right_y=or_p_right_y+cur_grid_space_y*spc(2);
        
            or_p_anter_x=x_index_tot_spc(labelstot(label_neighbor_or(9)));
            or_p_anter_y=y_index_tot_spc(labelstot(label_neighbor_or(9)));
            or_p_anter_z=z_index_tot_spc(labelstot(label_neighbor_or(9)));

            or_p_anter_x=or_p_anter_x-cur_grid_space_x*spc(1);
        
         

            or_p_poster_x=x_index_tot_spc(labelstot(label_neighbor_or(10)));
            or_p_poster_y=y_index_tot_spc(labelstot(label_neighbor_or(10)));
            or_p_poster_z=z_index_tot_spc(labelstot(label_neighbor_or(10)));
        

            or_p_poster_x=or_p_poster_x+cur_grid_space_x*spc(1);
         
         
       
        Jaco_1=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
        Jaco_1=Jaco_1*top_co;

        
        
        
        Jaco_2=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
        Jaco_2=Jaco_2*top_co;

        
        
        Jaco_3=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
        Jaco_3=Jaco_3*top_co;
        
        
       
        Jaco_4=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
        Jaco_4=Jaco_4*top_co; 
             
        
        
        Jaco_5=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
        Jaco_5=Jaco_5*top_co;

      
       
        Jaco_6=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
        Jaco_6=Jaco_6*top_co;
                     

        Jaco_7=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
        Jaco_7=Jaco_7*top_co;
        
                
        %Jaco8 point 1 3 6
        %Jaco_8=(orpx-p_left_x)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(p_up_x-orpx)-(orpz-p_left_z)*(p_poster_y-orpy)*(p_up_x-orpx)-(orpy-p_left_y)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(orpx-p_left_x);
        Jaco_8=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
        Jaco_8=Jaco_8*top_co;
        
        count_Jaco(index,:)=[Jaco_1,Jaco_2,Jaco_3,Jaco_4,Jaco_5,Jaco_6,Jaco_7,Jaco_8];
       
            
end
count1=length(find(count_Jaco<0));
end