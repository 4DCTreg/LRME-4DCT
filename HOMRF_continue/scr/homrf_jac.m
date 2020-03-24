function [tpenergy] = homrf_jac(tporpx,tporpy,tporpz,or_p_up_x,or_p_up_y,or_p_up_z,...
    or_p_down_x,or_p_down_y,or_p_down_z,or_p_left_x,or_p_left_y,or_p_left_z,...
    or_p_right_x,or_p_right_y,or_p_right_z,or_p_anter_x,or_p_anter_y,or_p_anter_z,...
    or_p_poster_x,or_p_poster_y,or_p_poster_z,Tcoe_n,top_co)
 %Jaco1 point 1 4 5
%   Jaco_1=(orpx-or_p_anter_x)*(or_p_right_y-orpy)*(or_p_up_z-orpz)+(orpy-or_p_anter_y)*(or_p_right_z-orpz)*(or_p_up_x-orpx)+(or_p_right_x-orpx)*(or_p_up_y-orpy)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(or_p_right_y-orpy)*(or_p_up_x-orpx)-(orpy-or_p_anter_y)*(or_p_right_x-orpx)*(or_p_up_z-orpz)-(or_p_right_z-orpz)*(or_p_up_y-orpy)*(orpx-or_p_anter_x);
Jaco_1=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
Jaco_1=Jaco_1*top_co;
if Jaco_1>=0
    orthird_1=0;
else
    orthird_1=log(-Jaco_1+1)*Tcoe_n;
end
        
%Jaco2 point 1 3 5
%                  Jaco_2=(orpx-or_p_anter_x)*(orpy-or_p_left_y)*(or_p_up_z-orpz)+(orpy-or_p_anter_y)*(orpz-or_p_left_z)*(or_p_up_x-orpx)+(orpx-or_p_left_x)*(or_p_up_y-orpy)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(orpy-or_p_left_y)*(or_p_up_x-orpx)-(orpy-or_p_anter_y)*(orpx-or_p_left_x)*(or_p_up_z-orpz)-(orpz-or_p_left_z)*(or_p_up_y-orpy)*(orpx-or_p_anter_x);
Jaco_2=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
Jaco_2=Jaco_2*top_co;
if Jaco_2>=0
    orthird_2=0;
else
    orthird_2=log(-Jaco_2+1)*Tcoe_n;
end
        
        
%Jaco3 point 2 4 5
%                 Jaco_3=(orpx-or_p_anter_x)*(or_p_right_y-orpy)*(orpz-or_p_down_z)+(orpy-or_p_anter_y)*(or_p_right_z-orpz)*(orpx-or_p_down_x)+(or_p_right_x-orpx)*(orpy-or_p_down_y)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(or_p_right_y-orpy)*(orpx-or_p_down_x)-(orpy-or_p_anter_y)*(or_p_right_x-orpx)*(orpz-or_p_down_z)-(or_p_right_z-orpz)*(orpy-or_p_down_y)*(orpx-or_p_anter_x); 
Jaco_3=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
Jaco_3=Jaco_3*top_co;
        
if Jaco_3>=0
    orthird_3=0;
else
    orthird_3=log(-Jaco_3+1)*Tcoe_n;
end
        
        
                %Jaco4 point 2 3 5
%                                Jaco_4=(orpx-or_p_anter_x)*(orpy-or_p_left_y)*(orpz-or_p_down_z)+(orpy-or_p_anter_y)*(orpz-or_p_left_z)*(orpx-or_p_down_x)+(orpx-or_p_left_x)*(orpy-or_p_down_y)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(orpy-or_p_left_y)*(orpx-or_p_down_x)-(orpy-or_p_anter_y)*(orpx-or_p_left_x)*(orpz-or_p_down_z)-(orpz-or_p_left_z)*(orpy-or_p_down_y)*(orpx-or_p_anter_x);
                Jaco_4=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
                Jaco_4=Jaco_4*top_co;
        
                if Jaco_4>=0
                    orthird_4=0;
                else
                    orthird_4=log(-Jaco_4+1)*Tcoe_n;
                end
        
        
                %Jaco5 point 1 4 6 fff
                %Jaco_5=(p_right_x-orpx)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(p_up_x-orpx)-(p_right_z-orpz)*(p_poster_y-orpy)*(p_up_x-orpx)-(p_right_y-orpy)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(p_right_x-orpx);
%                              Jaco_5=(or_p_poster_x-orpx)*(or_p_right_y-orpy)*(or_p_up_z-orpz)+(or_p_poster_y-orpy)*(or_p_right_z-orpz)*(or_p_up_x-orpx)+(or_p_right_x-orpx)*(or_p_up_y-orpy)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(or_p_right_y-orpy)*(or_p_up_x-orpx)-(or_p_poster_y-orpy)*(or_p_right_x-orpx)*(or_p_up_z-orpz)-(or_p_right_z-orpz)*(or_p_up_y-orpy)*(or_p_poster_x-orpx);
                Jaco_5=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
                Jaco_5=Jaco_5*top_co;
                if Jaco_5>=0
                   orthird_5=0;
                else
                   orthird_5=log(-Jaco_5+1)*Tcoe_n;
                end
        
        
        
                 %Jaco6 point 2 3 6
                 
%                   Jaco_6=(or_p_poster_x-orpx)*(orpy-or_p_left_y)*(orpz-or_p_down_z)+(or_p_poster_y-orpy)*(orpz-or_p_left_z)*(orpx-or_p_down_x)+(orpx-or_p_left_x)*(orpy-or_p_down_y)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(orpy-or_p_left_y)*(orpx-or_p_down_x)-(or_p_poster_y-orpy)*(orpx-or_p_left_x)*(orpz-or_p_down_z)-(orpz-or_p_left_z)*(orpy-or_p_down_y)*(or_p_poster_x-orpx);
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
%                      Jaco_7=(or_p_poster_x-orpx)*(or_p_right_y-orpy)*(orpz-or_p_down_z)+(or_p_poster_y-orpy)*(or_p_right_z-orpz)*(orpx-or_p_down_x)+(or_p_right_x-orpx)*(orpy-or_p_down_y)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(or_p_right_y-orpy)*(orpx-or_p_down_x)-(or_p_poster_y-orpy)*(or_p_right_x-orpx)*(orpz-or_p_down_z)-(or_p_right_z-orpz)*(orpy-or_p_down_y)*(or_p_poster_x-orpx);
               Jaco_7=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
               Jaco_7=Jaco_7*top_co;
               if Jaco_7>=0
                   orthird_7=0;
               else
                   orthird_7=log(-Jaco_7+1)*Tcoe_n;
               end
        
               %Jaco8 point 1 3 6
              
               %Jaco_8=(orpx-p_left_x)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(p_up_x-orpx)-(orpz-p_left_z)*(p_poster_y-orpy)*(p_up_x-orpx)-(orpy-p_left_y)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(orpx-p_left_x);
%           Jaco_8=(or_p_poster_x-orpx)*(orpy-or_p_left_y)*(or_p_up_z-orpz)+(or_p_poster_y-orpy)*(orpz-or_p_left_z)*(or_p_up_x-orpx)+(orpx-or_p_left_x)*(or_p_up_y-orpy)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(orpy-or_p_left_y)*(or_p_up_x-orpx)-(or_p_poster_y-orpy)*(orpx-or_p_left_x)*(or_p_up_z-orpz)-(orpz-or_p_left_z)*(or_p_up_y-orpy)*(or_p_poster_x-orpx);
               Jaco_8=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
               Jaco_8=Jaco_8*top_co;
          
                if Jaco_8>=0
                    orthird_8=0;
                else
                    orthird_8=log(-Jaco_8+1)*Tcoe_n;
                end
tpenergy = orthird_1+orthird_2+orthird_3+orthird_4+orthird_5+orthird_6+orthird_7+orthird_8;