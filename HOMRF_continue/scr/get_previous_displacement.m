function previous=get_previous_displacement(Knots_n,griddim)

dim_r=griddim(1);
dim_c=griddim(2);
dim_s=griddim(3);
previous=zeros(dim_r*dim_c*dim_s,3);
previous(:,1)=reshape(Knots_n(:,:,:,1),[dim_r*dim_c*dim_s,1]);
previous(:,2)=reshape(Knots_n(:,:,:,2),[dim_r*dim_c*dim_s,1]);
previous(:,3)=reshape(Knots_n(:,:,:,3),[dim_r*dim_c*dim_s,1]);
% for i=1:dim_r
%         for j=1:dim_c
%             for k=1:dim_s
%                 i_1=O_trans_X(i,j,k,1);
%                 i_2=O_trans_X(i,j,k,2);
%                 j_1=O_trans_Y(i,j,k,1);
%                 j_2=O_trans_Y(i,j,k,2);
%                 k_1=O_trans_Z(i,j,k,1);
%                 k_2=O_trans_Z(i,j,k,2);
%                 x_width=i_2-i_1+1;
%                 y_length=j_2-j_1+1;
%                 z_hight=k_2-k_1+1;
%                 v_zone=x_width*y_length*z_hight;
%                 tot_x=Knots_n(i_1:i_2,j_1:j_2,k_1:k_2,1);
%                 tot_y=Knots_n(i_1:i_2,j_1:j_2,k_1:k_2,2);
%                 tot_z=Knots_n(i_1:i_2,j_1:j_2,k_1:k_2,3);
%                 tot_x_avg=sum(tot_x(:))/v_zone;
%                 tot_y_avg=sum(tot_y(:))/v_zone;
%                 tot_z_avg=sum(tot_z(:))/v_zone;
%                 previous((k-1)*dim_r*dim_c+(j-1)*dim_r+i,1)=tot_x_avg;
%                 previous((k-1)*dim_r*dim_c+(j-1)*dim_r+i,2)=tot_y_avg;
%                 previous((k-1)*dim_r*dim_c+(j-1)*dim_r+i,3)=tot_z_avg;
% %                 sad_cg=sum(sad_cg(:))/v_zone;
% %                 unary((k-1)*dim_r*dim_c+(j-1)*dim_r+i,nlb)=sad_cg.*50;
%             end
%         end
% end

end