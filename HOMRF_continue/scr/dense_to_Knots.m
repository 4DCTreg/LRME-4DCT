function [Knots_n] = dense_to_Knots(Knots_n_tot,O_trans_X,O_trans_Y,O_trans_Z)

[r,c,s,N] = size(O_trans_X);
Knots_n = zeros(r,c,s,3);
i_start = O_trans_X(1,1,1,1);
i_end = O_trans_X(1,1,1,2);
j_start = O_trans_Y(1,1,1,1);
j_end = O_trans_Y(1,1,1,2);
k_start = O_trans_Z(1,1,1,1);
k_end = O_trans_Z(1,1,1,2);
x_width = i_end-i_start+1;
y_length = j_end-j_start+1;
z_height = k_end-k_start+1;

for xd = 1:r
    for yd = 1:c
        for zd = 1:s
            i_start = O_trans_X(xd,yd,zd,1);
            i_end = O_trans_X(xd,yd,zd,2);
            j_start = O_trans_Y(xd,yd,zd,1);
            j_end = O_trans_Y(xd,yd,zd,2);
            k_start = O_trans_Z(xd,yd,zd,1);
            k_end = O_trans_Z(xd,yd,zd,2);
            X = Knots_n_tot(i_start:i_end,j_start:j_end,k_start:k_end,1);
            Y = Knots_n_tot(i_start:i_end,j_start:j_end,k_start:k_end,2);
            Z = Knots_n_tot(i_start:i_end,j_start:j_end,k_start:k_end,3);
            Knots_n(xd,yd,zd,1) = mean(X(:));
            Knots_n(xd,yd,zd,2) = mean(Y(:));
            Knots_n(xd,yd,zd,3) = mean(Z(:));           
        end
    end
end
   
