function Knots_n=get_Knots_pre(Knots_pre,O_trans_X,O_trans_Y,O_trans_Z,griddim)

Knots_n=zeros(griddim(1),griddim(2),griddim(3),3);
for index=1:griddim(1)*griddim(2)*griddim(3)
    [idx,idy,idz]=ind2sub([griddim(1),griddim(2),griddim(3)],index);
    size_cur_x=O_trans_X(idx,idy,idz,2)-O_trans_X(idx,idy,idz,1)+1;
    size_cur_y=O_trans_Y(idx,idy,idz,2)-O_trans_Y(idx,idy,idz,1)+1;
    size_cur_z=O_trans_Z(idx,idy,idz,2)-O_trans_Z(idx,idy,idz,1)+1;
    vsize=size_cur_x*size_cur_y*size_cur_z;
    Knots_x=Knots_pre(O_trans_X(idx,idy,idz,1):O_trans_X(idx,idy,idz,2),O_trans_Y(idx,idy,idz,1):O_trans_Y(idx,idy,idz,2),O_trans_Z(idx,idy,idz,1):O_trans_Z(idx,idy,idz,2),1);
    Knots_n(idx,idy,idz,1)=sum(Knots_x(:))./vsize;
    Knots_y=Knots_pre(O_trans_X(idx,idy,idz,1):O_trans_X(idx,idy,idz,2),O_trans_Y(idx,idy,idz,1):O_trans_Y(idx,idy,idz,2),O_trans_Z(idx,idy,idz,1):O_trans_Z(idx,idy,idz,2),2);
    Knots_n(idx,idy,idz,2)=sum(Knots_y(:))./vsize;
    Knots_z=Knots_pre(O_trans_X(idx,idy,idz,1):O_trans_X(idx,idy,idz,2),O_trans_Y(idx,idy,idz,1):O_trans_Y(idx,idy,idz,2),O_trans_Z(idx,idy,idz,1):O_trans_Z(idx,idy,idz,2),3);
    Knots_n(idx,idy,idz,3)=sum(Knots_z(:))./vsize;
end
