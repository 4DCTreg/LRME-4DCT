function [x_index_tp,y_index_tp,z_index_tp]=Label_Coordinate_TP(labels,space_coefficient)
x_index_tp=zeros(labels.nlabels,1);
y_index_tp=zeros(labels.nlabels,1);
z_index_tp=zeros(labels.nlabels,1);
for num=1:labels.nlabels
    z_index_tp(num)=(fix((num-1)/labels.n_each_layer)-labels.hz)*space_coefficient;
    intermediate=mod(num,labels.n_each_layer);
      if intermediate==0
          
          intermediate=labels.n_each_layer;
      end
     y_index_tp(num)=(floor((intermediate-0.5)/labels.sx)-labels.hy)*space_coefficient;
     x_index_tp(num)=(mod((intermediate-1),labels.sx)-labels.hx)*space_coefficient;
end
