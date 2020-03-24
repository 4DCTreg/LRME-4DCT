function [x_index,y_index,z_index]=Label_Coordinate(labels)
x_index=zeros(labels.nlabels,1);
y_index=zeros(labels.nlabels,1);
z_index=zeros(labels.nlabels,1);
for num=1:labels.nlabels
    z_index(num)=fix((num-1)/labels.n_each_layer)-labels.hz;
    intermediate=mod(num,labels.n_each_layer);
      if intermediate==0
          
          intermediate=labels.n_each_layer;
      end
     y_index(num)=floor((intermediate-0.5)/labels.sx)-labels.hy;
     x_index(num)=mod((intermediate-1),labels.sx)-labels.hx;
end
