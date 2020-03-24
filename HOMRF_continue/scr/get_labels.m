function labelstot=get_labels(previous,x_tot,y_tot,z_tot)
labelstot=zeros(size(previous,1),1);
x_h=fix(x_tot/2)+1;
y_h=fix(y_tot/2)+1;
z_h=fix(z_tot/2)+1;
for n=1:size(previous,1)
    x_index=x_h+previous(n,1);
    y_index=y_h+previous(n,2);
    z_index=z_h+previous(n,3);
    labelstot(n)=sub2ind([x_tot,y_tot,z_tot],x_index,y_index,z_index);
end
end