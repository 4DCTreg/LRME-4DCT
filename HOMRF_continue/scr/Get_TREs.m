
function [TRE,STD]=Get_TREs(landregis,pts_mov,pts_fix,spc)

r_phy=spc(1);
c_phy=spc(2);
s_phy=spc(3);

T50=pts_fix;

T00=pts_mov;
Tregis=landregis;

distanceor=zeros(300,1);
distanceregis=zeros(300,1);
for i=1:300
    distanceor(i)=sqrt((T00(i,1)-T50(i,1))^2*r_phy^2+(T00(i,2)-T50(i,2))^2*c_phy^2+(T00(i,3)-T50(i,3))^2*s_phy^2);
    distanceregis(i)=sqrt((Tregis(i,1)-T00(i,1))^2*r_phy^2+(Tregis(i,2)-T00(i,2))^2*c_phy^2+(Tregis(i,3)-T00(i,3))^2*s_phy^2);
end
TRE=sum(distanceregis(:))/300;
STD=std(distanceregis);
end
