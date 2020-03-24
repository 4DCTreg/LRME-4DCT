function dist=Get_Smooth_tot(sx,sy,sz,spc)
% spc=[1,1,1];
dis_max=floor(sqrt(((sx-1)*spc(1))^2+((sy-1)*spc(2))^2+((sz-1)*spc(3))^2))-5;
dist=single(zeros(sx*sy*sz,sx*sy*sz));
totnumlabel=sx*sy*sz;
xc=spc(1);
yc=spc(2);
zc=spc(3);

parfor index=1:totnumlabel*totnumlabel
    r=mod(index,totnumlabel);
    if r==0
        r=totnumlabel;
    end
    z1=-fix((r-1)/(sx*sy))+floor(sz/2);
    intermediate=mod(r,sx*sy);
    if intermediate==0
        intermediate=sx*sy;
    end
    x1=mod((intermediate-1),sy)-floor(sx/2);
    y1=floor((intermediate-0.5)/sx)-floor(sy/2);
    
    c=fix((index-1)/totnumlabel)+1;
    z2=-fix((c-1)/(sx*sy))+floor(sz/2);
    intermediate=mod(c,sx*sy);
    if intermediate==0
       intermediate=sx*sy;
    end
    x2=mod((intermediate-1),sy)-floor(sx/2);
    y2=floor((intermediate-0.5)/sx)-floor(sy/2);
    dis_diff=sqrt(((x1-x2)*xc)^2+((y1-y2)*yc)^2+((z1-z2)*zc)^2);

        if dis_diff>=dis_max
            dist(index)=dis_max;
        else
        
        dist(index)=sqrt(((x1-x2)*xc)^2+((y1-y2)*yc)^2+((z1-z2)*zc)^2)
        end
end

end