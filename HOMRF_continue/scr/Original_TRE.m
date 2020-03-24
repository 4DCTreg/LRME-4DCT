function [TRE,TRE_std] = Original_TRE(pts_mov,pts_fix,spc)
koef = repmat(spc,[size(pts_mov,1),1]);
pt_errs_phys = sqrt( sum((  (pts_mov - pts_fix).*koef  ).^2, 2) );
 TRE = mean(pt_errs_phys);
 TRE_std = std(sqrt( sum((  (pts_mov - pts_fix).*koef  ).^2, 2) ));
