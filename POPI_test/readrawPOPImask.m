function G=readrawPOPImask(filename,init_size)
%readraw - read RAW format grey scale image file from Disk
%	Usuage  : G=readraw(filename)

	disp(['	Retrieving Image ' filename ' ...']);
	fid=fopen(filename,'rb');
	if (fid==-1)
	  	error('can not open imput image filem press CTRL-C to exit \n');
	  	pause
	end
	pixel=fread(fid,inf, 'unsigned char');
	fclose(fid);
    X= init_size(1);
    Y= init_size(2);
    Z= init_size(3);
%    [Y X]=size(pixel);
    G=reshape(pixel,[X,Y,Z]);

   end
