% temporarily unavailabel
function G=readvf(filename)

%readraw - read RAW format grey scale image file from Disk
%	Usuage  : G=readraw(filename)

	disp(['	Retrieving Image ' filename ' ...']);
	fid=fopen(filename,'r');
    if (fid==-1)
	  	error('can not open imput image filem press CTRL-C to exit \n');
	  	pause
    end
    head=fscanf(fid,'%20s',1);
    magic_number=fscanf(fid,'%10s',1);
    x=fscanf(fid,'%10s',1);
    y=fscanf(fid,'%10s',1);
    z=fscanf(fid,'%10s',1);
    size_x=str2double(x);
    size_y=str2double(y);
    size_z=str2double(z);
    grid_x=fscanf(fid,'%10s',1);
    grid_y=fscanf(fid,'%10s',1);
    grid_z=fscanf(fid,'%10s',1);
%     h=fscanf(fid,'%3f',1);
    fgetl(fid);
	pixel=fread(fid,inf,'float','l');
	fclose(fid);

%    [Y X]=size(pixel);
    G=reshape(pixel,[size_x,size_y,size_z,3]);

 end
