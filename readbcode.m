function [x,y,z,b] = readbcode(fname)
    M = dlmread(fname,'\t',1,0);
    x = M(:,1);
    y = M(:,2);
    z = M(:,3);
    b = M(:,4);
end