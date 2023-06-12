function jig2exo(fout, mesh, bcode, z)
% convert JIGSAW mesh to exo format
% fout :  output name
% mesh :  JIGSAW mesh
% bcode:  boundary code, 0-->interior node, 1-->boundary node, 2-->outlet node
% z    :  elevation [m]

if exist('/Users/xudo627/developments/petsc','dir')
    setenv('PETSC_DIR','/Users/xudo627/developments/petsc');
else
    prompt = "Enter your PETSc directory: ";
    txt = input(prompt,"s");
    if exist(txt,'dir')
        setenv('PETSC_DIR',txt);
    else
        error([txt ' does not exist!!!']);
    end
end

if any(isnan(z))
    error('NaNs exist in the elevation!');
end

coordx  = mesh.point.coord(:,1);
coordy  = mesh.point.coord(:,2);
connect = mesh.tria3.index(:,1:3);
num_of_nodes = length(coordx);

% Project the coordinate to [meter]
% proj = projcrs(32610);
% [coordx,coordy] = projfwd(proj,coordy,coordx);

coordx = round(coordx,1);
coordy = round(coordy,1);
z      = round(z,1);

inactive = false(size(connect,1),1);
for i = 1 : size(connect,1)
    if ~any(bcode(connect(i,:)) == 0)
        inactive(i) = true;
    end
end

if sum(inactive) > 0
    connect(inactive,:) = [];
    fprintf('\n****************************************\n');
    fprintf(['\n*****' num2str(sum(inactive)) ' inactive cells are removed!' '*****\n']);
    fprintf('\n****************************************\n');
end


fprintf('Writing to exoii file...');
write_exodus_file(coordx, coordy, [], connect, [fout '.exo']);
fprintf('Done!\n');

fprintf('\n Writing to .bcode file...');
fileID = fopen(strcat(fout,'.bcode'),'w');
fprintf(fileID, '%d \n', num_of_nodes);
fprintf(fileID, '%.2f \t %.2f \t %.2f \t %d \n', [coordx coordy z double(bcode)]');
fclose(fileID); 
fprintf('Done!\n');
    
end

