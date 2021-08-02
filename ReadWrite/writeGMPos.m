function writeGMPos(gmSurfaceMesh,pathnameHeadModel)

fid = fopen(fullfile(pathnameHeadModel,'/Utils/gmPos.txt'),'w');
for i = 1:size(gmSurfaceMesh.node,1)
    fprintf(fid,'%.6f %.6f %.6f \n',gmSurfaceMesh.node(i,:));
end
fclose(fid);  