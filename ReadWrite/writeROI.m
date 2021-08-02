function pathnameROI = writeROI(ROI,pathnameROI)

fid = fopen([pathnameROI(1:end-4) '.txt'],'w');
fprintf(fid,'%s\n',num2str(ROI.gmNodeList'));
fprintf(fid,'%s\n',num2str(ROI.solutionSpaceNodeList'));
fclose(fid);
