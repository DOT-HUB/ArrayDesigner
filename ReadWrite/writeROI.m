function pathnameROI = writeROI(ROI,pathnameROI)

fid = fopen(pathnameROI,'w');
fprintf(fid,'%s\n',num2str(ROI.gmNodeList'));
fprintf(fid,'%s\n',num2str(ROI.solutionSpaceNodeList'));
fclose(fid);

