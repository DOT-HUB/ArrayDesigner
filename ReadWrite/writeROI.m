function pathnameROI = writeROI(ROI)

ROI.gmNodeList;
ROI.scalpPosList;

pathnameROI = ['test.' 'txt'];
fid = fopen(outname,'w');
fprintf(fid,'%s\n',ROI.gmNodeList);
fprintf(fid,'%s\n',ROI.scalpPosList);
fclose(fid);

