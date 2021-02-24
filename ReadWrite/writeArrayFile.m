function writeArrayFile(filename,sources,detectors,refpts_10_2p5_label)

%REDUNDANT?
%Writes out array solution file

%Find 10-2p5 label:

for i = 1:length(sources)
    source_labels{i} = refpts_10_2p5_label{sources(i)};
end

for i = 1:length(sources)
    detector_labels{i} = refpts_10_2p5_label{detectors(i)};
end

save(filename,'sources','detectors','source_labels','detector_labels');
