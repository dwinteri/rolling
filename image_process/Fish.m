function [ frame, timestamp, binfish, theta ] = Fish( vid )
%Fish [ frame, timestamp, binfish, theta ] = Fish( vid )

trigger(vid)
[frame, timestamp]=getdata(vid);
binfish = ExtractTailH(frame);
t = regionprops(uint8(binfish),'Orientation');
theta = t.Orientation;


end