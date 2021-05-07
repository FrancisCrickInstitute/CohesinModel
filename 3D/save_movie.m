
if exist('VFN','var')
    global VFN
    VFN = VFN+1;
else
    global VFN
    VFN = 1;
end

v = VideoWriter(['newvideo' num2str(VFN) '.avi']);
v.FrameRate = 10;
v.Quality = 100;
open(v);

figure
hold on
view(az,el)
grid on
axis equal
S = dsave(1).Lp;
Nframe = length(dsave)+1;


cla
u.draw
for i=1:length(dsave)
    dsave(i).draw;
end
h=axis;
axis(h)
% 

cla
u.draw;
title(['Frame 1 / ' num2str(Nframe)])
axis(h)
grid on
writeVideo(v,getframe(gcf));

for i=1:length(dsave)
    cla
    dsave(i).draw;
    title(['Frame ' num2str(i+1) ' / ' num2str(Nframe)])
    pause(0.05)
    writeVideo(v,getframe(gcf));
end
close(v)