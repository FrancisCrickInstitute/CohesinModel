
if exist('VFN','var')
    global VFN
    VFN = VFN+1;
else
    global VFN
    VFN = 1;
end

v = VideoWriter(['testnewvideo' num2str(VFN) '.avi'],'MPEG-4');
v.FrameRate = 15;
v.Quality = 100;
open(v);

figure
hold on
view(az,el)
grid on
axis equal
S = cdsave(1).d.Lp;
Nframe = length(cdsave)+1;


cla
k = round(linspace(1,length(cdsave),200));
for i=1:length(cdsave)
%     cdsave(k(i)).draw(2);
    cdsave(i).draw;
end
h=axis;
axis(h)
% 

cla
u.draw(2);
title(['Frame 1 / ' num2str(Nframe)])
axis(h)
grid on
writeVideo(v,getframe(gcf));

for i=1:length(cdsave)
    cla
    cdsave(i).draw;
    title(['Frame ' num2str(i+1) ' / ' num2str(Nframe)])
    pause(0.05)
    writeVideo(v,getframe(gcf));
end
close(v)