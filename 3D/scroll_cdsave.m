function [] = scroll_cdsave(cdsave,nn)
global IC
% move figure with arrow keys.
S.fh = figure('units','pixels',...
              'name','move_fig',...
              'numbertitle','off',...
              'keypressfcn',@fh_kpfcn);
S.tx = uicontrol('style','text',...
                 'units','pixels',...
                 'position',[0 0 80 20],...
                 'fontweight','bold');
             
S.cdsave = cdsave;
if nargin>1
    S.nn = nn;
else
    S.nn = 0;
end


az = 32;
el = 22;

hold on
view(az,el)
grid on
axis equal
Nframe = length(cdsave)+1;

cla
k = round(linspace(1,length(cdsave),min(length(cdsave),50)));
for i=1:length(k)
    if nargin>1
        cdsave(k(i)).draw(nn);
    else
        cdsave(k(i)).draw;
    end
end
h=axis;
S.h = h;
axis(h)
% 

cla
    if nargin>1
        cdsave(1).draw(nn);
    else
        cdsave(1).draw;
    end
title(['Frame 1 / ' num2str(Nframe)])
axis(h)
grid on
IC = 1;

guidata(S.fh,S) 



function [] = fh_kpfcn(H,E)    
global IC
% Figure keypressfcn
S = guidata(H);
% P = get(S.fh,'position');
set(S.tx,'string',E.Key)
switch E.Key
    case 'rightarrow'
        if IC<length(S.cdsave)
            IC = IC + 1;
            cla
            if S.nn>0
               S.cdsave(IC).draw(S.nn);
            else
               S.cdsave(IC).draw;
            end
            title(['Frame ' num2str(IC) ' / ' num2str(length(S.cdsave))])
            axis(S.h)
            grid on
        end
    case 'leftarrow'
        if IC>1
            IC = IC - 1;
            cla
            if S.nn>0
               S.cdsave(IC).draw(S.nn);
            else
               S.cdsave(IC).draw;
            end
            title(['Frame ' num2str(IC) ' / ' num2str(length(S.cdsave))])
            axis(S.h)
            grid on
        end

    case 'uparrow'
%         set(S.fh,'pos',P+[0 5 0 0])
    case 'downarrow'
%         set(S.fh,'pos',P+[0 -5 0 0])
    otherwise  
end