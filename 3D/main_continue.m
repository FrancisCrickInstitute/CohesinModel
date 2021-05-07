% load longDNA_extrudes_02_super
% cd.Copy(cdsave(137));

cd.hindex = 0;
cd.RepKleisin = 5;


ID = load('C:\Users\molodtm\Documents\MATLAB\Cohesin_ID.txt');
ID = ID + 1;
save('C:\Users\molodtm\Documents\MATLAB\Cohesin_ID.txt','ID','-ASCII');
az = 32;
el = 22;

rng('shuffle')


view(az,el)
grid on
axis equal
cd.draw;
pause(0.2)

% This function creates Niter random rotation vectors with angle pi/2 and
% returns the hangle to this array
% You can pick next random vector using F.RandomSample
% When all samples are used the function generates new set of Niter random
% vectors. 

HingeDetach=zeros(1,1746);
astates=zeros(1,1746);
astates(1)=0;



hw = waitbar(0,'Preparing... ');
disp(['ID of the simulation: ' num2str(ID) ])

tic
for i=1:N/Niter
    if astates(i+1)~=astates(i)
        if astates(i+1)
            cd.c_state = 1;
            cd.c.ATPstate = 1;
            [cd,R,pos] = cd.smc3pos;
            cd.smc3index = pos;
            cd.smc1index = 0;
            cd.c.HHT = (4+8)/cd.c.Lp;
            cd.c.HingeState = -1;
        else
            cd.c_state = 0;
            cd.c.ATPstate = 0;
            cd.smc3index = 0;
%             [cd,R,pos] = cd.smc1pos;
            cd.smc1index = 0;
            cd.c.HHT = (14+8)/cd.c.Lp;
%             [cd,R,pos] = cd.hingepos;
            [cd,R,pos] = cd.smc3pos;
            cd.hindex = pos-1;
            cd.c.HingeState = 1;
        end
    end
    if i>3
        cd.RepKleisin = 4;
    end
    if HingeDetach(i)==1
        cd.hindex = 0;
    end
    [cd,Et]=cd.MMK(Niter,F);
    cla
    cd.draw;
    title(['Iteration ' num2str(i) ' of ' num2str(N/Niter) '; ID: ' num2str(ID) ])
    pause(0.05)
    % save copy of current state to cdsave
    cdsave(i) = CohesinDNA(SL,LL,[0 0 0]);
    cdsave(i).Copy(cd);
    
    E((i-1)*Niter+1:i*Niter) = Et;
    waitbar(i/N*Niter,hw,['ID: ' num2str(ID) '; Iteration ' num2str(i) ' of ' num2str(N/Niter) ' is done']);

end
toc
close(hw)

figure
hold on
w = zeros(length(cdsave),1);
e = zeros(length(cdsave),1);
for i=1:length(cdsave)
    w(i) = mean(E((i-1)*Niter+1:i*Niter));
    e(i) = std(E((i-1)*Niter+1:i*Niter));
    h = errorbar(i,w(i),e(i),'bo');
    if cdsave(i).c_state
        set(h,'Color',[1 0 0]);
    end
end





% DNA end-to-end distance
R = zeros(N/Niter,1);
for i=1:length(cdsave)
    r1 = cdsave(i).d.r(1,:);
    r2 = cdsave(i).d.r(cd.d.N,:);
    R(i) = cd.d.Lp*sqrt(sum((r2-r1).*(r2-r1)));
end
disp(['L ' num2str(cd.d.Lp*cd.d.N*cd.d.D) ' nm'])



% Extract cohesin parameters
alpha = zeros(length(cdsave),1);
beta = zeros(length(cdsave),1);
gamma = zeros(length(cdsave),1);
teta = zeros(length(cdsave),1);
head_head = zeros(length(cdsave),1);
s3h_d = zeros(length(cdsave),1);
s1h_d = zeros(length(cdsave),1);

hhxz = zeros(length(cdsave),1);

cohpos = zeros(length(cdsave),1);
astate = zeros(length(cdsave),1);

smc1 = zeros(length(cdsave),1);
smc3 = zeros(length(cdsave),1);
hinge = zeros(length(cdsave),1);
hingestate = zeros(length(cdsave),1);

for i=1:length(cdsave)
    v1 = cdsave(i).c.s3u;
    v2 = cdsave(i).c.e3u;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    alpha(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    v1 = cdsave(i).c.e3u;
    v2 = cdsave(i).c.hu;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    beta(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    v1 = cdsave(i).c.hu;
    v2 = cdsave(i).c.e1u;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    gamma(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    v1 = cdsave(i).c.e1u;
    v2 = cdsave(i).c.s1u;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    teta(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    R = cdsave(i).c.s3r - cdsave(i).c.s1r;
    Rperp = R - sum(R.*cdsave(i).c.s1u)*cdsave(i).c.s1u;
    head_head(i) = sqrt(sum(Rperp.*Rperp));
    R = cdsave(i).c.hr - cdsave(i).c.s3r;
    s3h_d(i) = sqrt(sum(R.*R));
    R = cdsave(i).c.hr - cdsave(i).c.s1r;
    s1h_d(i) = sqrt(sum(R.*R));
    
    r = [R(1) R(3)];
    hhxz(i) = sqrt(sum(r.*r));
    
    [cd,Rperp,pos,pos2]=cdsave(i).smc3pos;
    cohpos(i) = pos2;
    astate(i) = cdsave(i).c_state;
    
    smc1(i) = cdsave(i).smc1index;
    smc3(i) = cdsave(i).smc3index;
    hinge(i) = cdsave(i).hindex;
    
    hingestate(i) = cdsave(i).c.HingeState;
end


