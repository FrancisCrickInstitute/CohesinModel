

u = CohesinDNA(SL,LL,[70 0 -4]);


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
looplength = zeros(length(cdsave),1);

for i=1:length(cdsave)
    v1 = cdsave(i).c.s3u;
    v2 = cdsave(i).c.e3u;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    alpha(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    v1 = cdsave(i).c.hr-cdsave(i).c.e3r;
    v2 = cdsave(i).c.e3r-cdsave(i).c.s3r;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    beta(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    v1 = cdsave(i).c.hr-cdsave(i).c.e3r+cdsave(i).c.hr-cdsave(i).c.e1r;
    v2 = cdsave(i).c.e3r-cdsave(i).c.s3r+cdsave(i).c.e1r-cdsave(i).c.s1r;
%     CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    a1 = cdsave(i).c.e3r-cdsave(i).c.s3r;
    a2 = cdsave(i).c.s1r-cdsave(i).c.s3r;
    c = cross(a1,a2);
    CosTheta = dot(v1,c)/(norm(v1)*norm(c));
%     gamma(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    gamma(i) = (acosd(CosTheta));
    v1 = cdsave(i).c.hr-cdsave(i).c.e1r;
    v2 = cdsave(i).c.e1r-cdsave(i).c.s1r;
    CosTheta = dot(v1,v2)/(norm(v1)*norm(v2));
    teta(i) = sign(sum(cross(v1,v2)))*(acosd(CosTheta));
    R = cdsave(i).c.s3r - cdsave(i).c.s1r;
    Rperp = R - sum(R.*cdsave(i).c.s1u)*cdsave(i).c.s1u;
    head_head(i) = sqrt(sum(Rperp.*Rperp));
    R = cdsave(i).c.hr - cdsave(i).c.s3r;
    s3h_d(i) = sqrt(sum(R.*R));
    R = cdsave(i).c.hr - cdsave(i).c.s1r;
    s1h_d(i) = sqrt(sum(R.*R));
    
    r0 = cdsave(i).c.hr;
    r1 = cdsave(i).c.s3r;
    r2 = cdsave(i).c.s1r;
    hhxz(i) = norm(cross(r0-r1,r0-r2))/norm(r2-r1);
    
    [cd,Rperp,pos,pos2]=cdsave(i).smc3pos;
    cohpos(i) = pos2;
    astate(i) = cdsave(i).c_state;
    
    smc1(i) = cdsave(i).smc1index;
    smc3(i) = cdsave(i).smc3index;
    hinge(i) = cdsave(i).hindex;
    
    hingestate(i) = cdsave(i).c.HingeState;
    
    
    [cd,R,pos1] = cdsave(i).smc1pos;
    [cd,R,pos2] = cdsave(i).smc3pos;
    looplength(i) = abs(hinge(i) - pos2);
end

figure
plot(cohpos/0.1,'bo-')
hold on
plot(smc1,'o-')
plot(smc3,'o-');
plot(hinge,'o-');
plot(looplength,'o-');
legend('cohpos','smc1','smc3','hinge','loop')

% figure
% plot(astates,'o-')
% hold on
% plot(HingeDetach,'d-')
% legend('astates','HingeDetach')

