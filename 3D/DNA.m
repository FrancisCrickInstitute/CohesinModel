classdef DNA < handle
    properties
        D    % DNA Segment length
        N    % Number of segments
        Lp = 50;
        r    % positions
        u    % unit verctors
        E
        
        eb = 1.57;
        eta = 21;
        gamma = 0.98;
        eparal = 779;
        eperp = 1081;
        
    end
    methods
        function dna = DNA(D,N)
            dna.N = N;
            dna.D = D/dna.Lp;
            dna.r = zeros(N,3);
            dna.u = zeros(N,3);
            for i=1:N
                dna.r(i,1) = (i-1)*dna.D;
                dna.u(i,1) = 1;
            end
            dna.E = zeros(N,1);

            A=load('dssWLCparams.txt');
            
            dna.eb = interp1(A(:,1),A(:,2),dna.D);
            dna.gamma = interp1(A(:,1),A(:,3),dna.D);
            dna.eparal = interp1(A(:,1),A(:,4),dna.D);
            dna.eperp = interp1(A(:,1),A(:,5),dna.D);
            dna.eta = interp1(A(:,1),A(:,6),dna.D);

        end

        function dna = Energy(dna)
            for i = 2:dna.N
                e = 0;
                R = dna.r(i,:) - dna.r(i-1,:);
                Rperp = R - sum(R.*dna.u(i-1,:))*dna.u(i-1,:);
                var1 = dna.u(i,:) - dna.u(i-1,:) - dna.eta*Rperp;
                e = e + dna.eb*sum(var1.*var1)/2/dna.D;
                var2 = sum(R.*dna.u(i-1,:)) - dna.D*dna.gamma;
                e = e + dna.eparal*var2*var2/2/dna.D;
                e = e + dna.eperp*(sum(Rperp.*Rperp))/2/dna.D;
                dna.E(i) = e;
            end
        end
        function [dna,alpha] = Angles(dna)
            alpha = zeros(dna.N,1);
            for i=2:dna.N
                CosTheta = dot(dna.u(i,:),dna.u(i-1,:))/(norm(dna.u(i,:))*norm(dna.u(i-1,:)));
                alpha(i) = sign(sum(cross(dna.u(i,:),dna.u(i-1,:))))*(acosd(CosTheta));
            end
        end
        function dna = ReEnergy(dna,j)
            if j>1
                i = j;
                e = 0;
                R = dna.r(i,:) - dna.r(i-1,:);
                Rperp = R - sum(R.*dna.u(i-1,:))*dna.u(i-1,:);
                var1 = dna.u(i,:) - dna.u(i-1,:) - dna.eta*Rperp;
                e = e + dna.eb*sum(var1.*var1)/2/dna.D;
                var2 = sum(R.*dna.u(i-1,:)) - dna.D*dna.gamma;
                e = e + dna.eparal*var2*var2/2/dna.D;
                e = e + dna.eperp*(sum(Rperp.*Rperp))/2/dna.D;
                dna.E(i) = e;
            end
                if j<dna.N
                    i = j+1;
                    e = 0;
                    R = dna.r(i,:) - dna.r(i-1,:);
                    Rperp = R - sum(R.*dna.u(i-1,:))*dna.u(i-1,:);
                    var1 = dna.u(i,:) - dna.u(i-1,:) - dna.eta*Rperp;
                    e = e + dna.eb*sum(var1.*var1)/2/dna.D;
                    var2 = sum(R.*dna.u(i-1,:)) - dna.D*dna.gamma;
                    e = e + dna.eparal*var2*var2/2/dna.D;
                    e = e + dna.eperp*(sum(Rperp.*Rperp))/2/dna.D;
                    dna.E(i) = e;
                end
        end
        function dna = draw(dna,nn)
            hold on
            if nargin>1
                if sum(size(nn))>2
                    DU = nn(1:3,1:3);
                    for i=1:dna.N
%                         (DU*(dna.u(rn(i),:)'))'
                        x1 = dna.Lp*dna.r(i,:);
                        x1 = x1 + dna.Lp*nn(4,1:3);
                        x1 = (DU*(x1'))';
%                         x2 = dna.Lp*(dna.r(i,:) + dna.D*dna.u(i,:))+ dna.Lp*nn(4,1:3);
%                         x2 = (DU*(x2'))';
%                         plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'r');
                        if i==2
                            x2 = dna.Lp*dna.r(i-1,:);
                            x2 = x2 + dna.Lp*nn(4,1:3);
                            x2 = (DU*(x2'))';
                            h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'k');
                            set(h,'LineWidth',2)
                        elseif i>2
                            x2 = dna.Lp*dna.r(i-1,:);
                            x2 = x2 + dna.Lp*nn(4,1:3);
                            x2 = (DU*(x2'))';
                            h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'b');
                            set(h,'LineWidth',2)
                        end
                    end
                else
                    x1 = dna.Lp*dna.r(nn,:);
                    x2 = dna.Lp*dna.r(nn-1,:);
                    h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)]);
                    set(h,'LineWidth',2)
                end
            else
                for i=1:dna.N
                    x1 = dna.Lp*dna.r(i,:);
    %                 x2 = dna.Lp*(dna.r(i,:) + dna.D*dna.u(i,:));
    %                 plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'r');
                    if i==2
                        x2 = dna.Lp*dna.r(i-1,:);
                        h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'k');
                        set(h,'LineWidth',2)
                    elseif i>2
                        x2 = dna.Lp*dna.r(i-1,:);
                        h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'b');
                        set(h,'LineWidth',2)
                    end
                end
            end
        end
        function dna = makecircle(dna,phi)
            M = makehgtform('axisrotate',[0 1 0],phi*pi/180);
            for i=2:dna.N
                dna.r(i,:) = dna.r(i-1,:) + dna.D*dna.u(i-1,:);
                dna.u(i,:) = (M(1:3,1:3)*(dna.u(i-1,:)'))';
            end
            R = dna.D/2/tan(phi*pi/180/2);
            E1 = dna.D*dna.N/R/R/2;
            disp(['Energy is ' num2str(E1) ' kT'])
        end
        function [dna,EE] = MMK(dna,Niter,F)
            kT = 1;
            % F - random rotation matrix
            % choose whether to change position or orientation
            ri = randi(2,Niter,1)-1;
            % Chose element
            rn = randi(dna.N,Niter,1);
            % random addition to position
            rv = (rand(Niter,3)-0.5);
            % probability
            p = rand(Niter,1);
            
            EE = zeros(Niter,1);
            dna.Energy;
            E0Full = dna.E;
            E0 = sum(dna.E);
            
%             rt = randn(Niter,3);
%             rphi = (pi/2)*(rand(Niter,1)-0.5);
            
            for i=1:Niter
%                 rts = rt(i,:)./sqrt(sum(rt(i,:).*rt(i,:)));
%                 M = makehgtform('axisrotate',rts,rphi(i));
%                 DU = M(1:3,1:3);
%                 DU = F(:,:,i);
                if ri(i)==1
                    % change position of bead rn in this case
                    varsave = dna.r(rn(i),:);
                    dna.r(rn(i),:) = dna.r(rn(i),:) + 0.1*rv(i,:);
                else
                    % change orientation of bead rn in this case
                    [F,DU] = F.RandomSample;
                    varsave = dna.u(rn(i),:);
                    dna.u(rn(i),:) = (DU*(dna.u(rn(i),:)'))';
                end
                dna.ReEnergy(rn(i));
                E1 = sum(dna.E);

                if E1<E0
                    E0 = E1;
                    E0Full = dna.E;
                    update = 1;
                else
                    ro = exp(-(E1-E0)/kT);
                    if ro>p(i)
                        E0 = E1;
                        E0Full = dna.E;
                        update = 1;
                    else
                        if ri(i)==1
                            dna.r(rn(i),:) = varsave;
                        else
                            dna.u(rn(i),:) = varsave;
                        end
                        dna.E = E0Full;
%                         dna.ReEnergy(rn(i));
%                         qw = (dna.E==E0Full);
                        update = 0;
                    end
                end
                EE(i) = E0;
                at=1;
            end
        end
    end
end
