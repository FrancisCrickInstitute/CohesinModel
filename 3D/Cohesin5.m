classdef Cohesin5
    properties
        s3r
        s3u
        
        e3r
        e3u
        
        hr
        hu
        
        e1r
        e1u
        
        s1r
        s1u
        
        Lp = 50;
        LA = 30/50;  % Length of head to elbow
        UA = 20/50;  % (Upper Arm) length of elbow to hinge
        
        HHT = 4/50; % Head-head equilibrium distance in ATP state
        
        HingeState    % Hinge state. 1 - straight; -1 - folded
        HeadsState    % Head-to-head state. 1- engaged; 0 - apart
        ATPstate
    end
    methods
        function C = Cohesin5(R)
            if nargin<1
                R = [0 0 0];
            end
            C.s3r = R/C.Lp;
            C.s3u = [0 0 1];
            
            C.e3r = R/C.Lp + [0 0 C.LA];
            C.e3u = [0 0 1];
            
            C.hr = R/C.Lp + [C.HHT/2 0 C.LA+C.UA];
            C.hu = [0 0 1];
            
            C.e1r = R/C.Lp + [C.HHT 0 C.LA];
            C.e1u = [0 0 1];
            
            C.s1r = R/C.Lp + [C.HHT 0 0];
            C.s1u = [0 0 1];
            
            C.HingeState = 1;
            C.HeadsState = 1;
            C.ATPstate = 1;
        end
        function C = draw(C,nn)
            if nargin>1
                hold on
                DU = nn(1:3,1:3);
                A = C.Lp*[C.s3r;C.e3r;C.hr;C.e1r;C.s1r;C.s3r];
                sz = size(A);
                for i=1:sz(1)
                    x = A(i,:);
                    x = x + C.Lp*nn(4,1:3);
                    x = (DU*(x'))';
                    A(i,:) = x;
                end
                h = plot3(A(:,1),A(:,2),A(:,3));
                set(h,'LineWidth',2)
                set(h,'Color',[0 0.4470 0.7410])
                [x,y,z] = sphere(5);
                x1 = C.Lp*(0.08*x/2 + C.s3r(1));
                y1 = C.Lp*(0.08*y/2 + C.s3r(2));
                z1 = C.Lp*(0.08*z/2 + C.s3r(3));
                sz = size(x1);
                for i=1:sz(1)
                    for j=1:sz(2)
                        x = [x1(i,j) y1(i,j) z1(i,j)];
                        x = x + C.Lp*nn(4,1:3);
                        x = (DU*(x'))';
                        x1(i,j) = x(1);
                        y1(i,j) = x(2);
                        z1(i,j) = x(3);
                    end
                end
                if C.ATPstate
                    s = surf(x1,y1,z1, 'FaceColor', [1 0 0]);
                else
                    s = surf(x1,y1,z1, 'FaceColor', [0.8 0.4470 0.7410]);
                end
                s.EdgeColor = 'none';  s.EdgeColor = 'none';
                [x,y,z] = sphere(5);
                x1 = C.Lp*(0.08*x/2 + C.s1r(1));
                y1 = C.Lp*(0.08*y/2 + C.s1r(2));
                z1 = C.Lp*(0.08*z/2 + C.s1r(3));
                sz = size(x1);
                for i=1:sz(1)
                    for j=1:sz(2)
                        x = [x1(i,j) y1(i,j) z1(i,j)];
                        x = x + C.Lp*nn(4,1:3);
                        x = (DU*(x'))';
                        x1(i,j) = x(1);
                        y1(i,j) = x(2);
                        z1(i,j) = x(3);
                    end
                end
                if C.ATPstate
                    s = surf(x1,y1,z1, 'FaceColor', [1 0.4470 0.7410]);
                else
                    s = surf(x1,y1,z1, 'FaceColor', [0 0.4470 0.7410]);
                end
                s.EdgeColor = 'none';  s.EdgeColor = 'none';
            else
                hold on
                A = C.Lp*[C.s3r;C.e3r;C.hr;C.e1r;C.s1r;C.s3r];
                h = plot3(A(:,1),A(:,2),A(:,3));
                set(h,'LineWidth',2)
                set(h,'Color',[0 0.4470 0.7410])
                [x,y,z] = sphere(5);
                x1 = C.Lp*(0.08*x/2 + C.s3r(1));
                y1 = C.Lp*(0.08*y/2 + C.s3r(2));
                z1 = C.Lp*(0.08*z/2 + C.s3r(3));
                if C.ATPstate
                    s = surf(x1,y1,z1, 'FaceColor', [1 0 0]);
                else
                    s = surf(x1,y1,z1, 'FaceColor', [0.8 0.4470 0.7410]);
                end
                s.EdgeColor = 'none';  s.EdgeColor = 'none';
                x1 = C.Lp*(0.08*x/2 + C.s1r(1));
                y1 = C.Lp*(0.08*y/2 + C.s1r(2));
                z1 = C.Lp*(0.08*z/2 + C.s1r(3));
                if C.ATPstate
                    s = surf(x1,y1,z1, 'FaceColor', [1 0.4470 0.7410]);
                else
                    s = surf(x1,y1,z1, 'FaceColor', [0 0.4470 0.7410]);
                end
                s.EdgeColor = 'none';  s.EdgeColor = 'none';
            end
        end
        function E = energy(C)
            E = 0;
            % Smc3 - Elbow
            % this is strong coiled-coil interaction. Assume Lp = 150 nm.
            % Delta = Lp/LA = 150/25 = 6
            % Koslover params: 0.17433 1.5727 0.96886 344.63 392.95 13.816 0.00075187 0.0077757
            
            % Lp = 150 nm
            EA = 1.57;
            gamma = 0.97;
            eparall = 344;
            eperp = 392;
            eta = 13.8;
            
            % Lp = 50 nm;
            % 0.52983 1.5417 0.90171 77.102 36.91 5.0121 0.006012 0.031445
            EA = 1.54;
            gamma = 0.9;
            eparall = 77;
            eperp = 36;
            eta = 5;
            
            
           
            R = C.e3r - C.s3r;
            Rperp = R - sum(R.*C.s3u)*C.s3u;
            var1 = C.e3u - C.HingeState*C.s3u - eta*Rperp;  % This depends on the hinge state. It is either 0 or 180
            E = E + EA*sum(var1.*var1)/2/C.LA;
            var2 = sum(R.*C.s3u) - C.LA*gamma;
            E = E + eparall*var2*var2/2/C.LA;
            E = E + eperp*(sum(Rperp.*Rperp))/2/C.LA;
            
            % Smc3Elbow - Hinge
            
            R = C.hr - C.e3r;
            Rperp = R - sum(R.*C.e3u)*C.e3u;
            var1 = C.hu - C.e3u - eta*Rperp;
            E = E + EA*sum(var1.*var1)/2/C.UA;
            var2 = sum(R.*C.e3u) - C.UA*gamma;
            E = E + eparall*var2*var2/2/C.UA;
            E = E + eperp*(sum(Rperp.*Rperp))/2/C.UA;
            
            % Hinge-Smc1Elbow
            
            R = C.hr - C.e1r;
            Rperp = R - sum(R.*C.e1u)*C.e1u;
            var1 = C.hu - C.e1u - eta*Rperp;
            E = E + EA*sum(var1.*var1)/2/C.UA;
            var2 = sum(R.*C.e1u) - C.UA*gamma;
            E = E + eparall*var2*var2/2/C.UA;
            E = E + eperp*(sum(Rperp.*Rperp))/2/C.UA;
            
            % Smc1 head-elbow
            
            R = C.e1r - C.s1r;
            Rperp = R - sum(R.*C.s1u)*C.s1u;
            var1 = C.e1u - C.HingeState*C.s1u - eta*Rperp;  % This depends on the hinge state. It is either 0 or 180
            E = E + EA*sum(var1.*var1)/2/C.LA;
            var2 = sum(R.*C.s1u) - C.LA*gamma;
            E = E + eparall*var2*var2/2/C.LA;
            E = E + eperp*(sum(Rperp.*Rperp))/2/C.LA;
            
            % Smc1-Smc3 head
            
            MultFactor = 50;
            
            R = C.s3r - C.s1r;
            Rparall = sum(R.*C.s1u)*C.s1u;  % this is different because smc unity vectors are parallel
            var1 = C.s3u - C.s1u - MultFactor*eta*Rparall;  % This depends on the hinge state. It is either 0 or 180
            E = E + MultFactor*EA*sum(var1.*var1)/2/C.HHT;
            
            Rperp = R - sum(R.*C.s1u)*C.s1u;
            var2 = sqrt(sum(Rperp.*Rperp)) - C.HHT;
            E = E + MultFactor*eparall*var2*var2/2/C.HHT;

            E = E + MultFactor*eperp*(sum(Rparall.*Rparall))/2/C.HHT;
            
            % Repulsion
            v1 = C.hr-C.e3r+C.hr-C.e1r;
%             v2 = C.e3r-C.s3r+C.e1r-C.s1r;
            a1 = C.e3r-C.s3r;
            a2 = C.s1r-C.s3r;
            c = cross(a1,a2);
            CosTheta = dot(v1,c)/(norm(v1)*norm(c));
            Theta = acosd(CosTheta);
            if Theta>80
                E = E + 0.5*(Theta-80)*(Theta-80);
            end

             % Smc3-hinge binding
            if C.HingeState==-1
                R = C.hr - C.s3r;
                var1 = sqrt(sum(R.*R)) - 0.15;
                E = E + 100*var1*var1;
            end
            
       end
        function [C,varsave] = ChangeItem(C,i1,i2,DD,DU)
            % i1 = rn(i)
            % i2 = ri(i)
                if i1==1
                    % smc3 head
                    if i2==1
                        % change position 
                        varsave = C.s3r;
                        C.s3r = C.s3r + DD;
                    else
                        % change orientation
                        varsave = C.s3u;
                        C.s3u = (DU*(C.s3u'))';
                    end
                elseif i1==2
                    if i2==1
                        % change position 
                        varsave = C.e3r;
                        C.e3r = C.e3r + DD;
                    else
                        % change orientation
                        varsave = C.e3u;
                        C.e3u = (DU*(C.e3u'))';
                    end
                elseif i1==3
                    if i2==1
                        % change position 
                        varsave = C.hr;
                        C.hr = C.hr + DD;
                    else
                        % change orientation
                        varsave = C.hu;
                        C.hu = (DU*(C.hu'))';
                    end
                elseif i1==4
                    if i2==1
                        % change position 
                        varsave = C.e1r;
                        C.e1r = C.e1r + DD;
                    else
                        % change orientation
                        varsave = C.e1u;
                        C.e1u = (DU*(C.e1u'))';
                    end
                elseif i1==5
                    if i2==1
                        % change position 
                        varsave = C.s1r;
                        C.s1r = C.s1r + DD;
                    else
                        % change orientation
                        varsave = C.s1u;
                        C.s1u = (DU*(C.s1u'))';
                    end
                end
        end
        function C = ChangeItemBack(C,i1,i2,varsave)
            % i1 = rn(i)
            % i2 = ri(i)
                if i1==1
                    % smc3 head
                    if i2==1
                        % change position 
                        C.s3r = varsave;
                    else
                        C.s3u = varsave;
                    end
                elseif i1==2
                    if i2==1
                        C.e3r = varsave;
                    else
                        C.e3u = varsave;
                    end
                elseif i1==3
                    if i2==1
                        C.hr = varsave;
                    else
                        C.hu = varsave;
                    end
                elseif i1==4
                    if i2==1
                        C.e1r = varsave;
                    else
                        C.e1u = varsave;
                    end
                elseif i1==5
                    if i2==1
                        C.s1r = varsave;
                    else
                        C.s1u = varsave;
                    end
                end
        end
        function [C,EE] = MMK(C,Niter,F)
            kT = 1.4;
            % F - random rotation matrix
            % choose whether to change position or orientation
            ri = randi(2,Niter,1)-1;
            % Chose element
            rn = randi(5,Niter,1);
            % random addition to position
            rv = 0.1*(rand(Niter,3)-0.5);
            % probability
            p = rand(Niter,1);

            E0 = C.energy;
            for i = 1:Niter
                [C,varsave] = ChangeItem(C,rn(i),ri(i),rv(i,:),F(:,:,i));
                E1 = C.energy;
                if E1<E0
                    E0 = E1;
                    update = 1;
                else
                    ro = exp(-(E1-E0)/kT);
                    if ro>p(i)
                        E0 = E1;
                        update = 1;
                    else
                        C = ChangeItemBack(C,rn(i),ri(i),varsave);
                        update = 0;
                    end
                end
                EE(i) = E0;
            end
        end
    end
end
