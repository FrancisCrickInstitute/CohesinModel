classdef CohesinDNA < handle
    properties
        d   % DNA
        c   % Cohesin
        c_state    %  1 if ATP bound, 0 - no ATP
        int   % interaction with 0 - none, 1 - smc1, 2 - hinge, 3 - smc3
        s3n   % number of DNA unit sliding through smc3

        SMC3_Slide     % stiffness between DNA and Smc3 for sliding
        SMC1_Slide     % stiffness between DNA and Smc1 for perpendicular interaction
        SMC3   % Affinity
        SMC1   % Affinity
        Hinge  % Affinity
        
        smc3index   % dna segment index with which smc3 interacts
        smc1index   % dna segment index with which smc1 interacts
        hindex      % dna segment index with which hinge interacts
        
        K_bind          % Rate of ATP binding
        K_release       % Rate of ADP release
        
        
        K_h_bind        % rate of hinge binding to DNA given close proximity
        K_h_release     % rate of hinge release if bound
        
        K_hinge_bend      % rate of hinge bending
        K_hinge_unbend    % rate of hinge unbending
        
        Rep            % Repulsion potential
        RepKleisin     % This should be 5 for repulsiotn by all subunits of the ring and 4 for no open smc1-smc3

    end
    methods
        function cd = CohesinDNA(DNAD,DNAN,CohesinR)
            % Class constructor
            cd.d = DNA(DNAD,DNAN);
            cd.c = Cohesin5(CohesinR);
            cd.c_state = 1;
            cd.smc3index = 9;
            cd.smc1index = 0;
            cd.hindex = 0;
            
            cd.int = zeros(DNAN,1);
            cd.s3n = floor(DNAN/2);
            
            cd.SMC3_Slide = 1000;
            cd.SMC3 = 1000;
            cd.SMC1 = 1000;
            cd.Hinge = 1000;
            
            cd.SMC1_Slide = 1000;
            
            cd.K_bind = 1e-6;
            cd.K_release = 1e-6;

            cd.K_h_bind = 0*1e-2;
            cd.K_h_release = 1e99*1e-5;
            
            cd.K_hinge_bend = 0*1e-3;
            cd.K_hinge_unbend = 1e-3;
            
            cd.Rep = 30000;
            cd.RepKleisin = 5;
        end
        function cd = draw(cd,SW)
            % >> cd.draw(1) will draw the system with DNA start point at 0,0,0
            % and DNA end point on (1,0,0) axis
            % >> cd.draw(2) will draw the system with SMC3 head at 0,0,0
            % smc3-smc1 on (0,1,0) axis and smc3-elbow3 in Oyz plane
            if nargin>1
                if SW==1
                    R = cd.d.r(1,:);
                    U = cd.d.r(cd.d.N,:) - cd.d.r(1,:);
                    U=U/norm(U);
                    ax=cross([1 0 0],U);
                    CosTheta = dot([1 0 0],U);
                    alpha = sign(sum(ax))*(acos(CosTheta));
                    M = makehgtform('axisrotate',-sign(sum(ax))*ax,alpha);
                elseif SW==2
                    % Make first transform to bring smc3-smc1 vector to
                    % {0,1,0}
                    R = cd.c.s3r;
                    U = cd.c.s1r - cd.c.s3r;
                    U=U/norm(U);
                    ax=cross([0 1 0],U);
                    CosTheta = dot([0 1 0],U);
                    alpha = sign(sum(ax))*(acos(CosTheta));
                    if ~isnan(alpha)
                        M = makehgtform('axisrotate',-sign(sum(ax))*ax,alpha);
                    else
                        M = sign(sum([0 1 0].*U))*eye(4);
                    end
                    if alpha==0
                        M = eye(4);
                    end
                    % Find second transform to bring s3u vector to Ozy
                    % plane
                    DU = M(1:3,1:3);
                    u3 = cd.c.e3r - cd.c.s3r;
                    u3 = (DU*(u3'))';
                    u3 = [u3(1) 0 u3(3)];
                    CosTheta = dot([0 0 1],u3)/norm(u3);
                    alpha = (acos(CosTheta));
                    M2 = makehgtform('axisrotate',-[0 1 0],alpha);
%                     M(1:3,1:3) = DU;
                    M(1:3,1:3) = M2(1:3,1:3)*DU;
                    % check second transform and correct if needed
                    u3 = cd.c.e3r - cd.c.s3r;
                    u = ((M2(1:3,1:3)*DU)*(u3'))';
                    if (abs(u(1))) > 0.01
                        for i=0:0.01:2*pi
                            M2 = makehgtform('axisrotate',-[0 1 0],i);
                            u = ((M2(1:3,1:3)*DU)*(u3'))';
                            if (abs(u(1))<0.01)&&(u(3)>1e-10)
                                break
                            end
                        end
                        M(1:3,1:3) = M2(1:3,1:3)*DU;
                    end
                end
                M(4,1:3) = -R;
                cd.d.draw(M);
                cd.c.draw(M);
            else
                cd.d.draw;
                cd.c.draw;
                R = [0 0 0];
                M = eye(3,3);
            end
            if cd.hindex>0
                DU = M(1:3,1:3);
                x1 = cd.d.Lp*cd.c.hr;
                x1 = x1 - cd.d.Lp*R;
                x1 = (DU*(x1'))';
                x2 = cd.d.Lp*cd.d.r(cd.hindex,:);
                x2 = x2 - cd.d.Lp*R;
                x2 = (DU*(x2'))';
                h = plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'r');
                set(h,'LineWidth',3)
            end
        end
        function cd = InitLoop(cd,pos)
            % This function initializes loop of DNA and places cohesin at
            % 'pos' DNA segments from its tip
            % use pos = 0.5 for > turn
            if nargin<2
                pos = 1;
            end
            pos = sign(pos)*min(abs(pos),round(cd.d.N/2-2));
            n = round(cd.d.N/2);
            for i=1:n-1
                cd.d.r(i,:) = [-(i-1)*cd.d.D 0 0];
                cd.d.u(i,:) = [-1 0 0];
            end
            if round(pos)==pos
                cd.d.r(n,:) = [-(n-1)*cd.d.D 0 0];
                cd.d.u(n,:) = [0 -1 0];
                cd.d.r(n+1,:) = [-(n-1)*cd.d.D -cd.d.D 0];
                cd.d.u(n+1,:) = [1 0 0];
                ic = (n-1)*cd.d.D;
                for i=n+2:cd.d.N
                    ic = ic - cd.d.D;
                    cd.d.r(i,:) = [-ic -cd.d.D 0];
                    cd.d.u(i,:) = [1 0 0];
                end
            else
                cd.d.u(n-1,:) = [-1 -1 0]/sqrt(2);
                cd.d.r(n,:) = [-(n-2+0.5)*cd.d.D -0.5*cd.d.D 0];
                cd.d.u(n,:) = [1 -1 0]/sqrt(2);
                cd.d.r(n+1,:) = [-(n-2)*cd.d.D -cd.d.D 0];
                cd.d.u(n+1,:) = [1 0 0];
                ic = (n-2)*cd.d.D;
                for i=n+2:cd.d.N
                    ic = ic - cd.d.D;
                    cd.d.r(i,:) = [-ic -cd.d.D 0];
                    cd.d.u(i,:) = [1 0 0];
                end
            end
            pos = round(pos);
            
            for i=1:cd.d.N
                cd.d.r(i,:) = cd.d.r(i,:) + [0 0.05 0];
            end

%             if pos>0
%                 cd.c.s1r = [(n-1-pos)*cd.d.D -5/cd.d.Lp -5/cd.d.Lp];
%                 cd.c.e1r = cd.c.s1r + [0 0 cd.c.LA];
%                 cd.c.hr = [(n-1-pos)*cd.d.D cd.d.D/2 cd.c.LA+cd.c.UA];
%                 cd.c.s3r =[(n-1-pos)*cd.d.D 15/cd.d.Lp -5/cd.d.Lp];
%                 cd.c.e3r = cd.c.s3r + [0 0 cd.c.LA];
%             else
%                 pos =abs(pos);
%                 cd.c.s3r = [(n-1-pos)*cd.d.D -5/cd.d.Lp -5/cd.d.Lp];
%                 cd.c.e3r = cd.c.s3r + [0 0 cd.c.LA];
%                 cd.c.hr = [(n-1-pos)*cd.d.D cd.d.D/2 cd.c.LA+cd.c.UA];
%                 cd.c.s1r =[(n-1-pos)*cd.d.D 15/cd.d.Lp -5/cd.d.Lp];
%                 cd.c.e1r = cd.c.s1r + [0 0 cd.c.LA];
%             end
%             
%             cd.c.s3u = [0 0 1];
%             cd.c.e3u = [0 0 1];
%             cd.c.hu  = [0 0 1];
%             cd.c.s1u = [0 0 1];
%             cd.c.e1u = [0 0 1];
            FT = 3/8;
                pos = cd.d.N;
            a = 15 * pi/180;
            cd.c.s3r = [-round(FT*pos)*cd.d.D -2*cd.c.HHT/2 -2*cd.c.HHT/2];
            cd.c.e3r = cd.c.s3r + [-sin(a)*cd.c.LA -0.5 cos(a)*cd.c.LA];
            cd.c.hr = [-round(FT*pos)*cd.d.D-sin(a)*cd.c.LA -cd.c.HHT/2 cos(a)*cd.c.LA-cd.c.UA];
            cd.c.s1r = [-round(FT*pos)*cd.d.D 2*cd.c.HHT/2 -2*cd.c.HHT/2];
            cd.c.e1r = cd.c.s1r + [-sin(a)*cd.c.LA 0.5 cos(a)*cd.c.LA];

        end
        function cd = InitWalk(cd,pos)
            % This function initializes loop of DNA and places cohesin at
            % 'pos' DNA segments from its tip
            % use pos = 0.5 for > turn
            if nargin<2
                pos = 1;
            end
            pos = sign(pos)*min(abs(pos),round(cd.d.N/2-2));
            for i=1:cd.d.N
                cd.d.r(i,:) = [-(i-1)*cd.d.D+0.75 0 0];
                cd.d.u(i,:) = [-1 0 0];
            end
%             if round(pos)==pos
%                 cd.d.r(n,:) = [-(n-1)*cd.d.D 0 0];
%                 cd.d.u(n,:) = [0 -1 0];
%                 cd.d.r(n+1,:) = [-(n-1)*cd.d.D -cd.d.D 0];
%                 cd.d.u(n+1,:) = [1 0 0];
%                 ic = (n-1)*cd.d.D;
%                 for i=n+2:cd.d.N
%                     ic = ic - cd.d.D;
%                     cd.d.r(i,:) = [-ic -cd.d.D 0];
%                     cd.d.u(i,:) = [1 0 0];
%                 end
%             else
%                 cd.d.u(n-1,:) = [-1 -1 0]/sqrt(2);
%                 cd.d.r(n,:) = [-(n-2+0.5)*cd.d.D -0.5*cd.d.D 0];
%                 cd.d.u(n,:) = [1 -1 0]/sqrt(2);
%                 cd.d.r(n+1,:) = [-(n-2)*cd.d.D -cd.d.D 0];
%                 cd.d.u(n+1,:) = [1 0 0];
%                 ic = (n-2)*cd.d.D;
%                 for i=n+2:cd.d.N
%                     ic = ic - cd.d.D;
%                     cd.d.r(i,:) = [-ic -cd.d.D 0];
%                     cd.d.u(i,:) = [1 0 0];
%                 end
%             end
            pos = round(pos);
            
            for i=1:cd.d.N
                cd.d.r(i,:) = cd.d.r(i,:) + [0 0.05 0];
            end

            FT = 3/8;
                pos = cd.d.N;
            a = 15 * pi/180;
            cd.c.s3r = [-round(FT*pos)*cd.d.D -2*cd.c.HHT/2 -2*cd.c.HHT/2];
            cd.c.e3r = cd.c.s3r + [-sin(a)*cd.c.LA -0.5 cos(a)*cd.c.LA];
            cd.c.hr = [-round(FT*pos)*cd.d.D-sin(a)*cd.c.LA -cd.c.HHT/2 cos(a)*cd.c.LA-cd.c.UA];
            cd.c.s1r = [-round(FT*pos)*cd.d.D 2*cd.c.HHT/2 -2*cd.c.HHT/2];
            cd.c.e1r = cd.c.s1r + [-sin(a)*cd.c.LA 0.5 cos(a)*cd.c.LA];

        end
        function cd = InitLoopInCohesin(cd,pos)
            % This function initializes loop of DNA and places cohesin at
            % 'pos' DNA segments from its tip
            % use pos = 0.5 for > turn
            if nargin<2
                pos = cd.d.N/4;
            end
            pos = min(pos,round(cd.d.N/2-2));
            n = round(cd.d.N/2);
            for i=1:n-1
                cd.d.r(i,:) = [(i-1)*cd.d.D 0 0];
                cd.d.u(i,:) = [1 0 0];
            end
            if round(pos)==pos
                cd.d.r(n,:) = [(n-1)*cd.d.D 0 0];
                cd.d.u(n,:) = [0 1 0];
                cd.d.r(n+1,:) = [(n-1)*cd.d.D cd.d.D 0];
                cd.d.u(n+1,:) = [-1 0 0];
                ic = (n-1)*cd.d.D;
                for i=n+2:cd.d.N
                    ic = ic - cd.d.D;
                    cd.d.r(i,:) = [ic cd.d.D 0];
                    cd.d.u(i,:) = [-1 0 0];
                end
            else
                cd.d.u(n-1,:) = [1 1 0]/sqrt(2);
                cd.d.r(n,:) = [(n-2+0.5)*cd.d.D 0.5*cd.d.D 0];
                cd.d.u(n,:) = [-1 1 0]/sqrt(2);
                cd.d.r(n+1,:) = [(n-2)*cd.d.D cd.d.D 0];
                cd.d.u(n+1,:) = [-1 0 0];
                ic = (n-2)*cd.d.D;
                for i=n+2:cd.d.N
                    ic = ic - cd.d.D;
                    cd.d.r(i,:) = [ic cd.d.D 0];
                    cd.d.u(i,:) = [-1 0 0];
                end
            end
            pos = round(pos);
%             cd.c.s3r = [(n-1-pos)*cd.d.D 0 0];
%             cd.c.e3r = cd.c.s3r + [0 0 cd.c.LA];
%             cd.c.hr = cd.c.s3r + [0 cd.d.D/2 cd.c.LA+cd.c.UA];
%             cd.c.e1r = cd.c.s3r + [0 cd.d.D cd.c.LA];
%             cd.c.s1r = cd.c.s3r + [0 cd.d.D 0];
            cd.c.s1r = [pos*cd.d.D -10/cd.d.Lp 5/cd.d.Lp];
            cd.c.e1r = cd.c.s1r + [0 cd.c.LA 0];
            cd.c.hr = cd.c.s1r + [0 cd.c.LA+cd.c.UA -cd.d.D/2];
            cd.c.e3r = cd.c.s1r + [0 cd.c.LA -cd.d.D];
            cd.c.s3r = cd.c.s1r + [0 0 -cd.d.D];
            
            cd.c.s3u = [0 1 0];
            cd.c.e3u = [0 1 0];
            cd.c.hu  = [0 1 0];
            cd.c.s1u = [0 1 0];
            cd.c.e1u = [0 1 0];
        end
        function cd = PrepareATP(cd,pos)
            % This function prepares a state similar to Uhlmann's structure
            if nargin<2
                pos = cd.d.N;
            end
            for i=1:cd.d.N
                cd.d.r(i,:) = [(i-1)*cd.d.D 0 0];
                cd.d.u(i,:) = [1 0 0];
            end
            a = 30 * pi/180;
            cd.c.s3r = [round(pos/2)*cd.d.D -cd.c.HHT/2 -cd.c.HHT/2];
            cd.c.e3r = cd.c.s3r + [-sin(a)*cd.c.LA -0.5 cos(a)*cd.c.LA];
%             cd.c.hr = cd.c.e3r + [0 cd.c.HHT/2 -cd.c.UA];
            cd.c.hr = [round(pos/2)*cd.d.D-sin(a)*cd.c.LA 0 cos(a)*cd.c.LA-cd.c.UA];
            cd.c.s1r = [round(pos/2)*cd.d.D cd.c.HHT/2 -cd.c.HHT/2];
            cd.c.e1r = cd.c.s1r + [-sin(a)*cd.c.LA 0.5 cos(a)*cd.c.LA];

%             cd.c.s1r = [round(pos/2)*cd.d.D -cd.c.HHT/2 -cd.c.HHT/2];
%             cd.c.e1r = cd.c.s1r + [-sin(a)*cd.c.LA 0 cos(a)*cd.c.LA];
%             cd.c.hr = cd.c.e1r + [0 cd.c.HHT/2 -cd.c.UA];
%             cd.c.s3r = [round(pos/2)*cd.d.D cd.c.HHT/2 -cd.c.HHT/2];
%             cd.c.e3r = cd.c.s3r + [-sin(a)*cd.c.LA 0 cos(a)*cd.c.LA];
        
        end
        function cd = Perp(cd,pos)
            if nargin<2
                pos = cd.d.N/2;
            end
            pos = round(pos);
            for i=1:cd.d.N
                cd.d.r(i,:) = [(i-1)*cd.d.D 0 0];
                cd.d.u(i,:) = [1 0 0];
            end
%             cd.c.s1r = [pos*cd.d.D -cd.c.HHT/2 -cd.c.HHT/2];
%             cd.c.e1r = cd.c.s1r + [0 0 cd.c.LA];
%             cd.c.hr = [pos*cd.d.D 0 cd.c.LA+cd.c.UA];
%             cd.c.s3r = [pos*cd.d.D cd.c.HHT/2 -cd.c.HHT/2];
%             cd.c.e3r = cd.c.s3r + [0 0 cd.c.LA];

            cd.c.s3r = [pos*cd.d.D -cd.c.HHT/2 -cd.c.HHT/2];
            cd.c.e3r = cd.c.s3r + [0 0 cd.c.LA];
            cd.c.hr = [pos*cd.d.D 0 cd.c.LA+cd.c.UA];
            cd.c.s1r = [pos*cd.d.D cd.c.HHT/2 -cd.c.HHT/2];
            cd.c.e1r = cd.c.s1r + [0 0 cd.c.LA];
        end
        function [cd,Rperp,pos]=dist2(cd,v)
            % this function calculates shortest distance between vector v (smc3 if nothing provided) and dna
            % it finds three closest DNA segments and then splits distances in half 10 times to get precision of D/(2^10) 
            if nargin <2
                v = cd.c.s3r;
            end
            R = zeros(cd.d.N,1);
            for i=1:cd.d.N
                r = v - cd.d.r(i,:);
                R(i) = sum(r.*r);
            end
            [i,j]=min(R);
            pos = j(1);
            r = zeros(5,3);
            r(3,:) = cd.d.r(j(1),:);
            r(1,:) = cd.d.r(max(1,j(1)-1),:);
            r(5,:) = cd.d.r(min(j(1)+1,cd.d.N),:);
            R3 = i(1);
            R1 = R(max(1,j(1)-1));
            R5 = R(min(j(1)+1,cd.d.N));
            for k=1:10
                r(2,:) = (r(1,:)+r(3,:))/2;
                r(4,:) = (r(5,:)+r(3,:))/2;
                R = v - r(2,:);
                R2 = sum(R.*R);
                R = v - r(4,:);
                R4 = sum(R.*R);
                R = [R1 R2 R3 R4 R5];
                [i,j]=sort(R);
                R3 = i(1);
                r(3,:) = r(j(1),:);
                R1 = i(2);
                r(1,:) = r(j(2),:);
                R5 = i(3);
                r(5,:) = r(j(3),:);
            end
            Rperp = v - r(3,:);
        end
        function [cd,Rperp,pos,pos2]=smc3pos(cd)
            % shortest distance between smc3 and dna
            R = zeros(cd.d.N,1);
            S = zeros(cd.d.N,1);
            for i=1:cd.d.N
                r = cd.c.s3r - cd.d.r(i,:);
                R(i) = sum(r.*r);
                S(i) = sum(r.*cd.d.u(i,:));
%                 sum(R.*cd.d.u(i,:));
            end
            [i,j]=min(R);
            R = cd.c.s3r - cd.d.r(j,:);
%             plot3(50*[cd.c.s3r(1) cd.d.r(j,1)],50*[cd.c.s3r(2) cd.d.r(j,2)],50*[cd.c.s3r(3) cd.d.r(j,3)])
            Rperp = R - S(j)*cd.d.u(j,:);
%             pos = cd.d.D*j + S(j);
            pos = j;
            pos2 = cd.d.D*j + S(j);
        end
        function [cd,R,pos]=smc1pos(cd)
            % shortest distance between smc1 and dna
            R = zeros(cd.d.N,1);
            S = zeros(cd.d.N,1);
            for i=1:cd.d.N
                r = cd.c.s1r - cd.d.r(i,:);
                R(i) = sum(r.*r);
                S(i) = sum(r.*cd.d.u(i,:));
            end
            [i,j]=min(R);
%             for i=max(1,j-1):min(cd.d.N,j+1)
%                 R = cd.c.s3r - cd.d.r(i,:);
%                 sum(R.*cd.d.u(i,:));
%             end
            R = cd.c.s1r - cd.d.r(j,:);
%             Rperp = R - S(j)*cd.d.u(j,:);
%             pos = cd.d.D*j + S(j);
            pos = j;
        end
        function [cd,R,pos]=hingepos(cd)
            % shortest distance between hinge and dna segment
            R = zeros(cd.d.N,1);
            S = zeros(cd.d.N,1);
            for i=1:cd.d.N
                r = cd.c.hr - cd.d.r(i,:);
                R(i) = sum(r.*r);
                S(i) = sum(r.*cd.d.u(i,:));
            end
            [i,j]=min(R);
            R = cd.c.hr - cd.d.r(j,:);
%             Rperp = R - S(j)*cd.d.u(j,:);
            pos = j;
        end
        function [cd,Ec,Ed,Ecd]=energy(cd)
            % Stiffnesses
            
            % DNA and Cohesin energy
            Ec = cd.c.energy;
            cd.d.Energy;
            Ed = sum(cd.d.E);
            Ecd = 0;

            % Binding Point at Smc3
            u1 = (cd.c.s1r - cd.c.s3r);
            u1 = u1/norm(u1);
            u2 = (cd.c.e3r - cd.c.s3r);
            u2 = u2/norm(u2);
            S3 = cd.c.s3r + (u1*0.06 + u2*0.06)/2;
            
            % Sliding energy of Smc3 alng DNA
            [cd,R,pos] = cd.dist2(S3);   % Distance from SMC3 to DNA
            R = sum(R.*R);
%             if (R>0.0016-eps)&&(cd.SMC3_Slide>0)
%                 Ecd = 1e99;
%                 return
%             end
            Ecd = Ecd + cd.SMC3_Slide*R;
            % orientation of the SMC1-SMC3 vector with respect to DNA;
%             u = cd.c.s1r - cd.c.s3r;
            % For parallel vector use this
%             v1 = u/norm(u);
            
            % for perpendicular vector use this:
%             v1 = cross(cd.c.e3r-cd.c.s3r,u1);
            v1 = cross(cd.c.s3u,u1);
            v1 = v1/norm(v1);
            v1 = 2*v1-3*u1;         % Comment these for 90 degree angle
            v1 = v1/norm(v1);       % Comment these for 90 degree angle
            var1 = v1 - cd.d.u(pos,:);
            var2 = v1 + cd.d.u(pos,:);
            var3 = min(sum(var1.*var1),sum(var2.*var2));
            Ecd = Ecd + cd.SMC3_Slide*var3/25;
            
            
            % Binding energy
            if cd.smc3index > 0
                dr = cd.d.r(cd.smc3index,:);
                cr = cd.c.s3r;
                R = cr - dr;
                Ecd = Ecd + cd.SMC3*sum(R.*R);
            end
            if cd.smc1index > 0
                dr = cd.d.r(cd.smc1index,:);
                u1 = -u1;
                u2 = (cd.c.e1r - cd.c.s1r);
                u2 = u2/norm(u2);
                S3 = cd.c.s1r + (u1*0.06 + u2*0.06)/2;

                R = S3 - dr;
                Ecd = Ecd + cd.SMC1*sum(R.*R);
                
                 % for perpendicular to Smc1 vector use this:
                u = cd.c.s3r - cd.c.s1r;
%                 v1 = cross(cd.c.e1r-cd.c.s1r,u/norm(u));
                v1 = cross(cd.c.s1u,u/norm(u));
                v1 = v1/norm(v1);
                var1 = v1 + cd.d.u(cd.smc1index,:);
                var2 = v1 - cd.d.u(cd.smc1index,:);
                var3 = min(sum(var1.*var1),sum(var2.*var2));
                Ecd = Ecd + cd.SMC1_Slide*var3/25;                
            end
            if cd.hindex > 0
                dr = cd.d.r(cd.hindex,:);
                u1 = cd.c.e3r - cd.c.hr;
                u2 = cd.c.e1r - cd.c.hr;
                S3 = cd.c.hr + (u1*0.06 + u2*0.06)/2;
                R = S3 - dr;
                Ecd = Ecd + cd.Hinge*sum(R.*R);

                v1 = cross(u1,u2);
                v1 = v1/norm(v1);
                var1 = v1 + cd.d.u(cd.hindex,:);
                var2 = v1 - cd.d.u(cd.hindex,:);
                var3 = min(sum(var1.*var1),sum(var2.*var2));
                Ecd = Ecd + cd.Hinge*var3/25;                
            end
            
%             % Repulsion
            for i=2:cd.d.N
                for j=1:cd.RepKleisin
                    r = cd.DistBetweenSegments(i,j);
                    if r<0.04
                        Ecd = Ecd + cd.Rep*(r-0.04)*(r-0.04);
                    end
                end
            end
            
        end
        function [cd,Ec,Ed,Ecd]=ReEnergy(cd,j,E0c,E0d)
            if j>5
                cd.d.ReEnergy(j-5);
                Ec = E0c;
                Ed = sum(cd.d.E);
            else
                Ec = cd.c.energy;
                Ed = E0d;
            end
            Ecd = 0;

            % Binding Point at Smc3
            u1 = (cd.c.s1r - cd.c.s3r);
            u1 = u1/norm(u1);
            u2 = (cd.c.e3r - cd.c.s3r);
            u2 = u2/norm(u2);
            S3 = cd.c.s3r + (u1*0.06 + u2*0.06)/2;
            
            [cd,R,pos] = cd.dist2(S3);
            R = sum(R.*R);
%             if (R>0.0016-eps)&&(cd.SMC3_Slide>0)
%                 Ecd = 1e99;
%                 return
%             end
            Ecd = Ecd + cd.SMC3_Slide*R;
            % orientation of the SMC1-SMC3 vector with respect to DNA;
%             u = cd.c.s1r - cd.c.s3r;
            % For parallel vector use this
%             v1 = u/norm(u);
            % for perpendicular vector use this:
%             v1 = cross(cd.c.e3r-cd.c.s3r,u1);
            v1 = cross(cd.c.s3u,u1);
            v1 = v1/norm(v1);
            v1 = 2*v1-3*u1;         % Comment these for 90 degree angle
            v1 = v1/norm(v1);       % Comment these for 90 degree angle
            var1 = v1 + cd.d.u(pos,:);
            var2 = v1 - cd.d.u(pos,:);
            var3 = min(sum(var1.*var1),sum(var2.*var2));
            Ecd = Ecd + cd.SMC3_Slide*var3/25;
            
            if cd.smc3index > 0
                dr = cd.d.r(cd.smc3index,:);
                cr = cd.c.s3r;
                R = cr - dr;
                Ecd = Ecd + cd.SMC3*sum(R.*R);
            end
            if cd.smc1index > 0
                dr = cd.d.r(cd.smc1index,:);
                u1 = -u1;
                u2 = (cd.c.e1r - cd.c.s1r);
                u2 = u2/norm(u2);
                S3 = cd.c.s1r + (u1*0.06 + u2*0.06)/2;

                R = S3 - dr;
                Ecd = Ecd + cd.SMC1*sum(R.*R);
                 % for perpendicular to Smc1 vector use this:
                 
                u = cd.c.s3r - cd.c.s1r;
%                 v1 = cross(cd.c.e1r-cd.c.s1r,u/norm(u));
                v1 = cross(cd.c.s1u,u/norm(u));
                v1 = v1/norm(v1);
                var1 = v1 + cd.d.u(cd.smc1index,:);
                var2 = v1 - cd.d.u(cd.smc1index,:);
                var3 = min(sum(var1.*var1),sum(var2.*var2));
                Ecd = Ecd + cd.SMC1_Slide*var3/25;                
            end
            if cd.hindex > 0
                dr = cd.d.r(cd.hindex,:);
                u1 = cd.c.e3r - cd.c.hr;
                u2 = cd.c.e1r - cd.c.hr;
                S3 = cd.c.hr + (u1*0.06 + u2*0.06)/2;
                R = S3 - dr;
                Ecd = Ecd + cd.Hinge*sum(R.*R);

                v1 = cross(u1,u2);
                v1 = v1/norm(v1);
                var1 = v1 + cd.d.u(cd.hindex,:);
                var2 = v1 - cd.d.u(cd.hindex,:);
                var3 = min(sum(var1.*var1),sum(var2.*var2));
                Ecd = Ecd + cd.Hinge*var3/25;                
            end
            
%             % Repulsion
            for i=2:cd.d.N
                for k=1:cd.RepKleisin
                    r = cd.DistBetweenSegments(i,k);
                    if r<0.04
                        Ecd = Ecd + cd.Rep*(r-0.04)*(r-0.04);
                    end
                end
            end
        end
        function [intersect,pos] = CohesinDNASegmentIntersection(cd,n)
            % This function calcultes whether cohesin intersects DNA
            % segment n
            % Cohesin is split into 3 triangles:
            %           o hr
            %          / \  
            %         / 3 \  
            %        /_____\
            %   e3r  |\  2 |  e1r
            %        |  \  |  
            %        | 1  \|
            %        O-----O
            %      s3r     s1r
            %        
            % Algorythm is Möller and Trumbore (1997)
            % Modified from Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
            % intersect is 1 if segment intersects DNA

            intersect = 0;
            zero = eps;
            
            dir = cd.d.r(n,:)-cd.d.r(n-1,:);
            edge1 = cd.c.e3r-cd.c.s3r;          
            edge2 = cd.c.s1r-cd.c.s3r;
            tvec  = cd.d.r(n-1,:) - cd.c.s3r;          % vector from vert0 to ray origin
            pvec  = cross(dir,edge2);  % begin calculating determinant - also used to calculate U parameter
            det   = sum(edge1.*pvec);   % determinant of the matrix M = dot(edge1,pvec)
            angleOK = (abs(det)>eps);
            
            temp1 = cross(edge1,edge2);
            pos = sign(dot(temp1,tvec));
            
            if all(~angleOK)
                int1 = 0;
                intersect = intersect + int1;
            else
                det(~angleOK) = nan;              % change to avoid division by zero
                u    = sum(tvec.*pvec)./det;    % 1st barycentric coordinate
                v = nan+zeros(size(u));
                t=v;
                ok = (angleOK & u>=zero & u<=1.0+zero); % mask
                if ~any(ok)
                    int1 = ok;
                    intersect = intersect + int1;
                else
                    qvec = cross(tvec, edge1,2); % prepare to test V parameter
                    v(ok,:) = sum(dir.*qvec,2) ./ det; % 2nd barycentric coordinate
                    t(ok,:) = sum(edge2.*qvec,2)./det;
                    ok = (ok & v>=zero & u+v<=1.0+zero);
                    int1 = (ok & t>=zero & t<=1.0+zero);
                    intersect = intersect + int1;
                end
            end 
            
            edge1 = cd.c.e3r-cd.c.e1r;          
            edge2 = cd.c.s1r-cd.c.e1r;
            tvec  = cd.d.r(n-1,:) - cd.c.e1r;          
            pvec  = cross(dir,edge2);  
            det   = sum(edge1.*pvec);   
            angleOK = (abs(det)>eps);
            if all(~angleOK)
                int2 = 0;
                intersect = intersect + int2;
            else
                det(~angleOK) = nan;              
                u    = sum(tvec.*pvec)./det;    
                v = nan+zeros(size(u));
                t=v;
                ok = (angleOK & u>=zero & u<=1.0+zero); 
                if ~any(ok)
                    int2 = ok;
                    intersect = intersect + int2;
                else
                    qvec = cross(tvec, edge1,2); 
                    v(ok,:) = sum(dir.*qvec,2) ./ det; 
                    t(ok,:) = sum(edge2.*qvec,2)./det;
                    ok = (ok & v>=zero & u+v<=1.0+zero);
                    int2 = (ok & t>=zero & t<=1.0+zero);
                    intersect = intersect + int2;
                end
            end 

            
            edge1 = cd.c.e1r-cd.c.hr;          % find vectors for two edges sharing vert0
            edge2 = cd.c.e3r-cd.c.hr;
            tvec  = cd.d.r(n-1,:) - cd.c.hr;          % vector from vert0 to ray origin
            pvec  = cross(dir,edge2);  % begin calculating determinant - also used to calculate U parameter
            det   = sum(edge1.*pvec);   % determinant of the matrix M = dot(edge1,pvec)
            angleOK = (abs(det)>eps);
            if all(~angleOK)
                int3 = 0;
                intersect = intersect + int3;
            else
                det(~angleOK) = nan;              % change to avoid division by zero
                u    = sum(tvec.*pvec)./det;    % 1st barycentric coordinate
                v = nan+zeros(size(u));
                t=v;
                ok = (angleOK & u>=zero & u<=1.0+zero); % mask
                if ~any(ok)
                    int3 = ok;
                    intersect = intersect + int3;
                else
                    qvec = cross(tvec, edge1,2); % prepare to test V parameter
                    v(ok,:) = sum(dir.*qvec,2) ./ det; % 2nd barycentric coordinate
                    t(ok,:) = sum(edge2.*qvec,2)./det;
                    ok = (ok & v>=zero & u+v<=1.0+zero);
                    int3 = (ok & t>=zero & t<=1.0+zero);
                    intersect = intersect + int3;
                end
            end 
        end
        function [intersect,pos] = CohesinDNAIntersection(cd)
            % intersect is number of intersections between cohesin and DNA
            % pos calculates excess of DNA segments from one side of
            % cohesin
            intersect = 0;
            pos = 0;
            for i=2:cd.d.N
                [a,b] = cd.CohesinDNASegmentIntersection(i);
                intersect = intersect + a;
                pos = pos + b;
%                 if a~=0
%                     disp(i)
%                 end
            end
        end
        function [cd,varsave] = ChangeItem(cd,i1,i2,DD,F)
                % this function picks random bead from cohesin or DNA and
                % changes its position or orientation randomly for MMK
                % simulation
                % rn(1:5) -> change cohesin
                % rn(6,...) -> change DNA
                if i1>5
                    % Change one of DNA beads
%                     DD=0;
%                     DU=zeros(3,3);
                    if i2==1
                        % change position of bead rn in this case
                        varsave = cd.d.r(i1-5,:);
                        cd.d.r(i1-5,:) = cd.d.r(i1-5,:) + DD;
                    else
                        % change orientation of bead rn in this case
                        [F,DU] = F.RandomSample;
                        varsave = cd.d.u(i1-5,:);
                        cd.d.u(i1-5,:) = (DU*(cd.d.u(i1-5,:)'))';
                    end
                else
                    % Change one of cohesin beads
                    [F,DU] = F.RandomSample;
                    [cd.c,varsave] = cd.c.ChangeItem(i1,i2,DD,DU);
                end
        end
        function cd = ChangeItemBack(cd,i1,i2,varsave)
            % If new energy is not accepted this function returns variables
            % that were changed into its original value
                if i1>5
                    if i2==1
                        cd.d.r(i1-5,:) = varsave;
                    else
                        cd.d.u(i1-5,:) = varsave;
                    end
                else
                    cd.c = cd.c.ChangeItemBack(i1,i2,varsave);
                end
        end

        function [cd,EE] = MMK(cd,Niter,F)
            % Metropolis Monte-Carlo simulation of Niter iterations
            % F - handle to random rotation matrices
            kT = 1.1;
            % choose whether to change position or orientation
            ri = randi(2,Niter,1)-1;
            % Chose element
            rn = randi(cd.d.N+5,Niter,1);
            % random addition to position
            rv = 0.025*(rand(Niter,3)-0.5);
            % probability
            p = rand(Niter,1);

            EE = zeros(Niter,1);

            [cd,E0c,E0d,E0cd] = cd.energy;
            segment_parity0 = cd.CohesinDNAIntersection;
%             E0c = cd.c.energy;
%             cd.d.Energy;
%             E0d = sum(cd.d.E);
            E0dFull = cd.d.E;
            
            ratp = rand(Niter,1);
            rhingebind = rand(Niter,1);
            rhingefold = rand(Niter,1);
            
            for i=1:Niter
                % rn(1:5) -> change cohesin
                % rn(6,...) -> change DNA
                % ChangeItem randomly changes position or orientation of one
                % randomly chosed bead
                [cd,varsave] = cd.ChangeItem(rn(i),ri(i),rv(i,:),F);
                % Calculate new energy
                [cd,Ec,Ed,Ecd] = cd.ReEnergy(rn(i),E0c,E0d);
                % Check if number of Cohesin-DNA interactions is even
                % make istest=1 and allow change in energy only if either
                % parity did not change, or if it changed, then new number
                % is even and dna has large asymmetry
%                 [segment_parity,b] = cd.CohesinDNAIntersection;
%                 if segment_parity~=segment_parity0
%                     if (~mod(segment_parity,2))
% %                     if (~mod(segment_parity,2))&&(abs(b)>2*cd.d.N/3)
%                         istest = 1;
%                     else
%                         istest = 0;
%                     end
%                 else
%                     istest=1;
%                 end
%                 istest = ~mod(segment_parity,2);
                % Perform Metropolis Monte-Carlo
                istest = 1;
                if ((Ec+Ed+Ecd)<(E0c+E0d+E0cd))&&(istest)
                    E0c = Ec;
                    E0d = Ed;
                    E0cd = Ecd;
                    E0dFull = cd.d.E;
%                     segment_parity0=segment_parity;
                    update = 1;
                else
                    ro = exp(-((Ec+Ed+Ecd)-(E0c+E0d+E0cd))/kT);
                    if (ro>p(i))&&(istest)
                        E0c = Ec;
                        E0d = Ed;
                        E0cd = Ecd;
                        E0dFull = cd.d.E;
%                         segment_parity0=segment_parity;
                        update = 1;
                    else
                        cd = cd.ChangeItemBack(rn(i),ri(i),varsave);
                        cd.d.E = E0dFull;
                        update = 0;
                    end
                end
                need_update = 0;
%                 EE(i) = E0c+E0d+E0cd;
                % If molecule is in ATP state calculate probability that it
                % will change to ADP
%                 if cd.c_state
%                     pr = -log(1-ratp(i))/cd.K_release;
%                     if pr<1
%                         cd.c_state = 0;
%                         cd.c.ATPstate = 0;
%                         cd.smc3index = 0;
%                         % Change equilibrium distance between NBDs
%                         cd.c.HHT = 14/cd.c.Lp;
%                         % need to recalcualte E0 because it changes when
%                         % ATP state changes
%                         need_update = 1;    % energy needs to be updated before next iteration
%                     end
%                 else
%                     pb = -log(1-ratp(i))/cd.K_bind;
%                     if pb<1
%                         cd.c_state = 1;
%                         cd.c.ATPstate = 1;
%                         % Find DNA segment closest to Smc3 head and bind
%                         % them together
%                         [cd,R,pos] = cd.dist;
%                         cd.smc3index = pos;
%                         % Change equilibrium distance between NBDs
%                         cd.c.HHT = 4/cd.c.Lp;
%                         % Recalculate energy
%                         need_update = 1;
%                     end
%                 end

                % This model assumes that flipping of the hinge and its
                % binding to DNA are independent
%                 if cd.c_state == 0
%                     if cd.c.HingeState == 1
%                         pr = -log(1-rhingefold(i))/cd.K_hinge_bend;
%                         if pr<1
%                             cd.c.HingeState = -1;
%                             need_update = 1;
%                         end
%                     else
%                         pr = -log(1-rhingefold(i))/cd.K_hinge_unbend;
%                         if pr<1
%                             cd.c.HingeState = 1;
%                             need_update = 1;
%                         end
%                     end
%                 end
%                 if (cd.c_state)&&(cd.c.HingeState==-1)
%                     % if ATP is bound, then bring Hinge to straight state
%                     cd.c.HingeState = 1;
%                     need_update = 1;
%                 end
%                 
                % Binding of the hinge to DNA depends only on the distance
                % to DNA and does not depend on the hinge state
%                 if cd.hindex == 0
%                     % hinge is not bound
%                     % look whether it can bind DNA
%                     [cd,R,pos]=cd.hingepos;
%                     pr = -log(1-rhingebind(i))/cd.K_h_bind;
%                     R = sqrt(sum(R.*R));
%                     if (pr<1)&&(R<2*cd.d.D/3)
%                         cd.hindex = pos;
%                         need_update = 1;    % energy needs to be updated before next iteration
%                     end
%                 else
%                     pr = -log(1-rhingebind(i))/cd.K_h_release;
%                     if pr<1
%                         cd.hindex = 0;
%                         need_update = 1;   
%                     end
%                 end
                if need_update
                    % if configuration is changed due to chemical reactions
                    % (ATP binding/release, hinge folding), baseline energy
                    % needs to be recalculated before next iteration
                    [cd,E0c,E0d,E0cd] = cd.energy;
                    E0dFull = cd.d.E;
                end
%                 EE(i) = E0c+E0d+E0cd;
%                 EE(i) = b;
            end
        end
        function cd = Copy(cd,T)
            % create a copy of the current state
            cd.d.r = T.d.r;
            cd.d.u = T.d.u;
            cd.c = T.c;
            
            cd.c_state = T.c_state;
            cd.int = T.int;
            cd.s3n = T.s3n;
            cd.SMC3 = T.SMC3;
            cd.SMC1 = T.SMC1;
            cd.Hinge = T.Hinge;
            cd.SMC3_Slide = T.SMC3_Slide;
            cd.SMC1_Slide = T.SMC1_Slide;
            cd.smc3index = T.smc3index;
            cd.smc1index = T.smc1index;
            cd.hindex = T.hindex;
            
            cd.K_bind = T.K_bind;
            cd.K_release = T.K_release;
            cd.K_h_bind = T.K_h_bind;
            cd.K_h_release = T.K_h_release;

            cd.K_hinge_bend = T.K_hinge_bend;    
            cd.K_hinge_unbend = T.K_hinge_unbend;

            cd.Rep = T.Rep;
            cd.RepKleisin = T.RepKleisin;

        end
        function distance = DistBetweenSegments(cd,n,k)
            % This function calculates shortest distance between DNA
            % segment n and cohesin segment k:
            %           o
            %          / \  
            %    k=2  /   \  k=3
            %        /     \
            %        |     |
            %   k=1  |     |  k=4
            %        |     |
            %        O-----O
            %      Smc3    Smc1
            %          k=5
            % Algorythm is from DistBetween2Segments.m
            % http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
            
            u = cd.d.r(n,:) - cd.d.r(n-1,:);
            if k==1
                v = cd.c.e3r - cd.c.s3r;
                w = cd.d.r(n-1,:) - cd.c.s3r;
            elseif k==2
                v = cd.c.hr - cd.c.e3r;
                w = cd.d.r(n-1,:) - cd.c.e3r;
            elseif k==3
                v = cd.c.hr - cd.c.e1r;
                w = cd.d.r(n-1,:) - cd.c.e1r;
            elseif k==4
                v = cd.c.e1r - cd.c.s1r;
                w = cd.d.r(n-1,:) - cd.c.s1r;
            elseif k==5
                v = cd.c.s1r - cd.c.s3r;
                w = cd.d.r(n-1,:) - cd.c.s3r;
            end
            a = dot(u,u);
            b = dot(u,v);
            c = dot(v,v);
            d = dot(u,w);
            e = dot(v,w);
            D = a*c - b*b;
            sD = D;
            tD = D;
            SMALL_NUM = 0.00000001;
            % compute the line parameters of the two closest points
            if (D < SMALL_NUM)  % the lines are almost parallel
                sN = 0.0;       % force using point P0 on segment S1
                sD = 1.0;       % to prevent possible division by 0.0 later
                tN = e;
                tD = c;
            else                % get the closest points on the infinite lines
                sN = (b*e - c*d);
                tN = (a*e - b*d);
                if (sN < 0)   % sc < 0 => the s=0 edge is visible       
                    sN = 0.0;
                    tN = e;
                    tD = c;
                elseif (sN > sD)% sc > 1 => the s=1 edge is visible
                    sN = sD;
                    tN = e + b;
                    tD = c;
                end
            end

            if (tN < 0.0)            % tc < 0 => the t=0 edge is visible
                tN = 0.0;
                % recompute sc for this edge
                if (-d < 0.0)
                    sN = 0.0;
                elseif (-d > a)
                    sN = sD;
                else
                    sN = -d;
                    sD = a;
                end
            elseif (tN > tD)       % tc > 1 => the t=1 edge is visible
                tN = tD;
                % recompute sc for this edge
                if ((-d + b) < 0.0)
                    sN = 0;
                elseif ((-d + b) > a)
                    sN = sD;
                else 
                    sN = (-d + b);
                    sD = a;
                end
            end

            % finally do the division to get sc and tc
            if(abs(sN) < SMALL_NUM)
                sc = 0.0;
            else
                sc = sN / sD;
            end

            if(abs(tN) < SMALL_NUM)
                tc = 0.0;
            else
                tc = tN / tD;
            end

            dP = w + (sc * u) - (tc * v);  
            distance = norm(dP);
        
        end
    end
end
