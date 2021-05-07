classdef Loop5
    properties
        N    % Number of segments. Assume each 5 nm
        DS1 = 10;  % Diffusion coefficinet at Smc1 subunit
        DS3  % Diffusion coefficinet at Smc1 subunit
        HU  % Constant of hinge unfolding
        HD  % Constant of Hinge folding
        CB  % Constant of hinge binding DNA
        CU  % Constant of Hinge unbinding DNA
        HB = 0;  % Boolean reflecting whether hinge is bound (HB=1) or not (HB=0)
        HP = 0;  % Boolean definidng the hinge position at the Smc3 (HP=0) or far from it (HP=1)
        STEP = 6;  % Pulling step in number of 5 nm segments
        F
    end
    methods
        function L = Loop5(n)
            L.N = n;
            L.F = 0;
        end
        function [L,TS,NS,Hfolded,Hbound] = MMK(L,Niter)
            % Define constants
            global KD KD2 HFold Hunfold Force KHUglobal
            
            Kfactor = exp(-L.STEP*5*Force/4.14);

            
            K1U = KD2;  % free diffusion of DNA at Smc1 upwards . (0.5e6 -> 10 um2/s; )
            K1D = KD2;  % free diffusion of DNA at Smc1 downwards
            K3U = KD;  % diffusion of DNA at Smc3 upwards. 
            K3D = KD;  % diffusion of DNA at Smc3 downwards. (0.2e4 -> 0.05um2/s as in Koshland_Greene_2016. correspondence test >>test_rates.m)

            KHB = 1e10;   % Hinge binding event  (assume immediately in bend state and never in state far). 
            KHU = 1e-5;   
            KHU = KHUglobal; % Hinge unbinding event (only matters for the unfolded cohesin)
            
            HFLD = HFold;   % Hinge folding 
            HUNFLD = Hunfold;  % Hinge unfolding
            
            FR = 1;  % fraction to save
            
            TS = zeros(Niter/FR,1);
            NS = zeros(Niter/FR,1);
            Hfolded = zeros(Niter/FR,1);
            Hbound = zeros(Niter/FR,1);
            R = rand(Niter,6);
            
            % Do iterations
            Time = 0;
            for i=1:Niter
                if L.HB==0
                    % hinge is unbound
                    if L.HP==0
                        % hinge is unbound and at Smc3
                        % Diffusion is only next to unbound head
                        t1 = -log(1-R(i,1))/K1U; 
                        t2 = -log(1-R(i,2))/K1D; 
%                         t3 = -log(1-R(i,3))/K3U; 
%                         t4 = -log(1-R(i,4))/K3D; 
                        t5 = -log(1-R(i,5))/KHB; 
                        t6 = -log(1-R(i,6))/HUNFLD;
%                         t = [t1 t2 t3 t4 t5 t6];
                        t = [t1 t2 t5 t6];
                        [q,w]=min(t);
                        if w==1
                            L.N = L.N + 1;
                        elseif w==2
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==3
                            L.HB = 1;
                        elseif w==4
                            L.HP = 1;
                        end
                        Time = Time + q;
                    else
                        % hinge is unbound and far
                        t1 = -log(1-R(i,1))/K1U; 
                        t2 = -log(1-R(i,2))/K1D; 
                        t3 = -log(1-R(i,3))/K3U; 
                        t4 = -log(1-R(i,4))/K3D; 
                        t5 = -log(1-R(i,5))/(1e-20*KHB);   % assume never binds in this state
                        t6 = -log(1-R(i,6))/HFLD;
                        t = [t1 t2 t3 t4 t5 t6];
                        [q,w]=min(t);
                        if w==1
                            L.N = L.N + 1;
                        elseif w==2
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==3
                            L.N = L.N + 1;
                        elseif w==4
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==5
                            L.HB = 1;
                        elseif w==6
                            L.HP = 0;
                        end
                        Time = Time + q;
                    end
                else
                    % hinge is bound
                    if L.HP==0
                        % hinge is bound and at Smc3
                        t1 = -log(1-R(i,1))/K1U; 
                        t2 = -log(1-R(i,2))/K1D; 
                        t3 = -log(1-R(i,3))/KHU; 
                        t4 = -log(1-R(i,4))/(HUNFLD*Kfactor);
                        t = [t1 t2 t3 t4];
                        [q,w]=min(t);
                        if w==1
                            L.N = L.N + 1;
                        elseif w==2
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==3
                            L.HB = 0;
                        elseif w==4
                            L.HP = 1;
                            L.N = L.N + L.STEP;  % This size of DNA pulled by the hinge needs to be determined from 3D model
                        end
                        Time = Time + q;
                    else
                        % hinge is bound and far
                        t1 = -log(1-R(i,1))/K1U; 
                        t2 = -log(1-R(i,2))/K1D; 
                        t3 = -log(1-R(i,3))/K3U; 
                        t4 = -log(1-R(i,4))/K3D; 
                        t5 = -log(1-R(i,5))/KHU; 
                        t6 = -log(1-R(i,6))/HFLD;
                        t = [t1 t2 t3 t4 t5 t6];
                        [q,w]=min(t);
                        if w==1
                            L.N = L.N + 1;
                        elseif w==2
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==3
                            L.N = L.N + 1;
                        elseif w==4
                            if L.N>100
                                L.N = L.N - 1;
                            end
                        elseif w==5
                            L.HB = 0;
                        elseif w==6
                            % This assumes it cannot go and push back DNA
                            % because when ATP binds cohesin goes into
                            % gripping state. Therefore DNA needs to unbind
                            % first
                            L.HP = 0;
                            L.HB = 0;
%                             L.N = L.N - L.STEP;  % This size of DNA pulled by the hinge needs to be determined from 3D model
                        end
                        Time = Time + q;
                    end
                end
                NS(ceil(i/FR)) = L.N;
                Hfolded(ceil(i/FR)) = L.HP;
                Hbound(ceil(i/FR)) = L.HB;
                TS(ceil(i/FR)) = Time;
            end
        end
    end
end
