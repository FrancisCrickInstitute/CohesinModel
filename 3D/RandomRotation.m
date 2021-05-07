classdef RandomRotation < handle
% This class creates N random rotation vectors with angle phi and
% returns the hangle to this array
% You can pick next random vector using RandomSample function
% When all samples are used the function generates new set of N random
% vectors
    properties
        Niter % number of vectors to generate at one go
        phi   % random angle
        ic    % current number used
        F     % matrix containing vectors
        Rounds % number of times the matrix was reset
    end
    methods
        function R = RandomRotation(N,phi)
            R.Niter = N;
            R.phi = phi;
            R.ic = 1;
            
            rphi = phi*(rand(N,1)-0.5);
            rt = randn(N,3);
            R.F = zeros(3,3,N);
            for i=1:N
                % generate a random orientation vector
                rts = rt(i,:)./sqrt(sum(rt(i,:).*rt(i,:)));
                % rotate around it 
                M = makehgtform('axisrotate',rts,rphi(i));
                R.F(:,:,i) = M(1:3,1:3);
            end
            R.Rounds = 1;
        end
        function [R,DU] = RandomSample(R)
           if R.ic>R.Niter
                   % generate new matrix and reset
                rphi = R.phi*(rand(R.Niter,1)-0.5);
                rt = randn(R.Niter,3);
                for i=1:R.Niter
                    rts = rt(i,:)./sqrt(sum(rt(i,:).*rt(i,:)));
                    M = makehgtform('axisrotate',rts,rphi(i));
                    R.F(:,:,i) = M(1:3,1:3);
                end
                DU = R.F(:,:,1);
                R.ic = 2;
                R.Rounds = R.Rounds + 1;
           else
               DU = R.F(:,:,R.ic);
               R.ic = R.ic+1;
           end
        end
    end
end
