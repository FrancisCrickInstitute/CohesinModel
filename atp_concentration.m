% Selected ATP concentrations (mM)
atp=1e-3*[0.001 0.03 0.18 0.55 1];
% Same range with denser points for curve plotting
atpN = logspace(-6,-3,50);

x = [0;atp(end)];
y = [200;1.6];  % These are times between engagement at 0 and 1 mM ATP from experimental data
f = fit(x,y,'exp1');
t = f(atp);


% Number of points per each ATP concentration
N = 1500;

% Generate N exponentially distributed engagement times for each of the ATP
% concentrations 
times = zeros(N,length(t));
for i=1:length(t)
    q = exprnd(t(i),N,1);
    times(:,i) = q;
end

% Calculate mean and SEM
R = mean(times);
RS = std(times)/sqrt(N);

figure
hold on
errorbar(atp,R,RS,'o')
% Fit with single exponent
[f,g]=fit(atp',R','exp1');

plot(atpN,(f(atpN)))

