%% Facia Lata Model
% This code is associated with the article In vivo measurements of fascia lata effective mechanics combined to a memory
% fiber–recruitment–viscoelastic modeling approach" authored by F. Germain
% and T. Gibaud
% we have provided measurement to test the code



%function ramp_relax_main()

clc; clear all; close all;


%% ── PARAMETERS ──
% choose the ramp speed mm/s
P.v = 8;
% choose recrutement parameters 
% Left leg:   L0=55mm, sigma=0.40
% right leg:  L0=42mm, sigma=0.38
P.L0 = 42%55%42%55;       % recrutement parameters
P.sigma = .38%.4%.38%0.4; % recrutement parameters

% initialize the fit parameters
P.Ft = 150;
P.k_m = 0.22%0.15%0.22;
P.k_f = 1.7;
P.k1 = 0.5;
P.tau1 = 2.5;
P.k2 = 17%9%17;
P.tau2 =100% 120;

P0 = P; % store initial parameters
% constrain fit parameters if necessary
lb=[P0.Ft P0.k_m P.k_f P0.L0 P0.sigma P0.k1 P0.tau1 P0.k2 P0.tau2]*.1;
ub=[P0.Ft P0.k_m P.k_f P0.L0 P0.sigma P0.k1 P0.tau1 P0.k2 P0.tau2]*10;

% optimized recruitment parameters, fully constrained
lb(4)=P0.L0;ub(4)=P0.L0;
lb(5)=P0.sigma;ub(5)=P0.sigma;

% rop paramaters from calibration
KR = 35;
ALPHA = 0.031;
T_R0 = 0.16;
T_RELAX = 1000;

% ── LOAD DATA ──
[file,path] = uigetfile('*.txt','Select data file');
    if isequal(file,0)
        disp('No file selected'); return;
    end
[expRamp, expRelax] = parseFile(fullfile(path,file));

P.Ft=max(expRamp.F);
lb(1)=P.Ft;ub(1)=P.Ft;

%% Fit the experimental data
    % ── FIT (bounded, no toolbox) ──
    %tic
    names = {'Ft','k_m','k_f','L0','sigma','k1','tau1','k2','tau2'};
    
    
    x0 = cellfun(@(k) P.(k), names);
    z0 = invTransform(x0, lb, ub);
    
    costfun = @(z) costBounded(z, names, P, expRamp, expRelax, ...
                              KR, ALPHA, T_R0, T_RELAX, lb, ub);
    
    opts = optimset('Display','off','MaxIter',1000);
    zopt = fminsearch(costfun, z0, opts);
    
    xopt = transform(zopt, lb, ub);
    
    for i=1:length(names)
        P.(names{i}) = xopt(i);
    end
    
    % ── FINAL COST (MSR) ──
    finalCost = costFunction(xopt, names, P0, expRamp, expRelax, ...
                             KR, ALPHA, T_R0, T_RELAX);
    
    % ── DISPLAY RESULTS ──
    fprintf('\n===== FIT RESULTS =====\n');
    for i = 1:length(names)
        fprintf('%8s : %10.5f   (initial: %10.5f)\n', ...
            names{i}, P.(names{i}), P0.(names{i}));
    end
    fprintf('------------------------\n');
    fprintf('MSR = %.6f\n', finalCost);
    fprintf('========================\n\n');
    
    % ── SIMULATION ──
    ramp = simRamp(P);
    fgrid = buildFeqGrid(P, ramp.Lmax);
    % relax = simRelax(ramp, P, fgrid, KR, ALPHA, T_R0, T_RELAX);
    relax = simRelax(ramp, P, fgrid, KR, ALPHA, T_R0, T_RELAX, false);
    % ── PLOTS ──
    plotAll(ramp, relax, expRamp, expRelax, P);
    
    % ── EXPORT ──
    ii=1
    % exportResults(ramp, relax, P);
    ResultsFit(ii,:)=[P.Ft P.k_m P.k_f P.L0 P.sigma P.k1 P.tau1 P.k2 P.tau2 finalCost P.v];
    disp('Done.')
%toc


%% ================= FUNCTIONS =================

function [ramp, relax] = parseFile(filename)
fid = fopen(filename);
rL=[]; rF=[]; t=[]; f=[];
mode='ramp';

while ~feof(fid)
    line=strtrim(fgetl(fid));
    if isempty(line), continue; end

    if startsWith(line,'%')
        if contains(lower(line),'relax')
            mode='relax';
        end
        continue;
    end

    vals=sscanf(line,'%f');
    if numel(vals)<2, continue; end

    if strcmp(mode,'ramp')
        rL(end+1)=vals(1); rF(end+1)=vals(2);
    else
        t(end+1)=vals(1); f(end+1)=vals(2);
    end
end
fclose(fid);

if isempty(rL), ramp=[]; else ramp.L=rL; ramp.F=rF; end
if isempty(t), relax=[]; else relax.t=t; relax.F=f; end
end

function c = costBounded(z, names, P, expRamp, expRelax, KR, ALPHA, T_R0, T_RELAX, lb, ub)
x = transform(z, lb, ub);
c = costFunction(x, names, P, expRamp, expRelax, KR, ALPHA, T_R0, T_RELAX);
end

%%
function c = costFunction(x, names, P, expRamp, expRelax, KR, ALPHA, T_R0, T_RELAX)

% Update parameters
for i=1:length(names)
    P.(names{i}) = x(i);
end
try
    ramp = simRamp(P);
    fgrid = buildFeqGrid(P, ramp.Lmax);
catch
    c = 1e12;
    return;
end

sumErr = 0;
N = 0;

% ───── RAMP (vectorized) ─────
if ~isempty(expRamp)
    % Interpolate model at experimental L points
    Fm = interp1(ramp.pts(:,1), ramp.pts(:,5), ...
                 expRamp.L, 'linear', 'extrap');
    r = expRamp.F - Fm;
    sumErr = sumErr + sum(r.^2);
    N = N + numel(r);
end

% ───── RELAXATION (vectorized) ─────
if ~isempty(expRelax)
    try
        %relax = simRelax(ramp, P, fgrid, KR, ALPHA, T_R0, T_RELAX);
        relax = simRelax(ramp, P, fgrid, KR, ALPHA, T_R0, T_RELAX, true);
    catch
        c = 1e12;
        return;
    end
    % Interpolate model at experimental time points
    Fm = interp1(relax.pts(:,1), relax.pts(:,2), ...
                 expRelax.t, 'linear', 'extrap');
    r = expRelax.F - Fm;
    sumErr = sumErr + sum(r.^2);
    N = N + numel(r);
end

% ───── FINAL COST ─────
if N == 0
    c = 1e12;
else
    c = sumErr / N;   % Mean Squared Residual
end

end

%%
function ramp = simRamp(p)
dt=0.008; dL=p.v*dt;
e1=exp(-dt/p.tau1); e2=exp(-dt/p.tau2);

L=0.01; Feq=0; Fm1=0; Fm2=0; Fprev=0;
pts=[];

while (Feq+Fm1+Fm2)<p.Ft && L<800
    kA=kinst(L,p); kB=kinst(L+dL,p); G=gc(L,p);
    Feq=Feq+0.5*(kA+kB)*dL;
    Fm1=Fm1*e1+p.k1*G*L*(1-e1);
    Fm2=Fm2*e2+p.k2*G*L*(1-e2);

    Ftot=Feq+Fm1+Fm2;
    kapp=(Ftot-Fprev)/dL;

    Fprev=Ftot; L=L+dL;
    pts=[pts; L Feq Fm1 Fm2 Ftot kapp];
end

ramp.pts=pts;
ramp.Lmax=L;
ramp.Feq_end=Feq;
ramp.Fm1_end=Fm1;
ramp.Fm2_end=Fm2;
ramp.Fpeak=Feq+Fm1+Fm2;
end

%%
function relax = simRelax(ramp,p,fgrid,KR,ALPHA,T_R0,T_RELAX,fastMode)

Lmax = ramp.Lmax;
Feq_end = ramp.Feq_end;
Fm1_end = ramp.Fm1_end;
Fm2_end = ramp.Fm2_end;
Fpeak = ramp.Fpeak;

L_r0 = Fpeak / KR;

% ── number of points ──
if fastMode
    Npts = 40;    % fast fitting
else
    Npts = 350;   % smooth plotting
end

% ── log-spaced time vector ──
tvals = logspace(log10(0.01), log10(T_RELAX), Npts);

% preallocate (faster than growing array)
pts = zeros(Npts,5);

for i = 1:Npts

    t = tvals(i);

    % viscoelastic decay
    Fm1 = Fm1_end * exp(-t/p.tau1);
    Fm2 = Fm2_end * exp(-t/p.tau2);
    decF = Fm1 + Fm2;

    % rope creep
    rope = L_r0 * (1 - (T_R0/t)^ALPHA);

    % solve Feq iteratively
    Feq = Feq_end;
    for k = 1:10
        dLf = (Fpeak - (Feq + decF))/KR + rope;
        Fn = interpFeq(Lmax + dLf, fgrid);

        if abs(Fn - Feq) < 1e-7
            break;
        end

        Feq = 0.5 * (Feq + Fn);
    end

    % store results
    pts(i,:) = [t, Feq + decF, Feq, Fm1, Fm2];
end

relax.pts = pts;

end

%%
function g = buildFeqGrid(p,Lmax)
N=4000; Lext=max(Lmax*2.5,150); dL=Lext/N;
Fq=zeros(N+1,1);

for i=2:N+1
    La=(i-2)*dL; Lb=(i-1)*dL;
    Fq(i)=Fq(i-1)+0.5*(kinst(La,p)+kinst(Lb,p))*dL;
end

g.Fq=Fq; g.dL=dL; g.N=N; g.Lext=Lext;
end

function val = interpFeq(L,g)
if L<=0, val=0; return; end
if L>=g.Lext, val=g.Fq(end); return; end

idx=min(floor(L/g.dL)+1,g.N);
x0=(idx-1)*g.dL;
val=g.Fq(idx)+(L-x0)/g.dL*(g.Fq(idx+1)-g.Fq(idx));
end

function k = kinst(L,p)
k=p.k_m+p.k_f*gc(L,p);
end

function val = gc(L,p)
val=0.5*(1+erf((log(max(L,1e-6))-log(p.L0))/(p.sigma*sqrt(2))));
end

%%
function plotAll(ramp,relax,expRamp,expRelax,P)

figure(1);

% ── RAMP: F vs L ──
subplot(2,2,1); hold on;

% Model components
plot(ramp.pts(:,1), ramp.pts(:,5),'b','LineWidth',2)      % F total
plot(ramp.pts(:,1), ramp.pts(:,2),'--c','LineWidth',1.5)  % Finst (Feq)
plot(ramp.pts(:,1), ramp.pts(:,3),'--g','LineWidth',1.5)  % Fm1
plot(ramp.pts(:,1), ramp.pts(:,4),'--m','LineWidth',1.5)  % Fm2

% Experimental
if ~isempty(expRamp)
    scatter(expRamp.L,expRamp.F,20,'k','filled')
end

xlabel('L (mm)')
ylabel('F (N)')
title('Ramp: F vs L')
legend({'F total','F_{inst}','F_{m1}','F_{m2}','Exp'},'Location','northwest')
grid on


% ── RAMP: k vs L ──
subplot(2,2,3); hold on;

k = gradient(ramp.pts(:,5), ramp.pts(:,1));
plot(ramp.pts(:,1),k,'r','LineWidth',1.5)
if ~isempty(expRamp)
    scatter(expRamp.L(1:end-4-1),smooth(diff(expRamp.F(1:end-4))./diff(expRamp.L(1:end-4)),5),20,'k','filled')
end
xlabel('L (mm)')
ylabel('k (N/mm)')
title('k = dF/dL')
grid on


% ── RELAXATION: F vs t ──
subplot(2,2,2); hold on;

tlog = log10(relax.pts(:,1));

% Model components
plot(tlog, relax.pts(:,2),'b','LineWidth',2)      % total
plot(tlog, relax.pts(:,3),'--c','LineWidth',1.5)  % Feq
plot(tlog, relax.pts(:,4),'--g','LineWidth',1.5)  % Fm1
plot(tlog, relax.pts(:,5),'--m','LineWidth',1.5)  % Fm2

% Experimental
if ~isempty(expRelax)
    scatter(log10(expRelax.t),expRelax.F,20,'k','filled')
end

xlabel('log10(t)')
ylabel('F (N)')
title('Relaxation: F vs t')
legend({'F total','F_{inst}','F_{m1}','F_{m2}','Exp'},'Location','northeast')
grid on


% ── PARAMETER BOX ──
annotation('textbox',[0.65 0.05 0.3 0.3],...
    'String',sprintf(['Ft=%.1f\nk_m=%.2f\nk_f=%.2f\nL0=%.1f\nσ=%.2f\n',...
                      'k1=%.2f  τ1=%.1f\nk2=%.2f  τ2=%.1f'], ...
    P.Ft,P.k_m,P.k_f,P.L0,P.sigma,P.k1,P.tau1,P.k2,P.tau2),...
    'FitBoxToText','on');

end

%%

function exportResults(ramp,relax,P)

fid=fopen('result.txt','w');

fprintf(fid,'# RAMP\n');
for i=1:size(ramp.pts,1)
    fprintf(fid,'%f\t%f\n',ramp.pts(i,1),ramp.pts(i,5));
end

fprintf(fid,'# RELAX\n');
for i=1:size(relax.pts,1)
    fprintf(fid,'%f\t%f\n',relax.pts(i,1),relax.pts(i,2));
end

fclose(fid);
end

%%
function x = transform(z,lb,ub)
x = lb + (ub-lb)./(1+exp(-z));
end

%%
function z = invTransform(x,lb,ub)
x=(x-lb)./(ub-lb);
x=min(max(x,1e-6),1-1e-6);
z=log(x./(1-x));
end

