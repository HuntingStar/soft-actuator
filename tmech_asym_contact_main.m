%% tmech_asym_contact_main.m
% IEEE TMECH-style asymmetric adjacent-chamber contact model
% Units: mm, N, MPa (N/mm^2)
clear; clc; close all;

outDir = 'outputs'; if ~exist(outDir,'dir'), mkdir(outDir); end
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

par = defaultParameters();
p_sweep_kPa = linspace(1,120,180);
res = evalOnePressure(par, par.p_plot_kPa);
sweep = evalPressureSweep(par, p_sweep_kPa);

% parameter studies
studyDH = plotParametricStudy(par,'DeltaH',linspace(0,4,25));
studyHx = plotParametricStudy(par,'hx',linspace(0.3,1.5,25));
studyTm = plotParametricStudy(par,'tm',linspace(0.8,2.2,20));
studyMu = plotParametricStudy(par,'mu',linspace(0.08,0.30,20));

f1 = plotContactGeometry(par,res); saveFigurePair(f1,fullfile(outDir,'Fig1_contact_geometry'));
f2 = plotVariableSweep(sweep); saveFigurePair(f2,fullfile(outDir,'Fig2_pressure_sweep'));
f3 = plotStudyDeltaH(studyDH); saveFigurePair(f3,fullfile(outDir,'Fig3_height_difference'));
f4 = plotStudyHx(studyHx); saveFigurePair(f4,fullfile(outDir,'Fig4_gap_effect'));
f5 = plotEquivalentBending(par,sweep); saveFigurePair(f5,fullfile(outDir,'Fig5_bending_contours'));

save(fullfile(outDir,'results_all.mat'),'par','res','sweep','studyDH','studyHx','studyTm','studyMu');
writetable(struct2table(flattenSweep(sweep)),fullfile(outDir,'pressure_sweep.csv'));

fprintf('Done. Outputs saved to %s\n', outDir);

function par = defaultParameters()
par.a1=6; par.a2=5; par.b1=8; par.b2=8; par.hx=0.75;
par.tm1=1.5; par.tm2=1.5; par.mu1=0.198; par.mu2=0.198;
par.H=5; par.B=16; par.L=20; par.hL=4.2; par.Ec=0.408;
par.p_plot_kPa=40; par.Nbeam=300; par.rootTol=1e-9;
end

function theta_mf = solveThetaMF(p,a,mu,tm)
if p<=0, theta_mf=0; return; end
rhs = p*a/(2*mu*tm); f=@(t) sin(t)-sin(t).^7./t.^6-rhs;
g=linspace(1e-6,pi-1e-6,5000); v=f(g); k=find(v(1:end-1).*v(2:end)<=0,1,'first');
if isempty(k), theta_mf=NaN; else, theta_mf=fzero(f,[g(k),g(k+1)]); end
end

function u = uMem(Y,a,Rf,theta_mf)
inside=max(Rf.^2-(Y-a).^2,0); u=sqrt(inside)-Rf*cos(theta_mf);
end
function du = duMem(Y,a,Rf)
den=sqrt(max(Rf.^2-(Y-a).^2,1e-12)); du=-(Y-a)./den; end

function res=evalOnePressure(par,p_kPa)
p=p_kPa*1e-3; res=struct('p',p,'p_kPa',p_kPa);
res.theta_mf1=solveThetaMF(p,par.a1,par.mu1,par.tm1); res.theta_mf2=solveThetaMF(p,par.a2,par.mu2,par.tm2);
if any(isnan([res.theta_mf1,res.theta_mf2])), res=fillNoContact(res,par); return; end
res.Rf1=par.a1/sin(res.theta_mf1); res.Rf2=par.a2/sin(res.theta_mf2);
res.Ystar=(res.Rf2*par.a1+res.Rf1*par.a2)/(res.Rf1+res.Rf2);
res.rho1_star=res.Ystar-par.a1; res.rho2_star=res.Ystar-par.a2;
arg=(par.a1-par.a2)/(res.Rf1+res.Rf2); arg=max(min(arg,1),-1); res.theta_h_initial=asin(arg);
D=@(Y) uMem(Y,par.a1,res.Rf1,res.theta_mf1)+uMem(Y,par.a2,res.Rf2,res.theta_mf2)-2*par.hx;
res.Phi=D(res.Ystar); res.hasContact=res.Phi>0;
if ~res.hasContact, res=fillNoContact(res,par); return; end
Ymin=0; Ymax=min(2*par.a1,2*par.a2); Yg=linspace(Ymin,Ymax,2400); Dg=D(Yg); [~,im]=max(Dg);
li=find(Dg(1:im)<=0,1,'last'); ri=im-1+find(Dg(im:end)<=0,1,'first');
res.Yminus=iff(isempty(li),Ymin,fzero(D,[Yg(li),Yg(li+1)]));
res.Yplus=iff(isempty(ri)||ri<=im,Ymax,fzero(D,[Yg(ri-1),Yg(ri)]));
res.c1=max(res.Ystar-res.Yminus,0); res.c2=max(res.Yplus-res.Ystar,0); res.cy=max(res.Yplus-res.Yminus,0);
cz1=res.cy*par.b1*max(par.a1-abs(res.rho1_star),0)/max(par.a1^2,eps);
cz2=res.cy*par.b2*max(par.a2-abs(res.rho2_star),0)/max(par.a2^2,eps);
res.cz=max(min(cz1,cz2),0); res.Ac=pi/4*max(res.cy,0)*max(res.cz,0);
YcGrid=linspace(res.Yminus,res.Yplus,800); Dp=max(D(YcGrid),0); Aint=trapz(YcGrid,Dp);
res.Yc=iff(Aint<1e-12,res.Ystar,trapz(YcGrid,YcGrid.*Dp)/Aint);
sc=0.5*(duMem(res.Yc,par.a1,res.Rf1)-duMem(res.Yc,par.a2,res.Rf2)); res.theta_h=atan(sc);
res.Fp=p*res.Ac; res.Fx=res.Fp*cos(res.theta_h); res.FY=-res.Fp*sin(res.theta_h);
rx=par.hx; rY=par.H/2+res.Yc; res.Mm=abs(rx*res.FY-rY*res.Fx);
res.Mpi=res.Mm*sqrt(max(2*par.hx/(2*par.hx+par.hL),0)); Ic=par.B*par.H^3/12; res.kappa=res.Mpi/(par.Ec*Ic); res.thetaL=res.kappa*par.L;
[res.beamX,res.beamY]=beamCurve(par.L,res.kappa,par.Nbeam);
end

function out=evalPressureSweep(par,p_sweep_kPa)
for i=1:numel(p_sweep_kPa), rr(i)=evalOnePressure(par,p_sweep_kPa(i)); end %#ok<AGROW>
out.res=rr; out.p_kPa=p_sweep_kPa(:);
end
function [x,y]=beamCurve(L,k,N), s=linspace(0,L,N); if abs(k)<1e-12, x=s; y=zeros(size(s)); else, x=sin(k*s)/k; y=(1-cos(k*s))/k; end,end
function r=fillNoContact(r,par)
r.hasContact=false; r.Yminus=NaN; r.Yplus=NaN; r.Yc=r.Ystar; r.c1=0; r.c2=0; r.cy=0; r.cz=0; r.Ac=0; r.Fp=0; r.Fx=0; r.FY=0; r.theta_h=r.theta_h_initial; r.Mm=0; r.Mpi=0; r.kappa=0; r.thetaL=0; [r.beamX,r.beamY]=beamCurve(par.L,0,par.Nbeam); end
function v=iff(c,a,b), if c, v=a; else, v=b; end, end

function f=plotContactGeometry(par,res), f=figure('Color','w'); hold on; box on; axis equal; Y1=linspace(0,2*par.a1,600); Y2=linspace(0,2*par.a2,600);
plot(-par.hx*ones(size(Y1)),Y1,'k:'); plot(par.hx*ones(size(Y2)),Y2,'k:');
plot(-par.hx+uMem(Y1,par.a1,res.Rf1,res.theta_mf1),Y1,'b--'); plot(par.hx-uMem(Y2,par.a2,res.Rf2,res.theta_mf2),Y2,'r--');
title('Fig.1 Asymmetric chambers'); xlabel('x (mm)'); ylabel('Y (mm)'); end
function f=plotVariableSweep(sw), f=figure('Color','w'); t=flattenSweep(sw); subplot(2,2,1); plot(t.p_kPa,t.theta_h*180/pi); xlabel('p (kPa)'); ylabel('\theta_h (deg)'); subplot(2,2,2); plot(t.p_kPa,t.Ac); ylabel('A_c (mm^2)'); xlabel('p (kPa)'); subplot(2,2,3); plot(t.p_kPa,[t.Fp t.Fx t.FY]); legend('Fp','Fx','FY'); xlabel('p (kPa)'); subplot(2,2,4); plot(t.p_kPa,[t.Mm t.Mpi t.thetaL*180/pi]); legend('Mm','Mpi','thetaL'); xlabel('p (kPa)'); end
function f=plotEquivalentBending(par,sw), f=figure('Color','w'); hold on; idx=round(linspace(1,numel(sw.res),5)); for i=idx, plot(sw.res(i).beamX,sw.res(i).beamY,'LineWidth',1.4); end; xlabel('x (mm)'); ylabel('y (mm)'); title('Fig.5 equivalent bending'); legend(compose('p=%.0f kPa',sw.p_kPa(idx))); end
function st=plotParametricStudy(par,typ,vals)
for i=1:numel(vals), pp=par; switch typ, case 'DeltaH', pp.a2=pp.a1-vals(i)/2; case 'hx', pp.hx=vals(i); case 'tm', pp.tm1=vals(i); pp.tm2=vals(i); case 'mu', pp.mu1=vals(i); pp.mu2=vals(i); end, rr=evalOnePressure(pp,pp.p_plot_kPa); st.res(i)=rr; st.x(i)=vals(i); end
st.type=typ; end
function f=plotStudyDeltaH(st), f=figure('Color','w'); tt=[st.res]; subplot(2,2,1); plot(st.x,[tt.Ystar]); ylabel('Y*'); subplot(2,2,2); plot(st.x,[tt.theta_h]*180/pi); ylabel('theta_h (deg)'); subplot(2,2,3); plot(st.x,[tt.FY]./[tt.Fx]); ylabel('FY/Fx'); subplot(2,2,4); plot(st.x,[tt.thetaL]*180/pi); ylabel('thetaL (deg)'); xlabel('\DeltaH (mm)'); end
function f=plotStudyHx(st), f=figure('Color','w'); tt=[st.res]; subplot(2,2,1); plot(st.x,[tt.p]*1e3); ylabel('p (kPa)'); subplot(2,2,2); plot(st.x,[tt.Ac]); ylabel('Ac'); subplot(2,2,3); plot(st.x,[tt.Fp]); ylabel('Fp'); subplot(2,2,4); plot(st.x,[tt.Mpi]); ylabel('Mpi'); xlabel('h_x (mm)'); end
function t=flattenSweep(sw), rr=sw.res; fn={'p','p_kPa','theta_mf1','theta_mf2','Ystar','theta_h','Ac','Fp','Fx','FY','cy','cz','Mm','Mpi','kappa','thetaL'}; for k=1:numel(fn), t.(fn{k})=arrayfun(@(s)getfield(s,fn{k}),rr).'; end end
function saveFigurePair(f,base), exportgraphics(f,[base,'.png'],'Resolution',400); exportgraphics(f,[base,'.pdf']); end
