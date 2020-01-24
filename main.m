% Application of algorithms developed in IET Control to an 8th-order PEM Fuel Cell singularly perturbed linear system. Running main.m generates results.	

clear all

% Input Data
n=8; m=1; l=3; c=3;

A=[-6.30908 0 -10.9544 0 83.74458 0 0 24.05866;
    0 -161.083 0 0 51.52923 0 -18.0261 0;
    -18.7858 0 -46.3136 0 275.6592 0 0 158.3741;
    0 0 0 -17.3506 193.9373 0 0 0;
    1.299576 0 2.969317 0.3977 -38.7024 0.105748 0 0;
    16.64244 0 38.02522 5.066579 -479.384 0 0 0;
    0 -450.386 0 0 142.2084 0 -80.9472 0;
    2.02257 0 4.621237 0 0 0 0 -51.2108];

B=[0;0;0;3.946683;0;0;0;0];

C=[0 0 0 5.066579 -116.446 0 0 0;
    0 0 0 0 1 0 0 0;
    12.96889 10.32532 -0.56926 0 0 0 0 0];

D=[0;0;0];

% Transformation to put the system in explicit singularly perturbed form
I = eye(n);
P27 = swap_rows(I,2,7);
P36 = swap_rows(I,3,6);
P15 = swap_rows(I,1,5);
P24 = swap_rows(I,2,4);
P25 = swap_rows(I,2,5);

% Pij permutation matrices: switch i and j rows.
V=P25*P24*P27*P15*P36;
Asp=V*A*V';
Bsp=P25*P24*P27*P15*P36*B;
Csp=C*P36*P15*P27*P24*P25;
Dsp=D;
eps=0.157   % Calculated as ratio of eigenvalues with largest separation

A1=Asp(1:3,1:3);
A2=Asp(1:3,4:8);
A3=eps*Asp(4:8,1:3);
A4=eps*Asp(4:8,4:8);
det(A4);
A4inv=inv(A4);
B1=Bsp(1:3,1:1);
B2=eps*Bsp(4:8,1:1);
C1=Csp(1:3,1:3);
C2=Csp(1:3,4:8);
A0=A1-A2*A4inv*A3;
C0=C1-C2*A4inv*A3;
B0=B1-A2*A4inv*B2;
D0=Dsp-C2*A4inv*B2;

[Li,Hi] = calc_LH(A1,A2,A3,A4,eps);

% Slow fast matrices
As=A1-A2*Li;
Bs=B1-Hi*B2-eps*Hi*Li*B1;
Cs=C1-C2*Li;
Af=A4+eps*Li*A2;
Bf=B2+eps*Li*B1;
Cf=C2+eps*Cs*Hi;
%
% Chang Transformation
Tchang=[eye(3)-eps*Hi*Li -eps*Hi;Li eye(5)];
a=Tchang*Asp*inv(Tchang);
%
% Exact Slow-Fast Controllability Grammians, formula (52)
%
Wcs=lyap2(As,Bs*Bs');
Wcf=lyap2(Af,Bf*Bf');
Wcsf=lyap2(eps*As,Af',Bs*Bf');
Wc=inv(Tchang)*[Wcs Wcsf; Wcsf' Wcf/eps]*inv(Tchang'); % formula (54)
WcDirect=lyap2(Asp,Bsp*Bsp');
ERRcon=WcDirect-Wc;
%
% Exact Slow-Fast Observability Grammians, formula (53)
%
Wos=lyap2(As',Cs'*Cs);
Wof=lyap2(Af',Cf'*Cf);
Wosf=lyap2(eps*As',Af,Cs'*Cf);
WoSP=[Wos eps*Wosf; eps*Wosf' eps*Wof];
Wo=Tchang'*WoSP*Tchang; % formula (54)
WoDirect=lyap2(Asp',Csp'*Csp);
ERRobs=WoDirect-Wo;
lambdaobs=eig(WoDirect);
lambdaobs_sf=eig(Wo);
HSVslow=sqrt(eig(Wcs*Wos));
HSVfast=sqrt(eig(Wcf*Wof));
WconWobs=[Wcs*Wos+eps*Wcsf*Wosf' eps*(Wcs*Wosf+Wcsf*Wof);
                  Wcsf'*Wos+Wcf*Wosf' Wcf*Wof+eps*Wcsf'*Wosf];
sigma_i=sqrt(eig(WconWobs))

                  
%
% Conclusion: THE FAST SYSTEM IS STABILIZABLE, but  very weakly CONTROLLABLE
% 
% BALANCING
sysSP=ss(Asp,Bsp,Csp,Dsp);
Ds=zeros(3,1);
[sysb,Sigma,Tb,Tbinv]=balreal(sysSP);
Ab=Tb*Asp*Tbinv;
Bb=Tb*Bsp;
Cb=C*Tbinv;
Db=Dsp;
eig(Ab);
%
% Slow System Exact Balancing
%
sys_s=ss(As,Bs,Cs,Ds);
[sys_sb,Sigma_bs,Tbs,Tbinvs]=balreal(sys_s);
Abs=Tbs*As*Tbinvs;
Bbs=Tbs*Bs;
Cbs=Cs*Tbinvs;
Sigma_bs;
%
% Fast System Exact Balancing
%
Df=zeros(3,1);
sys_f=ss(Af/eps,Bf/eps,Cf,Df);
[sys_fb,Sigma_bf,Tbf,Tbinvf]=balreal(sys_f);
Abf=Tbf*Af*Tbinvf;
Bbf=Tbf*Bf;
Cbf=Cf*Tbinvf;
Sigma;
Sigma_bs;
Sigma_bf;
%
% O(e^2) approximation
Psf0=-Bs*Bf'*inv(Af');
Qsf0=-Cs'*Cf*inv(Af);
sigma2s=sqrt(eig(Wcs*Wos+eps*Psf0*Qsf0'));
sigma2f=sqrt(eig(Wcf*Wof+eps*Psf0'*Qsf0));
%
% NOT GOOD TO USE EXACT Wcs, Wos, Wcf, Wof, and approximate Psfo and Qsf0
% produces negative eigenvalues
%
% ApproximateSlow System Balancing
%
sys_sappr=ss(A0,B0,C0,D0);
[sys_sbappr,Sigma_bsappr,Tbsappr,Tbinvsappr]=balreal(sys_sappr);
Absappr=Tbsappr*A0*Tbinvsappr;
Bbsappr=Tbsappr*B0;
Cbsappr=C0*Tbinvsappr;
Sigma_bsappr;
%
% Appriximate Fast System Balancing
%
D2=zeros(3,1);
sys_fappr=ss(A4/eps,B2/eps,C2,D2);
[sys_fbappr,Sigma_bfappr,Tbfappr,Tbinvfappr]=balreal(sys_fappr);
Abfappr=Tbfappr*A4*Tbinvfappr;
Bbfappr=Tbfappr*B2;
Cbfappr=C2*Tbinvfappr;
Sigma_bfappr;
%
% Shahruz's Work
%
P1bar=lyap2(A0,B0*B0');
Q1bar=lyap2(A0',C0'*C0);
P3bar=lyap2(A4,B2*B2');
Q3bar=lyap2(A4',C2'*C2);
HSVslow_appr=sqrt(eig(P1bar*Q1bar));
HSVfast_appr=sqrt(eig(P3bar*Q3bar));
%
% O(e) Improvement of Shahruz's Work
%
% As1=(eye(3)-eps*A2*inv(A4)*A3)*A0
% Af1=A4+eps*inv(A4)*A3*A2
% Bs1=
% Bf1=
% Cs1=
% Cf1=
% P1bar=lyap2(A0,B0*B0');
% Q1bar=lyap2(A0',C0'*C0);
% P3bar=lyap2(A4,B2*B2');
% Q3bar=lyap2(A4',C2'*C2);
% HSVslow_appr=sqrt(eig(P1bar*Q1bar))
% HSVfast
%
% Summary of Hankel Singular Values
%
Sigma;
Sigma_bs;
%
Sigma_bf;
Sigma_bsappr;
Sigma_bfappr;
% Diagonalization is ok but MUST KEEP Cnew=[Cs Cf]
% to get the exact Hankel Singular Values, also B=[Bs;Bf/eps]
AAA=[As zeros(3,5);zeros(5,3) Af/eps];
BBB=[Bs;Bf/eps];
CCC=[Cs Cf];
DDD=Dsp;
sys_Diag=ss(AAA,BBB,CCC,DDD);
[sys_BiagBal,SIGMA]=balreal(sys_Diag);
SIGMA;

% Ok now, got the exact ones.
% Problems with Shahruz work for using the approximate DECOUPLED Systems
CMsappr=ctrb(A0,B0);
Crank_s_appr=rank(CMsappr);
CMfappr=ctrb(A4,B2);
Crank_f_appr=rank(CMfappr);
CMfexact=ctrb(Af,Bf);
Crank_fast_exact=rank(CMfexact);
delta_i=Sigma-sigma_i
percent=(delta_i./Sigma)*100


