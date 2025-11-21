%% modify parameters in this part
DG   = 4;                       % node = DG number
line = 3;                       % line number
load = 3;                       % load number

coef = 8;
RlineCoef = 0.2;
LlineCoef = 5;
RloadCoef = 0.5;
LloadCoef = 1.5;

% MG structure
% DG1 ->- line1 ->- DG2 ->- line2 ->- DG3 ->- line3 ->- DG4
%                    |                 |                 |
%                    V                 V                 V
%                  load1             load2             load3


% incidence ∇ for line
del = [ 1  -1   0   0 ;           
        0   1  -1   0 ;           
        0   0   1  -1];           
delT = del';      


% for load connect
phi = [ 0  1  0  0;
        0  0  1  0;
        0  0  0  1];
phiT  = phi';

F = [ -1  0  0  0;
      -1  0  0  0;
      -1  0  0  0];

% secondary connect
D = 2*eye(4);
A = [0  1  0  1;
     1  0  1  0;
     0  1  0  1;
     1  0  1  0];
L = D-A;
G = [1 0 0 0];
G = diag(G);
I=[1 1 1 1]';

K1=10*1e-4;
Kp_sec=[K1 0 0 0;
    0 K1 0 0;
    0 0 K1/2 0;
    0 0 0 K1/2];

k2=2e-4;
Kq_sec=[k2 0 0 0;
    0 k2 0 0;
    0 0 k2/2 0;
    0 0 0 k2/2];

c1 = 10;
c2 = 1;
c3 = 5;
c4 = 2;

% --------------------line / load impedance---------------------------------
R_load_1=14.1397*RloadCoef;
X_load_1=14.1397*LloadCoef;

R_load_2=11.1664*RloadCoef;
X_load_2=5.5832*LloadCoef;

R_load_3=8.3282*RloadCoef;
X_load_3=2.7761*LloadCoef;

R_line_12=0.1*RlineCoef;
X_line_12=(1e-4) *100*pi*LlineCoef;

R_line_23=0.2*RlineCoef;
X_line_23=(1e-4) *100*pi*LlineCoef;

R_line_34=0.1*RlineCoef;
X_line_34=(2e-4) *100*pi*LlineCoef;

rho_line = [R_line_12/X_line_12, R_line_23/X_line_23, R_line_34/X_line_34];
rho_line = diag(rho_line);

X_line_reci = [1/X_line_12, 1/X_line_23, 1/X_line_34];
X_line_reci = diag(X_line_reci);

rho_load = [R_load_1/X_load_1, R_load_2/X_load_2, R_load_3/X_load_3];
rho_load = diag(rho_load);

X_load_reci = [1/X_load_1, 1/X_load_2, 1/X_load_3];
X_load_reci = diag(X_load_reci);
% --------------------droop-----------------------------------
Kp1=1e-3;    % DG1
Kq1=2e-4;

Kp2=1e-3;    % DG2
Kq2=2e-4;

Kp3=5e-4;    % DG3
Kq3=1e-4;

Kp4=5e-4;    % DG4
Kq4=1e-4;

M = [Kp1,Kp2,Kp3,Kp4]./coef;
N = [Kq1,Kq2,Kq3,Kq4]./coef;

M = diag(M);
N = diag(N);

% --------------------system parameters----------------------------------------
w0 = 2*pi*50;    % rated angular speed           
wc = 37;    % filter frequency               
Vi = 311;   % rated voltage

% --------------------matrix index----------------------------------------
theta =   1                         : DG;
omega =   DG+1                      : 2*DG;
V     =   2*DG+1                    : 3*DG;
Id_load = 3*DG+ 1                   : 3*DG+ load; 
Iq_load = 3*DG+ load + 1            : 3*DG+ 2*load;
Id_line = 3*DG+ 2*load + 1          : 3*DG + 2*load + line;
Iq_line = 3*DG + 2*load + line + 1  : 3*DG + 2*load + 2*line; 
OMEGA   = 3*DG + 2*load + 2*line + 1: 4*DG + 2*load + 2*line;
LAMBDA  = 4*DG + 2*load + 2*line + 1: 5*DG + 2*load + 2*line;

Nx = 5*DG + 2*line+2*load;                    
A  = zeros(Nx);

%% -------------automatically calculate (no need adjustment)-------
% (1)  θ
A(theta,omega) = eye(DG);

% (2)  ω
A(omega,omega) = -wc*eye(DG)+c1*(-L-G);
A(omega,Id_line)    = -1.5*(wc*M+c2*L*Kp_sec)*Vi*delT;
A(omega,Id_load)    = -1.5*(wc*M+c2*Kp_sec)*Vi*phiT;
A(omega,OMEGA)= wc*eye(DG);

% (3)  V
A(V,V)  = -wc*eye(DG);
A(V,Iq_line) = 1.5*(wc*N+c4*L*Kq_sec)*Vi*delT;
A(V,Iq_load) = 1.5*(wc*N+c4*L*Kq_sec)*Vi*phiT;
A(V,LAMBDA)= wc*eye(DG)+c3*(-L);

% (4) load d 
A(Id_load,V)  =  w0 * X_load_reci * (phi+F);
A(Id_load,Id_load) = -w0 * rho_load;
A(Id_load,Iq_load) =  w0 * eye(load);

% (5) load q 
A(Iq_load,theta)  =  w0*Vi * X_load_reci *(phi+F);
A(Iq_load,Iq_load) = -w0 * rho_load;
A(Iq_load,Id_load) =  -w0 * eye(load);

% (6)  line d
A(Id_line,V)  =  w0 * X_line_reci * del;
A(Id_line,Id_line) = -w0 * rho_line;
A(Id_line,Iq_line) =  w0 * eye(line);

% (7)  line q
A(Iq_line,theta)=  w0 * Vi * X_line_reci * del;
A(Iq_line,Id_line)   = -w0 * eye(line);
A(Iq_line,Iq_line)   = -w0 * rho_line;

% (8) OMEGA
A(OMEGA,omega) = c1*(-L-G);
A(OMEGA,Id_line) = -1.5*c2*L*Kp_sec*Vi*delT;
A(OMEGA,Id_load) = -1.5*c2*L*Kp_sec*Vi*phiT;

% (9) LAMBDA
A(LAMBDA,LAMBDA) = c3*(-L);
A(LAMBDA,Iq_line) = 1.5*c4*L*Kq_sec*Vi*delT;
A(LAMBDA,Iq_load) = 1.5*c4*L*Kq_sec*Vi*phiT;

%% --------------------eigenvalue--------------------------------------
SecondaryModel = eig(A);
     
figure;  
plot(real(SecondaryModel),imag(SecondaryModel),'x'); grid on;
% xlim([-1000,0]);
hold on;
title('Simplified model Sec');
xlabel('Real'); 
ylabel('Imag');

