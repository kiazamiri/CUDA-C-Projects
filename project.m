clear all;
%
% parameters
N=30; % No. of snapshots
M=10;  % No. of sensors
tet_no=3; % nominal impinging angle
tet_ac=5; % actual impinging angle(In Example 2)
tet_mean = 3; % mean impinging angle for b1~b4 (in Example 3)
%
% interference generation
tet_I= [30, 50]; % impinging angles of interfering sources
J  = length(tet_I); % No. of interfering sources
a_I = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_I*pi/180.0)); % interfering steering matrix
INR  = 1000; % interference noise ratio
R_I_n=INR*a_I*a_I'+eye(M,M); % ideal interference-plus-noise covariance matrix % used to calculate SINR
%
SNR_dB = [-20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 ...
      -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];
for i=1:length(SNR_dB)
  SNR(i) = 10^(SNR_dB(i)/10);
end
%
SINR_SMI=zeros(1,length(SNR));
SINR_eig=zeros(1,length(SNR));
SINR_PRO=zeros(1,length(SNR));
SINR_LSMI=zeros(1,length(SNR));
%
% No. of simulations=200
for loop=1:20;
  a = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_no*pi/180.0));% nominal a
% % generation of actual signal steering vector a_tilde
% % Example 1: Exactly Known Signal Steering Vector
  a_tilde = a;
%  % Example 2: Signal Look Direction Mismatch
%     a_tilde = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_ac*pi/180.0));
%  % Example 3: Signal Signature Mismatch Due to Coherent Local Scattering
%     var=4;
%     b1 = exp(1i*2.0*[0:M-1]'*pi*0.5*sin((tet_mean+rand*sqrt(3*var))*pi/180.0)); % better to use randn for Gaussian noise
%     b2 = exp(1i*2.0*[0:M-1]'*pi*0.5*sin((tet_mean+rand*sqrt(3*var))*pi/180.0));
%     b3 = exp(1i*2.0*[0:M-1]'*pi*0.5*sin((tet_mean+rand*sqrt(3*var))*pi/180.0));
%     b4 = exp(1i*2.0*[0:M-1]'*pi*0.5*sin((tet_mean+rand*sqrt(3*var))*pi/180.0));
%     psi1 = 2*rand*pi; psi2 = 2*rand*pi; psi3 = 2*rand*pi; psi4 = 2*rand*pi;
%     a_tilde = a + exp(1i*psi1)*b1 + exp(1i*psi2)*b2 + exp(1i*psi3)*b3 + exp(1i*psi4)*b4; % actual a (tilde)
  a_tilde = sqrt(10)*a_tilde/norm(a_tilde)';
  %
  Ncount=0;
  %
  for ps=SNR;
      Ncount=Ncount+1;
      R_hat = zeros(M);
      for l=1:N
          s=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt(0.5);  % signal waveform
          Is = (randn(J,1) + 1i*randn(J,1)).*sqrt(INR)'*sqrt(0.5); % interference waveform
          n = (randn(M,1) + 1i*randn(M,1))*sqrt(0.5); % noise
          x = a_I*Is + s*a_tilde + n; % snapshot with signal
          R_hat = R_hat + x*x';
      end
      R_hat=R_hat/N;  % sample covariance matrix
      %
      % sample matrix inversion (SMI) algorithm
      w_SMI=inv(R_hat)*a;
      %
      % LSMI algorithm
      R_dl_hat = R_hat + eye(M,M)*10;
      w_LSMI=inv(R_dl_hat)*a;
      %
      % eigenspace-based robust beamformer
      [RR,D]=eig(R_hat);   % eigenvalue decomposition of covariance matrix
      RR_subs=RR(:,M-2:M);
      as_new=[RR_subs,zeros(M,M-3)]*a;  % projection of waveform to SIGNAL and INTERFERENCE space
      ws_EIG=inv(R_hat)*as_new;
      %
      % proposed robust beamformer using worst-case performance optimization
      R_hat = R_hat + eye(M,M)*0.01; % diagonal loading
      U=chol(R_hat);  % Cholesky factorization
      U_bowl=[real(U),-imag(U);imag(U),real(U)];  % 'U_bowl' defined based on 'U'
      M2 = 2*M;
      a_bowl=[real(a);imag(a)];  % 'a_bowl' difined based on 'a'
      a_bar = [imag(a);-real(a)];
      FT = [1,zeros(1,M2);zeros(M2,1),U_bowl;0,a_bowl';zeros(M2,1),3*eye(M2,M2);0,a_bar'];
      cvx_begin
          variable y(M2+1,1)
          minimize (y(1))
          c = FT*y;
          subject to
              c(1) >= norm(c(2:M2+1));
              c(M2+2)-1 >= norm(c(M2+3:2*M2+2));
              c(2*M2+3) == 0;
      cvx_end
      w_PRO_bowl = y(2:21);
      w_PRO = w_PRO_bowl(1:10) + 1i * w_PRO_bowl(11:20);
      %
      SINR_SMI_v(Ncount) = ps * (abs(w_SMI'*a_tilde) * abs(w_SMI'*a_tilde)) / abs(w_SMI'*R_I_n*w_SMI);  % SINR
      SINR_LSMI_v(Ncount) = ps * (abs(w_LSMI'*a_tilde) * abs(w_LSMI'*a_tilde)) / abs(w_LSMI'*R_I_n*w_LSMI);
      SINR_EIG_v(Ncount) = ps * (abs(ws_EIG'*a_tilde) * abs(ws_EIG'*a_tilde)) / abs(ws_EIG'*R_I_n*ws_EIG);
      SINR_PRO_v(Ncount) = ps * (abs(w_PRO'*a_tilde) * abs(w_PRO'*a_tilde)) / abs(w_PRO'*R_I_n*w_PRO); % SINR for robust beamformer
      SINR_opt(Ncount) = ps * a_tilde' * inv(R_I_n) * a_tilde;
  end
  SINR_SMI     = ((loop-1)/loop)*SINR_SMI + (1/loop)*SINR_SMI_v; % just average
  SINR_LSMI    = ((loop-1)/loop)*SINR_LSMI + (1/loop)*SINR_LSMI_v;
  SINR_eig     = ((loop-1)/loop)*SINR_eig + (1/loop)*SINR_EIG_v;
  SINR_PRO   = ((loop-1)/loop)*SINR_PRO + (1/loop)*SINR_PRO_v;
end;
save('Fig2_for_EX1','SINR_opt','SINR_SMI', 'SINR_LSMI','SINR_eig','SINR_PRO') % Ex1
% save('Fig6_for_EX2','SINR_opt','SINR_SMI', 'SINR_LSMI','SINR_eig','SINR_PRO') % Ex2
% save('Fig8_for_EX3','SINR_opt','SINR_SMI', 'SINR_LSMI','SINR_eig','SINR_PRO') % Ex3
