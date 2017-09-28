clear;clc;
% Start time
tic;
% Kernel Function
kernel_gaussian = @(u,theta)((sqrt(2*pi)*theta)^-1).*exp(-0.5.*u.^2/theta^2);
kernel_gaussian2 = @(u1,u2,theta)((2*pi*theta^2)^-1).*exp(-0.5.*(u1.^2+u2.^2)/theta^2);
kernel_gaussian3 = @(u1,u2,u3,theta)((2*pi*sqrt(2*pi)*theta^3)^-1).*exp(-0.5.*(u1.^2+u2.^2++u3.^2)/theta^2);
% Input path
InportPass='e:\MATLAB_WorkSpace\' ;
ExportPass='e:\MATLAB_WorkSpace\' ;
% Data name
DName=cell(4,1);
name=sprintf('V20SFL' ); DName{1}=name; % 20% vlotage dip, Symmetrical fault, Full load
name=sprintf('V20SPL' ); DName{2}=name; % 20% voltage dip, Symmetrical fault, Partial load
name=sprintf('V20UFL' ); DName{3}=name; % 20% voltage dip, Unsymmetrical fault, Full load
name=sprintf('V20UPL' ); DName{4}=name; % 20% voltage dip, Unsymmetrical fault, Partial load
clear name;
 
for SS=1
    % Data import
    fullname = [InportPass,DName{SS}, '.mat'];
    Input = importdata(fullname);
    %Imput (1k=50; 5k=10; 10k=5; 50k=1)
    t = downsample(double(Input.Data1_time_Ia), 50);
    Ia = downsample(double(Input.Data1_Ia), 50);
    Ib = downsample(double(Input.Data1_Ib), 50);
    Ic = downsample(double(Input.Data1_Ic), 50);
    % Step size
    tstep = t(2)-t(1);
    % Window size
    WINDOW_SIZE = floor((1/tstep)/50); 
    % Data length
    L = length(t);
    % Counter
    m=1;
    % TE Calculation
    for i = round(4/tstep:6/tstep) % from 4s to 6s
        % Data in one window
        j = round((i-WINDOW_SIZE):1:i);
        j=j';
        % Ia->Ib
        % Martix
        n=numel(j);
        I_mat1 = repmat(Ib(j),1,n); % In
        I_mat2 = repmat(Ib(j+1),1,n); % In+1
        I_mat3 = repmat(Ia(j),1,n); % Jn
        % P(In)
        I_K_In = kernel_gaussian((I_mat1 - I_mat1'),1); % theta=1
        I_P_In = sum(I_K_In')'/n;
        % P(In+1,In)
        I_K_In2 = kernel_gaussian2((I_mat1 - I_mat1'),(I_mat2 - I_mat2'),1); % theta=1
        I_P_In2 = sum(I_K_In2')'/n;
        % P(In,Jn)
        I_K_IJn = kernel_gaussian2((I_mat1 - I_mat1'),(I_mat3 - I_mat3'),1); % theta=1
        I_P_IJn = sum(I_K_IJn')'/n;
        % P(In+1,In,Jn)
        I_K_IJ3 = kernel_gaussian3((I_mat1 - I_mat1'),(I_mat2 - I_mat2'),(I_mat3 - I_mat3'),1); % theta=1
        I_P_IJ3 = sum(I_K_IJ3')'/n;
        % TE
        TE_AB(m,1) = t(i);
        TE_AB(m,2) = sum(I_P_IJ3.*log2(I_P_IJ3.*I_P_In./(I_P_In2.*I_P_IJn)));
        m=m+1;
    end
    %Data export
    fullname=[ExportPass, 'TE_I_',DName{SS}];
    save(fullname);
end
% End time
toc;
