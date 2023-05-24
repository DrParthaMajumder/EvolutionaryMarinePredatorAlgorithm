function [LB,UB,D,fobj] = Get_Functions_details(F)


switch F
    
    case 'Ackley_F1'
        fobj = @Ackley_F1;
        LB=-32.768;
        UB=32.768;
        D=30;
        
    case 'Beale_F2'
        fobj = @Beale_F2;
        LB=-4.5;
        UB=4.5;
        D=2;  
        
    case 'Bohachevsky_F3'
        fobj = @Bohachevsky_F3;
        LB=-100;
        UB=100;
        D=2; 
        
    case 'Booth_F4'
        fobj = @Booth_F4;
        LB=-10;
        UB=10;
        D=2; 
        
    case 'BUKINN6_F5'
        fobj = @BUKINN6_F5;
        LB=-10;
        UB=1;
        D=2; 
        
    case 'Colville_F6'
        fobj = @Colville_F6;
        LB=-10;
        UB=10;
        D=4;
        
    
    case 'Cross_In_Tray_F7'
        fobj = @Cross_In_Tray_F7;
        LB=-10;
        UB=10;
        D=2;
        
        
    case 'DejongN5_1_F8'
        fobj = @DejongN5_1_F8;
        LB=-65.536;
        UB=65.536;
        D=2; 
        
    case 'Dixonprice_F9'
        fobj = @Dixonprice_F9;
        LB=-10;
        UB=10;
        D=3;
        
    case 'Drop_Wave_F10'
        fobj = @Drop_Wave_F10;
        LB=-5.12;
        UB=5.12;
        D=2; 
        
        
     case 'EASOM1_F11'
        fobj = @EASOM1_F11;
        LB=-100;
        UB=100;
        D=2;
        
                
     case 'Eggholder_F12'
        fobj = @Eggholder_F12;
        LB=-512;
        UB=512;
        D=2;
        
                        
     case 'GoldsteinPrice_F13'
        fobj = @GoldsteinPrice_F13;
        LB=-2;
        UB=2;
        D=2;
        
     case 'GoldsteinPrice_Scaled_F14'
        fobj = @GoldsteinPrice_Scaled_F14;
        LB=-2;
        UB=2;
        D=2;
        
                   
     case 'Griewank_F15'
        fobj = @Griewank_F15;
        LB=-600;
        UB=600;
        D=30; 
        
                   
     case 'Hartmann_3D_F16'
        fobj = @Hartmann_3D_F16;
        LB=0;
        UB=1;
        D=3;
        
        
                           
     case 'Hartmann_4D_F17'
        fobj = @Hartmann_4D_F17;
        LB=0;
        UB=1;
        D=4;
    
                                          
     case 'Hartmann_6D_F18'
        fobj = @Hartmann_6D_F18;
        LB=0;
        UB=1;
        D=6;
        
                       
     case 'Holder_Table_F19'
        fobj = @Holder_Table_F19;
        LB=-10;
        UB=10;
        D=2;
        
     
     case 'Langermann_F20'
        fobj = @Langermann_F20;
        LB=0;
        UB=10;
        D=2;
        
             
     case 'Levy_F21'
        fobj = @Levy_F21;
        LB=-10;
        UB=10;
        D=5;
        
                     
     case 'LevyN13_F22'
        fobj = @LevyN13_F22;
        LB=-10;
        UB=10;
        D=2;
        
                             
     case 'Matyas_F23'
        fobj = @Matyas_F23;
        LB=-10;
        UB=10;
        D=2;
        
        
                                     
     case 'Mccormick_F24'  
        fobj = @Mccormick_F24;
        LB=[-1.5,-3];
        UB=[4,4];
        D=2;
     
     case 'Michalewicz1_F25'  
        fobj = @Michalewicz1_F25;
        LB=0;
        UB=pi;
        D=10;  % D=2; D=5; D=10
          
     case 'EggcratefcnF26'  
        fobj = @EggcratefcnF26;
        LB=-5;
        UB=5;
        D=2;     
                   
     case 'Permdb_F27'  
        fobj = @Permdb_F27;
        LB=1;
        UB=5;
        D=5; 
        
        
     case 'Powell_F28'      
        fobj = @Powell_F28;
        LB=-4;
        UB=5;
        D=10;
        
                     
     case 'Sum_Power_F29'      
        fobj = @Sum_Power_F29;
        LB=-1;
        UB=1;
        D=10;
        
                             
     case 'Rastrigin1_F30'      
        fobj = @Rastrigin1_F30;
        LB=-5.12;
        UB=5.12;
        D=30;
                        
     case 'Rosenbrock1_F31'      
        fobj = @Rosenbrock1_F31;
        LB=-2.048;
        UB=2.048;
        D=30;
        
                                
     case 'Rotted_hyper_ellipsoid_F32'      
        fobj = @Rotted_hyper_ellipsoid_F32;
        LB=-65.536;
        UB=65.536;
        D=10;
        
                                        
     case 'SchafferN2_F33'      
        fobj = @SchafferN2_F33;
        LB=-100;
        UB=100;
        D=2;
        
                
                                        
     case 'SchafferN4_F34'      
        fobj = @SchafferN4_F34;
        LB=-100;
        UB=100;
        D=2; 
        
                                                
     case 'Schwef1_F35'      
        fobj = @Schwef1_F35;
        LB=410;
        UB=430;
        D=10;
        
                
                                                
     case 'Shekel_F36'      
        fobj = @Shekel_F36;
        LB=0;
        UB=10;
        D=4; 
                
                                                
     case 'Shubert_F37'      
        fobj = @Shubert_F37;
        LB=-5.12;
        UB=5.12;
        D=2;
        
        
        
                                                        
     case 'Six_Hump_Camel_F38'      
        fobj = @Six_Hump_Camel_F38;
        LB=[-3,-2];
        UB=[3,2];
        D=2;
        
  
                                                                
     case 'Sphere1_F39'      
        fobj = @Sphere1_F39;
        LB=-5.12;
        UB=5.12;
        D=30;
        
        
                                                                        
     case 'Styblinski_Tang_F40'      
        fobj = @Styblinski_Tang_F40;
        LB=-5;
        UB=5;
        D=30;
        
                                                                                
     case 'Sum_Square_Function_F41'      
        fobj = @Sum_Square_Function_F41;
        LB=-5.12;
        UB=5.12;
        D=10;
        
                                                                                        
     case 'Three_Hump_Camel_F42'      
        fobj = @Three_Hump_Camel_F42;
        LB=-5;
        UB=5;
        D=2;   
        
        
        
                                                                                                
     case 'Trid_F43'      
        fobj = @Trid_F43;
        LB=-100;
        UB=100;
%         D=6; 
        D=10;
        
                                                                                                
     case 'HimmelblaufcnF44'      
        fobj = @HimmelblaufcnF44;
        LB=-6;
        UB=6;
        D=2; 
        
                                                                                                        
     case 'Shubert4fcn45'      
        fobj = @Shubert4fcn45;
        LB=-10;
        UB=10;
%          D=2; 
        D=30;
        
        
        
                                                                                                                
     case 'alpinen2fcn46'      
        fobj = @alpinen2fcn46;
%         LB=7.916;
%         UB=7.918;
        LB=0;
        UB=10;
        D=10;
        
        
        
end

end

%%
function F = Ackley_F1(x,D)
% D=D
% xg=[0,0,0,0,0....0]
% Fg=0
sum1 = 0;
sum2 = 0; 
for kk=1:D
    sum1 = sum1+x(kk)^2;
    sum2 = sum2+cos((2*pi)*x(kk));
end
F = 20 + exp(1)-20*exp(-0.2*sqrt(1/D*sum1))-exp(1/D*sum2);
end

%%
function [ F ] = Beale_F2( x,D )
% D=2 (A fixed Value)
% xg=[3,0.5]
% Fg=0
x1 = x(1);
x2 = x(2);
term1 = (1.5 - x1 + x1*x2)^2;
term2 = (2.25 - x1 + x1*x2^2)^2;
term3 = (2.625 - x1 + x1*x2^3)^2;

F = term1 + term2 + term3;
end

%%
function [ F ] = Bohachevsky_F3( x,D )
% D=2;
% xg=[0,0]
% Fg=0

x1=x(1);
x2=x(2);
term1 = x1^2;
term2 = 2*x2^2;
term3 = -0.3 * cos(3*pi*x1);
term4 = -0.4 * cos(4*pi*x2);

F = term1 + term2 + term3 + term4 + 0.7;
end

%%
function [ F ] = Booth_F4(x,D)
% D=2;
% xg=[1,3];
% Fg=0;
x1=x(1);
x2=x(2);
term1 = (x1 + 2*x2 - 7)^2;
term2 = (2*x1 + x2 - 5)^2;
F = term1 + term2;
end

%%
function [ F ] = BUKINN6_F5( x,D )
% D=2;
% xg=[-10,1]
% Fg=0
x1=x(:,1);
x2=x(:,2);
s1 = 100 * sqrt(abs(x2 - 0.01*x1^2));
s2 = 0.01 * abs(x1+10);
F = s1 + s2;
end

%% 
function [ F ] = Colville_F6(x,D)
% D=4;
% xg=[1,1,1,1];
% Fg=0;

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

term1 = 100 * (x1^2-x2)^2;
term2 = (x1-1)^2;
term3 = (x3-1)^2;
term4 = 90 * (x3^2-x4)^2;
term5 = 10.1 * ((x2-1)^2 + (x4-1)^2);
term6 = 19.8*(x2-1)*(x4-1);

F = term1 + term2 + term3 + term4 + term5 + term6;

end


%%
function [ F ] = Cross_In_Tray_F7(x,D)
% D=2;
% xg=[1.3491,1.3491];
% Fg=-2.06261;

x1=x(:,1);
x2=x(:,2);

fact1 = sin(x1)*sin(x2);
fact2 = exp(abs(100 - sqrt(x1^2+x2^2)/pi));

F = -0.0001 * (abs(fact1*fact2)+1)^0.1;
end


function [ F ] = DejongN5_1_F8(x,D)
% D=2
% xg=[-31.9769  -31.9804]
% Fg=0
% @dejong5fcn : Matlab inbuilt function

x1=x(:,1);
x2=x(:,2);
sum = 0;

AA = zeros(2, 25);
aa = [-32, -16, 0, 16, 32];
AA(1, :) = repmat(aa, 1, 5);
ar = repmat(aa, 5, 1);
ar = ar(:)';
AA(2, :) = ar;

for ii = 1:25
    a1i = AA(1, ii);
    a2i = AA(2, ii);
    term1 = ii;
    term2 = (x1 - a1i)^6;
    term3 = (x2 - a2i)^6;
    new = 1 / (term1+term2+term3);
    sum = sum + new;
end

F = 1 / (0.002 + sum);
end



function [ F ] = Dixonprice_F9( x ,D)
% D=D;
% xg=[1.00024173484239         0.707172787008294         0.594605228353372];
% Fg=0;
%xi=2^(-((2^i)-2)/(2^i))
term1 = (x(1)-1)^2;
sum = 0;
for kk = 2:1:D
	xin = x(kk);
	xold = x(kk-1);
	new = kk * (2*xin^2 - xold)^2;
	sum = sum + new;
end
F = term1 + sum;
end


function [ F ] = Drop_Wave_F10(x,D)
% D=2
% xg=[0, 0];
% Fg=-1

x1=x(:,1);
x2=x(:,2);

frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
frac2 = 0.5*(x1^2+x2^2) + 2;

F = -frac1/frac2;
end


function [ F ] = EASOM1_F11(x,D)
% D=2
% xg=[pi, pi];
% Fg=-1;

x1=x(:,1);
x2=x(:,2);
fact1 = -cos(x1)*cos(x2);
fact2 = exp(-(x1-pi)^2-(x2-pi)^2);
F = fact1*fact2;
end


function [ F ] = Eggholder_F12( x,D )
D=2;
% xg=[512, 404.2319];
% Fg=-959.6407;
x1=x(1);
x2=x(2);
term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
term2 = -x1 * sin(sqrt(abs(x1-(x2+47))));
F = term1 + term2;
end


function [F] = GoldsteinPrice_F13(x,D)
% D=2
% xg=[0,-1];
% Fg=3;
x1 = x(1);
x2 = x(2);

fact1a = (x1 + x2 + 1)^2;
fact1b = 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2;
fact1 = 1 + fact1a*fact1b;

fact2a = (2*x1 - 3*x2)^2;
fact2b = 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2;
fact2 = 30 + fact2a*fact2b;

F = fact1*fact2;
end


function [F] = GoldsteinPrice_Scaled_F14(x,D)
% D=2
% xg=[0.5,0.25];
% Fg=3;
x1 = 4*x(1) - 2;
x2 = 4*x(2) - 2;

fact1a = (x1 + x2 + 1)^2;
fact1b = 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2;
fact1 = 1 + fact1a*fact1b;

fact2a = (2*x1 - 3*x2)^2;
fact2b = 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2;
fact2 = 30 + fact2a*fact2b;

prod = fact1*fact2;

F = (log(prod) - 8.693) / 2.427;
end


function [ F ] = Griewank_F15(x,D)
% D=D
% xg=[0,0,0,0,0....0]
% Fg=0

 Sum=0;
 pp=1;
 for kk=1:1:D
     Sum=Sum+(x(kk)^2)/4000;
     pp = pp*cos(x(kk)/sqrt(kk));
 end
 F=Sum-pp+1;
end


function [F] = Hartmann_3D_F16(x,D)
% D=3
% xg=[0.114614, 0.555649, 0.852547]
% Fg=-3.86278
% xg=[0.721704002493248         0.542379320681405         0.854264412555292]
% Fg=-3.67597355769227

Alpha_V = [1.0, 1.2, 3.0, 3.2]';
AA = [3.0, 10, 30;
     0.1, 10, 35;
     3.0, 10, 30;
     0.1, 10, 35];
PP = 10^(-4) * [3689, 1170, 2673;
               4699, 4387, 7470;
               1091, 8732, 5547;
               381, 5743, 8828];

outer = 0;
for kk = 1:1:4
	inner = 0;
	for kkk = 1:1:D
		x_kk = x(kkk);
		A_kkk = AA(kk, kkk);
		P_kkk = PP(kk, kkk);
		inner = inner + A_kkk*(x_kk-P_kkk)^2;
	end
	new = Alpha_V(kk) * exp(-inner);
	outer = outer + new;
end

F = -outer;
end


function [F] = Hartmann_4D_F17(x,D)
% D=4
% xg = c(0.1873, 0.1906, 0.5566, 0.2647)
% Fg = -3.135474

Alpha_V = [1.0, 1.2, 3.0, 3.2]';
AA = [10, 3, 17, 3.5, 1.7, 8;
     0.05, 10, 17, 0.1, 8, 14;
     3, 3.5, 1.7, 10, 17, 8;
     17, 8, 0.05, 10, 0.1, 14];
PP = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 381];
outer = 0;
for kk = 1:1:4
	inner = 0;
	for kkk = 1:1:D
		x_kk = x(kkk);
		A_kkk = AA(kk, kkk);
		P_kkk = PP(kk, kkk);
		inner = inner + A_kkk*(x_kk-P_kkk)^2;
	end
	new = Alpha_V(kk) * exp(-inner);
	outer = outer + new;
end

F = (1.1 - outer) / 0.839;
end



function [F] = Hartmann_6D_F18(x,D)
% D=6
% xg=[0.20169, 0.150011, 0.476874, 0.275332,0.311652, 0.6573]\
% Fg=-3.04245773783059
Alpha_V = [1.0, 1.2, 3.0, 3.2]';
AA = [10, 3, 17, 3.5, 1.7, 8;
     0.05, 10, 17, 0.1, 8, 14;
     3, 3.5, 1.7, 10, 17, 8;
     17, 8, 0.05, 10, 0.1, 14];
PP = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 381];
outer = 0;
for kk = 1:1:4
	inner = 0;
	for kkk = 1:1:D
		x_kk = x(kkk);
		A_kkk = AA(kk, kkk);
		P_kkk = PP(kk, kkk);
		inner = inner + A_kkk*(x_kk-P_kkk)^2;
	end
	new = Alpha_V(kk) * exp(-inner);
	outer = outer + new;
end
F = -(2.58 + outer) / 1.94;
end


function [ F ] = Holder_Table_F19(x,D )
% D=2;
% xg=[-8.05503496507683          -9.6645717318531]
% Fg=-19.2085
x1=x(1);
x2=x(2);
fact1 = sin(x1)*cos(x2);
fact2 = exp(abs(1 - sqrt(x1^2+x2^2)/pi));
F = -abs(fact1*fact2);
end


function [ F ] = Langermann_F20( x ,D)
% D=2
% xg=[2.91333111584265, 1.31556959589056];
% Fg=-4.05404569816266;
% Accuracy is a concern

M_V = 5;
C_V = [1, 2, 5, 2, 3];
A_V = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9];
  
outer = 0;
for kk = 1:1:M_V
    inner = 0;
    for kkk = 1:1:D
        x_kkk = x(kkk);
        AkK = A_V(kk,kkk);
        inner = inner + (x_kkk-AkK)^2;
    end
    new = C_V(kk) * exp(-inner/pi) * cos(pi*inner);
    outer = outer + new;
end

F = outer;
end


function [ F ] = Levy_F21(x,D)
%D=D
%xg=[1,1,1,....1];
%Fg=0; 

for kk = 1:1:D
	w(kk) = 1 + (x(kk) - 1)/4;
end
term1 = (sin(pi*w(1)))^2;
term3 = (w(D)-1)^2 * (1+(sin(2*pi*w(D)))^2);
sum = 0;
for kk = 1:1:(D-1)
	w_kk = w(kk);
        new = (w_kk-1)^2 * (1+10*(sin(pi*w_kk+1))^2);
	sum = sum + new;
end
F = term1 + sum + term3;  
end


function [ F ] = LevyN13_F22(x,D)
% D=2
% xg=[1,1];
% Fg=0;

x1 = x(1);
x2 = x(2);

term1 = (sin(3*pi*x1))^2;
term2 = (x1-1)^2 * (1+(sin(3*pi*x2))^2);
term3 = (x2-1)^2 * (1+(sin(2*pi*x2))^2);
F = term1 + term2 + term3;
end

function [ F ] = Matyas_F23(x,D)
% D=2;
% xg=[0,0];
% Fg=0
x1=x(1);
x2=x(2);
term1 = 0.26 * (x1^2 + x2^2);
term2 = -0.48*x1*x2;
F = term1 + term2;
end

function [ F ] = Mccormick_F24(x,D)
% D=2;
% xg=[-0.54719, -1.54719];
% Fg=-1.9133;
x1=x(1);
x2=x(2);
term1 = sin(x1 + x2);
term2 = (x1 - x2)^2;
term3 = -1.5*x1;
term4 = 2.5*x2;
F = term1 + term2 + term3 + term4 + 1;
end



function [ F ] = Michalewicz1_F25(x,D)
% D=2; 
% xg=[2.20, 1.57]
% Fg=-1.8013
% 
% D=5;
% Fg=-4.687658;
% 
% D=10;
% Fg=-9.66015;


mv=10; % Constant Parameter
sum = 0;

for kk = 1:1:D
	x_k = x(kk);
	new = sin(x_k) * (sin(kk*x_k^2/pi))^(2*mv);
	sum  = sum + new;
end

F = -sum;
end

function F = EggcratefcnF26(x,D)
    % D=2
    % xg=[0,0]
    % Fg=0;
        
    x1 = x(1);
    x2 = x(2);    
    F = x1^2 + x2^2 + (25 * (sin(x1)^2 + sin(x2)^2));
end

function [ F ] = Permdb_F27( x, D)
% D=D
% xg=[1,2,3,4,5] 
% Fg=0
bb = 0.5;
outer = 0;
for kk = 1:1:D
	inner = 0;
	for mm = 1:1:D
		xkm = x(mm);
        inner = inner + (mm^kk+bb)*((xkm/mm)^kk-1);
    end
	outer = outer + inner^2;
end

F = outer;
end




function [F] = Powell_F28(x,D)
% D=D
% xg=[0,0,0,0.....0];
% Fg=0
sum = 0;
for kk = 1:1:(D/4)
	term1 = (x(4*kk-3) + 10*x(4*kk-2))^2;
	term2 = 5 * (x(4*kk-1) - x(4*kk))^2;
	term3 = (x(4*kk-2) - 2*x(4*kk-1))^4;
	term4 = 10 * (x(4*kk-3) - x(4*kk))^4;
	sum = sum + term1 + term2 + term3 + term4;
end

F = sum;
end


function [ F ] = Sum_Power_F29(x,D)
% D=D
% xg=[0,0,0,....0];
% Fg=0;

sum = 0;
for kk = 1:1:D
    xk = x(kk);
    new = (abs(xk))^(kk+1);
    sum = sum + new;
end

F = sum;
end


function [ F ] = Rastrigin1_F30( x ,D)
% D=D
% xg=[0,0,....0];
% Fg=0

Sum=0;
    for kk=1:1:D
    Sum=Sum+(x(kk)^2-10*cos(2*pi*x(kk))+10);    
    end
F=Sum;
end



function [ F ] = Rosenbrock1_F31( x ,D)
% D=D
% xg=[1,1,1,...1]
% Fg=0
sum = 0;
for kk = 1:1:(D-1)
	x_kk = x(kk);
	xnext = x(kk+1);
	new = 100*(xnext-x_kk^2)^2 + (x_kk-1)^2;
	sum = sum + new;
end
F = sum;
end


function [ F ] = Rotted_hyper_ellipsoid_F32(x,D)
% D=D
% xg=[0,0,0,....0]
%Fg=0
outer = 0;
for kk = 1:1:D
    inner = 0;
    for kkk = 1:kk
        x_kkk = x(kkk);
        inner = inner + x_kkk^2;
    end
    outer = outer + inner;
end
F = outer;  
end


function [F] = SchafferN2_F33(x,D)
% D=2;
% xg=[0,0]
% Fg=0
x1 = x(1);
x2 = x(2);

fact1 = (sin(x1^2-x2^2))^2 - 0.5;
fact2 = (1 + 0.001*(x1^2+x2^2))^2;
F = 0.5 + fact1/fact2;
end

function [F] = SchafferN4_F34(x,D)
% D=2;
% xg=[0,1.253115]
% Fg=0.292579
x1 = x(1);
x2 = x(2);
numeratorcomp = (cos(sin(abs(x1 .^ 2 - x2 .^ 2))) .^ 2) - 0.5; 
denominatorcomp = (1 + 0.001 * (x1 .^2 + x2 .^2)) .^2 ;
F = 0.5 + numeratorcomp ./ denominatorcomp;
end


function [F] = Schwef1_F35(x,D)
% D=D
% xg=[420.9687,420.9687,..... 420.9687]
% Fg=0
sum = 0;
for kk = 1:1:D
	x_kk = x(kk);
	sum = sum + x_kk*sin(sqrt(abs(x_kk)));
end

F = 418.9829*D - sum;
end

function [ F ] = Shekel_F36( x,D )
% D=4;
% xg=[4,4,4,4]
% Fg=-10.53628
mm = 10;
bb = 0.1 * [1, 2, 2, 4, 4, 6, 3, 7, 5, 5]';
cc = [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0;
     4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6;
     4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0;
     4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6];

outer = 0;
for kk = 1:1:mm
	b_kk = bb(kk);
	inner = 0;
	for mm = 1:1:D
		xkkmm = x(mm);
		cckkmm = cc(mm, kk);
		inner = inner + (xkkmm-cckkmm)^2;
	end
	outer = outer + 1/(inner+b_kk);
end

F = -outer;
end




function [ F ] = Shubert_F37( x,D )
% D=2;
% xg=[-0.800247198673221,        -1.42516777282714]
% Fg=-186.7309
x1=x(1);
x2=x(2);

sum1 = 0;
sum2 = 0;

for kk = 1:5
	new1 = kk * cos((kk+1)*x1+kk);
	new2 = kk * cos((kk+1)*x2+kk);
	sum1 = sum1 + new1;
	sum2 = sum2 + new2;
end

F = sum1 * sum2;
end


function [ F ] = Six_Hump_Camel_F38(x,D)
% D=2;
% xg=[-0.0898410125071629         0.712656335200722]
% Fg=-1.031628453486
x1=x(:,1);
x2=x(:,2);
F=(4*x1.^2)-(2.1*x1.^4)+((x1.^6)/3)+(x1.*x2)-(4*x2.^2)+(4*x2.^4);
end


function [ F ] = Sphere1_F39( x,D )
% D=D
% xg=[0,0,0,....,0]
% Fg=0;
 Sum=0;
 for kk=1:1:D
     Sum=Sum+x(kk)^2;      
 end
F=Sum; 
end


function [F] = Styblinski_Tang_F40(x,D)
% D=D
% xg=[-2.90359711129225         -2.90347094574996,............, -2.90359711129225         -2.90347094574996]
% Fg=-39.16599*D
sum = 0;
for kk = 1:1:D
	x_k = x(kk);
	new = x_k^4 - 16*x_k^2 + 5*x_k;
	sum = sum + new;
end
F = sum/2;
end

function [ F ] = Sum_Square_Function_F41(x,D)
% D=D
% xg=[0,0,....0];
% Fg=0
sum = 0;
for kk = 1:1:D
	x_kk = x(kk);
	sum = sum + kk*x_kk^2;
end
F = sum; 
end


function [ F ] = Three_Hump_Camel_F42(x,D)
% D=2;
% xg=[0,0]
% Fg=0;
x1=x(:,1);
x2=x(:,2);
term1 = 2*x1^2;
term2 = -1.05*x1^4;
term3 = x1^6 / 6;
term4 = x1*x2;
term5 = x2^2;

F = term1 + term2 + term3 + term4 + term5;
end


function [ F ] = Trid_F43(x,D)
% D=D
% Fg=-50; D=6;
% Fg=-200; D=10;
sum1 = 0;
sum2 = 0;
for kk = 1:1:D;
    sum1 = sum1+(x(kk)-1)^2;    
end
for kk = 2:1:D;
    sum2 = sum2+x(kk)*x(kk-1);    
end
F = sum1-sum2;
end


function [ F ] = HimmelblaufcnF44(x,D)
% D=2
% xg=[(3,2);(2.805118, 3.283186); (3.779310, 3.283186); (3.58445,1.848126)]
% Fg=0
x1 = x(:, 1);
x2 = x(:, 2); 
F = ((x1.^ 2 + x2 - 11).^2) + ((x1 + (x2.^ 2) - 7) .^ 2);
end
   

function F = Shubert4fcn45(x,D)
    % D=D
    % Fg=-25.70458 D=2
    F = 0;
    for kk = 1:1:D
        for kkk = 1:5
            F = F + kkk * cos(((kkk + 1) * x(:, kk)) + kkk);
        end
    end
end


function [ F ] = alpinen2fcn46(x,D)
% x=[ 7.917, 7.917, 7.917,........7.917]
% D=D
% Fg=2.808^n
F= prod(sqrt(x) .* sin(x), 2);
end





