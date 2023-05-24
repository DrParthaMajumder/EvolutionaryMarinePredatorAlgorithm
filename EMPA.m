clc
clear all
close all
format long g



%% Get Details of objective function
%% Unimodal function
% [LB,UB,D,fobj] = Get_Functions_details('Sphere1_F39'); % Fg=0;  D=D
% [LB,UB,D,fobj] = Get_Functions_details('Rosenbrock1_F31'); % Fg=0; D=D
% [LB,UB,D,fobj] = Get_Functions_details('Griewank_F15');  % Fg=0; D=D

%% Multimodal function
% [LB,UB,D,fobj] = Get_Functions_details('Ackley_F1');  % Fg=0; D=D
% [LB,UB,D,fobj] = Get_Functions_details('Rastrigin1_F30'); % Fg=0; D=D


% [LB,UB,D,fobj] = Get_Functions_details('HimmelblaufcnF44'); % Fg=0; D=2
% [LB,UB,D,fobj] = Get_Functions_details('Styblinski_Tang_F40');  % Fg=-39.16599*D 
% [LB,UB,D,fobj] = Get_Functions_details('Shubert4fcn45');  % Fg=-386.12
[LB,UB,D,fobj] = Get_Functions_details('alpinen2fcn46');  % Fg=2.808^n

%% Other Functions
% [LB,UB,D,fobj] = Get_Functions_details('Six_Hump_Camel_F38');  % Fg=-1.031628453486

%% DE PARAMETERS (1)
beta_min=0.2;     % Lower Bound of Scaling Factor
beta_max=0.8;     % Upper Bound of Scaling Factor
pCR=0.2;          % Crossover Probability
VarSize=[1 D];
%% DE PARAMETERS (1)

N=80; 
itmax=400;



if length(LB)==1
for kk=1:1:D
    lb(1:N,kk)=LB; 
    ub(1:N,kk)=UB;
end
end

if length(LB)~=1
for kk=1:1:D
    lb(1:N,kk)=LB(kk); 
    ub(1:N,kk)=UB(kk);
    
end
end


x=lb+(ub-lb).*rand(N,D);
  


Top_predator_pos=zeros(1,D);
Top_predator_fit=inf;
Fgbest_vect=zeros(1,itmax);
stepsize=zeros(N,D);
F=inf(N,1);

Xmin=lb;
Xmax=ub;
it=1;
FADs=0.2;
P=0.5;

%% DE initial computation
for ii=1:1:N
    F_DE(ii) = fobj(x(ii,:),D);
end
[F_g_bestDE,pp]=min(F_DE);
g_best_DE=x(pp,:);
%% DE initial computation




for it=1:1: itmax
    for kk=1:1:D
        x(:,kk)=min(x(:,kk),ub(:,kk));
        x(:,kk)=max(x(:,kk),lb(:,kk));
    end
    
    for ii=1:1:N
        ii=ii;
        %% Elite Opposition Based Learning
        x_O=lb(ii,:)+ub(ii,:)-x(ii,:);  % Opposite number
        F_O(ii,1)=fobj(x_O,D);
        %% Elite Opposition Based Learning        
        F(ii,1)= fobj(x(ii,:),D);
        
        if F_O(ii,1)<F(ii,1)
            disp('improvement from opposition based learning:')
            F(ii,1)=F_O(ii,1);
            x(ii,:)=x_O;
        end
        
        
        if F(ii,1)<Top_predator_fit
            Top_predator_fit=F(ii,1); 
            Top_predator_pos=x(ii,:);
        end
    end
    
    if it==1
        fit_old=F;
        x_old=x;
    end
    
    Inx=(fit_old<F);
    Indx=repmat(Inx,1,D);
    x=Indx.*x_old+~Indx.*x;
    F=Inx.*fit_old+~Inx.*F;
    fit_old=F;
    x_old=x;
    
    Elite=repmat(Top_predator_pos,N,1); 
    CF=(1-it/itmax)^(2*it/itmax);    
    RL=0.05*levyMPA(N,D,1.5);  
    RB=randn(N,D);
    
    for ii=1:1:N
        for jj=1:1:D 
            R=rand();
            if it<itmax/3
                stepsize(ii,jj)=RB(ii,jj)*(Elite(ii,jj)-RB(ii,jj)*x(ii,jj));
                x(ii,jj)=x(ii,jj)+P*R*stepsize(ii,jj);
            elseif it>itmax/3 && it<2*itmax/3
                if ii>N/2
                    stepsize(ii,jj)=RB(ii,jj)*(Elite(ii,jj)-RB(ii,jj)*x(ii,jj));
                    x(ii,jj)=Elite(ii,jj)+P*CF*stepsize(ii,jj);
                else
                    stepsize(ii,jj)=RL(ii,jj)*(Elite(ii,jj)-RL(ii,jj)*x(ii,jj));
                    x(ii,jj)=x(ii,jj)+P*R*stepsize(ii,jj);
                end
            else
                stepsize(ii,jj)=RL(ii,jj)*(Elite(ii,jj)-RL(ii,jj)*x(ii,jj));
                x(ii,jj)=Elite(ii,jj)+P*CF*stepsize(ii,jj);
            end
        end
    end
    for kk=1:1:D
        x(:,kk)=min(x(:,kk),ub(:,kk));
        x(:,kk)=max(x(:,kk),lb(:,kk));
    end
    
    for ii=1:1:N
        F(ii,1)= fobj(x(ii,:),D);
        if F(ii,1)<Top_predator_fit
            Top_predator_fit=F(ii,1);
            Top_predator_pos=x(ii,:);
        end
    end
    
    if it==1
        fit_old=F;
        x_old=x;
    end    
    Inx=(fit_old<F);
    Indx=repmat(Inx,1,D);
    x=Indx.*x_old+~Indx.*x;
    F=Inx.*fit_old+~Inx.*F;
    fit_old=F;
    x_old=x;
    
    if rand()<FADs
        U=rand(N,D)<FADs;
        x=x+CF*((Xmin+rand(N,D).*(Xmax-Xmin)).*U);
    else
        
        r=rand();
        Rs=size(x,1);
        stepsize=(FADs*(1-r)+r)*(x(randperm(Rs),:)-x(randperm(Rs),:));
        x=x+stepsize;
    end
    
     %% DE
    for ii=1:1:N
        x_vecT=x(ii,:);
        vect=randperm(N);
        vect(vect==ii)=[];
        ri1=vect(1);
        ri2=vect(2);
        ri3=vect(3);
        
        %% Mutation
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=x(ri1,:)+beta.*(x(ri2,:)-x(ri3,:));
        y = max(y, lb(1,:));
        y = min(y, ub(1,:));
        %% Crossover
        z_c=zeros(size(x_vecT));       %Crossover vector
        J0=randi([1 length(x_vecT)]);
        for jj=1:length(x_vecT)
            if jj==J0 || rand<=pCR
                z(jj)=y(jj);
            else
                z(jj)=x_vecT(jj);
            end
        end
        
        Fz= fobj(z,D);
        if  Fz<F_DE(ii)
            x(ii,:)=z;
            F_DE(ii)=Fz;
            if F_DE(ii)<F_g_bestDE
                F_g_bestDE=F_DE(ii);
                g_best_DE=x(ii,:);
            end
        end
    end
    
    for kk=1:1:D
        x(:,kk)=min(x(:,kk),ub(:,kk));
        x(:,kk)=max(x(:,kk),lb(:,kk));
    end
    
    %% DE

    %% Elimination Mechanism
    percent=0.2;
    [value, INDEX] = sort(F_DE,'ascend');
    for iiu=1:1:length(INDEX)
        x_sorted(iiu,:)=x(INDEX(iiu),:);
    end
    
    x=x_sorted;
    N_sort=round(N*rand*(percent));  % Number of population those will be eliminated from the group
    N_sort_Index=N-N_sort+1;     % Staring index number for elimination 
    
    if length(LB)==1
        for kk=1:1:D
            lb_Ele(1:N_sort,kk)=LB;
            ub_Ele(1:N_sort,kk)=UB;
        end
    end
    
    if length(LB)~=1
        for kk=1:1:D
            lb_Ele(1:N_sort,kk)=LB(kk);
            ub_Ele(1:N_sort,kk)=UB(kk);
        end
    end
        
    x_born=lb_Ele+(ub_Ele-lb_Ele).*rand(N_sort,D);
    
     
    niiv=1;
    for iiv=N_sort_Index:1:N
        x(iiv,:)=x_born(niiv,:);
        iiv=iiv;
        niiv=niiv;        
        niiv=niiv+1;
    end
    
    clear lb_Ele;
    clear ub_Ele;
    clear x_born; 
    clear iiv; 
    
    for kk=1:1:D
        x(:,kk)=min(x(:,kk),ub(:,kk));
        x(:,kk)=max(x(:,kk),lb(:,kk));
    end
       
    %% Elimination Mechanism 
    
    
 Fgbest_vect(it)=Top_predator_fit;
 F_g_bestDE_vect(it)=F_g_bestDE;
 
       
end


gbest=Top_predator_pos
Fg_best=Top_predator_fit
F_g_bestDE=F_g_bestDE
plot (Fgbest_vect)
hold on
%plot(F_g_bestDE_vect)
%grid on



break_point=1;













