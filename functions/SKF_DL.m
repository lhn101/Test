function [x, V, VV, S, loglik,U,D] = SKF_DL(data,model,option)
% INPUTS:
% y(:,t)   - the observation at time t
% A - the system matrix
% C - the observation matrix
% Q - the system covariance
% R - the observation covariance
% Z - Transition probabilities for model classes
% init_x - the initial state (column) vector
% init_V - the initial state covariance
% init_S - the initial model class probabilities
%
% OUTPUTS (where X is the hidden state being estimated)
% x(:,t) = E[X(:,t) | y(:,1:t)]
% V(:,:,t) = Cov[X(:,t) | y(:,1:t)]
% loglik = sum{t=1}^T log P(y(:,t))
%

T = data.nb_steps;            %Nunber of time stamps
M = model.nb_class;           %Number of model classes

if strcmp(option.method,'UD')
    error('UD method not supported')
elseif strcmp(option.method,'kalman')
    U=[];
    D=[];
else
    error('estimation method not recognized | SKF.m')
end

%% Initialization
x = cell(1,M);
V = cell(1,M);
VV = cell(1,M);
x_ij=cell(1,M);
V_ij=cell(1,M);
VV_ij=cell(1,M);

%% Preallocate the matrices A, C, R, Z, Q
A=cell(M,T);
C=cell(M,T);
R=cell(M,T);
Z=cell(1,T);
Q=cell(M,M,T);

for j=1:M
    ss=size(model.hidden_states_names{1},1);
    x{j} = zeros(ss,T);
    V{j} = zeros(ss, ss,T);
    VV{j} = zeros(ss, ss,T);
    
    x_ij{j} = zeros(ss,M);
    V_ij{j} = zeros(ss,ss,M);
    VV_ij{j} = zeros(ss,ss,M);
end

S = zeros(T,M);
W=zeros(M,M);
LL=zeros(M,M);

loglik = 0;
%% Estimate state for each t
for t=1:T
    log_S_marginal = zeros(M,M);
    lik_merge=0;
    for j=1:M       % transition model
        if (t==1|(data.dt_steps(t)~=data.dt_steps(1:t-1)))
            A{j,t} = model.A{j}(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
            C{j,t} = model.C{j}(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
            R{j,t} = model.R{j}(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
            Z{t} = model.Z(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
        else
            idx=find(data.dt_steps(t)==data.dt_steps(1:t-1),1,'first');
            A{j,t}=A{j,idx};
            C{j,t} = model.C{j}(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
            R{j,t}=R{j,idx};
            Z{t}=Z{idx};
            
        end
        for i=1:M   % starting model
            if (t==1|(data.dt_steps(t)~=data.dt_steps(1:t-1)))
                Q{j,i,t} = model.Q{j}{i}(model.parameter(model.p_ref),data.timestamps(t),data.dt_steps(t));
            else
                idx=find(data.dt_steps(t)==data.dt_steps(1:t-1),1,'first');
                Q{j,i,t}=Q{j,i,idx};
            end
            
            if t==1
                prevX = model.initX{i};
                prevV = model.initV{i};
                prevS = model.initS{i};
            else
                prevX = x{i}(:,t-1);
                prevV = V{i}(:,:,t-1);
                prevS = S(t-1,i);
            end
            
            %% LSTM correction term
            nb_DL_samples=100;
            [x_DL,V_DL]=LSTM_DL_pred(prevX(model.DL.x_idx),prevV(model.DL.x_idx,model.DL.x_idx),option.name,model.DL.norm,nb_DL_samples);
            B_DL=zeros(length(prevX),1);
            B_DL(model.DL.x_idx)=x_DL;
            
            W_DL=zeros(length(prevX),length(prevX));
            W_DL(model.DL.x_idx,model.DL.x_idx)=x_DL;
            
            %% Kalman filter
            warning('off','all')
            [x_ij{j}(:,i), V_ij{j}(:,:,i), VV_ij{j}(:,:,i), LL(i,j,t)] = KF(A{j,t},C{j,t},Q{j,i,t},R{j,t}, data.values(t,:)', prevX, prevV,'B',B_DL,'W',W_DL);
            warning('on','all')
            
            if isnan(LL(i,j,t))
                LL(i,j,t)=1;
                loglik=-inf;
                disp(['warning: LL{' num2str(i) ',' num2str(j) '}(t=' num2str(t) ')=nan | SKF.m'])
                break
            end
            
            log_S_marginal(i,j)=LL(i,j,t) + log(Z{t}(i,j)) + log(prevS);
            lik_merge=lik_merge + exp(log_S_marginal(i,j));
            
            if any(any(isnan(VV_ij{j}(:,:,i))))
                disp(['warning: VV_ij{' num2str(i) ',' num2str(j) '}(t=' num2str(t) ')=nan | SKF.m'])
            end
        end
        if loglik==-inf
            disp(['warning: LL{' num2str(i) ',' num2str(j) '}(t=' num2str(t) ')=-inf | SKF.m'])
            break
        end
    end
    if loglik==-inf
        break
    end
    
    loglik=loglik+log(lik_merge);
    if loglik==-inf
        disp('warning: loglik=-inf | SKF.m')
        break
    elseif isnan(loglik)
        disp('warning: loglik=-inf | SKF.m')
    end
    
    S_marginal=exp(log_S_marginal);
    if any(any(S_marginal==0))
        S_marginal=exp(log_S_marginal+(299-max(max(log_S_marginal))));
    end
    if any(any(S_marginal==0))
        S_marginal(S_marginal==0)=1E-100;
    end
    S_norm = sum(sum(S_marginal));
    
    %% posterior for state j at time t
    S_marginal = S_marginal/S_norm;
    for j=1:M
        S(t,j) = sum(S_marginal(:,j));
    end
    S(t,S(t,:)==0)=1E-99;
    
    %% weights of state components
    for j=1:M
        for i=1:M
            W(i,j) = S_marginal(i,j)/S(t,j);
        end
    end
    W(W==0)=1E-99;
    
    if isnan(W(i,j))
        disp(['warning: W(i,j)(t=' num2str(t) ')=nan | SKF.m'])
    end
    
    %% approximate new continuous state
    for j=1:M
        x{j}(:,t) = x_ij{j}(:,:) * W(:,j);
        for i=1:M
            m = x_ij{j}(:,i) - x{j}(:,t);
            V{j}(:,:,t) = V{j}(:,:,t) + W(i,j)*(V_ij{j}(:,:,i) + m*m');
            VV{j}(:,:,t) = VV{j}(:,:,t) + W(i,j)*(VV_ij{j}(:,:,i) + m*m');
        end
    end
end



