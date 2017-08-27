%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program : State_estimation
% Author : James-A. Goulet
% Date : Feb. 3th 2015
% Last update : May. 18th 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function estim=state_estimation(data,model,option,varargin)

estim.smooth=0;
disp_flag=1;

args = varargin;
for t=1:2:length(args)
    switch args{t}
        case 'smooth', estim.smooth = args{t+1};
        case 'disp_flag', disp_flag = args{t+1};
        otherwise, error(['unrecognized argument ' args{t}])
    end
end
if disp_flag==1
    disp('    ...in progress')
    disp(' ')
end

%% Kalman/UD filter
[estim.x_M, estim.V_M, estim.VV_M, estim.S, estim.LL,estim.U,estim.D]=SKF(data,model,option);

if disp_flag==1
    disp(['    -> log-likelihood:    ' num2str(estim.LL)])
end

%% Kalman smoother
if estim.smooth==1
    [estim.x_M, estim.V_M, estim.VV_M, estim.S]=RTS_SKS(estim,data,model);
end

%% Collapse multiple model classes in a single one
estim.x=zeros(size(estim.x_M{1}));
estim.V=zeros(size(estim.x_M{1}));
estim.y=zeros(model.nb_obs,data.nb_steps);
estim.Vy=zeros(model.nb_obs,data.nb_steps);

mx=zeros(size(estim.x_M{1},1),model.nb_class);
for t=1:data.nb_steps
    for j=1:model.nb_class
        estim.x(:,t)=estim.x(:,t)+estim.S(t,j)*estim.x_M{j}(:,t);
        mx(:,j)=estim.x_M{j}(:,t);
        estim.V(:,t)=estim.V(:,t)+estim.S(t,j)*diag(estim.V_M{j}(:,:,t));
        
        C_j=model.C{j}(model.parameter,data.timestamps(t),data.dt_steps(t));
        estim.y(:,t)=estim.y(:,t)+estim.S(t,j)*(C_j*estim.x_M{j}(:,t));
        my(:,j)=(C_j*estim.x_M{j}(:,t));
        estim.Vy(:,t)=estim.Vy(:,t)+estim.S(t,j)*diag(C_j*estim.V_M{j}(:,:,t)*C_j');
    end
    for j=1:model.nb_class
        estim.V(:,t)=estim.V(:,t)+estim.S(t,j)*diag((mx(:,j)-estim.x(:,t))*(mx(:,j)-estim.x(:,t))');
        estim.Vy(:,t)=estim.Vy(:,t)+estim.S(t,j)*diag((my(:,j)-estim.y(:,t))*(my(:,j)-estim.y(:,t))');
    end
end

if disp_flag==1
    disp('    \\done')
end
