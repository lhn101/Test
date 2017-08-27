function model=matrices_construction(model,data)

dt_ref=data.dt_ref; %Refence time step length

%% Block Component definition

%#11 Local level
model.components.idx{11}='LL';
LL.A=@(p,t,dt) 1;
LL.pA=[];
LL.pA0=[];
LL.C=@(p,t,dt) 1;
LL.pC=[];
LL.pC0=[];
LL.I=@(p,t,dt) 0;
LL.pI=[];
LL.pI0=[];
LL.Q=@(p,t,dt) p^2*dt/dt_ref;
LL.pQ={'\sigma_w','LL',[],[],[nan,nan]};
LL.pQ0={'0*nanstd(data.values(:,obs))'};
LL.x={'x^{LL}',[],[]};
LL.init={'nanmean(data.values(:,obs))','(2*nanstd(data.values(:,obs)))^2'};

%#12 Local trend
model.components.idx{12}='LT';
LT.A=@(p,t,dt) [1 dt;0 1];
LT.pA=[];
LT.pA0=[];
LT.C=@(p,t,dt) [1 0];
LT.pC=[];
LT.pC0=[];
LT.I=@(p,t,dt) [0 0];
LT.pI=[];
LT.pI0=[];
LT.Q=@(p,t,dt) p^2*dt/dt_ref*[dt^3/3 dt^2/2;dt^2/2 dt];
LT.pQ={'\sigma_w','LT',[],[],[0,inf]};
LT.pQ0={'1E-7*nanstd(data.values(:,obs))'};
LT.x={'x^{LL}',[],[];'x^{LT}',[],[]};
LT.init={'[nanmean(data.values(:,obs)) 0]','[(2*nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2]'};

%#13 Local acceleration
model.components.idx{13}='LA';
LA.A=@(p,t,dt)[1 dt dt^2;0 1 dt;0 0 1];
LA.pA=[];
LA.pA0=[];
LA.C=@(p,t,dt)[1 0 0];
LA.pC=[];
LA.pC0=[];
LA.I=@(p,t,dt) [0 0 0];
LA.pI=[];
LA.pI0=[];
LA.Q=@(p,t,dt) p^2*dt/dt_ref*[dt^5/20 dt^4/8 dt^3/6;dt^4/8 dt^3/3 dt^2/2;dt^3/6 dt^2/2 dt];
LA.pQ={'\sigma_w','LA',[],[],[0,inf]};
LA.pQ0={'1E-8*nanstd(data.values(:,obs))'};
LA.x={'x^{LL}',[],[];'x^{LT}',[],[];'x^{LA}',[],[]};
LA.init={'[nanmean(data.values(:,obs)) 0 0]','[(2*nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2]'};

%#21 Local level compatible with LT
model.components.idx{21}='LcT';
LcT.A=@(p,t,dt) [1 0;0 0];
LcT.pA=[];
LcT.pA0=[];
LcT.C=@(p,t,dt) [1 0];
LcT.pC=[];
LcT.pC0=[];
LcT.I=@(p,t,dt) [0 0];
LcT.pI=[];
LcT.pI0=[];
LcT.Q=@(p,t,dt) p^2*dt/dt_ref*[1 0;0 1E-15/(p^2*dt/dt_ref)]; % the last diagonal term needs to be set to 1E-15 for UD-filter
LcT.pQ={'\sigma_w','LcT',[],[],[0,inf]};
LcT.pQ0={'1E-7*nanstd(data.values(:,obs))'};
LcT.x={'x^{LL}',[],[];'x^{LTc}',[],[]};
LcT.init={'[nanmean(data.values(:,obs)) 0]','[(2*nanstd(data.values(:,obs)))^2 (1E-6*nanstd(data.values(:,obs)))^2]'};

%#22 Local level compatible with LA
model.components.idx{22}='LA';
LcA.A=@(p,t,dt) [1 0 0;0 0 0; 0 0 0]';
LcA.pA=[];
LcA.pA0=[];
LcA.C=@(p,t,dt) [1 0 0];
LcA.pC0=[];
LcA.pC=[];
LcA.I=@(p,t,dt) [0 0 0];
LcA.pI=[];
LcA.pI0=[];
LcA.pA0=[];
LcA.Q=@(p,t,dt) p^2*dt/dt_ref*[1 0 0;0 0 0; 0 0 1E-15/(p^2*dt/dt_ref)];
LcA.pQ={'\sigma_w','LcA',[],[],[0,inf]};
LcA.pQ0={'1E-7*nanstd(data.values(:,obs))'};
LcA.x={'x^{L}',[],[];'x^{LTc}',[],[];'x^{LAc}',[],[]};
LcA.init={'[nanmean(data.values(:,obs)) 0 0]','[(2*nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2]'};

%#23 Local trend compatible with LA
model.components.idx{23}='TcA';
TcA.A=@(p,t,dt) [1 dt 0;0 1 0; 0 0 0];
TcA.pA=[];
TcA.pA0=[];
TcA.C=@(p,t,dt) [1 0 0];
TcA.pC=[];
TcA.pC0=[];
TcA.I=@(p,t,dt) [0 0 0];
TcA.pI=[];
TcA.pI0=[];
TcA.pC0=[];
TcA.Q=@(p,t,dt) p^2*dt/dt_ref*[dt^3/3 dt^2/2 0;dt^2/2 dt 0; 0 0 1E-15/(p^2*dt/dt_ref)];
TcA.pQ={'\sigma_w','TcA',[],[],[0,inf]};
TcA.pQ0={'1E-8*nanstd(data.values(:,obs))'};
TcA.x={'x^{LL}',[],[];'x^{LT}',[],[];'x^{LAc}',[],[]};
TcA.init={'[nanmean(data.values(:,obs)) 0 0]','[(2*nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2 (nanstd(data.values(:,obs)))^2]'};

%#31 Periodic component
model.components.idx{31}='PD';
PD.A=@(p,t,dt) [cos(2*pi*dt/p(1)) sin(2*pi*dt/p(1));-sin(2*pi*dt/p(1)) cos(2*pi*dt/p(1))];
PD.pA={'p','PD',[],[nan,nan]};
PD.pA0={'365.2422','1','365.2422/2','7'};
PD.C=@(p,t,dt) [1 0];
PD.pC=[];
PD.pC0=[];
PD.I=@(p,t,dt) [p 0];
PD.pI={'\phi','PD,I',[],[],[-inf,inf]};
PD.pI0={'0.01'};
PD.Q=@(p,t,dt) p^2*dt/dt_ref*[1 0;0 1];
PD.pQ={'\sigma_w','PD',[],[],[nan,nan]};
PD.pQ0={'0*nanstd(data.values(:,obs))'};
PD.x={'x^{S1}',[],[];'x^{S2}',[],[]};
PD.init={'[0,0]','[(2*nanstd(data.values(:,obs)))^2,(2*nanstd(data.values(:,obs)))^2]'};

%#41 Autoregressive component
model.components.idx{41}='AR';
AR.A=@(p,t,dt) p(1)^(dt/dt_ref);
AR.pA={'\phi','AR',[],[],[0,1]};
AR.pA0={'0.75'};
AR.C=@(p,t,dt) 1;
AR.pC=[];
AR.pC0=[];
AR.I=@(p,t,dt) p;
AR.pI={'\phi','AR,I',[],[],[-inf,inf]};
AR.pI0={'0.01'};
AR.Q=@(p,t,dt) p^2*dt/dt_ref;
AR.pQ={'\sigma_w','AR',[],[],[0,inf]};
AR.pQ0={'1E-1*nanstd(data.values(:,obs))'};
AR.x={'x^{AR}',[],[]};
AR.init={'0','(nanstd(data.values(:,obs)))^2'};


%#51 Dynamic regression with hidden component
if ~isfield(model.components,'nb_DH_p')
    model.components.nb_DH_p=5;
end
model.components.idx{51}='DH';
DH.A=@(p,t,dt) 1;
DH.pA=[];
DH.pA0=[];
DH.C=@(p,t,dt) hidden_component_WT(p,t);
DH.pC=[];
DH.pC0=[];
for i=1:model.components.nb_DH_p
    DH.pC=[DH.pC;{['ycp' num2str(i)],'DH' ,[],[],[-1,1]}];
    DH.pC0=[DH.pC0,{'0.5'}];
end
DH.pC=[DH.pC;{['x_0' num2str(i)],'DH' ,[],[],[-inf,inf]};{['d_x' num2str(i)],'DH' ,[],[],[nan,nan]}];
DH.pC0=[DH.pC0,{num2str(datenum('2017-1-1 00:00'))},{num2str(365.2224)}];
clear i
DH.I=@(p,t,dt) p;
DH.pI={'\phi','DH,I',[],[],[-inf,inf]};
DH.pI0={'0.01'};
DH.Q=@(p,t,dt) p^2*dt/dt_ref;
DH.pQ={'\sigma_w','DH',[],[],[nan,nan]};
DH.pQ0={'0'};
DH.x={'x^{DH}',[],[]};
DH.init={'0','(nanstd(data.values(:,obs)))^2'};

%#61 LSTM component
model.components.idx{61}='DL';
DL.A=@(p,t,dt) p(1)^(dt/dt_ref);
DL.pA={'\phi','DL',[],[],[0,1]};
DL.pA0={'0.75'};
DL.C=@(p,t,dt) 1;
DL.pC=[];
DL.pC0=[];
DL.I=@(p,t,dt) p;
DL.pI={'\phi','DL,I',[],[],[-inf,inf]};
DL.pI0={'0.01'};
DL.Q=@(p,t,dt) p^2*dt/dt_ref;
DL.pQ={'\sigma_w','DL',[],[],[0,inf]};
DL.pQ0={'1E-1*nanstd(data.values(:,obs))'};
DL.x={'x^{DL}',[],[]};
DL.init={'0','(nanstd(data.values(:,obs)))^2'};

%% Block initialization
param_properties={};        %Reference vector for parameter names in A,C,Q,R matrices
param_idx=1;                %Index for parameter count in A,C,Q,R matrices
x_count=zeros(model.nb_obs,100);
parameter=[];               %Default values for parameters
nb_param_class=1;           %number of parameter per model class
for class_from=1:model.nb_class  %Loop over each model class
    h_count=0;                   %Count of the number of periodic component
    A{class_from}=[];            %Matrices initialization
    R{class_from}=[];
    initX{class_from}=[];        %Initial state variable expected value
    initV{class_from}=[];        %Initial state variable variance
    hidden_states_names{class_from}={};     %Reference vector for hidden state variables names
    for class_to=1:model.nb_class%Loop over each model class
        Q{class_from}{class_to}=[];
        CI_block=cell(model.nb_obs);
        for obs=1:model.nb_obs   %Loop over each observation
            p_count=0;                   %Count of the number of periodic component
            nb_hidden_states{class_from}{obs}=0;
            C{class_from}{obs,obs}=[];
            C_block=[];
            for block=1:length(model.components.block{class_from}{obs}) %loop over each model block
                block_idx=model.components.block{class_from}{obs}(block);
                block_name=model.components.idx{block_idx};
                if class_from==class_to
                    x_count(obs,block_idx)=x_count(obs,block_idx)+1;    %count the number of components
                    if x_count(obs,block_idx)>1
                        count_label=num2str(x_count(obs,block_idx));
                    else
                        count_label=num2str([]);
                    end
                    name=eval([block_name '.x']);
                    for i=1:size(name,1)
                        name{i,2}=num2str(class_from);
                        name{i,3}=data.labels{obs};
                        hidden_states_names{class_from}=[hidden_states_names{class_from};name(i,:)];
                    end
                    % initial hidden state values
                    initX{class_from}=[initX{class_from},eval(eval([block_name '.init{1}']))];
                    initV{class_from}=[initV{class_from},eval(eval([block_name '.init{2}']))];
                    
                    % A - Transition matrix
                    if class_from>1&&model.components.const{class_from}{obs}(block)==1
                        A{class_from}=[A{class_from},',',[block_name '.A(p([' num2str(A_param_ref{1}{obs}{block}) ']),t,dt)']];
                    else
                        nb_param=eval(['size(' block_name '.pA,1)']);
                        if nb_param>0
                            p_idx=param_idx:param_idx-1+nb_param;
                        else
                            p_idx=[];
                        end
                        A{class_from}=[A{class_from},',',[block_name '.A(p([' num2str(p_idx) ']),t,dt)']];
                        if any(strcmp(eval([block_name '.pA']),'p'))
                            p_count=p_count+1;
                            name={'p' ,['PD' num2str(p_count)],[],[],eval([block_name '.pQ{5}'])};
                            p_value=eval([block_name '.pA0{' num2str(p_count) '}']);
                            if ~isempty(p_value)
                                parameter=[parameter;eval(p_value)];
                            end
                        else
                            name=eval([block_name '.pA']);
                            p_value=eval([block_name '.pA0']);
                            if ~isempty(p_value)
                                parameter=[parameter;eval(p_value{:})];
                            end
                        end
                        if ~isempty(name)
                            name{3}=[num2str(class_from)];
                            name{4}=data.labels{obs};
                        end
                        param_properties=[param_properties;name];
                        A_param_ref{class_from}{obs}{block}=p_idx;
                        param_idx=param_idx+nb_param;
                    end
                    
                    % Q - Model transition error covariance
                    if class_from>1&&model.components.const{class_from}{obs}(block)==1
                        Q{class_from}{class_to}=[Q{class_from}{class_to},',',[block_name '.Q(p([' num2str(Q_param_ref{1}{obs}{block}) ']),t,dt)']];
                        Q_param_ref{class_from}{obs}{block}=Q_param_ref{1}{obs}{block};
                    else
                        nb_param=eval(['size(' block_name '.pQ,1)']);
                        if nb_param>0
                            p_idx=param_idx:param_idx-1+nb_param;
                        else
                            p_idx=[];
                        end
                        Q{class_from}{class_to}=[Q{class_from}{class_to},',',[block_name '.Q(p([' num2str(p_idx) ']),t,dt)']];
                        if any(strcmp(eval([block_name '.pA']),'p'))
                            name={'\sigma_w',['PD' num2str(p_count)],[],[],eval([block_name '.pQ{5}'])};
                        else
                            name=eval([block_name '.pQ']);
                        end
                        if ~isempty(name)
                            p_value=eval(eval([block_name '.pQ0{:}']));
                            parameter=[parameter;p_value];
                            name{3}=[num2str(class_from)];
                            name{4}=data.labels{obs};
                        end
                        param_properties=[param_properties;name];
                        Q_param_ref{class_from}{obs}{block}=p_idx;
                        param_idx=param_idx+nb_param;
                    end
                    
                    % C - Observation matrix
                    if class_from>1&&model.components.const{class_from}{obs}(block)==1
                        C_block=[C_block,',',[block_name '.C(p([' num2str(C_param_ref{1}{obs}{block}) ']),t,dt)']];
                    else
                        nb_param=eval(['size(' block_name '.pC,1)']);
                        if nb_param>0
                            p_idx=param_idx:param_idx-1+nb_param;
                        else
                            p_idx=[];
                        end
                        C_block=[C_block,',',[block_name '.C(p([' num2str(p_idx) ']),t,dt)']];
                        name=eval([block_name '.pC']);
                        if ~isempty(name)
                            if any(strcmp(eval([block_name '.pC{1,1}']),'ycp1'))
                                for i=1:length(eval([block_name '.pC']))
                                    h_count=h_count+1;
                                    p_value=eval([block_name '.pC0{' num2str(h_count) '}']);
                                    parameter=[parameter;eval(p_value)];
                                    name{i,3}=[num2str(class_from)];
                                    name{i,4}=data.labels{obs};
                                end
                            else
                                p_value=eval([block_name '.pC0{' num2str(p_count) '}']);
                                parameter=[parameter;eval(h_value)];
                                name{3}=[num2str(class_from)];
                                name{4}=data.labels{obs};
                            end
                        end
                        param_properties=[param_properties;name];
                        C_param_ref{class_from}{obs}{block}=p_idx;
                        param_idx=param_idx+nb_param;
                    end
                    nb_hidden_states{class_from}{obs}=nb_hidden_states{class_from}{obs}+length(eval([block_name '.A(ones(1,1000),1,1)']));
                end
            end
            if class_from==class_to
                
                R{class_from}=[R{class_from},',p([' num2str(param_idx) '])^2'];
                C{class_from}{obs,obs}=[C{class_from}{obs,obs},',[' C_block(2:end) ']'];
                param_properties=[param_properties;{['\sigma_v'],[],[num2str(class_from)],[data.labels{obs}],[0,inf]}];
                parameter=[parameter;0.05*nanstd(data.values(:,obs))];
                param_idx=param_idx+1;
                
                % PCA inter-component regression
                if isfield(model.components,'PCA')
                    %PCA_coeff=model.components.PCA{obs};
                    for pca_idx=1:length(model.components.ic{obs})
                        PCA_reg_param{obs}(pca_idx)=param_idx;
                        param_properties=[param_properties;{['\phi'],[ data.labels{obs}(1) '|' 'PC' num2str(pca_idx) ],[num2str(class_from)],[data.labels{obs}],[-inf,inf]}];
                        parameter=[parameter;1];
                        param_idx=param_idx+1;
                    end
                    
                end
                
                % I - Inter-component regression
                for obs_idx=1:model.nb_obs
                    ic_idx=find([model.components.ic{obs_idx}]==obs);
                    if ~isempty(ic_idx)
                        %param_reg=0;
                        for block=1:length(model.components.block{class_from}{obs}) %loop over each model block
                            block_idx=model.components.block{class_from}{obs}(block);
                            block_name=model.components.idx{block_idx};
                            nb_param=eval(['size(' block_name '.pI,1)']);
                            if nb_param>0
                                %if param_reg==0
                                p_idx=param_idx:param_idx-1+nb_param;
                                param_reg=p_idx;
                                param_properties=[param_properties;{'\phi', [data.labels{obs_idx}(1) '|' data.labels{obs}(1) '(' block_name ')' ],[num2str(class_from)],[data.labels{obs}],[-inf,inf]}];
                                parameter=[parameter;eval(eval([block_name '.pI0{:}']))];
                                param_idx=param_idx+1;
                                %else
                                %    p_idx=param_reg;
                                %end
                            else
                                p_idx=[];
                            end
                            
                            if class_from>1&&model.components.const{class_from}{obs}(block)==1
                                C{class_from}{obs_idx,obs}=[C{class_from}{obs_idx,obs},',',[block_name '.I(p([' num2str(IC_param_ref{1}{obs}{block}) ']),t,dt)']];
                            else
                                C{class_from}{obs_idx,obs}=[C{class_from}{obs_idx,obs},',',[block_name '.I(p([' num2str(p_idx) ']),t,dt)']];
                            end
                            if isfield(model.components,'PCA')
                                C{class_from}{obs_idx,obs}=[C{class_from}{obs_idx,obs},'*','sum(p([' num2str(PCA_reg_param{obs_idx}) ']).*model.components.PCA{' num2str(obs_idx) '}(' num2str(obs_idx) ',:)'')' ];
                            end
                            IC_param_ref{class_from}{obs}{block}=p_idx;
                            
                            %if eval([block_name '.I(parameter([' num2str(p_idx-1) ']))'])~=0&&class_from==1%||model.components.const{class_from}{obs}(block)==0
                            
                            %end
                        end
                    elseif obs_idx~=obs
                            C{class_from}{obs_idx,obs}=[',' num2str(zeros(1,nb_hidden_states{class_from}{obs}))];
                    end
                end
            end
        end
        if class_from==class_to
            for obs=1:model.nb_obs  %Loop over each observation
                Cc{class_from}{obs,1}=[',[' cell2mat(C{class_from}(obs,:)) ']'];
            end
            model.A{class_from}=eval(['@(p,t,dt) blkdiag(' A{class_from}(2:end) ')']);
            model.C{class_from}=eval(['@(p,t,dt) reshape([' ([Cc{class_from}{:}]) '],[' num2str(sum([nb_hidden_states{class_from}{:}])) ','  num2str(model.nb_obs) '])''']);
            model.R{class_from}=eval(['@(p,t,dt) blkdiag(' R{class_from}(2:end) ')']);
            model.Q{class_from}{class_to}=eval(['@(p,t,dt) blkdiag(' Q{class_from}{class_to}(2:end) ')']);
        end
    end
end

%% Q_ij model process noise transition errors
for class_from=1:model.nb_class  %Loop over each model class
    for class_to=setdiff([1:model.nb_class],class_from)%Loop over each model class
        Q{class_from}{class_to}=[];
        if class_from~=class_to
            for obs=1:model.nb_obs   %Loop over each observation
                for block=1:length(model.components.block{class_from}{obs}) %loop over each model block
                    block_idx=model.components.block{class_from}{obs}(block);
                    block_name=model.components.idx{block_idx};
                    % Q - Model transition error covariance
                    if block_idx>10&&block_idx<30
                        nb_param=size(eval([block_name '.Q(1,1,1)']),1);
                        if nb_param>0
                            p_idx=param_idx:param_idx-1+nb_param;
                        else
                            p_idx=[];
                        end
                        if class_from>1&&model.components.const{class_from}{obs}(block)==1
                            Q{class_from}{class_to}=[Q{class_from}{class_to},',',[block_name '.Q(p([' num2str(Q_param_ref{1}{obs}{block}) ']).^2,t,dt)']];
                        else
                            Q{class_from}{class_to}=[Q{class_from}{class_to},',',['diag(p([' num2str(p_idx) ']).^2)']];
                            
                            for i=1:nb_param
                                name=[eval([block_name '.pQ'])];
                                if ~isempty(name)
                                    if nb_param>1
                                        name{1}=[name{1},'(' num2str(i),num2str(i) ')'];
                                    end
                                    parameter=[parameter;1E3*eval(eval([block_name '.pQ0{:}']))];
                                    name{3}=[num2str(class_from) num2str(class_to)];
                                    name{4}=data.labels{obs};
                                end
                                param_properties=[param_properties;name];
                            end
                            param_idx=param_idx+nb_param;
                        end
                    else
                        Q{class_from}{class_to}=[Q{class_from}{class_to},',',[block_name '.Q(p(' num2str(Q_param_ref{class_from}{obs}{block}) '),t,dt)']];
                    end
                end
            end
        end
        if class_from~=class_to
            model.Q{class_from}{class_to}=eval(['@(p,t,dt) blkdiag(' Q{class_from}{class_to}(2:end) ')']);
        end
    end
end

%% Initial nanmean vector and covariance matrix

if ~isfield(model,'initX')
    for class_from=1:model.nb_class
        model.initX{class_from}=initX{class_from}';
    end
end
if ~isfield(model,'initV')
    for class_from=1:model.nb_class
        model.initV{class_from}=diag(initV{class_from});
    end
end
if ~isfield(model,'initS')
    for class_from=1:model.nb_class
        model.initS{class_from}=1/model.nb_class;
    end
end

%% Transition matrix
model.Z=[];
if model.nb_class==1
    model.Z=@(p,t,dt) 1;
elseif model.nb_class==2
    model.Z=@(p,t,dt) eval(['[p(' num2str(param_idx) ') 1-p(' num2str(param_idx) ') ;1-p(' num2str(param_idx) '+1) p(' num2str(param_idx) '+1) ]^(dt/' num2str(dt_ref) ')']);
    param_properties=[param_properties;{'Z(11)',[],'11',[],[0,1]};{'Z(22)',[],'22',[],[0,1]}];
    parameter=[parameter;1-1/(data.dt_ref*365*10);1-1/(data.dt_ref*365*10)];
    param_idx=param_idx+2;
    
else
    for i=1:model.nb_class
        for j=i:model.nb_class
            param_idx=param_idx+1;
            model.Z=[model.Z;{ 'p(param_idx)'}];
            param_properties=[param_properties;{'\pi_z',[],[num2str(i),num2str(j)], [] ,[nan,nan]}];
            parameter=[parameter;0.5];
        end
    end
    model.Z=@(p,t,dt) [reshape(model.Z),[model.nb_class,model.nb_class]]^(dt/dt_ref);
end
if ~isfield(model,'parameter')
    model.parameter=parameter;
end
if ~isfield(model,'p_ref')
    model.p_ref=1:length(parameter);
end
if ~isfield(model,'param_properties')
    model.param_properties=param_properties;
end
if ~isfield(model.components,'PCA')
    model.components.PCA=cell(1,model.nb_obs);
end
model.hidden_states_names=hidden_states_names;
