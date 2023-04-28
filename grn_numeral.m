clear;clc;
%% loading data
load('GRN_EXP_unhealthy.mat');load('GRN.mat');load('GRN_name.mat');
interaction = GRN; gene_expre = GRN_EXP_unhealthy; Name = GRN_name;

% ######## you need to specify the columns of microRNA ########
micro = 13314:13482;
% #############################################################

%% parameters for ID
patients = size(gene_expre,2);  % n times
gene_num = size(gene_expre,1);  %number of gene
Final_reg_ability = double(zeros(size(interaction),'like',interaction));
basal = zeros(size(interaction,1),1);
[BB,index] = sort(sum(interaction,2));
id = find(BB~=0,1,'first');   %you can use for the start point

fprintf('start GRN ID:\n');
acc = 200;
accelerate = 1:acc;  %used when forward model selection
All = 1:size(interaction,2);  % all binds (used for diff to get the genes that don't bind with ith gene)
Ones = ones(patients,1);  %used for phi in linear regression
cons = gene_num+1;

%----parameters for lsqlin (algorithm : 'interior-point'(default) )----
options = optimset('Display','off');
lb_tmp = -1*Inf(gene_num+2,1);    %lower bound
ub_tmp = Inf(gene_num+2,1);    %upper bound
ub_tmp(micro) = zeros(length(micro),1);

% default : 'MaxIter',200,'TolCon',1e-8,'TolFun',1e-8,'TolX',1e-12


%% Start (or continue. To continue, you need to load the saved variables and enter the lastest j.)
start = id; % start order  (or enter the last j to continue)
total_time = tic;
for j = start:gene_num
    tic
    i = index(j); % index of the instant protein
    remaining = gene_num-j;  % the amount of protein left without processing (only be used to print information)
    bind1 = find(interaction(i,:));  % ppi = 1 with the instant (ith) protein
    % initial AIC_value
    AIC_value_initial = 90000;
	
    %----linear regression (with all known genes interacting with ith gene)----
    X = gene_expre(i,:)';
    if ~any(X) || isempty(bind1)
        continue
    end
    pi = gene_expre(i,:)';
    phi_temp = gene_expre';
	phi_temp(:,micro) = phi_temp(:,micro).*pi;
	phi_temp = cat(2,phi_temp,Ones);  % [PPI Pi(t) Pi(t) 1], Pi(t) are use to multiply the translation and degration rate
    % above variable would skip in the following linear regression
    
    
    % Foward
    Bind1 = {}; Theta = {}; AIC_value_pack = {};
    binding = bind1;
    if length(bind1) < acc  %to accelerate, here skip some combinations of parameter
        tmp = 1:length(bind1);
    else
        tmp = accelerate;
    end
    for m = tmp  
        AIC_value_temp = AIC_value_initial;
        for u = 1:length(binding)
            if m == 1
                bindtemp = binding(u);
            else
                bindtemp = cat(2,Bind1{m-1},binding(u));
            end
            
            %----linear regression forward----
            Sort = [sort(bindtemp),cons];
            phi = phi_temp(:,Sort);
            lb = lb_tmp(Sort);
            ub = ub_tmp(Sort);
            phi_len = size(phi,2);
            
            [theta_temp,resnorm_temp] = lsqlin(phi,X,[],[],[],[],lb,ub,[],options); % return theata that minimize the function and get the value(resnorm)
            %----linear regression end and calculate AIC ----
            pre_AIC_value = log( resnorm_temp/(patients-1) ) + 2*phi_len/(patients-1);  %AIC
            
            if u == 1 && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                AIC_value_temp = pre_AIC_value;
                bindinfo = bindtemp;
            elseif pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                AIC_value_temp = pre_AIC_value;
                bindinfo = bindtemp;
            end
        end
        Bind1{m} = bindinfo; Theta{m} = theta_temp; AIC_value_pack{m} = AIC_value_temp;  %record variable
        binding = setdiff(bind1,Bind1{m});
    end
    [~,min_index] = min(cell2mat(AIC_value_pack));
    
    %Forward for the Nodes that the number of interactions may be larger
    if min_index >= acc-30
        if length(bind1) < 170  %to accelerate, here skip some combinations of parameter
            tmp = m+1:length(bind1);
        else
            tmp = m+1:170;
        end
        for m = tmp
            AIC_value_temp = AIC_value_initial;
            for u = 1:length(binding)
                if m == 1
                    bindtemp = binding(u);
                else
                    bindtemp = cat(2,Bind1{m-1},binding(u));
                end

                %----linear regression forward----
                Sort = [sort(bindtemp),cons];
                phi = phi_temp(:,Sort);
                lb = lb_tmp(Sort);
                ub = ub_tmp(Sort);
                phi_len = size(phi,2);

                [theta_temp,resnorm_temp] = lsqlin(phi,X,[],[],[],[],lb,ub,[],options); % return theata that minimize the function and get the value(resnorm)
                %----linear regression end and calculate AIC ----
                pre_AIC_value = log( resnorm_temp/(patients-1) ) + 2*phi_len/(patients-1);  %AIC

                if u == 1 && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                elseif pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                end
            end
            Bind1{m} = bindinfo; Theta{m} = theta_temp; AIC_value_pack{m} = AIC_value_temp;  %record variable
            binding = setdiff(bind1,Bind1{m});
        end
        [~,min_index] = min(cell2mat(AIC_value_pack));
    end
    
    bindout = Bind1{min_index};
    index1 = min_index;
    Bind1 = {bindout}; Theta = {Theta{min_index}}; AIC_value_pack = {AIC_value_pack{min_index}};
    

    % Backward
    if length(bindout) > 1
        for y = 1:length(bindout)-2  
            AIC_value_temp = AIC_value_initial;
            for z = 1:length(bindout)
                bindtemp = bindout;
                bindtemp(z) = [];
                
                %----linear regression backward----
                Sort = [sort(bindtemp),cons];
                phi = phi_temp(:,Sort);
                lb = lb_tmp(Sort);
                ub = ub_tmp(Sort);
                phi_len = size(phi,2);
                
                [theta_temp,resnorm_temp] = lsqlin(phi,X,[],[],[],[],lb,ub,[],options);  % return theata that minimize the function and get the value(resnorm)
                %----linear regression end  and calculate AIC----
                pre_AIC_value = log( resnorm_temp/(patients-1) ) + 2*phi_len/(patients-1);  %AIC
                
                if z == 1 && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                elseif pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                end
            end
            Bind1{1+y} = bindinfo; Theta{1+y} = theta_temp; AIC_value_pack{1+y} = AIC_value_temp;
            bindout = bindinfo;
        end
    end
    
    [~,min_index] = min(cell2mat(AIC_value_pack));
    thetaout = Theta{min_index};
    bindout = Bind1{min_index};
	
    
    basal(i) = thetaout(end); % basal level for epigenetic regulation
    Final_reg_ability(i,bindout) = thetaout(1:end-1);
    
    
%     if j < 15000
%         if mod(j,3000) == 0
%             save GRN_ID_TEMP_;
%         end
%     elseif j < 23500
%         if mod(j,5) == 0
%             save GRN_ID_TEMP_;
%         end
%     else
%         save GRN_ID_TEMP_;
%     end

    if min_index == 1
        flag2 = 'F';
    else
        flag2 = 'B';
    end
    t = toc; % start time
    fprintf('j= %d(%d)(%d -> %d nodes)(index: %d,%d,%c), time: %.0f(s) \n', j,remaining,length(bind1),length(bindout),index1,min_index,flag2,t);
end
fprintf('Done\n');
save GRN_ID_patient_finish;
toc(total_time)  % elapsed time
A = full(Final_reg_ability);
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))

GRN_PNP = A;
GRN_Edge = A;
Basal_GRN = table(Name,basal);
save('Result/GRN_unhealthy/GRN_PNP.mat','GRN_PNP', '-v7.3')
save('Result/GRN_unhealthy/GRN_Edge.mat','GRN_Edge', '-v7.3')
writetable(Basal_GRN,'Result/GRN_healthy/Basal_GRN.txt')