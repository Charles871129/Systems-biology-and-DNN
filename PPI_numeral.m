clear,clc;
%% loading data
load('PPI_EXP_healthy.mat'); load('PPI.mat'); load('PPI_name.mat');
interaction = PPI; gene_expre = PPI_EXP_healthy; Name = PPI_name;
for i = 1 : size(PPI,1)
    interaction(i,i) = 0;
end

%% parameters for ID
patients = size(PPI_EXP_healthy,2);  % number of patients
Final_reg_ability = double(zeros(size(interaction),'like',interaction));
basal = zeros(size(interaction,1),1);
[BB,index] = sort(sum(interaction,2));
id = find(BB~=0,1,'first');   %you can use for the start point

fprintf('start PPI ID:\n');
accelerate = 1:250;  %used when forward model selection
All = 1:size(interaction,2);  % all binds (used for diff to get the genes that don't bind with ith gene)
Ones = ones(patients,1);  %used for phi in linear regression

%----parameters for lsqlin (algorithm : 'interior-point'(default) )----
options = optimset('Display','off');
% default : 'MaxIter',200,'TolCon',1e-8,'TolFun',1e-8,'TolX',1e-12

%% Start (or continue. To continue, you need to load the saved variables and enter the lastest j.)
total_time = tic;
start = 1; % start order  (or enter the last j to continue)
for j = start:8800
    tic % start time
    flag = 0;
    i = index(j);  % index of the instant protein
    remaining = size(gene_expre,1)-j;  % the amount of protein left without processing (only be used to print information)
    bind1 = find(interaction(i,:));  % ppi = 1 with the instant (ith) protein
    % initial AIC_value
    AIC_value_initial = 90000;
    
    %----linear regression (with all known genes interacting with ith gene)----
    X = gene_expre(i,:)';
    if ~any(X) || isempty(bind1)
        continue
    end
    phi_temp = gene_expre';
    phi_temp = cat(2,phi_temp.*X,Ones);  % [PPI 1]
    % above variables would skip in the following linear regression
    
    % Foward
    Bind1 = {}; Theta = {}; AIC_value_pack = {};
    binding = bind1;
    if length(bind1) < 250  %to accelerate, here skip some combinations of parameter
        tmp = 1:length(bind1);
    else
        tmp = accelerate;
    end
    for m = tmp  %Originally, it should be 1:length(bind1), but it's too long to calculate the result.
        AIC_value_temp = AIC_value_initial;
        for u = 1:length(binding)
            if m == 1
                bindtemp = binding(u);
            else
                bindtemp = cat(2,Bind1{m-1},binding(u));
            end

            %----linear regression forward----
            phi = phi_temp(:,[sort(bindtemp), 20039]);
            phi_len = size(phi,2);

            %----parameters for lsqlin----
            theta_temp = mldivide(phi,X);
            % theta_temp = matlab.internal.math.nowarn.mldivide(phi,X);
            resnorm_temp = sum((phi*theta_temp-X).^2);
            
            %----linear regression end and calculate AIC ----
            pre_AIC_value = log( resnorm_temp/patients ) + 2*phi_len/patients;  %AIC

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
                phi = phi_temp(:,[sort(bindtemp), 20039]);
                phi_len = size(phi,2);

                %----parameters for lsqlin----
                theta_temp = mldivide(phi,X);
                % theta_temp = matlab.internal.math.nowarn.mldivide(phi,X);
                resnorm_temp = sum((phi*theta_temp-X).^2);
                
                %----linear regression end  and calculate AIC----
                pre_AIC_value = log( resnorm_temp/patients ) + 2*phi_len/patients;  %AIC

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
    
    if j < 11000
        if mod(j,100) == 0
            save PPI_ID_TEMP;
        end
    elseif j < 17330
         if mod(j,30) == 0
            save PPI_ID_TEMP;
         end
    else
         if mod(j,5) == 0
            save PPI_ID_TEMP;
         end
    end
             
    
    if min_index == 1
        flag2 = 'F';
    else
        flag2 = 'B';
    end
    t = toc; % start time
    fprintf('j= %d(%d)(%d -> %d nodes)(index: %d,%d,%c), time: %.0f(s)\n', j,remaining,length(bind1),length(bindout),index1,min_index,flag2,t)
end
toc(total_time);
A = full(Final_reg_ability);
fprintf('Done\n')
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))
toc % elapsed time
%%
PPI_PNP = A;
PPI_Edge = (A+A')./2;
Basal_PPI = table(Name,basal);
save('Result/PPI_healthy_start/PPI_PNP_start.mat','PPI_PNP', '-v7.3')
save('Result/PPI_healthy_start/PPI_Edge_start.mat','PPI_Edge', '-v7.3')
writetable(Basal_PPI,'Result/PPI_healthy_start/Basal_PPI_start.txt')