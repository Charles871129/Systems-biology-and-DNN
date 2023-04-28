clear all, clc;
load data: PPI and GRNname and regulation ability from PPI and GRN
load('PPI_name.mat');  load('GRN_name.mat');
load('PPI_PNP');  load('GRN_PNP');
%% 

ppi = PPI_PNP;
grn = GRN_PNP;


PNP_up_name = GRN_name;
PNP_down_name = [PPI_name;GRN_name];
save('PNP_name.mat','PNP_up_name', 'PNP_down_name');


PPI = sparse(size(ppi,1),size(grn,2));
PPI(1:size(ppi,1),1:size(ppi,2)) = ppi;
GRN = grn;

%%  Construct the network matrix
% M = sparse(length(ppi)+18263+length(miRNA)+length(lncRNA)+11,length(ppi)+length(miRNA)+length(lncRNA));
% M(1:size(PPI,1),1:size(PPI,2)) = PPI;

M = cat(1,PPI,GRN);

A = M;
lower=-5*10^(2); upper=5*10^(2);
M_low = find(M<=lower);
M_up = find(M>=upper);
M(M_low) = lower;  M(M_up) = upper;
tic
[U, S, V]=svd(full(M)); % U1=mxm S1=mxn V1=nxn
toc
Sdiag = diag(S);
Sdiag = Sdiag.^2;

e = 0;
i = 1;
e85 = 0.85*sum(Sdiag);
while(e < e85)
    e = e+Sdiag(i);
    i = i+1;
end

% Rank
V = V(:,1:i);  downstream = A*V;
U = U(:,1:i);  upstream = A'*U;      
downstream(isnan(downstream))=0;
down = sum(downstream.^2,2);
upstream(isnan(upstream))=0;
up = sum(upstream.^2,2);

save('stream_85.mat','downstream','upstream', 'down','up','-v7.3')
% save('Result/BB1.mat','BB1', '-v7.3')
% save('Result/BB2.mat','BB2', '-v7.3')
% save('Result/Gup.mat','Gup', '-v7.3')
% save('Result/Gdn.mat','Gdn', '-v7.3')
