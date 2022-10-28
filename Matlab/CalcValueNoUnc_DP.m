function vo = CalcValueNoUnc_DP(maxp,minp,c, P, eta, vi, ed, iC, iD)
%%
% Title: Calculate Risk-Neutral value function using deterministic price
% Inputs:
%   maxp - discharge price
%   minp - charge price
%   c - marginal discharge cost
%   P - power rating w.r.t to energy rating and sampling time, i.e., a
%   1MW/2MWh battery with 5min resolution -> P = 0.5/12 
%   eta - efficiency
%   vi - input value function for the next time period, which equals to
%   v_t(e) where e is sampled from 0 to 1 at the granularity e
%   ed - granularity at which vi is sampled, in p.u. to energy rating
% 
% Outputs:
%   vo - value function for the current time period sampled at ed

%% 

% create samples of e
% es = (0:ed:1)'; 
% Ne = numel(es); % number of samples

% add a large number of upper and lower v, where the first point is
% v_t(0-) = +infty, and the second point is v_t(0), the second largest one is
% v_t(1), and the largest one is v_t(1+) = -infty
lNum = 1e5;
v_foo = [lNum;vi;-lNum];

% % calculate soc after charge vC = (v_t(e+P*eta))
% eC = es + P*eta; 
% % round to the nearest sample 
% iC = ceil(eC/ed)+1;
% iC(iC > (Ne+1)) = Ne + 2;
% iC(iC < 2) = 1;
vC = v_foo(iC);


% % calculate soc after discharge vC = (v_t(e-P/eta))
% eD = es - P/eta; 
% % round to the nearest sample 
% iC = floor(eD/ed)+1;
% iC(iC > (Ne+1)) = Ne + 2;
% iC(iC < 2) = 1;
vD = v_foo(iD);


% calculate CDF and PDF
FtCC_min = (vC*eta > minp); % F_t(v_t(e+P*eta)*eta)
FtEC_min = (vi*eta > minp); % F_t(v_t(e)*eta)
% FtED_min = ((vi/eta + c).*((vi/eta + c) > 0) > minp); % F_t(v_t(e)/eta + c) 
% FtDD_min = ((vD/eta + c).*((vD/eta + c) > 0) > minp); % F_t(v_t(e-P/eta)/eta + c)

% FtCC_max = (vC*eta > maxp); % F_t(v_t(e+P*eta)*eta)
% FtEC_max = (vi*eta > maxp); % F_t(v_t(e)*eta)
FtED_max = ((vi/eta + c).*((vi/eta + c) > 0) > maxp); % F_t(v_t(e)/eta + c) 
FtDD_max = ((vD/eta + c).*((vD/eta + c) > 0) > maxp); % F_t(v_t(e-P/eta)/eta + c) 

% IFtCCmin = find(FtCC_min==1,1,'last'); % find the last SoC index 
% IFtECmin = find(FtEC_min==1,1,'last'); % find the last SoC index
% IFtEDmax = find(FtED_max==1,1,'last');
% IFtDDmax = find(FtDD_max==1,1,'last');

IFtCCmin = nnz(FtCC_min); % find the last SoC index 
IFtECmin = nnz(FtEC_min); % find the last SoC index
IFtEDmax = nnz(FtED_max);
IFtDDmax = nnz(FtDD_max);

%% Method 2
% Term1 = vC .* FtCC_min;
% Term2 = (minp/eta) .* (FtEC_min - FtCC_min);
% Term3 = vi .* (FtED_max - FtEC_min);
% Term4 = (maxp-c) * eta .* (FtDD_max - FtED_max);
% Term6 = vD .* (1-FtDD_max);

%% Method 1
if IFtECmin <= IFtEDmax
    Term1 = vC .* FtCC_min;
    Term2 = (minp/eta) .* (FtEC_min - FtCC_min);
    Term3 = vi .* (FtED_max - FtEC_min);
    Term4 = (maxp-c) * eta .* (FtDD_max - FtED_max);
    Term6 = vD .* (1-FtDD_max);
elseif (IFtECmin <= IFtDDmax) && (IFtECmin > IFtEDmax)
    if IFtCCmin <= IFtEDmax
        Term1 = vC .* FtCC_min;
        Term2 = (minp/eta) .* (FtED_max - FtCC_min);
        Term3 = vi .* (FtEC_min - FtED_max);
        Term4 = (maxp-c) * eta .* (FtDD_max - FtEC_min);
        Term6 = vD .* (1-FtDD_max);
    else
        Term1 = vC .* FtED_max;
        Term2 = vi .* (FtCC_min - FtED_max);
        Term3 = vi .* (FtEC_min - FtCC_min);
        Term4 = (maxp-c) * eta .* (FtDD_max - FtEC_min);
        Term6 = vD .* (1-FtDD_max);
    end
elseif IFtECmin > IFtDDmax
    if IFtCCmin > IFtDDmax
        Term1 = vC .* FtED_max;
        Term2 = vi .* (FtDD_max - FtED_max);
        Term3 = vi .* (FtCC_min - FtDD_max);
        Term4 = vi .* (FtEC_min - FtCC_min);
        Term6 = vD .* (1-FtEC_min);
    elseif (IFtCCmin > IFtEDmax) && (IFtCCmin <= IFtDDmax)
        Term1 = vC .* FtED_max;
        Term2 = vi .* (FtCC_min - FtED_max);
        Term3 = vi .* (FtDD_max - FtCC_min);
        Term4 = vi .* (FtEC_min - FtDD_max);
        Term6 = vD .* (1-FtEC_min);
    elseif IFtCCmin <= IFtEDmax
        Term1 = vC .* FtCC_min;
        Term2 = (minp/eta) .* (FtED_max - FtCC_min);
        Term3 = vi .* (FtDD_max - FtED_max);
        Term4 = vi .* (FtEC_min - FtDD_max);
        Term6 = vD .* (1-FtEC_min);
    end
end

%% Method 3
% if IFtECmin <= IFtEDmax
%     Term1 = vC .* FtCC_min;
%     Term2 = (minp/eta) .* (FtEC_min - FtCC_min);
%     Term3 = vi .* (FtED_max - FtEC_min);
%     Term4 = (maxp-c) * eta .* (FtDD_max - FtED_max);
%     Term6 = vD .* (1-FtDD_max);
% elseif (IFtECmin <= IFtDDmax) && (IFtECmin > IFtEDmax)
%     if IFtCCmin <= IFtEDmax
%         Term1 = vC .* FtCC_min;
%         Term2 = (minp/eta) .* (FtED_max - FtCC_min);
%         Term3 = vi .* (FtEC_min - FtED_max);
%         Term4 = (maxp-c) * eta .* (FtDD_max - FtEC_min);
%         Term6 = vD .* (1-FtDD_max);
%     else
%         Term1 = vC .* FtED_max;
%         Term2 = (minp/eta) .* (FtCC_min - FtED_max);
%         Term3 = vi .* (FtEC_min - FtCC_min);
%         Term4 = (maxp-c) * eta * (FtDD_max - FtEC_min);
%         Term6 = vD .* (1-FtDD_max);
%     end
% elseif IFtECmin > IFtDDmax
%     if IFtCCmin > IFtDDmax
%         Term1 = vC .* FtED_max;
%         Term2 = (minp/eta) .* (FtDD_max - FtED_max);
%         Term3 = vi .* (FtCC_min - FtDD_max);
%         Term4 = (maxp-c) * eta .* (FtEC_min - FtCC_min);
%         Term6 = vD .* (1-FtEC_min);
%     elseif (IFtCCmin > IFtEDmax) && (IFtCCmin <= IFtDDmax)
%         Term1 = vC .* FtED_max;
%         Term2 = (minp/eta) .* (FtCC_min - FtED_max);
%         Term3 = vi .* (FtDD_max - FtCC_min);
%         Term4 = (maxp-c) * eta .* (FtEC_min - FtDD_max);
%         Term6 = vD .* (1-FtEC_min);
%     elseif IFtCCmin <= IFtEDmax
%         Term1 = vC .* FtCC_min;
%         Term2 = (minp/eta) .* (FtED_max - FtCC_min);
%         Term3 = vi .* (FtDD_max - FtED_max);
%         Term4 = (maxp-c) * eta .* (FtEC_min - FtDD_max);
%         Term6 = vD .* (1-FtEC_min);
%     end
% end

% output new value function samped at ed
vo = Term1 + Term2 + Term3 + Term4 + Term6;
% vo(vo >= 1e6) = lNum;
% vo(vo <= -1e6) = -lNum;