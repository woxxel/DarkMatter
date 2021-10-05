
% Copr: (C) 2017 by Mantas Gabrielaitis, MPI ds Goettingen, mantas@nld.ds.mpg.de

function [X, Y, Yrec, F, ext, dN, RAO, qRangeY, qRangeYrec, qRangeX] = ...
    f_1d_noise_deconvolution_qrange_adjusted (qval, y, F0, TikhFact, Nboots, optionsQP)
% >
% F_1D_NOISE_DECONVOLUTION_QRANGE_ADJUSTED performs noise deconvolution of a given
% probability distribution, with a specified noise kernel and Tikhonov factor, such
% that the reconstructed distribution has the same quantile range as the original
% one.
%
% -----------------------------------------------------------------------------------
% i n p u t
% [qval]-> Array of two elements: qval(1) - left quantile, qval(2) - right quantile.
% [y]-> A column vector corresponding to the probability distribution of the original
% data.
% [F0]->  1D cell array with information about the noise kernel: first size(Y,1)
% elements of the array are column vectors of the noise kernel (of different
% lengths in general), the last two elements are a vector of length size(Y,1) with
% the center indexes of each of the kernel columns and a vector of length size(Y,1)
% with the lengths of each of the kernel columns.
% [TikhFact]-> Tikhonov factor.
% [Nboots]-> size of the bootstrapping sample (Nboots < 0 -> no bootstrapping).
% [optionsQP]-> options structure for the quadratic programming minimizer.
%
% o u t p u t
% [X]-> Column vector of the noise-deconvolved probability distribution on the
% extended grid.
% [Y]-> Column vector of the original probability distribution on the extended grid.
% [Yrec]-> Column vector of the reconstructed probability distribution y on the
% extended grid.
% [F]-> Square matrix of the noise kernel on the extended grid.
% [ext]-> Array of two elements: ext(1) - the number of points to be added to the
% left of the original grid when creating the extended grid, ext(1) - the number of
% points to be added to the right to the original grid when creating the extended
% grid.
% [dN]-> Array of two elements: dN(1) - the number of the leftmost entries of the
% original distribution that will be set to 0 in the case of the noise-deconvolved
% distribution, dN(2) - the number of the rightmost entries of the original
% distribution that will be set to 0 in the case of the noise-deconvolved
% distribution.
% [RAO]-> Relative area overlap of the original and reconstructed distributions y on
% the extended grid.
% [qRangeY]-> The quantile range of the probability distribution of the original data
% corresponding to the quantile values defined in the input variable qval.
% [qRangeYrec]-> The quantile range of the reconstructed probability distribution of
% the original data corresponding to the quantile values defined in the input
% variable qval.
% [qRangeX]-> The quantile range of the noise deconvolved probability distribution
% of the original data corresponding to the quantile values defined in the input
% variable qval.
%
% d e p e n d e n c i e s
% One external function is used: f_1d_noise_deconvolution.m.
% -----------------------------------------------------------------------------------


% Quantile range of the original probability distribution Y
[~, Y, ~, ~] = f_1d_noise_deconvolution (y, F0, [0 0], TikhFact, Nboots, optionsQP);
aux = cumsum(Y);
qRangeY = [find(aux>=qval(1),1,'first'); find(aux>=qval(2),1,'first')];


le = size(y,1);

dN = Inf*ones(2,le,le);
qRangeYrec = Inf*ones(2,le,le);
qRangeX = Inf*ones(2,le,le);
qRangeYdiff = Inf*ones(le,le);
RAO = Inf*ones(le,le);


for i=0:le-1
    for j=0:le-1-i
        
        dN(:,i+1,j+1) = [i j];
        
        [X, Y, F, ~] = f_1d_noise_deconvolution (y, F0, dN(:,i+1,j+1), TikhFact, Nboots, optionsQP);
        Yrec = F*X;

        % Quantile range of the reconstructed probability distribution Y
        aux = cumsum(Yrec);
        qRangeYrec(:,i+1,j+1) = [find(aux>=qval(1),1,'first'); find(aux>=qval(2),1,'first')];
        
        % Quantile range of the noise-deconvolved probability distribution X
        aux = cumsum(X);
        qRangeX(:,i+1,j+1) = [find(aux>=qval(1),1,'first'); find(aux>=qval(2),1,'first')];
        
        
%         qRangeYdiff(i+1,j+1) = sum(abs(qRangeY - qRangeYrec(:,i+1,j+1)));
        qRangeYdiff(i+1,j+1) = max(abs(qRangeY - qRangeYrec(:,i+1,j+1)));
        
        RAO(i+1,j+1) = 1 - sum(abs(Y-Yrec)) / 2;
        
    end
end

% finds the element corresponding to the maximum RAO of those which minimize the
% difference in quantile ranges between the original and reconstructed probability
% distribution Y.
minval = min(qRangeYdiff(:));
ind = find(qRangeYdiff == minval);
[~,ind2] = max(RAO(ind));
indg = ind(ind2);

% selects values of dN, qRangeYrec, and qRangeX corresponding to the best element
dN = dN(:,indg);
qRangeYrec = qRangeYrec(:,indg);
qRangeX = qRangeX(:,indg);
RAO = RAO(indg);

% calculates, X, Yrec, and ext corresponding to the best element
[X, ~, F, ext] = f_1d_noise_deconvolution (y, F0, dN, TikhFact, Nboots, optionsQP);
Yrec = F*X;