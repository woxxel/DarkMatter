
% Copr: (C) 2017 by Mantas Gabrielaitis, MPI ds Goettingen, mantas@nld.ds.mpg.de

function [X, Y, F, ext] = f_1d_noise_deconvolution (y, F0, dN, TikhFact, Nboots, optionsQP)
% >
% F_1D_NOISE_DECONVOLUTION performs noise deconvolution of a given probability
% distribution, with a specified noise kernel, non-zero-range, and Tikhonov factor.
%
% -----------------------------------------------------------------------------------
% i n p u t
% [y]-> A column vector corresponding to the probability distribution of the original
% data. Can also be concatenation of (bootstrapping) samples
% [F0]-> 1D cell array with information about the noise kernel: first size(Y,1)
% elements of the array are column vectors of the noise kernel (of different
% lengths, in general), the last two elements are a vector of length size(Y,1) with
% the center indexes of each of the kernel columns and a vector of length size(Y,1)
% with the lengths of each of the kernel columns.
% [dN]-> Array of two elements: dN(1) - the number of the leftmost entries of the
% original distribution that will be set to 0 in the case of the noise-deconvolved
% distribution, dN(2) - the number of the rightmost entries of the original
% distribution that will be set to 0 in the case of the noise-deconvolved
% distribution.
% [TikhFact]-> Tikhonov factor.
% [Nboots]-> size of the bootstrapping sample (Nboots < 0 -> no bootstrapping).
% [optionsQP]-> Options structure for the quadratic programming minimizer.
%
% o u t p u t
% [X]-> Column vector of the noise-deconvolved probability distribution on the
% extended grid. (if Nboots>1 this is average)
% [Y]-> Column vector of the original probability distribution on the extended grid.
% [F]-> Square matrix of the noise kernel on the extended grid.
% [ext]-> Array of two elements: ext(1) - the number of points to be added to the
% left of the original grid when creating the extended grid, ext(1) - the number of
% points to be added to the right to the original grid when creating the extended
% grid.
%
% d e p e n d e n c i e s
% One external function is used: quadprog (Optimization Toolbox).
% -----------------------------------------------------------------------------------

szF = size(F0);
szY = size(y);
Nboots = Nboots * (Nboots > 0);

if length(szF) > 2 || (szF(1) ~= 1 && szF(2) ~= 1)
    error('Cell array F0 has has to be a vector!');
elseif length(szY) > 2 || szY(2) ~= Nboots+1
    error('The probability distribution vector y has has to be a column vector!');
elseif (dN(1)+dN(2) >= szY(1))
    error('The sum dN(1)+dN(2) has to be smaller than the length of the probability distribution vector y!')
elseif szY(1) ~= length(F0)-2;
    error('The number of rows in the cell array F0 has to be larger by two than the number of rows  of the probability distribution vector y!');
elseif any(y(:)<-5*eps) || any(abs(sum(y,1)-1)>5*eps)
    error('The sum of elements of the probability distribution vector y has to be equal to 1! All elements of y have to be nonnegative!');
else
    for i=1:szY(1)
        if any(F0{i}<-5*eps) || abs(sum(F0{i})-1)>10*eps
            error('The sum of elements in each column vector of F0 has to be equal to 1 and all its elements have to be nonnegative!');
        end
    end
end

%%% NOISE KERNEL F

% % center indexes of each column of F0
% index position of the leftmost element with respect to the central element
iL = F0{end-1};
% % length of each column of F0
% index position of the rightmost element with respect to the central element
iR = F0{end};

% the number of entries we have to extend the grid by to the left (1) and to the
% right (2)
ext = [max([0,-iL-(0:szY(1)-1)]), max([0, iR-(szY(1)-1:-1:0)])];

% the number of entries from the left and from the right which are set to zero for
% the noise-deconvolved distribution
extL2 = ext(1) + dN(1);
extR2 = ext(2) + dN(2);

% the length of the new grid
LE = szY(1) + ext(1) + ext(2);

% extended vector of the original distribution of the data
Y = [zeros(ext(1),Nboots+1); y; zeros(ext(2),Nboots+1)];

% noise matrix
F = zeros(LE);
for i=1:LE-extL2-extR2
    F(extL2+i+iL(i+dN(1)):extL2+i+iR(i+dN(1)),extL2+i) = F0{i+dN(1)};
end


%%% MATRIX H and VECTOR F for QUADPROG
aux = TikhFact * diag(ones(1,LE));
aux(1:extL2,1:extL2) = 0;
aux(end-extR2+1:end,end-extR2+1:end) = 0;

H = F'*F + aux;

f = -F'*Y;

%%% OPTIMIZATION CONSTRAINTS
lb = zeros(LE,1);                           % lower bound for prob
ub = ones(LE,1);                            % upper bound for prob

Aeq = ones(1,LE);                           % sum over columns has to be = 1 (normalization constraint)
beq = 1;                                    % normalization

% forces non-wanted elements of reconstructed distribution to have prob = 0
Aeq = [Aeq; zeros(extL2+extR2,LE)];
Aeq(2:1+extL2,1:extL2) = eye(extL2);
Aeq(2+extL2:1+extL2+extR2,end-extR2+1:end) = eye(extR2);
beq = [beq; zeros(extL2+extR2,1)];

% initial conditions
x0 = ones(LE,1)/LE;
x0(1:extL2) = 0;
x0(end-extR2+1:end) = 0;

%%% NOISE DECONVOLUTION
X = zeros(size(Y,1),1);
for i=1:Nboots+1
    X = X + quadprog(H, f(:,i), [], [], Aeq, beq, lb, ub, x0, optionsQP);
end

X = X / (Nboots+1);
Y = mean(Y,2);