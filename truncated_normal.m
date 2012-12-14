function [out] = truncated_normal(a, b, mu, sigma, n)

% [OUT] = TRUNCATED_NORMAL(A, B, MU, SIGMA, N)
% 
% Purpose
% 
% To generate N samples from a truncated normal distribution with mean=mu, sigma=sigma and with bounds A and B
%  
% Input
%
% --A: lower bound
% --B: upper bound
% --mu: mean
% --sigma: sigma, standard deviation
% --N: number of samples
%
% Output
%
% --out: mean + truncated Gaussian noise with mean=mu, sigma=sigma, between bounds A and B 
%
% Example usage:
%
% truncated_normal(-1, 1, 0.1, 0.001, 1000)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the P-CIT toolbox released under the BSD license.
% Copyright (c) 2012, Princeton University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in %the documentation and/or other materials provided with the distribution.
% Neither the name of the Princeton University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT %NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, error('Missing input aruguments!'); end
if (b - a) < 0, error('Lower bound is greater then upper bound!'); end
if sigma <= 0, error('Sigma is <= 0!'); end

PHIl = 0.5 .* erfc(-(((a - mu) ./ sigma) ./ sqrt(2)));
PHIr = 0.5 .* erfc(-(((b - mu) ./ sigma) ./ sqrt(2)));

% Refer to http://www.maths.uq.edu.au/~chancc/10Fstat3001/ass4sol.pdf for truncated normal dist sampling below -- If this source does not exist then refer to code in the link below,
% http://www.wiley.com/legacy/wileychi/koopbayesian/supp/normt_rnd.m
out = mu + sigma .* ( sqrt(2) .* erfinv(2 .* (PHIl + (PHIr - PHIl) .* rand(n, 1)) - 1) );

