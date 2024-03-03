function emb_vec = embedding(rr_ints,tau,m)

% ----------------------------------------------------------------------- %
% phase-space reconstruction using time-delay embedding, according to
% Taken's theorem

% input RR = input vector of RR intervals
% input tau = lag used for time-delay embedding
% input m = embedding dimension
% output emb_vec = matrix of vectors (in rows) of the reconstructed
% attractor
% ----------------------------------------------------------------------- %

% ensure rr_ints are column vector
if size(rr_ints,2) > size(rr_ints,1)
    rr_ints = rr_ints';
end
num_ints = length(rr_ints);

% pre-allocate matrix of embedded vectors
emb_vec = [];
for k = 1:num_ints-(m-1)*tau
    emb_vec = [emb_vec; rr_ints(k:tau:k+(m-1)*tau)'];
end

end