function out = midpoint_op(N)
bands = ones(N,1);
out = 0.5.*(spdiags(bands,1,N-1,N) + spdiags(bands,0,N-1,N) );
