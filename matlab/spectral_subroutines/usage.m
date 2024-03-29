fprintf('Usage : demo_sbm(N,q,cin,cout,seed,method) where:\n');
fprintf('   N is the number of nodes (will be rounded to a multiple of the number of communities)\n');
fprintf('   q is the number of communities\n');
fprintf('   cin is the average intra-community connectivity\n');
fprintf('   cout is the average extra-community connectivity\n');
fprintf('   seed is the seed of the random number generator\n');
fprintf('   method is a string that can be either :\n');
fprintf('     - ''BH'' : clustering with the Bethe Hessian\n');
fprintf('     - ''Bprime'' : clustering with the reduced form of the non-backtracking matrix\n');
fprintf('     - ''A'' : clustering with the adjacency matrix\n');
fprintf('     - ''Lap'' : clustering with the symmetrically normalized laplacian\n');
