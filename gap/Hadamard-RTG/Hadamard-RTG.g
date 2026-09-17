# Purpose: construct RTG on n^2 vertices, and descendant graphs, from an n x n Hadamard (symmetric)
# matrix H.
# Method: construct n of these, removing a different row of H each time
# Print automorphism group of the graph to file ("hadamard-aut.txt")
# as well as aut group of the subconstituent graph for each v. 
# file will have n aut groups, (one for each row of H), then n^2 subconstituents of each (one for each vertex)
# Repeating this for each vertex is relevant when the graph is not vertex-transitive

# see Fickus et all for the method involved in constructing SRGs this way
#
# usage: enter the file with H at command line, or add a line: Read(<myfile>)
#  e.g. gap hadamard8.txt Hadamard-RTG-construction.g
#
###############################################################################

LoadPackage("AssociationSchemes");;
Read("eye-J.txt");;

n := Length(H);;

pairs := [];;
for i in [1..n/2] do
    Append( pairs, [[2*i-1,2*i]] );;
od;

P := Combinations( [1..n] , 2);;
SubtractSet(P, pairs);
# now construct incidence matrix with first rows from pairs
# rest of the rows from P
nc2 := n*(n-1)/2;;

# matrix B is incidence (n choose 2) x n
B := NullMat(nc2,n);;
for row in [1..n/2] do
    B[ row, pairs[row][1]] := 1;;
    B[ row, pairs[row][2]] := 1;;
od;
for i in [n/2+1..nc2] do
    B[ i, P[i-n/2][1]] := 1;;
    B[ i, P[i-n/2][2]] := 1;;
od;


#  B has row sum 2 and col sum n-1

sub1 := OutputTextFile("subconstituents.txt", true);;
autgp := OutputTextFile("hadamard-aut.txt", true);;
graphs := OutputTextFile("hadamard-graphs.txt", true);;
subgraphs := OutputTextFile("hadamard-subgraphs.txt", true);;

for k in [1..n] do
    	rows := [1..n];;
    	Remove(rows, k);  # H minus the kth row will be Hm
    	Hm := ExtractSubMatrix(H, rows, [1..n]);;
	phi := NullMat(nc2 , n^2);;
	# count maintains position in Hm; loops 2-8 for each column of B
	for j in [1..n] do
    		count := 1;;
    		colstart :=n*(j-1)+1;;
    		colsend := n*j;;
    		for i in [1..nc2] do
        		if B[i,j] = 1 then
            		CopySubMatrix( Hm, phi, [count], [i], [1..n], [colstart..colsend] );;
            		count := count+1;;
        		fi;
    		od;
	od;
	
	# test whether we have phi correct
	G := TransposedMat(phi)*phi/(n-1);;
	if not( G*G = 2*n/(n-1)*G ) then
		Print("mistake in k=", k);
		k := n+1;;
	fi;

	Seidelmat := (n-1)*(eye(n^2) - G);;
	A := ( (n-1)*G - n*eye(n^2) + J(n^2) )/2;;
	#  Acomp := J - A - IdentityMat(n^2);;
	gamma := StronglyRegularGraphScheme( A );;

	AppendTo( graphs, "mat", k, " := ", A, ";;\n\n");;

        autgamma := AutomorphismGroup(gamma);;
	# write aut group to file	
	AppendTo( autgp, "\n", "group when row ", k, " removed is: ", StructureDescription( autgamma ), "\n", "parameters ", StronglyRegularGraphParameters(gamma), "\n", " Orbits: ", Orbits( autgamma ), "\n" );;


	# test nbhds for srg (each vertex)

	for v in [1..n^2] do
    		nbrs_v := [];;  
    		for i in [1..n^2] do
        		if A[v,i] = 1 then
            		Add(nbrs_v,i);;
       		fi;
    		od;
    		A1 := ExtractSubMatrix(A, nbrs_v, nbrs_v );;
		# gamma1 := StronglyRegularGraphScheme(A1);;  # fails, perhaps
    		if not(StronglyRegularGraphScheme(A1) = fail) then
			# write A1 to file but not expecting any
			gamma1 := StronglyRegularGraphScheme(A1);; 
			AppendTo( sub1, "\n", "Subconstituent DRG k = ", k, " vertex ", v, "\n" , StronglyRegularGraphParameters(gamma1), "\n\n");;
		fi;


		# construct descendant by isolating v
		# for switching, need a diagonal +-1 matrix with -1 on neighbours of v
		switchingvec := allones(n^2);;
		for i in [1..n^2] do
    			if i in nbrs_v then
        			switchingvec[i] := -1;;
    			fi;
		od;
		switchingmat := DiagonalMat(switchingvec);;

		# switch the Seidel matrix 
		Seidel_isolated := switchingmat*Seidelmat*switchingmat;;
		vertices := [1..n^2];;
		Remove( vertices, v );;
		Seidel0 := ExtractSubMatrix(Seidel_isolated, vertices, vertices);;
		N := n^2-1;;

		A0 := (J(N) - Seidel0 - eye(N))/2;;
		# A0comp := 1/2*(J(N) + Seidel0 - IdentityMat(N));;
		gamma0 := StronglyRegularGraphScheme(A0);;
		
		AppendTo( autgp, "descendant group for vertex ", v, " is: ", StructureDescription( AutomorphismGroup(gamma0) ), "\n", "parameters ", StronglyRegularGraphParameters(gamma0), "\n");;
		AppendTo( subgraphs, "mat", k, "v", v, " := ", A0, ";;\n\n");;

	od;  # end loop on vertices
od;  # end loop on rows of H

CloseStream(sub1);;
CloseStream(autgp);;
CloseStream(graphs);;
CloseStream(subgraphs);;

#  All of the above repeated with each matrix H, and with each of the n-1 rows removed.
#  Possibly: permute which row of H' replaces which 1 in each col of the incidence matrix.
