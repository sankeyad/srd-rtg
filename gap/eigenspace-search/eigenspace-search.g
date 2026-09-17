# This code looks for specific eigenvectors that delineate a strongly regular decomposition
# Input required is a graph (gamma) adjacency matrix along with parameters for that graph and the two expected
#  strongly regular subgraphs (format specified below)
# WARNING: This runs through 2^f cases where f is the dimension of the eigenspace of gamma. 
# The reasonable limit on a laptop is f=35.
# enter the graph file on command line, as in
    # gap M22.txt eigenspace-search.g
    
# usage: gap [<mygraph>] eigenspace-search.g
# for example:
    # gap M22.txt eigenspace-search.g
# <mygraph> must supply adjacency matrix, parameters v(_i), k(_i), lambda(_i), mu(_i), r, s, f, g
# format is:
#
#		gamma :=  [[ ...adjacency matrix... ]];;
#		v := 100;;
#		k := 22;;
#		lambda := 0;;
#		mu := 6;;
#		r := 2;;
#		s := -8;;
#		f := 77;;
#		g := 22;;
#		v1 := 50;;
#		k1 := 7;;
#		lambda1 := 0;;
#		mu1 := 1;;
#		v2 := 50;;
#		k2 := 7;;
#		lambda2 := 0;;
#		mu2 := 1;;
# see sample M22.txt if this is unclear


LoadPackage("AssociationSchemes");;   # can also put this on command line

Read("eye-J.txt");;

S := StronglyRegularGraphScheme(gamma);

S := StronglyRegularGraphScheme(gamma);;
M := IntersectionMatrices(S);;
P := MatrixOfEigenvalues(S);;
Q := P^(-1)*v;;
A := AdjacencyMatrices(S);;

# define the idempotents
# requires functions J(n) and eye(n)

eid := [];;
eid[1] := 1/v*( Q[1,2]*eye(v) + Q[2,2]*A[2] + Q[3,2]*A[3] );;
eid[2] := 1/v*( Q[1,3]*eye(v) + Q[2,3]*A[2] + Q[3,3]*A[3] );;

# determine which eigenspace to use
if k1+k2-k = r then
dim := f;;
num := 1;; # which idempotent
else
dim := g;;
num := 2;;
fi;

# can add verification for debugging:
# if not( A[2]*eid[1] = r*eid[1] ) then  # write some bail-out code
# A[2]*eid[2] = s*eid[2];

# get nice basis for eigenspace
mat := TransposedMat( eid[num] );;
bm := BaseMat( mat );;
temp := TriangulizedMat(bm);;
matt := TransposedMat( temp );;
# mat is vxv, bm is fxv, matt is vxf

# computations over the rationals given integer eigenvalues

# testing, for debugging
# A[2]*matt = s*matt;

# set up good/bad columns
goodcols := [];;
for i in [1..dim] do
    # find a 1 in each row
    j := 1;;
    while (temp[i][j] = 0 ) do
        j := j+1;;
    od;;
    # j is now col where a leading 1 is
    Append(goodcols, [j] );;
od;;

badcols := [1..v];;
SubtractSet(badcols, goodcols);;

# testing, for debugging
# Length(goodcols)=g;
# Length(badcols)=f+1;

id := ExtractSubMatrix( matt, goodcols, [1..dim] );;
id = eye(dim);

# set up the search
L := AsList( badcols );;
# looking for S2^v1, (-S1)^v2 in eigenspace k1+k2-k,
# where S2 = k-k1, S1 = k-k2
S2 := k-k1;;
S1 := k-k2;;
yes := 0;;
yesvecs := [];;
#
count := 0;;
for vec in IteratorOfTuples( [S2, -S1] , dim ) do   
    i := 1;;
    while i<(v-dim) and (vec*matt[ L[i] ] in AsSet([S2, -S1]) ) do
        i := i+1;;
    od;;
    if i=(v-dim) then   # not f or g here, v-dim is number of entries to check
        Append(yesvecs, [vec] );;
        yes := yes+1;;
    fi;
od;

splits := OutputTextFile("possiblesplits.txt", true);;
AppendTo(splits, "splits := ", yesvecs, ";;", "\n" , "# Number found: ", yes, "\n");;
CloseStream(splits);;

# the 'possible splits' define eigenvectors that indicate a partition
# into two *regular* graphs. They are not *strongly regular* necessarily, so
# what follows will locate the actual strongly regular decompositions

winners := [];;
evecs := [];;
for i in [1..Length(yesvecs)] do
    split := matt*yesvecs[i];;  
    # what is stored in yesvecs is the coefficients of a particular linear combination 
    # of the eigenspace basis vectors.
    # to expand it to a full eigenvector we multiply by the vxf matrix (matt) whose columns
    # form the basis. We then get 'split' as a vector of length v with entries S2, -S1.
    V1 := [];;
    V2 := [];;
    for j in [1..v] do
        if split[j] = S2 then
            Add(V1, j);;
        else
            Add(V2, j);;
        fi;
    od;
    # test this partition
    gamma1 := ExtractSubMatrix( A[2], V1, V1);;
    gamma2 := ExtractSubMatrix( A[2], V2, V2);;
    if gamma1^2 = k1*eye(v1) + lambda1*gamma1 + mu1*(J(v1)-eye(v1)-gamma1 ) and
gamma2^2 = k2*eye(v2) + lambda2*gamma2 + mu2*(J(v2)-eye(v2)-gamma2 ) then
       Append( evecs, [split] );; 	   
       Append(winners, [i]);;
    fi;
od;

foundSRDs := OutputTextFile("foundSRDs.txt", true);;
AppendTo( foundSRDs, "# Eigenvectors permitting a strongly regular decomposition", "\n", "splits := ", evecs, ";;", "\n");;
AppendTo( foundSRDs, "# Number found: ", Length(winners), "\n\n" );;
AppendTo( foundSRDs, "# Successful vectors from the 2^f combinations are: ", winners, "\n" );;
CloseStream(foundSRDs);;




