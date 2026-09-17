# Strongly regular decompositions in descendants of regular twographs

## Description  
This site houses GAP code and text files for the constructions and analysis of strongly regular graphs in the family  
$\Gamma^{\pm}=(4m^2, 2m^2\pm m, m^2\pm m, m^2\pm m)$.  
These are representatives of *regular twographs* and therefore have strongly regular descendant graphs in the family  
$\Gamma=(4m^2-1, 2m^2, m^2, m^2)$.  

## Available files
1. GAP code and sample input
   1. [construction of $\Gamma^{\pm}$](gap/Hadamard-RTG) from a symmetric Hadamard matrix H of order $2m$, and automorphism group. This is
done once for each row of H in .
   2. [eigenspace search](gap/eigenspace-search) for strongly regular decompositions within a candidate graph. This requires an input file with
the adjacency matrix and candidate parameters, as in the sample input file.
2. adjacency matrices
   1. [m=2](matrices/m=2)
   2. [m=4](matrices/m=4)
   3. [m=8](matrices/m=8)
   4. [m=10](matrices/m=10)
