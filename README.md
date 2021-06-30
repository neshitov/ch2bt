# ch2bt
Calculations of torsion in codimension 2 cycles on classifying spaces of tori in GAP computer algebra system [GAP21]. Calculates are done using the group Phi(G,M) introduced in [Merkurjev19]

## Requirements:
Calculations require additional GAP packages:
* Gauss, https://www.gap-system.org/Packages/gauss.html
* LocalSNF, https://github.com/neshitov/localsnf

### Acknowledgements:
Computations here use the following results which are included in this repository for imports:

* CARAT tables (file "carat_tables.tar.gz"), enumeration of finite groups of unimodular matrices, made by A. Hoshi and A. Yamasaki [HY12] using the CARAT package [OPS19]

* Algorithm for finding a flabby resolution (code in folder "./RatProbAlgTori") by A. Hoshi and A. Yamasaki [HY12]

### Installation:
Run the ./configure script to unpack CARAT tables.<br>

### Files:
* Example.g provides example of the group Phi(G,M) calculations
* ComputePhi.g contains the necessary code for calculation of Phi(G,M)
* CoflasqueCover.g: code to compute coflasque resolutions, based in algorithm of [HY12]
* filter_groups.g: code to find CARAT ids of all groups possibly having nontrivial Phi(G,M)
* check.g: code to check all the groups in dimensions 3,4,5,6
* ./dim_../carat_ids.txt: CARAT ids of all 2-groups in corresponding dimension that can possibly have nontrivial Phi(G,M), is produced by filter_groups.g
* ./dim_../result.txt: Result of computations for groups listed in ./dim_../carat_ids.txt
* ./dim_6/result_table.tex: LaTeX table containing the list of all 2-subgroups G of GL(6,Z) with nontrivial Phi(G,M)
### References:
* [GAP21] The GAP Group, GAP -- Groups, Algorithms, and Programming, Version 4.11.1; 2021. (https://www.gap-system.org)
* [HY12] Hoshi, A., Yamasaki, A., Rationality Problem for Algebraic Tori. Memoirs of the American Mathematical Society. 248. 10.1090/memo/1176 (2012). https://www.math.kyoto-u.ac.jp/~yamasaki/Algorithm/RatProbAlgTori/
* [Merkurjev19] A. Merkurjev "On a pairing for algebraic tori
Mathematische Nachrichten, 292 (2019), 2283-2293."
https://www.math.ucla.edu/~merkurev/papers/pairing3.pdf
* [OPS19] Opgenorth, J., Plesken, W. and Schulz, T., Carat, Crystallographic AlgoRithms And Tables, Version 2.1b1 (2008)
https://github.com/lbfm-rwth/carat/.
