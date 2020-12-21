# ch2bt
GAP calculations of torsion in codimension 2 cycles on classifying spaces of tori.


## Requirements:
Calculations require additional GAP packages:
* LocalSNF, https://github.com/neshitov/localsnf
* Gauss, https://www.gap-system.org/Packages/gauss.html

### Acknowledgements:
Computations here use the following results which are included in this repository for imports:

* CARAT tables (file "carat_tables.tar.gz"), part of the CARAT package:<br> [OPS19] Opgenorth, J., Plesken, W. and Schulz, T., Carat, Crystallographic AlgoRithms And Tables, Version 2.1b1 (2008)
https://github.com/lbfm-rwth/carat/.
* Algorithm for finding a flabby resolution (folder "./RatProbAlgTori") by A. Hoshi and A. Yamasaki:<br>
[HY12] Hoshi, A., Yamasaki, A., Rationality Problem for Algebraic Tori. Memoirs of the American Mathematical Society. 248. 10.1090/memo/1176 (2012). https://www.math.kyoto-u.ac.jp/~yamasaki/Algorithm/RatProbAlgTori/

### Installation and usage:
Run the ./configure script to unpack CARAT tables.<br>
File Example.g provides an example of calculation for the case
G = CARAT(6, 2377, 10), M=Z^6.<br> It computes the group Phi(G, M) of the paper
"On a pairing for algebraic tori
Mathematische Nachrichten, 292 (2019), 2283-2293." by A. Merkurjev. <br>
In this case Phi(G, M) = Z/2, thus CH^2(BT)_tors = Z/2 when T is versal with splitting Galois group G.
