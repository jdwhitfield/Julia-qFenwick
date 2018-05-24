# Julia-qFenwick

Quantum simulation of fermionic systems using qubits will play an important role in current and
emerging quantum computing architectures. Limitations of small scale quantum computers will
help inform the roadmap to the future and here we present new tools for creating qubit operators
for fermionic systems. 

The Fenwick tree approach to fermion-to-qubit mappings is implemented in the Julia Language. A similar Python implementation found in https://github.com/quantumlib/OpenFermion as contributed by @VojtaHavlicek.


In the technical note, we give an introduction and a summary of our generalized
approach to Jordan-Wigner and Bravyi-Kitaev transforms and beyond. In particular, the
parity tree forms the crux of our general approach. The primary reference for this work is:

 > Operator Locality in Quantum Simulation of Fermionic Models  
 > Vojtěch Havlíček, Matthias Troyer, James D. Whitfield  
 > Phys. Rev. A 95, 032332 (2017)  
 > https://arxiv.org/pdf/1701.07072.pdf  


James D. Whitfield @ Dartmouth.Edu
