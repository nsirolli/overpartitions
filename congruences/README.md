- In the folder _interesting_ we store the data which supports the congruences
  stated in Theorems 9.3 and 9.5. Each file references in its name to a pair
  $(l,Q)$, and contains in its lines the pairs $(n,\overline p(n) \pmod{l})$
  required to prove the congruences in Algorithm 9.1.

- On the other hand, for pairs $(l,Q)$ for which Algorithm 9.1 returns False, in
  the folder _uninteresting_ we store data explaining the failure: each file has
  the format:

  $n$

  $(n, \overline p(n) \pmod{l})$

  $(nQ^2, \overline p(nQ^2) \pmod{l})$

  $(n/Q^2, \overline p(n/Q^2) \mod{l})$

  $res$

  where $res$ is the nonzero value corresponding to $n$ according to line 12 of
  the algorithm.
