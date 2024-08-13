# Numerical methods in practice

## Linear algebra

### Cost scaling

For each of the following linear algebra operations, confirm the computational cost scaling with matrix size ($N \times N$) and verify it numerically by timing the same operation with matrices of increasing size.

1. Matrix condensation $O(N^2)$: `np.sum`, `np.max`, `np.linalg.norm`
2. Matrix multiplication $O(N^3)$: `A = np.dot(B,C)`
3. Linear solve $O(N^3)$: `A = np.linalg.solve(A,b)`
4. Eigenvalue decomposition of symmetric matrix $O(N^3)$: `e, U = np.linalg.eigh(A)`
