# 微分可能なNeedleman-Wunschアルゴリズム
import jax
import jax.numpy as jnp
from jax.nn import logsumexp

def differentiable_nw(
    seq1, seq2, match_score=1.0, mismatch_score=-1.0, gap_penalty=-1.0, temperature=0.01):
    
    len1, len2 = len(seq1), len(seq2)
    
    # initialize score matrix
    S = jnp.zeros((len1 + 1, len2 + 1)) - jnp.inf
    S = S.at[0, 0].set(0.0)
    
    # create substitution matrix
    def substitution(a, b):
        if (a == b):
            return match_score
        else:
            return mismatch_score
    
    subs_matrix = jnp.array([[substitution(a, b) for b in seq2] for a in seq1])
    
    # compute scores
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = S[i - 1, j - 1] + subs_matrix[i - 1, j - 1]
            delete = S[i - 1, j] + gap_penalty
            insert = S[i, j - 1] + gap_penalty

            S = S.at[i, j].set(logsumexp(jnp.array([match, delete, insert]) / temperature) * temperature)
    
    return S[len1, len2]

# example
seq1 = "GATTACA"
seq2 = "GCTTGCA"
score = differentiable_nw(seq1, seq2, temperature=0.1)
print("Alignment score:", score)
