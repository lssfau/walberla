def _factors(n):
    while n > 1:
        for i in range(2, n + 1):
            if n % i == 0:
                n //= i
                yield i
                break


def block_decomposition(processes):
    """Compute a 3D factorization of 'processes' in a 3 factors returned as tuple"""
    result = [1, 1, 1]
    for factor in _factors(processes):
        min_idx = result.index(min(result))
        result[min_idx] *= factor
    assert result[0] * result[1] * result[2] == processes
    return tuple(result)
