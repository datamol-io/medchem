def _in_range(x, min_val: float = -float("inf"), max_val: float = float("inf")):
    """Check if a value is in a range
    Args:
        x: value to check
        min_val: minimum value
        max_val: maximum value
    """
    return min_val <= x <= max_val
