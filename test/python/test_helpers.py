import pytest
from assembly_theory import _validate_bounds

def test_validate_bounds():
    # Case 1: bounds is None and no_bounds is True → should return an empty set
    assert _validate_bounds(None, True) == set()

    # Case 2: bounds is None and no_bounds is False → should return {"intchain", "vecchain"}
    assert _validate_bounds(None, False) == {"intchain", "vecchain"}

    # Case 3: bounds is a set and no_bounds is True → should raise ValueError
    with pytest.raises(ValueError, match="bounds specified but `no_bounds` is True."):
        _validate_bounds({"log"}, True)

    # Case 4: bounds is a set and no_bounds is False → should return the given set
    assert _validate_bounds({"log"}, False) == {"log"}
