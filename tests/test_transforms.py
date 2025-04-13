import pytest
from pymath import dft, idft


def test_dft_idft():
    signal = [0.0, 1.0, 0.0, -1.0]
    transformed = dft(signal)
    restored = idft(transformed)
    assert restored == pytest.approx(signal, abs=1e-3)
