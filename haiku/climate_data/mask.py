from __future__ import annotations
from dataclasses import dataclass
import pickle
import numpy.typing as npt


@dataclass
class Mask:

    indices: npt.NDArray

    def serialize(self, filepath: str) -> None:
        if not filepath.endswith(".pkl"):
            filepath += ".pkl"
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)

    def deserialize(filepath: str) -> Mask:
        with open(filepath, 'rb') as f:
            return pickle.load(f)
