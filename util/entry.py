from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np

@dataclass
class Entry:
    arr: np.ndarray
    header: Optional[Tuple[str, ...]] = None
    axis: Optional[Tuple[np.ndarray, ...]] = None