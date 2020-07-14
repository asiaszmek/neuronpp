import numpy as np
from typing import Optional


class Dist:
    seed = None

    def __init__(self, dtype="float"):
        self.dtype = dtype

    @classmethod
    def set_seed(cls, seed):
        np.random.seed(seed)
        cls.seed = seed


class UniformDist(Dist):
    def __init__(self, low=0, high=1, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        self.low = low
        self.high = high


class NormalDist(Dist):
    def __init__(self, mean, std, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        self.mean = mean
        self.std = std


class NormalTruncatedDist(NormalDist):
    def __init__(self, mean, std, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        if mean <= 0 or std <= 0:
            raise ValueError("mean and std cannot be <= 0 for Truncated Normal Distribution.")
        NormalDist.__init__(self, mean=mean, std=std)


class UniformSegDist(UniformDist):
    def __init__(self):
        """
        Uniform Segment Distribution, mean is 0 high is 1.
        """
        Dist.__init__(self, dtype="float")
        self.low = 0
        self.high = 1


class NormalTruncatedSegDist(NormalTruncatedDist):
    def __init__(self, std, mean: Optional[float] = None):
        """
        Normal distribution for segment choosing.
        Normal means that choosing segments are clustered around mean with standard deviation std.
        :param mean:
            Provided in normalized arbitrary unit between 0-1.

            Default is None, which means - mean will be selected with random uniform distribution
            from all provided segments.

            It is normalized mean (between 0-1), where all provided segments are organized as list
            and first element has location=0 and the last location=1
            During computation this number will be change for appropriate mean in um.
        :param std:
            Provided in um.
            standard deviation of the cluster of distribution.
        """
        Dist.__init__(self, dtype="float")

        if mean > 1:
            raise ValueError("Param mean cannot be > 1")
        if mean < 0 or std <= 0:
            raise ValueError("mean and std cannot be <= 0 for Truncated Normal Distribution.")
        NormalTruncatedDist.__init__(self, mean=mean, std=std)


class Proba(Dist):
    def __init__(self, expected: float = 0.5, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        self.expected = expected


class UniformProba(Proba, UniformDist):
    def __init__(self, expected: float = 0.5, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        Proba.__init__(self, expected=expected)
        UniformDist.__init__(self, low=0, high=1)


class NormalProba(Proba, NormalTruncatedDist):
    def __init__(self, mean: float = 0.5, std: float = 0.1, dtype="float"):
        Dist.__init__(self, dtype=dtype)
        Proba.__init__(self, expected=mean)
        NormalTruncatedDist.__init__(self, mean=mean, std=std)
