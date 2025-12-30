import numpy as np
from dataclasses import dataclass

@dataclass
class DistributionModelParams:
    gamma: float
    delta: float
    nu_max: float

    ## careful: kappa should already enter in calculation of all parameters above! and not necessarily linearly
    kappa: float = 1.

    def are_values_ok(self):
        for field in self.__dataclass_fields__:
            # print(field)
            if np.isnan(getattr(self,field)) or getattr(self,field)<0:
                return False
        return True