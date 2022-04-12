import math as m
import numpy as np

class Transformacje:
    
    def __init__(self, model:str = "wgs84"):
        
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518
            
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
            
        else:
            raise NotImplementedError(f"{model} model not implemented")
        
        self.flattening = (self.a-self.b)/ self.a
        self.ecc = 2*self.flattening - self.flattening**2
