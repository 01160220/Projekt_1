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

    def xyz_2_blh(X, Y, Z, self):
        """
        Funkcja przelicza współrzędne geocentryczne (ECEF) 
        na współrzędne geodezyjne (Algorytm Hirvonena).

        Parameters
        ----------
        X   : [float] : współrzędna geocentryczna (ECEF) [m]
        Y   : [float] : współrzędna geocentryczna (ECEF) [m]
        Z   : [float] : współrzędna geocentryczna (ECEF) [m]
        a   : [float] : dłuższa półoś elipsoidy [m]
        e2  : [float] : mimośrod elipsoidy [niemianowana]
        Returns
        -------
        fi  : [float] : szerokość geodezyjna [rad]
        lam : [float] : długość geodezyjna [rad]
        h   : [float] : wysokość elipsoidalna [m]

        """
        r = np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(r*(1-self.ecc)))
        eps = 0.000001/3600*m.pi/180
        fi0 = fi*2

        while abs((fi - fi0).all()) > eps:
            fi0 = fi
            N = self.a/np.sqrt(1-self.ecc*np.sin(fi)**2)
            h = r/np.cos(fi)-N
            fi = np.arctan(Z/(r*(1-self.ecc*(N/(N+h)))))

        lam = np.arctan(Y/X)
        N = self.a/np.sqrt(1-self.ecc*np.sin(fi)**2)
        h = r/np.cos(fi)-N
        return fi, lam, h 