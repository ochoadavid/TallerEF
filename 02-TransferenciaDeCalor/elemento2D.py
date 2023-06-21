import numpy as np
import matplotlib.pyplot as plt


class ElTriTC:
    def __init__(self, x, y, n, k, h12=0, h23=0, h31=0, q12=0, q23=0, q31=0, Q=0):
        '''Elemento triangular para transferencia de calor.'''
        # Esta función simplemente guarda toda la información, calcula el área y las longitudes de los lados
        # Los parámetros igualados a cero son opicionales y su valor por defecto es cero
        self.x = np.array(x)
        self.y = np.array(y)
        self.n = n
        self.k = float(k)
        self.A = np.abs(self.x[0]*(self.y[1]-self.y[2])+
                        self.x[1]*(self.y[2]-self.y[0])+
                        self.x[2]*(self.y[0]-self.y[1]))/2
        self.l12 = np.sqrt((self.y[1]-self.y[0])**2 + (self.x[1]-self.x[0])**2)
        self.l23 = np.sqrt((self.y[2]-self.y[1])**2 + (self.x[2]-self.x[1])**2)    
        self.l31 = np.sqrt((self.y[0]-self.y[2])**2 + (self.x[0]-self.x[2])**2)
        self.h12 = h12
        self.h23 = h23
        self.h31 = h31
        self.q12 = q12
        self.q23 = q23
        self.q31 = q31
        self.Q = Q
        
    def MatrizK(self, expNodos=False):
        '''Calcula la matríz K del elemento, incluyendo la convección.'''
        # Cálculo de b y c
        b1=self.y[1]-self.y[2]
        c1=self.x[2]-self.x[1]
        b2=self.y[2]-self.y[0]
        c2=self.x[0]-self.x[2]
        b3=self.y[0]-self.y[1]
        c3=self.x[1]-self.x[0]
        Kb = np.array([[b1*b1, b1*b2, b1*b3],
                       [b2*b1, b2*b2, b2*b3],
                       [b3*b1, b3*b2, b3*b3]])
        Kc = np.array([[c1*c1, c1*c2, c1*c3],
                       [c2*c1, c2*c2, c2*c3],
                       [c3*c1, c3*c2, c3*c3]])
        # Matriz de conductividad
        K = self.k / (4*self.A) * (Kb + Kc)
        # Para convección
        H12 = (self.h12 * self.l12 / 6) * np.array([[2,1,0],[1,2,0],[0,0,0]]) 
        H23 = (self.h23 * self.l23 / 6) * np.array([[0,0,0],[0,2,1],[0,1,2]]) 
        H31 = (self.h31 * self.l31 / 6) * np.array([[2,0,1],[0,0,0],[1,0,2]])
        K = K + H12 + H23 + H31
        # Para expandir la matriz
        if expNodos:
            Kexp = np.zeros((expNodos,expNodos))
            for i, j in enumerate(self.n):
                for k, m in enumerate(self.n):
                    Kexp[j - 1,m - 1] = Kexp[j - 1,m - 1] + K[i,k]
            K = Kexp
        return K

    def VectorF(self, expNodos=False):
        '''Calcula el vector F del elemento, incluyendo la generación.'''
        # Vectores de q y hTinf
        F12 = (self.q12 * self.l12 / 2) * np.array([1,1,0]) 
        F23 = (self.q23 * self.l23 / 2) * np.array([0,1,1]) 
        F31 = (self.q31 * self.l31 / 2) * np.array([1,0,1])
        # Generacion
        FQ = (self.Q * self.A / 3) * np.array([1,1,1])
        F = FQ + F12 + F23 + F31
        # Para expandir el vector
        if expNodos:
            Fexp = np.zeros((expNodos))
            for i, j in enumerate(self.n):
                    Fexp[j - 1] = Fexp[j - 1] + F[i]
            F = Fexp
        return F
    
    def Grafica(self):
        '''Función para graficar el elemento '''
        try:
            fig = plt.gcf()
            endplot = False
        except:
            fig = plt.figure(figsize=(5,5))
            endplot = True
        plt.plot(self.x.tolist()+[self.x[0]],self.y.tolist()+[self.y[0]])
        for x, y, n in zip(self.x,self.y,self.n):
             plt.text(x,y,n)
        if endplot:
            plt.grid()
            plt.axis('equal')
            plt.show()

    def out_xlsx(self, archivo, nodos):
        ''' Graba la matriz del elemento en un archivo de excel'''
        import pandas as pd
        T = ['T{}'.format(i+1) for i in range(nodos)]
        empty = np.zeros(nodos)
        empty[:] = np.nan
        out = np.column_stack((self.MatrizK(expNodos=nodos),empty,self.VectorF(expNodos=nodos)))
        df = pd.DataFrame(out, index = T, columns = T + ['', 'F'])
        df.to_excel(archivo)
        