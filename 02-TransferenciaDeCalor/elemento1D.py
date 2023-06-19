import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
sp.init_printing()

class PropiedadesTC:
    def __init__(self, grado, **kargs):
        if 'K' in kargs:
            if isinstance(kargs['K'], (int, float)):
                self.K = kargs['K']
            elif isinstance(kargs['K'], (list, tuple)):
                if grado == 1:
                    self.K = np.mean(np.array(kargs['K']))
            else:
                raise TypeError('El tipo de la variable K no fue reconocido')
        for var in ['Q', 'q1', 'q2', 'T1', 'T2', 'h1', 'h2', 'Tinf1', 'Tinf2']:
            if var in kargs:
                if isinstance(kargs[var], (int, float)):
                    setattr(self, var, kargs[var])
                else:
                    setattr(self, var, None)
            else:
                setattr(self, var, None)

    def to_dict(self, solocero=False):
        out = {}
        for var in ['K','Q','q1','q2','T1','T2','h1','h2','Tinf1','Tinf2']:
            if getattr(self, var) is not None:
                if not solocero:
                    out.update({var:getattr(self, var)})
            else:
                out.update({var:0.0})
        return out
            

def ModeloTC(grado, axisimetrico=False):
    if grado not in (1,):
        raise NotImplementedError('El grado {} no ha sido implementado'.format(grado))
    K, Q, q1, q2, T1, T2, h1, h2, Tinf1, Tinf2 = sp.symbols('K, Q, q1, q2, T1, T2, h1, h2, Tinf1, Tinf2')
    h = sp.symbols('h')
    if not axisimetrico:
        A = K/h*sp.Matrix([[1, -1],[-1, 1]]) + sp.Matrix([[h1,0],[0,0]]) + sp.Matrix([[0,0],[0,h2]])
        b = sp.Matrix([Tinf1*h1+q1,0]) + sp.Matrix([0,Tinf2*h2+q2]) + Q*h*sp.Matrix([1,1])
    else:
        r1, r2 = sp.symbols('r1, r2')
        h = r2 - r1
        A = K*(r1+r2)/(2*h)*sp.Matrix([[1, -1],[-1, 1]]) + sp.Matrix([[r1*h1,0],[0,0]]) + sp.Matrix([[0,0],[0,r2*h2]])
        b = sp.Matrix([(Tinf1*h1+q1)*r1,0]) + sp.Matrix([0,(Tinf2*h2+q2)*r2]) # + Q*h*sp.Matrix([1,1])
    return A, b



class Elemento1D:
    def __init__(self, grado, fisica, nodos_globales=None, numero_nodos=None, x=np.array([0,1]), axisimetrico=False, **kargs):
        if fisica not in ('TC',):
            raise NotImplementedError('La física {} no ha sido implementada'.format(fisica))
        else:
            self.fisica = fisica
        if grado not in (1,):
            raise NotImplementedError('El grado {} no ha sido implementado'.format(grado))
        else:
            self.grado = grado
        self.x = np.linspace(x[0], x[1], num=grado+1)
        self.h = np.abs(float(x[1]-x[0]))
        self.n = nodos_globales if nodos_globales is not None and len(nodos_globales) == len(self.x) else None
        self.num_nodos = numero_nodos
        self.axisimetrico = axisimetrico
        if len(kargs) != 0 and fisica == 'TC':
            self.props = PropiedadesTC(grado, **kargs)

    def out_model(self, simplificado=True):
        if self.fisica == 'TC' and self.grado == 1:
            A, b = ModeloTC(1, axisimetrico=self.axisimetrico)
            return A.subs(self.props.to_dict(solocero=True)), b.subs(self.props.to_dict(solocero=True))
        else:
            raise NotImplementedError

    def out_matrix(self, extend=False):
        if self.fisica == 'TC' and self.grado == 1:
            lmatrix = ModeloTC(1, axisimetrico=self.axisimetrico)[0].subs(self.props.to_dict()).subs({'h': self.h})
            if self.axisimetrico:
                lmatrix = lmatrix.subs({'r1': self.x[0], 'r2': self.x[1]})
            if not extend:
                return lmatrix
            else:
                if self.n is None:
                    raise TypeError('No se especificaron los nodos')
                if self.num_nodos is None:
                    raise TypeError('No se especifico el número de nodos')
                gmatrix = np.zeros([self.num_nodos, self.num_nodos])
                for i in range(2):
                    for j in range(2):
                        gmatrix[self.n[i]-1, self.n[j]-1] = gmatrix[self.n[i]-1, self.n[j]-1] + lmatrix[i,j]
                return gmatrix
        else:
            raise NotImplementedError

    def out_force(self, extend=False):
        if self.fisica == 'TC' and self.grado == 1:
            lmatrix = ModeloTC(1, axisimetrico=self.axisimetrico)[1].subs(self.props.to_dict()).subs({'h': self.h})
            if self.axisimetrico:
                lmatrix = lmatrix.subs({'r1': self.x[0], 'r2': self.x[1]})
            if not extend:
                return lmatrix
            else:
                if self.n is None:
                    raise TypeError('No se especificaron los nodos')
                if self.num_nodos is None:
                    raise TypeError('No se especifico el número de nodos')
                gmatrix = np.zeros([self.num_nodos])
                for i in range(2):
                    gmatrix[self.n[i]-1] = gmatrix[self.n[i]-1] + lmatrix[i]
                return gmatrix
        else:
            raise NotImplementedError

    def Grafica(self):
        '''Función para graficar el elemento '''
        try:
            fig = plt.gcf()
            endplot = False
        except:
            fig = plt.figure(figsize=(6,2))
            endplot = True
        var = 'r' if self.axisimetrico else 'x'
        plt.plot([self.x[0], self.x[-1]],[0, 0],'k')
        tam = self.x[-1] - self.x[0]
        for x, n in zip(self.x, self.n):
            plt.text(x,tam/4,"${}_{{{}}} = {}$".format(var, n,x), ha='center')
            plt.plot((x,x),(-tam/20,tam/20),'k')
        if endplot:
            plt.axis('off')
            plt.axis('equal')
            plt.show()
        
