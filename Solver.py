import numpy as np
import sympy as sym
from sympy.core.rules import Transform
from sympy import *


class Solver:

    def __init__(self, matrix, vector_column, T, start = 0):
        self.m = matrix.shape[0]  # height
        self.n = matrix.shape[1]  # width
        self.A = matrix
        self.b = vector_column
        self.start = start
        self.T = T  # upper bound of time variable
        self.t_symbol = sym.Symbol('t')
        self.t_i_symbol = sym.Symbol('t_i')
        self.t_j_symbol = sym.Symbol('t_j')
        self.P_1 = None
        self.P_1_inv = None
        self.A_nu_values = []
        self.det = None # value of a unity check
        self.epsilon = np.array([[0]])  # accuracy
        self._unity_flag = None  # flag to indicate unity of the solution
        self.n_steps = 1 # number of steps to split an interval [0;T]
        self.ndigits=3

        self._nu_set_length = 5
        self._nu_low = -10
        self._nu_high = 10
        self.nu_values = []
        self.solution = dict()

    def _init_nu(self):
        vals = [-2,-1,0,1,2]
        for i in range(self._nu_set_length):
            self.nu_values.append(
                sym.Matrix([[str(vals[i]) + "*t"] for _ in range(self.n)]))
        #    self.nu_values.append(sym.Matrix([[str(np.random.randint(self._nu_low,self._nu_high))+"*t"] for _ in range(self.n)]))

    @staticmethod
    def pinv(matrix: sym.Matrix):
        matrix_np = np.array(matrix,dtype=np.float64)
        matrix_np_inv = np.linalg.pinv(matrix_np)
        return sym.Matrix(matrix_np_inv)

    def _compute_p_1(self):
        self.P_1 = sym.integrate(self.A@self.A.T,(self.t_symbol,self.start,self.T)).evalf()

    def _compute_p_1_inv(self):
        self.P_1_inv = self.pinv(self.P_1)

    def _compute_a_nu(self):
        for nu in self.nu_values:
            self.A_nu_values.append(sym.integrate(self.A@nu,(self.t_symbol,self.start,self.T)))

    def _check_for_unity(self, N):
        """Checks a solution for unity"""

        def split_interval(start,end,steps):
            return [start+(i+1)*(end-start)/steps for i in range(steps)]


        N_values = split_interval(self.start,self.T,N)
        big_matrix = sym.zeros(N*self.n,N*self.n)
        for i in range(N):
            for j in range(N):
                temp = self.A.T.evalf(subs={self.t_symbol:N_values[i]}) @ self.A.evalf(subs={self.t_symbol:N_values[j] })
                big_matrix[i*self.n:self.n*(i+1) , j*self.n:self.n*(j+1)] = temp
        np_matrix = np.array(big_matrix,dtype=np.float64)
        self.det = np.linalg.det(np_matrix)
        print(np_matrix)
        print(f"Det: {self.det}")
        return self.det>0

    def _calculate_accuracy(self):
        if not self._unity_flag:
            pass
        self.epsilon = self.b.T @ self.b - self.b.T @ self.P_1 @ self.P_1_inv @ self.b

    def _solve(self):
        for nu in range(self._nu_set_length):
            self.solution[tuple(self.nu_values[nu])] = sym.simplify(self.A.T@self.P_1_inv@self.b + self.nu_values[nu] - self.A.T@self.P_1_inv@self.A_nu_values[nu]).evalf()
            self.solution[tuple(self.nu_values[nu])] = (self.solution[tuple(self.nu_values[nu])]).xreplace(Transform(lambda x: x.round(self.ndigits), lambda x: isinstance(x, Float)))
            # to round?

    def main(self):
        self._init_nu()
        self._compute_p_1()
        self._compute_p_1_inv()
        self._compute_a_nu()
        self._unity_flag = self._check_for_unity(self.n_steps)
        self._calculate_accuracy()
        self._solve()
        return self.solution, self.epsilon, self.det, self._unity_flag


if __name__ == "__main__":
    t = sym.Symbol('t')
    matrix = sym.Matrix([[t,2*t],[sin(t),cos(t)]])
    b = np.array([[1],[16]])
    T = 1
    solver = Solver(matrix,b,T)
    res, eps, det, unity_flag = solver.main()
    print(eps)
    for v in res.values():
        for i in range(matrix.shape[1]):
            print(v[i])
        print()

