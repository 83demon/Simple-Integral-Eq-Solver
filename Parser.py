import numpy as np
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr

class Parser:
    def __init__(self,m,n,matrix,b):
        self.m = int(m)
        self.n = int(n)
        self.raw_matrix = matrix
        self.raw_b = b
        self.matrix = sym.zeros(m,n)
        self.b = sym.zeros(m,1)

    def _parse_matrix(self):
        for i in range(self.m*self.n):
            self.matrix[i//self.n,i%self.n] = parse_expr(self.raw_matrix[i])

    def _parse_b_vector(self):
        for i in range(self.m):
            self.b[i] = parse_expr(self.raw_b[i])

    def main(self):
        self._parse_matrix()
        self._parse_b_vector()
        return self.matrix, self.b