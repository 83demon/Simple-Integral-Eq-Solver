import PySimpleGUI as sg
from Parser import Parser
from Solver import Solver
import numpy as np
import matplotlib.pyplot as plt


class GUI:
    def __init__(self):
        self._default_dim_value = 2
        self.title = "Simple Matrix Equation Solver"
        self._greeting_text = """Hi! This is a simple matrix equation solver. \n
                                    Please write down height and width of your matrix:"""
        self._input_matrix_text = "Please, write down your matrix data."
        self._input_b_vector_text = "Please, write down your b vector data."
        self._input_T_text = "Pleasem write down your T value."
        self._result_text = "Your result is:"
        self.layout = [[sg.Text(self._greeting_text)],
                       [sg.Input(default_text=self._default_dim_value),sg.Input(default_text=self._default_dim_value)],
                       [sg.Submit()]]
        self._m_n = None
        self._m = None
        self._n = None
        self._marker_width = 3
        self._input_box_size = (6,1)
        self._output_box_size = (35,1)
        self._default_value = '0'
        self.color = 'forestgreen'
        self.handler()

    def _first_run(self):
        self.window = sg.Window(self.title, self.layout)
        event, self._m_n = self.window.read()
        self._m = int(self._m_n[0])
        self._n = int(self._m_n[1])
        sg.popup(f'You entered: \nHeight: {self._m}\t Width: {self._n}')
        self.window.close()

    def _construct_matrix_input(self):
        self.layout = [[sg.Text(self._input_matrix_text)]]
        for i in range(int(self._m)):
            temp = [sg.Input(size=self._input_box_size,default_text=self._default_value) for _ in range(int(self._n))]
            self.layout.append(temp)
        self.layout.append([sg.Submit()])
        self.window = sg.Window(self.title, self.layout)
        event, self.raw_matrix = self.window.read()
        #print("You entered:", self.raw_matrix)
        self.window.close()

    def _construct_b_vector_input(self):
        self.layout = [[sg.Text(self._input_b_vector_text)]]
        for i in range(int(self._m)):
            temp = [sg.Input(size=self._input_box_size,default_text=self._default_value)]
            self.layout.append(temp)
        self.layout.append([sg.Submit()])
        self.window = sg.Window(self.title, self.layout)
        event, self.raw_b_vector = self.window.read()
        #print("You entered:", self.raw_b_vector)
        self.window.close()

    def _construct_T_input(self):
        self.layout = [[sg.Text(self._input_T_text)],
                       [sg.Input(size=self._input_box_size,default_text=self._default_value)],
                       [sg.Submit()]]
        self.window = sg.Window(self.title, self.layout)
        event, self.T = self.window.read()
        #print("You entered:", self.raw_b_vector)
        self.window.close()

    def _show_result(self, solution, det, epsilon, unity_check):
        self.layout = [[sg.Text(self._result_text)]]
        for i in range(int(self._n)):
            temp = [sg.Text(size=self._output_box_size,text=f"{solution[key][i]}") for key in solution.keys()]
            self.layout.append(temp)

        self.layout.append([sg.Text("")])
        k=1
        for v in solution.keys():
            self.layout.append([sg.Text("v"+str(k)+": " + " ".join([str(v[i]) for i in range(self._n)]))])
            k+=1
        #self.layout.append([sg.Text(f"Pseudo incerse matrix is:")])

        #pseudo inverse matrix has shape (n,m) while original one has shape of (m,n)
        #for i in range(self._n):
        #    self.layout.append([sg.Text(" ".join([str(pinv_matrix[i][j]) for j in range(self._m)]))])
        self.layout.append([sg.Text(f"Epsilon^2: {epsilon[0,0]}")])
        self.layout.append([sg.Text(f"There is only one solution: {unity_check}")])
        self.layout.append([sg.Text(f"Det: {det}")])
        self.layout.append([sg.CloseButton("Bye!")])
        self.window = sg.Window(self.title, self.layout)

        while True:

            event, values = self.window.read()

            if event == sg.WIN_CLOSED:
                break

        self.window.close()

    def _show_graphs(self,solution,unity_flag):
        raise NotImplementedError()


    def handler(self):
        self._first_run()
        if self._n<=0 or self._m <=0:
            raise ValueError(f'Cannot construct matrix with shape {self._m,self._n}')
        else:
            self._construct_matrix_input()
            self._construct_b_vector_input()
            self._construct_T_input()
            parser = Parser(self._m,self._n,self.raw_matrix,self.raw_b_vector,self.T)
            solver = Solver(*parser.main())
            solution, eps, det, unity_check = solver.main()
            console_output(solution, det, eps, unity_check)
            self._show_result(solution, det, eps, unity_check)
        """else:
            matrix, b , T = read_from_console(self._m, self._n)
            parser = Parser(self._m,self._n,matrix,b,T)
            solver = Solver(*parser.main())
            solution, det, eps, unity_check = solver.main()
            console_output(solution, det, eps, unity_check)
            self._show_graphs(solution,unity_check)"""




def read_from_console(m,n):
    print(f"Enter a {(m,n)} matrix data separated by whitespace:\n")
    matrix_data = []
    for i in range(m):
        temp = input().split(" ")
        if temp[-1] == "":
            temp.pop()
        if len(temp)!=n:
            raise ValueError(f"You entered {len(temp)} values")
        matrix_data.append(temp)


    print(f"Enter a {(m,1)} b vector data separated by whitespace:\n")
    b = input().split(" ")
    print("Enter your T value:\n")
    T = {0:input()}

    return matrix_data, b, T


def console_output(solution, det, eps, unity_check):
    print("Your solution:\n")
    for k,v in solution.items():
        print(f"{k} :\n{v}")
    print()
    #print("Pseudo inverted matrix is:\n")
    #for i in range(pinv_matrix.shape[0]):
    #   print(f"{' '.join(str(pinv_matrix[i][j]) for j in range(pinv_matrix.shape[1]))}")
    print()
    print(f'Epsilon^2 is: {eps}')
    print(f'There is only one solution: {unity_check}')
    print()