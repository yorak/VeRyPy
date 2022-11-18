#!/usr/bin/python

try:
    import Tkinter.filedialog as filedialog
    import Tkinter as tk
except ImportError:
    import tkinter.filedialog as filedialog
    import tkinter as tk
from os.path import exists, splitext

from verypy import get_algorithms
from verypy.shared_cli import read_and_solve_a_problem
from verypy.cvrp_io import as_OPT_solution, as_VRPH_solution

algo_list = get_algorithms()

UP_SPACING = 10

class VeRyPyGuiApp:
    def __init__(self):
        self.create_gui()

        self.solution = None
        self.solution_cost = None

    def on_algorithm_information(self, algo_name_desc, algo_f):
        # TODO: This is a hack, better way would be to have a custom info window
        #  with 80 column scrollable text field
        tk.messagebox.showinfo(title=algo_name_desc, message=algo_f.__doc__)

    def on_select_vrp_file(self, store_in_entry):
        filetypes = (
            ('TSPLIB95 VRP files', '*.vrp'),
            ('TSPLIB95 TSP files', '*.tsp'),
            ('Pickled data', '*.pickle'),
            ('text files', '*.txt'),
            ('All files', '*.*')
        )

        filename = filedialog.askopenfilename(
            title='Open a problem file',
            filetypes=filetypes)

        if exists(filename):
            store_in_entry.delete(0, tk.END)
            store_in_entry.insert(0, filename)


    def on_solve_vrp_file(self, problem_filename, algorithm_f):
        try:
            sol, _, objf_c, elapsed = \
                read_and_solve_a_problem(problem_filename, algorithm_f, minimize_K=False)

            self.result_label_var.set("Solved in %.1f s to value %.1f (K=%d)"%
                                (elapsed, objf_c, sol.count(0)-1))
            self.solution = sol
            self.solution_cost = objf_c
        except Exception as e:
            tk.messagebox.showerror(title="VeRyPy solve error", message=str(e))

    def on_export_solution(self):
        if not self.solution:
            tk.messagebox.showerror(title="VeRyPy export error",
                message="There is no solution to export.\nPress solve first.")

        filetypes = (
            ('OPT file', '*.opt'),
            ('SOL (VRPH) file', '*.sol'),
        )
        filename = filedialog.asksaveasfilename(
            title='Save a solution file',
            filetypes=filetypes)
        with open(filename, 'w') as wf:
            if splitext(filename)[1].lower()==".opt":
                wf.write(as_OPT_solution(self.solution_cost, self.solution))
            else:
                wf.write(as_VRPH_solution(self.solution))

    def create_gui(self):
        self.root = tk.Tk()
        self.result_label_var = tk.StringVar()
        self.result_label_var.set('Choose problem and press solve!')

        self.root.title('VeRyPy GUI')
        self.root.resizable(False, False)
        self.root.geometry('320x440')

        header = tk.Frame(self.root, width=400)#, bg='grey')
        header.grid(row=0, column=0, padx=UP_SPACING, pady=UP_SPACING/2)
        vrp_file_field = tk.Entry(header, width=32)
        vrp_file_field.insert(0, "<your problem file here>")
        vrp_file_field.pack(expand=True)
        vrp_file_field.grid(row=0, column=0, sticky="nsew")
        open_button = tk.Button(header, text='...', 
            command=lambda:self.on_select_vrp_file(vrp_file_field))
        open_button.grid(row=0, column=1)

        algof_from_option = {name+":"+desc:f for _, name, desc, f in algo_list}
        algo_options = list(algof_from_option.keys())
        active_algo = tk.StringVar(self.root, algo_options[0], 'PY_ALGO')
        algo_option_menu = tk.OptionMenu(header, active_algo, *algo_options)
        algo_option_menu.config(width=32)
        algo_option_menu.configure(anchor='w')
        algo_option_menu.grid(row=1, column=0)
        info_button = tk.Button(header, text='🛈',
            command=lambda:self.on_algorithm_information(active_algo.get(),
                algof_from_option[active_algo.get()]))
        info_button.grid(row=1, column=1)

        solve_button = tk.Button(header, text='solve',
            command=lambda:self.on_solve_vrp_file(
                vrp_file_field.get(),
                algof_from_option[active_algo.get()]))
        solve_button.grid(row=3, column=0, columnspan=2, sticky="news")

        canvas = tk.Canvas(self.root, bg="white", height=300, width=300)
        canvas.grid(row=1, column=0, padx=UP_SPACING, pady=UP_SPACING/2)
        canvas.create_text(150, 50, text="Placeholder, TODO:\nimplement solution visualization",
            fill="black", font=('Helvetica 12 bold'))
        
        footer = tk.Frame(self.root, width=400)#, bg='grey')
        footer.grid(row=2, column=0)
        result_label = tk.Label(footer, textvariable=self.result_label_var, borderwidth=1, relief="solid")
        result_label.grid(row=0, column=0, sticky="news", padx=UP_SPACING, pady=UP_SPACING/2)
        export_button = tk.Button(footer, text='export',
            command=lambda:self.on_export_solution())
        export_button.grid(row=0, column=1)


        self.root.iconbitmap("doc/logo_nodes.ico")
    
    def run(self):
        self.root.mainloop()

if __name__=="__main__":
    app = VeRyPyGuiApp()
    app.run()