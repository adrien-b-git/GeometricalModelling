"""
author: Adrien Bouchet
description: this script is a tkinter GUI curve editor which allows user to:
-add/remove point on  the canvas
-drag the points on the canvas
-select and edit coordinates of a point on th canvas (validate with enter)
-choose the algorithm used to compute the curve
-change the resolution of final linear interpolation of the curve
-reset button clear all points

To add algorithm, see the __init__ doc of CurveEditorWindow
"""
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import math

class CurveEditorWindow(Tk):
    def __init__(self, compute_algorithms, derivation_methods) -> None:
        """compute_algorithms is a list of tuple {"name" : string, "algo" : function) which 
        -string is the text in button selection of algorithm
        -function has the signature: def <name>(points, T) 
            where points is the list of points drawn on the canvas and T the array of t in [0; 1] with resolution given by the slider"""
        super().__init__()

        self.title("Curvy editor")
        self.geometry("750x520")

        self._selected_data = {"x": 0, "y": 0, 'item': None}
        self.radius = 5
        self.curve = None
        self.compute_algorithms = compute_algorithms
        self.derivation_methods = derivation_methods

        self.rowconfigure([0, 1], weight=1, minsize=200)
        self.columnconfigure([0, 1], weight=1, minsize=200)

        self.setup_canvas()
        self.setup_panel()

    # Setup
    def setup_canvas(self):
        self.graph = Canvas(self, bd=2, cursor="plus", bg="#fbfbfb")
        self.graph.grid(column=0,
                        padx=2,
                        pady=2,
                        rowspan=2,
                        columnspan=2,
                        sticky="nsew")

        self.graph.bind('<Button-1>', self.handle_canvas_click)
        self.graph.tag_bind("control_points", "<ButtonRelease-1>",
                            self.handle_drag_stop)
        self.graph.bind("<B1-Motion>", self.handle_drag)

    def setup_panel(self):
        # Right panel for options
        self.frame_pannel = Frame(self, relief=RAISED, bg="#e1e1e1")
        self.frame_curve_type = Frame(self.frame_pannel)
        self.frame_edit_type = Frame(self.frame_pannel)
        self.frame_edit_position = Frame(self.frame_pannel)
        self.frame_sliders = Frame(self.frame_pannel)
        self.frame_derivation_modes = Frame(self.frame_pannel)

        # Selection of curve type
        curve_types = [algo['name'] for algo in self.compute_algorithms]
        curve_types_val = list(range(len(self.compute_algorithms)))
        self.curve_type = IntVar()
        self.curve_type.set(curve_types_val[0])

        self.radio_curve_buttons = [None] * len(self.compute_algorithms)
        for i in range(len(self.compute_algorithms)):
            self.radio_curve_buttons[i] = Radiobutton(self.frame_curve_type,
                                                      variable=self.curve_type,
                                                      text=curve_types[i],
                                                      value=curve_types_val[i],
                                                      bg="#e1e1e1")
            self.radio_curve_buttons[i].pack(side='left', expand=1)
            self.radio_curve_buttons[i].bind(
                "<ButtonRelease-1>",
                lambda event: self.graph.after(100, lambda: self.draw_curve()))

        # Selection of edit mode
        edit_types = ['Add', 'Remove', 'Drag', 'Select']
        edit_types_val = ["add", "remove", "drag", "select"]
        self.edit_types = StringVar()
        self.edit_types.set(edit_types_val[0])

        self.radio_edit_buttons = [None] * 4
        for i in range(4):
            self.radio_edit_buttons[i] = Radiobutton(self.frame_edit_type,
                                                     variable=self.edit_types,
                                                     text=edit_types[i],
                                                     value=edit_types_val[i],
                                                     bg="#e1e1e1")
            self.radio_edit_buttons[i].pack(side='left', expand=1)
            self.radio_edit_buttons[i].bind(
                "<ButtonRelease-1>", lambda event: self.reset_selection())

        # Selection of derivation mode
        derivation_modes = [method['name'] for method in self.derivation_methods]
        derivation_modes_val = list(range(len(self.derivation_methods)))
        self.derivation_modes = IntVar()
        self.derivation_modes.set(derivation_modes_val[0])

        self.radio_deriv_buttons = [None] * len(self.derivation_methods)
        for i in range(len(self.derivation_methods)):
            self.radio_deriv_buttons[i] = Radiobutton(self.frame_derivation_modes,
                                                     variable=self.derivation_modes,
                                                     text=derivation_modes[i],
                                                     value=derivation_modes_val[i],
                                                     bg="#e1e1e1")
            self.radio_deriv_buttons[i].pack(side='left', expand=1)
            self.radio_deriv_buttons[i].bind(
                "<ButtonRelease-1>", lambda event: self.graph.after(100, lambda: self.draw_curve()))

        # Edit position of selected point widget
        self.label_pos_x = Label(self.frame_edit_position, text='x: ')
        self.label_pos_y = Label(self.frame_edit_position, text='y: ')
        self.pos_x = StringVar()
        self.pos_y = StringVar()
        self.entry_position_x = Entry(self.frame_edit_position,
                                      textvariable=self.pos_x)
        self.entry_position_y = Entry(self.frame_edit_position,
                                      textvariable=self.pos_y)
        self.label_pos_x.grid(row=0, column=0)
        self.entry_position_x.grid(row=0, column=1)
        self.label_pos_y.grid(row=1, column=0)
        self.entry_position_y.grid(row=1, column=1)

        self.entry_position_x.bind("<FocusOut>", self.update_pos)
        self.entry_position_x.bind("<KeyPress-Return>", self.update_pos)
        self.entry_position_x.bind("<KeyPress-KP_Enter>", self.update_pos)

        self.entry_position_y.bind("<FocusOut>", self.update_pos)
        self.entry_position_y.bind("<KeyPress-Return>", self.update_pos)
        self.entry_position_y.bind("<KeyPress-KP_Enter>", self.update_pos)

        # Slider for parameter update
        # resolution
        self.label_resolution = Label(self.frame_sliders, text="Resolution: ")
        self.slider_resolution = Scale(self.frame_sliders,
                                       from_=5,
                                       to=500,
                                       orient=HORIZONTAL,
                                       bg="#e1e1e1")
        self.slider_resolution.set(50)
        self.label_resolution.grid(row=0, column=0)
        self.slider_resolution.grid(row=0, column=1)
        self.slider_resolution.bind("<ButtonRelease-1>",
                                    lambda event: self.draw_curve())
        # tension parameter
        self.label_tension = Label(self.frame_sliders, text="Tension: ")
        self.slider_tension = Scale(self.frame_sliders,
                                       from_=0,
                                       to=100,
                                       orient=HORIZONTAL,
                                       bg="#e1e1e1")
        self.slider_tension.set(1/2)
        self.label_tension.grid(row=1, column=0)
        self.slider_tension.grid(row=1, column=1)
        self.slider_tension.bind("<ButtonRelease-1>",
                                    lambda event: self.draw_curve())
        
        # tension parameter for Kochanek–Bartels spline
        self.label_tensionKB = Label(self.frame_sliders, text="TensionKB: ")
        self.slider_tensionKB = Scale(self.frame_sliders,
                                       from_=0,
                                       to=200,
                                       orient=HORIZONTAL,
                                       bg="#e1e1e1")
        self.slider_tensionKB.set(1/2)
        self.label_tensionKB.grid(row=2, column=0)
        self.slider_tensionKB.grid(row=2, column=1)
        self.slider_tensionKB.bind("<ButtonRelease-1>",
                                    lambda event: self.draw_curve())

        # bias parameter for Kochanek–Bartels spline
        self.label_biasKB = Label(self.frame_sliders, text="biasKB: ")
        self.slider_biasKB = Scale(self.frame_sliders,
                                       from_=0,
                                       to=200,
                                       orient=HORIZONTAL,
                                       bg="#e1e1e1")
        self.slider_biasKB.set(1/2)
        self.label_biasKB.grid(row=3, column=0)
        self.slider_biasKB.grid(row=3, column=1)
        self.slider_biasKB.bind("<ButtonRelease-1>",
                                    lambda event: self.draw_curve())
                                
        # continuity parameter for Kochanek–Bartels spline
        self.label_continuityKB = Label(self.frame_sliders, text="continuityKB: ")
        self.slider_continuityKB = Scale(self.frame_sliders,
                                       from_=0,
                                       to=200,
                                       orient=HORIZONTAL,
                                       bg="#e1e1e1")
        self.slider_continuityKB.set(1/2)
        self.label_continuityKB.grid(row=4, column=0)
        self.slider_continuityKB.grid(row=4, column=1)
        self.slider_continuityKB.bind("<ButtonRelease-1>",
                                    lambda event: self.draw_curve())
        

        self.frame_pannel.grid(row=0,
                               column=2,
                               padx=2,
                               pady=2,
                               rowspan=2,
                               sticky="nswe")
        self.frame_curve_type.pack()
        self.frame_edit_type.pack()
        self.frame_derivation_modes.pack()
        self.frame_edit_position.pack()
        self.frame_sliders.pack()

        self.button_reset = Button(self.frame_pannel, text="Reset")
        self.button_reset.pack(side=BOTTOM, fill="x")
        self.button_reset.bind("<ButtonRelease-1>",
                               lambda event: self.graph.delete("all"))
        self.button_all = Button(self.frame_pannel, text="Draw all")
        self.button_all.pack(side=TOP, fill="x")
        self.button_all.bind("<ButtonRelease-1>",
                               lambda event: self.draw_curve(True)) 
        self.button_courbure = Button(self.frame_pannel, text="Courbure")
        self.button_courbure.pack(side=TOP, fill="x")
        self.button_courbure.bind("<ButtonRelease-1>",
                               lambda event: self.draw_courbure())
        self.all_courbure = Button(self.frame_pannel, text="All courbure")
        self.all_courbure.pack(side=TOP, fill="x")
        self.all_courbure.bind("<ButtonRelease-1>",
                               lambda event: self.draw_courbure(True))                    

    # Drawing
    def get_points(self):
        points = []
        for item in self.graph.find_withtag("control_points"):
            coords = self.graph.coords(item)
            points.append([
                float(coords[0] + self.radius),
                float(coords[1] + self.radius)
            ])  # Ensure curve accuracy
        return points

    def create_point(self, x, y, color):
        """Create a token at the given coordinate in the given color"""
        item = self.graph.create_oval(x - self.radius,
                                      y - self.radius,
                                      x + self.radius,
                                      y + self.radius,
                                      outline=color,
                                      fill=color,
                                      tags="control_points")
        return item

    def draw_polygon(self):
        self.graph.delete("control_polygon")
        points = self.get_points()
        for i in range(0, len(points) - 1):
            self.graph.create_line(points[i][0],
                                   points[i][1],
                                   points[i + 1][0],
                                   points[i + 1][1],
                                   fill="blue",
                                   tags="control_polygon")
    

    def draw_curve(self, all=False):
        self.graph.delete("curve")
        couleurs = ["green", "blue", "red", "black", "yellow"]
        if (all == True):
            for j in range(len(self.compute_algorithms)):
                points = self.get_points()
                if len(points) <= 1:
                    return
                T = np.linspace(0, 1, self.slider_resolution.get())
                result = self.compute_algorithms[j]['algo'](
                np.array(points), T, self.slider_tension.get()/100, -1+self.slider_tensionKB.get()/100, -1+self.slider_biasKB.get()/100, -1+self.slider_continuityKB.get()/100, self.derivation_methods[self.derivation_modes.get()]['method'])

                self.curve = np.array(result)

                for i in range(0, self.curve.shape[0] - 1):
                    self.graph.create_line(self.curve[i, 0],
                                   self.curve[i, 1],
                                   self.curve[i + 1, 0],
                                   self.curve[i + 1, 1],
                                   fill= couleurs[j],
                                   width=3,
                                   tags="curve")
        else :     
            points = self.get_points()
            if len(points) <= 1:
                return
            T = np.linspace(0, 1, self.slider_resolution.get())
            result = self.compute_algorithms[self.curve_type.get()]['algo'](
            np.array(points), T, self.slider_tension.get()/100, -1+self.slider_tensionKB.get()/100, -1+self.slider_biasKB.get()/100, -1+self.slider_continuityKB.get()/100, self.derivation_methods[self.derivation_modes.get()]['method'])

            self.curve = np.array(result)

            for i in range(0, self.curve.shape[0] - 1):
                self.graph.create_line(self.curve[i, 0],
                                   self.curve[i, 1],
                                   self.curve[i + 1, 0],
                                   self.curve[i + 1, 1],
                                   fill=couleurs[self.curve_type.get()],
                                   width=3,
                                   tags="curve")

    def draw_courbure(self, all=False):
        couleurs = ['g', 'b', 'r', 'k', 'y']
        if (all==False):
            courbe = self.curve.tolist()
            T = np.linspace(0, 1, self.slider_resolution.get())
            self.courbure(self.derivate1(courbe, T)[0], self.derivate1(courbe, T)[1], self.derivate2(courbe, T)[0], self.derivate2(courbe, T)[1], T, couleurs[self.curve_type.get()])
            plt.show()  
        else:
            for j in range(len(self.compute_algorithms)):
                points = self.get_points()
                if len(points) <= 1:
                    return
                T = np.linspace(0, 1, self.slider_resolution.get())
                result = self.compute_algorithms[j]['algo'](
                np.array(points), T, self.slider_tension.get()/100, -1+self.slider_tensionKB.get()/100, -1+self.slider_biasKB.get()/100, -1+self.slider_continuityKB.get()/100, self.derivation_methods[self.derivation_modes.get()]['method'])
                self.courbure(self.derivate1(result, T)[0], self.derivate1(result, T)[1], self.derivate2(result, T)[0], self.derivate2(result, T)[1], T, couleurs[j])
            plt.show()
                    


    # Event handling
    def find_closest_with_tag(self, x, y, radius, tag):
        distances = []
        for item in self.graph.find_withtag(tag):
            c = self.graph.coords(item)
            d = (x - c[0])**2 + (y - c[1])**2
            if d <= radius**2:
                distances.append((item, c, d))

        return min(distances,
                   default=(None, [0, 0], float("inf")),
                   key=lambda p: p[2])

    def reset_selection(self):
        if self._selected_data['item'] is not None:
            self.graph.itemconfig(self._selected_data['item'], fill='red')

        self._selected_data['item'] = None
        self._selected_data["x"] = 0
        self._selected_data["y"] = 0

    def handle_canvas_click(self, event):
        self.reset_selection()

        if self.edit_types.get() == "add":
            item = self.create_point(event.x, event.y, "red")
            self.update_pos_entry(item)

            points = self.get_points()

            if len(points) > 1:
                self.graph.create_line(points[-2][0],
                                       points[-2][1],
                                       points[-1][0],
                                       points[-1][1],
                                       fill="blue",
                                       tag="control_polygon")
                self.draw_curve()

        elif self.edit_types.get() == "remove":
            self._selected_data[
                'item'], coords, _ = self.find_closest_with_tag(
                    event.x, event.y, 3 * self.radius, "control_points")
            if self._selected_data['item'] is not None:
                self.graph.delete(self._selected_data['item'])

                self.draw_polygon()
                self.draw_curve()

        elif self.edit_types.get() == "drag":
            self._selected_data[
                'item'], coords, _ = self.find_closest_with_tag(
                    event.x, event.y, 3 * self.radius, "control_points")

            if self._selected_data['item'] is not None:
                self._selected_data["x"] = event.x
                self._selected_data["y"] = event.y
                self.graph.move(self._selected_data['item'],
                                event.x - coords[0] - self.radius,
                                event.y - coords[1] - self.radius)

        else:
            self._selected_data[
                'item'], coords, _ = self.find_closest_with_tag(
                    event.x, event.y, 3 * self.radius, "control_points")
            if self._selected_data['item'] is not None:
                self.graph.itemconfig(self._selected_data['item'],
                                      fill='orange')  # Mark as selected
                self.update_pos_entry(self._selected_data['item'])

    def handle_drag_stop(self, event):
        """End drag of an object"""
        if self.edit_types.get() != "drag":
            return
        self.reset_selection()

    def handle_drag(self, event):
        """Handle dragging of an object"""
        if self.edit_types.get() != "drag" or self._selected_data[
                'item'] is None or "control_points" not in self.graph.gettags(
                    self._selected_data['item']):
            return

        # compute how much the mouse has moved
        delta_x = event.x - self._selected_data["x"]
        delta_y = event.y - self._selected_data["y"]
        # move the object the appropriate amount
        self.graph.move(self._selected_data['item'], delta_x, delta_y)
        # record the new position
        self._selected_data["x"] = event.x
        self._selected_data["y"] = event.y

        self.update_pos_entry(self._selected_data['item'])
        self.draw_polygon()
        self.draw_curve()

    def update_pos_entry(self, item):
        coords = self.graph.coords(item)
        self.entry_position_x.delete(0, END)
        self.entry_position_x.insert(0, int(coords[0]))
        self.entry_position_y.delete(0, END)
        self.entry_position_y.insert(0, int(coords[1]))

    def update_pos(self, event):
        if self.edit_types.get(
        ) != "select" or self._selected_data['item'] is None:
            return

        coords = self.graph.coords(self._selected_data['item'])
        self.graph.move(self._selected_data['item'],
                        float(self.pos_x.get()) - coords[0],
                        float(self.pos_y.get()) - coords[1])

        self.draw_polygon()
        self.draw_curve()

    def derivate1(self, vecteur, T):
        """
        Calcul de la dérivée du vecteur
        - vecteur est une liste de points en  [:][0] pour avoir tous les x, [:][1] pour avoir tous les y
        - T est une liste avec tous les points
        On prend en entrée le résultat d'un des deux algos Hermite ou Casteljau
        Formule pour (u(x+h)-u(x-h))/2h pour u'(x) à l'ordre 2 sachant que la première et la dernière dérivée sont données par le m
        """
        derivee1_x = []
        derivee1_y = []
        pas = T[1] - T[0]   # pas est la valeur du h
        n = len(vecteur)
        # gérer la première dérivée avec (u(x+h)-u(x))/2
        derivee1_x.append((vecteur[1][0] - vecteur[0][0])/pas) 
        derivee1_y.append((vecteur[1][1] - vecteur[0][1])/pas) 
        # gérer les points du milieu
        for i in range(1, n-1):
            derivee1_x.append(
                (vecteur[i+1][0] - vecteur[i-1][0])/2*pas
            )  # Ensure curve accuracy
            derivee1_y.append(
                (vecteur[i+1][1] - vecteur[i-1][1])/2*pas
            )  # Ensure curve accuracy
        # gérer la dernière dérivée avec (u(x+h)-u(x))/2
        derivee1_x.append((vecteur[n-1][0] - vecteur[n-2][0])/pas) 
        derivee1_y.append((vecteur[n-1][1] - vecteur[n-2][1])/pas) 
        #plt.plot(derivee1_x, derivee1_y)
        #plt.show()
        return derivee1_x, derivee1_y

    def derivate2(self, vecteur, T):
        """
        calcul de la dérivée seconde du vecteur
        pas est la valeur du h
        - vecteur est une liste de points en  [:][0] pour avoir tous les x, [:][1] pour avoir tous les y
        - T est une liste avec tous les points
        On prend en entrée le résultat d'un des deux algos Hermite ou Casteljau
        Formule pour (u(x+2h)-2u(x+h)+u(x))/h**2 pour u''(x) à l'ordre 2 toujours.
        Pour les premières et dernières dérivées on change la formule pour prendre en compte les bons points
        """
        derivee2_x = []
        derivee2_y = []
        pas = T[1] - T[0]   # pas est la valeur du h
        n = len(vecteur)
        # gérer les premiers points
        for i in range(0, n-2):
            derivee2_x.append(
                (vecteur[i+2][0] - 2*vecteur[i+1][0] + vecteur[i][0])/pas**2
            )  # Ensure curve accuracy
            derivee2_y.append(
                (vecteur[i+2][1] - 2*vecteur[i+1][1] + vecteur[i][1])/pas**2
            )  # Ensure curve accuracy
        # gérer l'avant dernière dérivée
        derivee2_x.append((vecteur[n-2][0] - 2*vecteur[n-3][0] + vecteur[n-4][0])/pas**2) 
        derivee2_y.append((vecteur[n-2][1] - 2*vecteur[n-3][1] + vecteur[n-4][1])/pas**2)

        # gérer la dernière dérivée
        derivee2_x.append((vecteur[n-1][0] - 2*vecteur[n-2][0] + vecteur[n-3][0])/pas**2) 
        derivee2_y.append((vecteur[n-1][1] - 2*vecteur[n-2][1] + vecteur[n-3][1])/pas**2)
        #plt.plot(derivee2_x, derivee2_y)
        #plt.show()
        return  derivee2_x, derivee2_y

    def courbure(self, x_prime_x, x_prime_y, x_seconde_x, x_seconde_y, T, couleur = 'b'):
        # Ÿ(t) = det( [xÕ (t),xÕÕ (t)] ) / ÎxÕ (t)Î3
        n = len(x_prime_x)
        courbure = []
        nb_morceaux_courbes = n // len(T)
        # liste_intermédiaire = T # cette liste de points permet d'avoir le bon domaine de définition pour la courbure k(t)
        pas = T[1] - T[0]   # pas est la valeur du h
        domaine_definition = np.linspace(0, nb_morceaux_courbes, len(T)*nb_morceaux_courbes)
        for i in range(n):
            courbure.append((x_prime_x[i]*x_seconde_y[i] - x_prime_y[i]*x_seconde_x[i])/np.linalg.norm((x_prime_x[i], x_prime_y[i]))**3)
        plt.plot(domaine_definition, courbure, couleur)
        return
    

# ---------- Here add algorithms ----------


def DeCasteljau(points, T, c, tension, bias, continuity, tangente_calc): 
    n = points.shape[0] - 1
    result = []
    for t in T:
        r = points.copy()
        for k in range(0, n):
            for i in range(0, n - k):
                r[i, :] = (1 - t) * r[i, :] + t * r[i + 1, :]

        result.append(r[0, :])
    return np.array(result)

def Aitken_Neville(points, T, c, tension, bias, continuity, tangente_calc): 
    n = points.shape[0] - 1
    nb_morceaux_courbes = n
    result = []
    definition_domain = np.linspace(0, nb_morceaux_courbes, len(T)*nb_morceaux_courbes)
    for t in definition_domain:
        r = points.copy()
        for k in range(1, n+1):
            for i in range(0, n+1 - k):
                r[i, :] = (i+k - t)/k * r[i, :] + (t-i)/k * r[i + 1, :]

        result.append(r[0, :])
    return np.array(result)

def HermiteSpline(points, T, c, tension, bias, continuity, tangente_calc): 
    n = points.shape[0] - 1
    result = []
    if n < 1:
        # Cas avec zéro point
        return np.array(result)
    if n == 1: 
        # Cas avec seulement deux points
        # Gestion de la première courbe
        y_1 = 2*points[0, :] - points[1, :] # Calcul du point y_{-1} 
        y_2 = 2*points[1, :] - points[0, :] # Calcul du point y_{2} 
        y_3 = 2*y_2 - points[1, :] # Calcul du point y_{3} 
        m0 = tangente_calc(points[1, :], y_1, c, points[0, :], y_2, tension, bias, continuity, 0)
        m1 = tangente_calc(y_2, points[0, :], c, points[1, :], y_3, tension, bias, continuity, 1)
        for t in T:
            t2 = t * t
            t3 = t2 * t
            # Calcul des fonctions de base d'une spline de Hermite
            H0 = 1 - 3 * t2 + 2 * t3
            H1 = 3 * t2 - 2 * t3
            H2 = t - 2 * t2 + t3
            H3 = - t2+ t3
            result.append(H0 * points[0, :] + H1 * points[1, :] + H2 * m0 + H3 * m1)
        return np.array(result)
    
    if n == 2: 
        # Cas avec seulement trois points
        # Gestion de la première courbe
        y_1 = 2*points[0, :] - points[1, :] # Calcul du point y_{-1} 
        y_3 = 2*points[2, :] - points[1, :] # Calcul du point y_{3} 
        m0 = tangente_calc(points[1, :], y_1, c, points[0, :], points[2, :], tension, bias, continuity, 0)
        m1 = tangente_calc(points[2, :], points[0, :], c, points[1, :], y_3, tension, bias, continuity, 1)
        for t in T:
            t2 = t * t
            t3 = t2 * t
            # Calcul des fonctions de base d'une spline de Hermite
            H0 = 1 - 3 * t2 + 2 * t3
            H1 = 3 * t2 - 2 * t3
            H2 = t - 2 * t2 + t3
            H3 = - t2+ t3
            result.append(H0 * points[0, :] + H1 * points[1, :] + H2 * m0 + H3 * m1)
        # Gestion de la deuxième courbe
        y_4 = 2*y_3 - points[2, :] # Calcul du point y_{4} 
        m0 = tangente_calc(points[2, :], points[0, :], c, points[1, :], y_3, tension, bias, continuity, 0)
        m1 = tangente_calc(y_3, points[1, :], c, points[2, :], y_4, tension, bias, continuity, 1)
        for t in T:
            t2 = t * t
            t3 = t2 * t
            # Calcul des fonctions de base d'une spline de Hermite
            H0 = 1 - 3 * t2 + 2 * t3
            H1 = 3 * t2 - 2 * t3
            H2 = t - 2 * t2 + t3
            H3 = - t2+ t3
            result.append(H0 * points[1, :] + H1 * points[2, :] + H2 * m0 + H3 * m1)
        
        return np.array(result)
    
    # cas avec plus de 3 points
    # Gestion de la première courbe
    y_1 = 2*points[0, :] - points[1, :] # Calcul du point y_{-1} 
    m0 = tangente_calc(points[1, :], y_1, c, points[0, :], points[2, :], tension, bias, continuity, 0)
    m1 = tangente_calc(points[2, :], points[0, :], c, points[1, :], points[3, :], tension, bias, continuity, 1)
    for t in T:
        t2 = t * t
        t3 = t2 * t
        # Calcul des fonctions de base d'une spline de Hermite
        H0 = 1 - 3 * t2 + 2 * t3
        H1 = 3 * t2 - 2 * t3
        H2 = t - 2 * t2 + t3
        H3 = - t2+ t3
        result.append(H0 * points[0, :] + H1 * points[1, :] + H2 * m0 + H3 * m1)
    
    # Gestion des points intermédiaires
    for i in range(1, n-2):
        m0 = tangente_calc(points[i+1, :], points[i-1, :], c, points[i, :], points[i+2, :], tension, bias, continuity, 0)
        m1 = tangente_calc(points[i+2, :], points[i, :], c, points[i+1, :], points[i+3, :], tension, bias, continuity, 1)
        for t in T:
            t2 = t * t
            t3 = t2 * t

            # Calcul des fonctions de base d'une spline de Hermite
            H0 = 1 - 3 * t2 + 2 * t3
            H1 = 3 * t2 - 2 * t3
            H2 = t - 2 * t2 + t3
            H3 = - t2+ t3

            result.append(H0 * points[i, :] + H1 * points[i + 1, :] + H2 * m0 + H3 * m1)
    
    # Gestion de l'avant dernière courbe
    y_n1 = 2*points[n, :] - points[n-1, :] # Calcul du point y_{n+1} 
    m0 = tangente_calc(points[n-1, :], points[n-3, :], c, points[n-2, :], points[n, :], tension, bias, continuity, 0)
    m1 = tangente_calc(points[n, :], points[n-2, :], c, points[n-1, :], y_n1, tension, bias, continuity, 1)
    for t in T:
        t2 = t * t
        t3 = t2 * t
        # Calcul des fonctions de base d'une spline de Hermite
        H0 = 1 - 3 * t2 + 2 * t3
        H1 = 3 * t2 - 2 * t3
        H2 = t - 2 * t2 + t3
        H3 = - t2+ t3
        result.append(H0 * points[n-2, :] + H1 * points[n-1, :] + H2 * m0 + H3 * m1)



    # Gestion de la dernière courbe
    y_n1 = 2*points[n, :] - points[n-1, :] # Calcul du point y_{n+1} 
    y_n2 = 2*y_n1 - points[n, :] # Calcul du point y_{n+2} 
    m0 = tangente_calc(points[n, :], points[n-2, :], c, points[n-1, :], y_n1, tension, bias, continuity, 0)
    m1 = tangente_calc(y_n1, points[n-1, :], c, points[n, :], y_n2, tension, bias, continuity, 1)
    for t in T:
        t2 = t * t
        t3 = t2 * t
        # Calcul des fonctions de base d'une spline de Hermite
        H0 = 1 - 3 * t2 + 2 * t3
        H1 = 3 * t2 - 2 * t3
        H2 = t - 2 * t2 + t3
        H3 = - t2+ t3
        result.append(H0 * points[n-1, :] + H1 * points[n, :] + H2 * m0 + H3 * m1)
    return result

def cubic_spline(x, y, tol = 1e-100):
    """
    Interpolate using natural cubic splines.
    
    Generates a strictly diagonal dominant matrix then applies Jacobi's method.
    
    Returns coefficients:
    b, coefficient of x of degree 1
    c, coefficient of x of degree 2
    d, coefficient of x of degree 3
    """ 
    x = np.array(x)
    size = len(y)
    delta_x = np.diff(x)
    delta_y = np.diff(y)
    ### Get matrix A
    A = np.zeros(shape = (size,size))
    b = []
    A[0,0] = 2
    A[-1,-1] = 2
    if (size > 1):
        A[0, 1] = 1
        A[-1, -2] = 1
    for i in range(1,size-1):
        A[i, i-1] = 1
        A[i, i+1] = 1
        A[i,i] = 4

    ### Get matrix b    
    b.append([3*(y[1][0] - y[0][0]), 3*(y[1][1] - y[0][1])])
    for j in range(1, size-1):
        b.append([3*(y[j+1][0] - y[j-1][0]), 3*(y[j+1][1] - y[j-1][1])])
    b.append([3*(y[-1][0] - y[-2][0]), 3*(y[-1][1] - y[-2][1])])
    ### Solves for c in Ac = b
    # Si vous souhaitez exécuter la résolution avec la méthode de jacobi, décommentez
    # c = jacobi(A, b, np.zeros((len(A), 2)), tol = tol, n_iterations=1000)
    c = solve_matrix(A, b) 

    return c

def jacobi(A, b, x0, tol, n_iterations=300):
    """
    Performs Jacobi iterations to solve the line system of
    equations, Ax=b, starting from an initial guess, ``x0``.
    
    Returns:
    x, the estimated solution
    """
    
    n = A.shape[0]
    b = np.array(b)
    x = x0.copy()
    x_prev = x0.copy()
    counter = 0
    x_diff = tol+1
    
    while (x_diff > tol) and (counter < n_iterations): #iteration level 
        
        for i in range(0, n): #element wise level for x
            s_x = 0
            s_y = 0
            for j in range(0,n): #summation for i !=j
                if i != j:
                    s_x += A[i,j] * x_prev[j][0] 
                    s_y += A[i,j] * x_prev[j][1] 

            x[i][0] = (b[i][0] - s_x) / A[i,i]
            x[i][1] = (b[i][1] - s_y) / A[i,i]
        #update values
        counter += 1
        x_diff = (np.sum((x-x_prev)**2))**0.5 
        x_prev = x.copy() #use new x for next iteration
        
        
    
    print("Number of Iterations: ", counter)
    print("Norm of Difference: ", x_diff)
    return x

def solve_matrix(a, b):
    return np.linalg.solve(a, b)

def Natural_Splines(points, T, c, tension, bias, continuity, tangente_calc): 
    """
    Perform natural splines computation. 
    C2 condition is respected. 
    """
    n = points.shape[0] - 1
    x = list(range(n+1))
    x = np.array(x)
    y = np.array(points)
    result = []
    d = cubic_spline(x, y, 10**(-2))
    for i in range(n):
        for t in T:
            coeff = Bernstein3_poly(t)
            result.append(points[i]*coeff[0] + (points[i] + 1/3*d[i])*coeff[1] + (points[i+1]-1/3*d[i+1])*coeff[2] + points[i+1]*coeff[3])
    return result

def Bernstein3_poly(t):
    b_0_3 = math.comb(3, 0) * (1-t)**3 * t**0 
    b_1_3 = math.comb(3, 1) * (1-t)**2 * t**1 
    b_2_3 = math.comb(3, 2) * (1-t)**1 * t**2 
    b_3_3 = math.comb(3, 3) * (1-t)**0 * t**3 
    return b_0_3, b_1_3, b_2_3, b_3_3




# ---------- Here add derivation functions ----------

def Cardinal_Splines(psup, pinf, c, p, psupsup, tension, bias, continuity, mode):
    return (1-c)*(psup - pinf)/2

def Kochanek_Bartels(psup, pinf, c, p, psupsup, tension, bias, continuity, mode):
    if (mode == 0):
        return (1-tension)*(1+bias)*(1+continuity)/2*(p-pinf) + (1-tension)*(1-bias)*(1-continuity)/2*(psup-p)
    if (mode == 1):
        return (1-tension)*(1+bias)*(1-continuity)/2*(psup-p) + (1-tension)*(1-bias)*(1+continuity)/2*(psupsup-psup)
    
    

if __name__ == "__main__":
    window = CurveEditorWindow([{"name": "Bezier", "algo": DeCasteljau}, {"name": "Hermite", "algo": HermiteSpline}, {"name": "Lagrange", "algo": Aitken_Neville}, {"name": "Natural Splines", "algo": Natural_Splines}], [{"name": "Cardinal Splines", "method": Cardinal_Splines}, {"name": "Kochanek-Bartels", "method": Kochanek_Bartels}])
    window.mainloop()
