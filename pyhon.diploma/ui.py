import tkinter as tk
from tkinter import messagebox
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import cmath

from logic import (
    calculate_nums, 
    complex_roots, 
    koshi_integral,
    koshi_integral_with_steps,
    mnatsqneri_integral,
    mnatsqneri_integral_with_steps, 
    complex_cucchayin, 
    complex_graphic
)

class ComplexCalcApp:
    def __init__(self, root): #ստեղծվում է կոնստրուկտր, որն ավտոմատ կանչվում է ծրագիրն աշխատեցնելիս
        self.root = root
        self.root.title("Complex Analysis Calculator") #վերևի փախ անկյունում գրվում է ծրագրի անունը
        self.root.geometry("1000x950") 
        self.root.configure(bg="#1e1e2f")
        self.z_sym = sp.Symbol('z')
        self.fig, self.ax = plt.subplots(figsize=(5, 5), dpi=100)
        self.fig.patch.set_facecolor('#1e1e2f')
        self.canvas = None #գրաֆիկի համար
        
        self.main_frame = tk.Frame(self.root, bg="#1e1e2f")
        self.main_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        self.main_menu()

    def clear_screen(self): #մաքրում է էկրանը, ցիկլով անցնելով main.frame-ի միջի բոլոր էլեմենտերի վրայով
        for widget in self.main_frame.winfo_children():
            widget.destroy()
        self.canvas = None

    def back_button(self, callback): #ավելացնում ենք back button 
        btn = tk.Button(self.main_frame, text="← ՀԵՏ", font=("Segoe UI", 12, "bold"), 
                        bg="#5a5a8a", fg="white", bd=0, padx=20, pady=10, 
                        command=callback)
        btn.pack(anchor="nw", pady=10)

    def main_menu(self): #գլխավոր մենյու
        self.clear_screen()
        tk.Label(self.main_frame, text="ԿՈՄՊԼԵՔՍ ՀԱՇՎԻՉ", font=("Segoe UI", 32, "bold"), 
                 fg="#00ffcc", bg="#1e1e2f").pack(pady=(40, 20))

        theory_text = (
            "Կոմպլեքս անալիզը մաթեմատիկական անալիզի այն բաժինն է, "
            "որն ուսումնասիրում է կոմպլեքս թվերը և կոմպլեքս փոփոխականի ֆունկցիաները։\n"
            "Այս հավելվածը թույլ է տալիս կատարել գործողություններ կոմպլեքս թվերի հետ, հաշվել հաշվել ինտեգրալներ։"
        )
        tk.Label(self.main_frame, text=theory_text, font=("Segoe UI", 11), 
                 fg="#cfcfda", bg="#1e1e2f", justify="center", wraplength=650).pack(pady=20, padx=40)

        btn_style = {"font": ("Segoe UI", 16, "bold"), "width": 25, "height": 2, "bg": "#3a3a5a", "fg": "white", "bd": 0}
        tk.Button(self.main_frame, text="ԿՈՄՊԼԵՔՍ ԹՎԵՐ", command=self.complex_numbers_menu, **btn_style).pack(pady=15)
        tk.Button(self.main_frame, text="ԿՈՄՊԼԵՔՍ ՖՈՒՆԿՑԻԱՆԵՐ", command=self.functions_menu, **btn_style).pack(pady=15)

    def complex_numbers_menu(self): #կոմպլեքս թվերի մենյու, +,*,արմատ, ցուցչային, գրաֆիկ
        self.clear_screen()
        self.back_button(self.main_menu)
        tk.Label(self.main_frame, text="Կոմպլեքս թվեր", font=("Segoe UI", 24, "bold"), fg="white", bg="#1e1e2f").pack(pady=20)
        
        btn_style = {"font": ("Segoe UI", 14), "width": 35, "height": 2, "bg": "#404070", "fg": "white", "bd": 0}
        tk.Button(self.main_frame, text="Գործողություններ", command=self.show_operations, **btn_style).pack(pady=10)
        tk.Button(self.main_frame, text="Արմատների հաշվում", command=self.show_demoivre, **btn_style).pack(pady=10)
        tk.Button(self.main_frame, text="Ցուցչային ներկայացում", command=self.show_conversion, **btn_style).pack(pady=10)
        tk.Button(self.main_frame, text="Գրաֆիկական պատկերում", command=self.show_graphic_ui, **btn_style).pack(pady=20)

    def functions_menu(self): #ֆունկցիաների մենյու, կոշի, մնացնքերի մեթոդով ինտգեգրալ
        self.clear_screen()
        self.back_button(self.main_menu)
        tk.Label(self.main_frame, text="Կոմպլեքս ֆունկցիաներ", font=("Segoe UI", 24, "bold"), fg="white", bg="#1e1e2f").pack(pady=20)
        
        btn_style = {"font": ("Segoe UI", 14), "width": 35, "height": 2, "bg": "#404070", "fg": "white", "bd": 0}
        tk.Button(self.main_frame, text="Ինտեգրալ (Կոշի մեթոդ)", command=self.show_cauchy, **btn_style).pack(pady=10)
        tk.Button(self.main_frame, text="Ինտեգրալ (Մնացքների մեթոդ)", command=self.show_residue, **btn_style).pack(pady=10)

    def prepare_ax(self): #ստեղծում ենք Matplotlib-ի գրաֆիկական առանցքները
        self.ax.cla()
        self.ax.set_facecolor('#2b2b4a')
        self.ax.axhline(0, color="white", lw=1, alpha=0.5)
        self.ax.axvline(0, color="white", lw=1, alpha=0.5)
        self.ax.grid(True, alpha=0.2, linestyle='--')
        self.ax.tick_params(colors='white')

    def draw_contour(self, contour_str): #շրջանագծի պատկերում
        self.prepare_ax()
        try:
            contour_str = contour_str.replace(" ", "").lower()
            if "|z|" in contour_str:
                radius = float(contour_str.split("=")[1])
                center = 0j
            elif "|z-" in contour_str or "|z+" in contour_str:
                parts = contour_str.split("=")
                radius = float(parts[1])
                inner = parts[0].replace("|", "").replace("z", "")
                center = -complex(inner.replace("i", "j"))
            else:
                radius = 2
                center = 0j

            theta = np.linspace(0, 2*np.pi, 100)
            x = center.real + radius * np.cos(theta)
            y = center.imag + radius * np.sin(theta)
            
            self.ax.plot(x, y, color="#00ffcc", lw=3, label="Contour")
            self.ax.fill(x, y, color="#00ffcc", alpha=0.1)
            
            limit = radius + abs(center) + 1
            self.ax.set_xlim(-limit, limit)
            self.ax.set_ylim(-limit, limit)
            self.ax.set_aspect('equal')
            self.canvas.draw()
        except Exception as e:
            print(f"Contour drawing error: {e}")

    def update_plot(self): #կոմպլեքս թիվը որպես վկտոր ներկայացնելու համար
        try:
            z = complex_graphic(self.plot_ent.get())
            x, y = float(z.real), float(z.imag)
            
            self.prepare_ax()
            self.ax.quiver(0, 0, x, y, angles='xy', scale_units='xy', scale=1, color='#00ffcc', width=0.015)
            self.ax.text(x, y, f' ({x:g}, {y:g})', color='white', fontsize=10, fontweight='bold')

            limit = max(abs(x), abs(y), 2) + 1
            self.ax.set_xlim(-limit, limit)
            self.ax.set_ylim(-limit, limit)
            self.ax.set_aspect('equal')
            self.canvas.draw()
        except Exception as e:
            messagebox.showerror("Սխալ", f"Անվավեր արտահայտություն: {e}")

    def update_operation_plot(self, z_res): #2 կոմպլեքս թվեր որոնք գումարվել են, դա որպես վեկտոր ներկայացնելու համար
        try:
            res_val = complex(z_res)
            x, y = res_val.real, res_val.imag
            self.prepare_ax()
            
            self.ax.quiver(0, 0, x, y, angles='xy', scale_units='xy', scale=1, color='#ff00ff', width=0.015)
            self.ax.text(x, y, f' ({x:g}, {y:g})', color='#00ffcc', fontsize=12, fontweight='bold')
            
            limit = max(abs(x), abs(y), 2) + 1
            self.ax.set_xlim(-limit, limit)
            self.ax.set_ylim(-limit, limit)
            self.ax.set_aspect('equal')
            self.canvas.draw()
        except Exception as e:
            print(f"Plot error: {e}")

    def show_operations(self): #գումարում, բազմապատկում հրամանների ֆունկցիայի տեսք, կառուցում matplotlib-ով գրաֆիկը
        self.clear_screen()
        self.back_button(self.complex_numbers_menu)
        f = tk.Frame(self.main_frame, bg="#1e1e2f")
        f.pack(pady=10)
        self.ent1 = self.create_input(f, "z1 =", "1+i")
        self.ent2 = self.create_input(f, "z2 =", "2-3i")
        btn_frame = tk.Frame(self.main_frame, bg="#1e1e2f")
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Գումարել", command=lambda: self.calc_op("add"), bg="#5a5a8a", fg="white", font=("Segoe UI", 12, "bold"), width=15).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Բազմապատկել", command=lambda: self.calc_op("mul"), bg="#5a5a8a", fg="white", font=("Segoe UI", 12, "bold"), width=15).pack(side="left", padx=5)
        self.res_txt = tk.Text(self.main_frame, height=3, width=50, bg="#2b2b4a", fg="#00ffcc", font=("Consolas", 14, "bold"), padx=10, pady=10)
        self.res_txt.pack(pady=10)
        self.canvas_frame = tk.Frame(self.main_frame, bg="#1e1e2f")
        self.canvas_frame.pack(pady=5, fill="both", expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.prepare_ax()

    def calc_op(self, op): #վերցնում է տեքստը փոխանցում logic-ի calculate_nums, արդյունքը տպում էկրանին և գծում գրաֆիկ, հակառակ դեպքում error
        try:
            res = calculate_nums(self.ent1.get(), self.ent2.get(), op)
            self.res_txt.delete("1.0", tk.END)
            self.res_txt.insert(tk.END, f"Արդյունք: {res['formatted']}\nՄոդուլ: {res['modulus']:.4f}")
            self.update_operation_plot(res['result'])
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def show_graphic_ui(self):  #գրաֆիկի պատկերում
        self.clear_screen()
        self.back_button(self.complex_numbers_menu)
        f = tk.Frame(self.main_frame, bg="#1e1e2f")
        f.pack(pady=10)
        self.plot_ent = tk.Entry(f, width=15, font=("Arial", 20), justify="center")
        self.plot_ent.insert(0, "1+2i")
        self.plot_ent.pack(side="left", padx=10)
        tk.Button(f, text="ՊԱՏԿԵՐԵԼ", font=("Segoe UI", 14, "bold"), command=self.update_plot, bg="#00ffcc", fg="black").pack(side="left")
        self.canvas_frame = tk.Frame(self.main_frame, bg="#1e1e2f")
        self.canvas_frame.pack(pady=20, fill="both", expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.update_plot()

    def show_cauchy(self):  # Կոշի ինտեգրալի էկրան  
        self.clear_screen()
        self.back_button(self.functions_menu)
        f = tk.Frame(self.main_frame, bg="#1e1e2f")
        f.pack(pady=10)
        self.ent_f = self.create_input(f, "f(z) =", "1/(z-1)")
        self.ent_c = self.create_input(f, "Կոնտուր =", "|z|=2")
        tk.Button(self.main_frame, text="ՀԱՇՎԵԼ", bg="#5a5a8a", fg="white", font=("Segoe UI", 14, "bold"), command=self.calc_cauchy).pack(pady=10)
        self.res_txt = tk.Text(self.main_frame, height=20, width=100, bg="#2b2b4a", fg="#00ffcc", font=("Consolas", 11), padx=10, pady=10)
        self.res_txt.pack(pady=5)
        self.canvas_frame = tk.Frame(self.main_frame, bg="#1e1e2f")
        self.canvas_frame.pack(pady=5, fill="both", expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.prepare_ax()

    def calc_cauchy(self): # քայլերը հերթով կոշի լուծման էկրանին
        try:
            c = self.ent_c.get()
            result, steps = koshi_integral_with_steps(self.ent_f.get(), c, self.z_sym)
            
            self.res_txt.delete("1.0", tk.END)
            
            for step in steps:
                self.res_txt.insert(tk.END, step + "\n")
            self.res_txt.insert(tk.END, f"∮ f(z) dz = {result}\n")
            
            self.draw_contour(c)
        except Exception as e:
            messagebox.showerror("Սխալ", str(e))

    def show_residue(self): # մնացքների մեթոդի էկրան
        self.clear_screen()
        self.back_button(self.functions_menu)
        f = tk.Frame(self.main_frame, bg="#1e1e2f")
        f.pack(pady=10)
        self.ent_f_res = self.create_input(f, "f(z) =", "1/(z**2+1)")
        self.ent_r_res = self.create_input(f, "R =", "2")
        tk.Button(self.main_frame, text="ՀԱՇՎԵԼ", bg="#5a5a8a", fg="white", font=("Segoe UI", 14, "bold"), command=self.calc_residue).pack(pady=10)
        self.res_txt = tk.Text(self.main_frame, height=20, width=100, bg="#2b2b4a", fg="#00ffcc", font=("Consolas", 11), padx=10, pady=10)
        self.res_txt.pack(pady=5)
        self.canvas_frame = tk.Frame(self.main_frame, bg="#1e1e2f")
        self.canvas_frame.pack(pady=5, fill="both", expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.prepare_ax()

    def calc_residue(self):  # մնացքներովի հերթականությամբ լուծումը էկրանին
        try:
            r = self.ent_r_res.get()
            result, steps = mnatsqneri_integral_with_steps(self.ent_f_res.get(), r, self.z_sym)
            
            self.res_txt.delete("1.0", tk.END)
            for step in steps:
                self.res_txt.insert(tk.END, step + "\n")
            self.res_txt.insert(tk.END, f"∮ f(z) dz = {result}\n")
            
            self.draw_contour(f"|z|={r}")
        except Exception as e:
            messagebox.showerror("Սխալ", str(e))

    def show_demoivre(self): # ցուցչայինի էկրան
        self.clear_screen()
        self.back_button(self.complex_numbers_menu)
        f = tk.Frame(self.main_frame, bg="#1e1e2f")
        f.pack(pady=20)
        self.ent_z = self.create_input(f, "z =", "16")
        self.ent_n = self.create_input(f, "n =", "4")
        tk.Button(self.main_frame, text="Հաշվել արմատները", command=self.calc_roots, bg="#00ffcc", fg="black", font=("Segoe UI", 16, "bold")).pack(pady=20)
        self.res_txt = tk.Text(self.main_frame, height=10, width=50, bg="#2b2b4a", fg="#00ffcc", font=("Consolas", 14, "bold"), padx=15)
        self.res_txt.pack(pady=10, fill="both", expand=True)

    def calc_roots(self): # արմատների էկրան
        try:
            n = int(self.ent_n.get())
            roots = complex_roots(self.ent_z.get(), n)
            self.res_txt.delete("1.0", tk.END)
            for k, root in enumerate(roots):
                self.res_txt.insert(tk.END, f"w{k} = {root}\n")
        except:
            messagebox.showerror("Սխալ", "Ստուգեք տվյալները")

    def show_conversion(self): # փախակերպման էկրան
        self.clear_screen()
        self.back_button(self.complex_numbers_menu)
        self.ent_conv = self.create_input(self.main_frame, "z =", "1+i")
        tk.Button(self.main_frame, text="ՎԵՐԼՈՒԾԵԼ", font=("Segoe UI", 16, "bold"), command=self.calc_bever, bg="#5a5a8a", fg="white").pack(pady=20)
        self.res_txt = tk.Text(self.main_frame, height=10, width=50, bg="#2b2b4a", fg="#00ffcc", font=("Consolas", 18, "bold"), padx=15)
        self.res_txt.pack(pady=10, fill="both", expand=True)

    def calc_bever(self): # բևեռի էկրան
        try:
            res = complex_cucchayin(self.ent_conv.get())
            self.res_txt.delete("1.0", tk.END)
            self.res_txt.insert(tk.END, f"Մոդուլ: {res['modulus']:.4f}\nԱրգումենտ: {res['argument']:.4f}\nԲևեռային: {res['polar_form']}")
        except:
            messagebox.showerror("Սխալ", "Ստուգեք տվյալները")

    def create_input(self, parent, label, default=""): #օժանդակ մեթոդ, մուտքագրման համար
        row = tk.Frame(parent, bg="#1e1e2f")
        row.pack(pady=10)
        tk.Label(row, text=label, fg="white", bg="#1e1e2f", font=("Segoe UI", 16, "bold"), width=12).pack(side="left")
        e = tk.Entry(row, width=15, font=("Arial", 18), justify="center")
        e.insert(0, default)
        e.pack(side="left", padx=10)
        return e