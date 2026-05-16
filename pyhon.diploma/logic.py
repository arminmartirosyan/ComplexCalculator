from matplotlib import contour
import sympy as sp
import cmath
import re
import numpy as np

def complex_nums(input_str: str):         # բացատների հեռացում, i-ն j-ով փոխարինում
    s = input_str.replace(" ", "").replace("i", "j")
    try:
        if s == 'j': return 0, 1
        if s == '-j': return 0, -1
        z_val = complex(s)
        return z_val.real, z_val.imag
    except:
        matches = re.findall(r'[+-]?\d*\.?\d*j?', s)
        real, imag = 0, 0
        for m in matches:
            if not m or m in ('+', '-'): continue
            if 'j' in m:
                if m in ('j', '+j'): imag += 1
                elif m == '-j': imag -= 1
                else: imag += float(m.replace('j', ''))
            else: real += float(m)
        return real, imag

def complex_nums_form(z):   #իրական կեղծ մասեր
    real = round(z.real, 3) #թվերի կլորացում մինչև 3 նիշ
    imag = round(z.imag, 3)

    if imag == 0:
        return f"{real:g}"
    if real == 0:
        if imag == 1: return "i"
        if imag == -1: return "-i"
        return f"{imag:g}i"
    if imag > 0:
        imag_part = "i" if imag == 1 else f"{imag:g}i"
        return f"{real:g} + {imag_part}"
    else:
        imag_part = "i" if imag == -1 else f"{abs(imag):g}i"
        return f"{real:g} - {imag_part}"

def contour(expr: str, z_sym):
    expr = expr.replace(" ", "")  #հանում ենք բացատները
    match = re.match(r"\|(.+)\|=(.+)", expr)  #regex-ով առանձնացնում ենք կոնտուրի ներսի արտահայտությունն ու շառավիղը
    if not match:
        raise ValueError("Սխալ կոնտուր: Օրինակ՝ |z|=2")
    inside_str = match.group(1).replace("^", "**").replace("i", "I")  #^, i-ի փախակերպում **, I sympy-ի համար
    radius = float(sp.sympify(match.group(2)))
    g = sp.sympify(inside_str)
    sol = sp.solve(g, z_sym)
    return sol[0], radius  

def calculate_nums(z1_str, z2_str, operation):  #կոմպլեքս թվերի գումարում և բազմապատկում
    r1, i1 = complex_nums(z1_str) 
    r2, i2 = complex_nums(z2_str)
    z1, z2 = complex(r1, i1), complex(r2, i2)
    res = z1 + z2 if operation == "add" else z1 * z2
    return {
        'result': res,
        'formatted': complex_nums_form(res),
        'modulus': abs(res), #abs() ֆունկցիան հաշվում է ավտոմատ կոմպլեքս թվի մոդուլը
        'argument': cmath.phase(res) #արգումետը հաշվելու համար
    }
    
def complex_roots(z_str, n):  #կոմպլեքս թվի n aստիճանի արմատ
    r, i = complex_nums(z_str) 
    z = complex(r, i)
    m, a = abs(z), cmath.phase(z)
    roots = [cmath.rect(m**(1/n), (a + 2*cmath.pi*k)/n) for k in range(n)]
    results = []
    for root in roots:
        root_str = complex_nums_form(root)
        root_str = root_str.replace('sqrt', '√')
        results.append(root_str)
    return results

def koshi_integral_with_steps(f_str, contour_str, z_sym): #Կոշի ինտեգրալ
    steps = []
    
    steps.append("ՔԱՅԼ 1: Ֆունկցիայի վերլուծություն")
    f_expr = sp.sympify(f_str.replace("i", "I"))
    steps.append(f"   Տրված ֆունկցիան՝ f(z) = {f_expr}")
    
    steps.append("\nՔԱՅԼ 2: Կոնտուրի վերլուծություն")
    z0, R = contour(contour_str, z_sym)
    steps.append(f"   Կոնտուր՝ |z - {z0}| = {R}")
    steps.append(f"   Կենտրոնը՝ z0 = {z0}")
    steps.append(f"   Շառավիղը՝ R = {R}")
    
    steps.append("\nՔԱՅԼ 3: Բևեռների որոշում")
    denominator = sp.denom(f_expr)
    steps.append(f"   Հայտարարը՝ {denominator}")
    
    poles = sp.solve(denominator, z_sym)
    steps.append(f"   Բևեռները՝ {poles}")
    
    steps.append("\nՔԱՅԼ 4: Կոնտուրի ներսում գտնվող բևեռները")
    inside = []
    for p in poles:
        p_complex = complex(sp.N(p))
        z0_complex = complex(sp.N(z0))
        distance = abs(p_complex - z0_complex)
        if distance < R:
            steps.append(f"   {p} < {R} => Գտնվում է կոնտուրի ներսում")
            inside.append(p)
        else:
            steps.append(f"   {distance} >= {R} => Գտնվում է կոնտուրից դուրս")
    
    if not inside:
        steps.append("\nԿոնտուրի ներսում բևեռներ չկան, ինտեգրալը = 0")
        return 0, steps
    
    steps.append(f"\n   Կոնտուրի ներսում գտնվող բևեռներ՝ {inside}")
    
    steps.append("\nՔԱՅԼ 5: Մնացքների հաշվում")
    residues = []
    for p in inside:
        residue = sp.residue(f_expr, z_sym, p)
        steps.append(f" Res(f, {p}) = {residue}")
        residues.append(residue)
    
    steps.append("\nՔԱՅԼ 6: Մնացքների գումարում")
    residue_sum = sum(residues)
    sum_expression = " + ".join([str(r) for r in residues])
    steps.append(f"   Σ Res(f, zₖ) = {sum_expression} = {residue_sum}")
    
    steps.append("\nՔԱՅԼ 7: Կոշիի մնացքների թեորեմի կիրառում")
    steps.append("   ∮ f(z) dz = 2πi · Σ Res(f, zₖ)")
    steps.append(f"   = 2πi · ({residue_sum})")
    
    result = 2 * sp.pi * sp.I * residue_sum
    simplified = sp.simplify(result)
    steps.append(f"   = {simplified}")
    
    return simplified, steps


def mnatsqneri_integral_with_steps(f_str, R_str, z_sym): #մնացքների մեթոդով ինտեգրալ
    steps = []
    
    steps.append("ՔԱՅԼ 1: Ֆունկցիայի վերլուծություն")
    f_expr = sp.sympify(f_str.replace("i", "I").replace("^", "**"))
    steps.append(f"   Տրված ֆունկցիան՝ f(z) = {f_expr}")
    
    steps.append("\nՔԱՅԼ 2: Կոնտուրի շառավիղ")
    R = float(R_str)
    steps.append(f"   |z| = {R}")
    
    steps.append("\nՔԱՅԼ 3: Բևեռների որոշում")
    denominator = sp.denom(f_expr)
    steps.append(f"   Հայտարարը՝ {denominator}")
    
    poles = sp.solve(denominator, z_sym)
    steps.append(f"   Բևեռները՝ {poles}")
    
    steps.append(f"\nՔԱՅԼ 4: |z| < {R} պայմանին բավարարող բևեռների ընտրություն")
    inside = []
    for p in poles:
        p_complex = complex(sp.N(p))
        distance = abs(p_complex)
        steps.append(f"   |{p}| = {distance:.4f}")
        if distance < R:
            steps.append(f"   {distance:.4f} < {R} => Գտնվում է շրջանի ներսում")
            inside.append(p)
        else:
            steps.append(f"   {distance:.4f} >= {R} => Գտնվում է շրջանից դուրս")
    
    if not inside:
        steps.append("\nՇրջանի ներսում բևեռներ չկան, ինտեգրալը = 0")
        return 0, steps
    
    steps.append(f"\n   Շրջանի ներսում գտնվող բևեռներ՝ {inside}")
    
    steps.append("\nՔԱՅԼ 5: Մնացքների հաշվում")
    residues = []
    for p in inside:
        residue = sp.residue(f_expr, z_sym, p)
        steps.append(f"\n  Res(f, {p}) = {residue}")
        residues.append(residue)
    
    steps.append("\nՔԱՅԼ 6: Մնացքների գումարում")
    residue_sum = sum(residues)
    sum_expression = " + ".join([str(r) for r in residues])
    steps.append(f"   Σ Res(f, zₖ) = {sum_expression} = {residue_sum}")
    
    steps.append("\nՔԱՅԼ 7: Մնացքների թեորեմի կիրառում")
    steps.append("   ∮ f(z) dz = 2πi · Σ Res(f, zₖ)")
    steps.append(f"   = 2πi · ({residue_sum})")
    
    result = 2 * sp.pi * sp.I * residue_sum
    simplified = sp.simplify(result)
    steps.append(f"   = {simplified}")
    
    return simplified, steps

def koshi_integral(f_str, contour_str, z_sym):
    result, _ = koshi_integral_with_steps(f_str, contour_str, z_sym)
    return result

def mnatsqneri_integral(f_str, R_str, z_sym):
    result, _ = mnatsqneri_integral_with_steps(f_str, R_str, z_sym)
    return result

def complex_cucchayin(z_str): #ցուցչայինի ֆունկցիա
    r, i = complex_nums(z_str)
    z = complex(r, i)
    modulus = abs(z)
    argument = cmath.phase(z)
    return {
        'modulus': modulus,
        'argument': argument,
        'polar_form': f"{modulus:.2f} * e^(i*{argument:.2f})"
    }

def complex_graphic(z_str): #գրաֆիկի պատկերման համար
    r, i = complex_nums(z_str)
    return complex(r, i)