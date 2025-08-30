import sympy as sp
import time

# شروع اندازه‌گیری زمان
start_time = time.time()

# مقادیر نمونه برای N و K
N = 5  # مقدار کوچک برای تست
K = 10 # مقدار کوچک برای تست
L = 6
M = 60.21
b = 0.6

# تعریف پارامترهای a_{i, 2m-2} برای i از 1 تا 2N-1 و m
a_ij = sp.Matrix([[sp.symbols(f'a_{i}_{2*m-2}') for m in range(1, K + 1)] for i in range(1, 2 * N)])

# تنظیم مقادیر a_{i,0} برابر صفر
for i in range(1, 2 * N):
    a_ij[i - 1, 0] = 0

# تابعی برای محاسبه C_n^{2m-2} به صورت نمادین
def compute_C_n_2m_2(n, m, M, L, N, b, i):
    C_n_2m_2 = sp.Sum((2 / (M * L)) * a_ij[sp.symbols('i')-1, m-1] * sp.sin(n * sp.symbols('i') * b * sp.pi / L),
                      (sp.symbols('i'), 1, 2 * N - 1)).doit()
    return C_n_2m_2

# متغیر برای ذخیره معادلات
equations = []

# حلقه برای تغییر مقدار k از 1 تا maxk
for k in range(1, K+1):
    # حلقه برای تغییر مقدار l از 1 تا 2N-1
    for l in range(1, 2 * N):
        result_l = 0  # معادله برای هر مقدار l
        for n in range(1, 4):  # n به عنوان شمارنده
            # محاسبه omega به عنوان تابعی از n
            omega_n = 24.62 * n ** 2 * 3.14**2
            for m in range(1, k+1):  # m به عنوان شمارنده
                # محاسبه C_n^{2m-2} از معادله دوم
                C_n_2m_2 = compute_C_n_2m_2(n, m, M, L, N, b, l)
                # اضافه کردن به معادله اول
                term = (-1)**(k + m) * C_n_2m_2 * omega_n ** (2*(k - m)) * sp.factorial(2 * m - 2)
                term *= sp.sin((n * sp.pi * l) / (2 * N))
                result_l += term
        # ساخت معادله برای هر l به فرم معادله = 0
        equations.append(sp.Eq(result_l, 0))

# حذف معادلاتی که به True یا False تبدیل شده‌اند
valid_equations = [eq for eq in equations if isinstance(eq, sp.Equality)]

# ضرایب a_ij را تا دو رقم اعشار رند کنیم
rounded_equations = [sp.Eq(sp.N(eq.lhs, 2), sp.N(eq.rhs, 2)) for eq in valid_equations]

# حل دستگاه معادلات نسبت به a_ij غیر از a_{i,0}
remaining_vars = a_ij[:, 1:]  # حذف a_{i,0} ها برای حل دستگاه
solution = sp.solve(rounded_equations, remaining_vars)

# چاپ حل دستگاه
sp.pprint(solution)

# چاپ حل دستگاه خط به خط در فایل
with open('aij_solutions.txt', 'w') as file:
    file.write("Solutions for a_ij (with coefficients rounded to 2 decimals):\n")
    for sol in solution:
        file.write(f'{sol} = {solution[sol]}\n')

# پایان زمان اجرا
end_time = time.time()
execution_time = end_time - start_time

# ذخیره زمان اجرا در فایل
with open('aij_solutions.txt', 'a') as file:  # 'a' برای الحاق به فایل
    file.write("\nExecution Time: {:.4f} seconds\n".format(execution_time))

print("مقادیر a_ij، نتایج و زمان اجرا در فایل 'aij_solutions.txt' ذخیره شدند.")
