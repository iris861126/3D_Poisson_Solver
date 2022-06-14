# 3D_Poisson_Solver
使用successive over-relaxation (SOR)迭代求解三維Poisson Equation。

執行以下以產生.so檔便於python呼叫函式:
```
f2py -m poisson_ter -c poisson_ter.f95 --opt=-O3
```
