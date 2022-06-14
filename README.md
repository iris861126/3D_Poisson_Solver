# 3D_Poisson_Solver
求解三維Poisson Equation

執行以下以產生.so檔便於Python呼叫函式:
```
f2py -m poisson_ter -c poisson_ter.f95 --opt=-O3
```


於Python程式碼中以
```
from poisson_ter import poisson3d_solver_sor
```
導入Fortran函式。
