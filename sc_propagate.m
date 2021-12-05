Thi_TFI = Thi*TFI;
Tho_TFI = Tho*TFI;
Thi_TFO = Thi*TFO;
Tho_TFO = Tho*TFO;

Mcmb{1,1} = Thi*TOhx - Thi_TFO*TOhx;
Mcmb{2,1} = Tho*TIhx - Tho_TFI*TIhx;
Mcmb{3,1} = Thi*TIhx - Thi_TFO*TIhx;
Mcmb{4,1} = Tho*TOhx - Tho_TFI*TOhx;
Mcmb{5,1} = p0kI;
Mcmb{6,1} = p0kO;

Tri_TFI = Tri*TFI;
Tro_TFI = Tro*TFI;
Tri_TFO = Tri*TFO;
Tro_TFO = Tro*TFO;  

Mcmb{1,2} = Tri*TOhx - Tri_TFO*TOhx;
Mcmb{2,2} = Tro*TIhx - Tro_TFI*TIhx;
Mcmb{3,2} = Tri*TIhx - Tri_TFO*TIhx;
Mcmb{4,2} = Tro*TOhx - Tro_TFI*TOhx;
Mcmb{5,2} = p0mI;
Mcmb{6,2} = p0mO;