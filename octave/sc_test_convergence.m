for ii = 1:(length(nR)/2-1)

display([ii,length(nR)/2-1])  
ndx_filter = [ii:(length(nR)-ii),(length(nR)+ii):(2*length(nR)-ii)];

Tri_filter = Tri(:,ndx_filter);
Tro_filter = Tro(:,ndx_filter);
TOhx_filter = TOhx(ndx_filter,:);
TIhx_filter = TIhx(ndx_filter,:);

Tri_TOhx_filter = Tri_filter*TOhx_filter;
Tro_TIhx_filter = Tro_filter*TIhx_filter;
Tri_TIhx_filter = Tri_filter*TIhx_filter;
Tro_TOhx_filter = Tro_filter*TOhx_filter;

TFI_filter1 = TFI(ndx_filter,ndx_filter);
TFO_filter1 = TFO(ndx_filter,ndx_filter);

TFIq1_filter2 = toeplitz(TFI_filter1(1,1:end/2),TFI_filter1(1,1:end/2));
TFIq2_filter2 = toeplitz(TFI_filter1(1,(end/2+1):end),TFI_filter1(1,(end/2+1):end));
TFIq3_filter2 = toeplitz(TFI_filter1((end/2+1),1:end/2),TFI_filter1((end/2+1),1:end/2));
TFIq4_filter2 = toeplitz(TFI_filter1((end/2+1),(end/2+1):end),TFI_filter1((end/2+1),(end/2+1):end));

TFOq1_filter2 = toeplitz(TFO_filter1(1,1:end/2),TFO_filter1(1,1:end/2));
TFOq2_filter2 = toeplitz(TFO_filter1(1,(end/2+1):end),TFO_filter1(1,(end/2+1):end));
TFOq3_filter2 = toeplitz(TFO_filter1((end/2+1),1:end/2),TFO_filter1((end/2+1),1:end/2));
TFOq4_filter2 = toeplitz(TFO_filter1((end/2+1),(end/2+1):end),TFO_filter1((end/2+1),(end/2+1):end));

TFI_filter2 = cat(2,cat(1,TFIq1_filter2,TFIq3_filter2),cat(1,TFIq2_filter2,TFIq4_filter2));
TFO_filter2 = cat(2,cat(1,TFOq1_filter2,TFOq3_filter2),cat(1,TFOq2_filter2,TFOq4_filter2));

Tri_TFI_filter1 = Tri_filter*TFI_filter1;
Tro_TFI_filter1 = Tro_filter*TFI_filter1;
Tri_TFO_filter1 = Tri_filter*TFO_filter1;
Tro_TFO_filter1 = Tro_filter*TFO_filter1; 

Tri_TFI_filter2 = Tri_filter*TFI_filter2;
Tro_TFI_filter2 = Tro_filter*TFI_filter2;
Tri_TFO_filter2 = Tri_filter*TFO_filter2;
Tro_TFO_filter2 = Tro_filter*TFO_filter2;

error_filter{1,1}(ii) = sum(sum(abs(Mcmb{1,2}-(Tri_TOhx_filter - Tri_TFO_filter1*TOhx_filter))/max(abs(Mcmb{1,2}(:)))));
error_filter{2,1}(ii) = sum(sum(abs(Mcmb{2,2}-(Tro_TIhx_filter - Tro_TFI_filter1*TIhx_filter))/max(abs(Mcmb{2,2}(:)))));
error_filter{3,1}(ii) = sum(sum(abs(Mcmb{3,2}-(Tri_TIhx_filter - Tri_TFO_filter1*TIhx_filter))/max(abs(Mcmb{3,2}(:)))));
error_filter{4,1}(ii) = sum(sum(abs(Mcmb{4,2}-(Tro_TOhx_filter - Tro_TFI_filter1*TOhx_filter))/max(abs(Mcmb{4,2}(:)))));

error_filter{1,2}(ii) = sum(sum(abs(Mcmb{1,2}-(Tri_TOhx_filter - Tri_TFO_filter2*TOhx_filter))/max(abs(Mcmb{1,2}(:)))));
error_filter{2,2}(ii) = sum(sum(abs(Mcmb{2,2}-(Tro_TIhx_filter - Tro_TFI_filter2*TIhx_filter))/max(abs(Mcmb{2,2}(:)))));
error_filter{3,2}(ii) = sum(sum(abs(Mcmb{3,2}-(Tri_TIhx_filter - Tri_TFO_filter2*TIhx_filter))/max(abs(Mcmb{3,2}(:)))));
error_filter{4,2}(ii) = sum(sum(abs(Mcmb{4,2}-(Tro_TOhx_filter - Tro_TFI_filter2*TOhx_filter))/max(abs(Mcmb{4,2}(:)))));


end