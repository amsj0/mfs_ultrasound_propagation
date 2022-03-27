resp = [];
sc_prepare
for unit = Mcmb_l
  Mcmb = unit{1};
  sc_integrate
  rr = diag(response.range.pid/g.model_scale+(response.range.prr+response.range.prl));
  resp = [resp,rr];
%  pause
end