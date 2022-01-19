function rd = analyze(Mcmb);
     sc_integrate
     rr0 = diag(response.range.pid);
     rr1 = diag(response.range.prl);
     rr2 = diag(response.range.prr);
     rd.rr = [rr0,rr1,rr2];
     rd.dd = field.range.prr+field.range.prl;
endfunction
