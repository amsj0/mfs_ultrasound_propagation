function dThi = fn_propagator_inc_ref(p1,v1,p2,v2)

v1 = inv(v1);
p1 = p1*v1;
v1 = inv(p2-p1*v2);
p1 = inv(v2-inv(p1)*p2);
dThi = cat(1,cat(2,p2*v1,p2*p1),cat(2,v2*v1,v2*p1));

end