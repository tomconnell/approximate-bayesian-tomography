function [loglikelihood] = abcfun(parameters,data,edge_effect,gtolerance,ttolerance,toprint)

g_data = data(1,1:20);
t_data = data(1,21:140);

vp_parameters = brocherizegrid(parameters);

[g_ll,g_e] = abcgravfun(parameters,g_data,edge_effect,gtolerance);

[t_ll,t_e] = abctomfun(vp_parameters,t_data,ttolerance);

g_ll = (g_ll*10^(-1.5));

if toprint
    
    dlmwrite('abc_misfit.csv',t_e+g_e,'-append')
    
    dlmwrite('abc_gravity_misfit.csv',g_e,'-append')
    
    dlmwrite('abc_tom_misfit.csv',t_e,'-append')
    
    dlmwrite('abc_gravity_likelihood.csv',g_ll,'-append')
    
    dlmwrite('abc_tom_likelihood.csv',t_ll,'-append')
    
end

loglikelihood = g_ll + t_ll;