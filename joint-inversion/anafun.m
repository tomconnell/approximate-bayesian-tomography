function [loglikelihood] = anafun(parameters,data,edge_effect,toprint)

g_data = data(1,1:20);
t_data = data(1,21:140);

vp_parameters = brocherizegrid(parameters);

[g_ll,g_e] = gravfun(parameters,g_data,edge_effect);

[t_ll,t_e] = tomfun(vp_parameters,t_data);

g_ll = g_ll*10^(-1.5);

if toprint
    
    dlmwrite('ana_misfit.csv',t_e+g_e,'-append')
    
    dlmwrite('gravity_misfit.csv',g_e,'-append')
    
    dlmwrite('tom_misfit.csv',t_e,'-append')
    
    dlmwrite('gravity_likelihood.csv',g_ll,'-append')
    
    dlmwrite('tom_likelihood.csv',t_ll,'-append')
    
end

loglikelihood = g_ll + t_ll;