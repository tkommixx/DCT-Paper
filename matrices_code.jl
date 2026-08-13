function _get_prob_test(x::Human,test::Int64)
    M = get_matrix(x,test)
    asymp_red = 0.5 # This is the reduction rate that determines the likelihood of an asymptomatic individual testing positive.

    if !x.got_inf
        prob = 0
    elseif x.daysinf+1 > size(M,2)
        prob = 0
    else
        d = x.daysinf+1 #first row is 0
        prob = M[x.incubationp,d]
        prob = x.wentto > 0 ? prob*(asymp_red^(x.wentto-1)) : prob
    end
    return prob
end