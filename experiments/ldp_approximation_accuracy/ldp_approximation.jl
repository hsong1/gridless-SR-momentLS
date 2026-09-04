using Random, Distributions
using Polynomials

function coeffs_random_poly(n, seed)
  Random.seed!(Int(seed))
  n = Int(n)
  coeff = rand(Normal(0,1/sqrt(n)), n+1)
  return(coeff)
end

function coeffs_truncated_rational(n, a)
    return [a^j for j in 0:n]
end




function ldp_test(cf, tol, verbose = true)
# test
    n = length(cf)-1 # highest degree in p(theta)
  	# number of sub-intervals in [0,1]
  	J = Int(ceil(log2(n)))
  	J0 = J+1 
    id_int = [1:J0;-J0:-1]
    
    # exact polynomial ~ approx polynomial??
    ldp_dict = Dict(id_int .=> [estimate_ldp(i,cf,tol,abs(i)==J+1) for i in id_int])
    
    # test range
    xs = range(0.0, 1.0; length=1001)
    # initialization
    diff_vec = zeros(length(id_int))
    normalized_diff = zeros(length(id_int))
    deg_vec = zeros(length(id_int))
    
    # l1 norm of coefficients
    cf_1 = sum(abs,cf)
    #@polyvar x
    
    if verbose println("comparision:") end
    
    for (j,i) in enumerate(id_int)
      tx = [transform(xi,i,abs(i)==J+1) for xi in xs]
      v0 = [evalpoly(tx_i, cf) for tx_i in tx]
      #println(v0)
      v1 = [evalpoly(xi,ldp_dict[i][1]) for xi in xs]
      #println(v1)
      #v0 = [p_exact(x => xi) for xi in xs]
      #v1 = [p_approx(x => xi) for xi in xs]
      diff = v0 .- v1
      max_abs_diff = maximum(abs, diff)
      imax = argmax(abs.(diff))
      if verbose println("id: $i (endpoint = $(ldp_dict[i][2])), max |p0 - p1| = $max_abs_diff at x = $(xs[imax])") end
      
      # output
      diff_vec[j] = max_abs_diff
      normalized_diff[j] = max_abs_diff / (tol*cf_1)
      deg_vec[j] = length(ldp_dict[i][1])-1
      end
    
    return(diff = diff_vec,normalized_diff=normalized_diff,maxdegree=deg_vec)
end