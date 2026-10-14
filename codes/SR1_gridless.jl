import Clarabel
# import SumOfSquares as SOS
using SumOfSquares
using Polynomials
import DynamicPolynomials as DP
using PolynomialRoots
# using MosekTools
using RCall
using Random
using Distributions
using LinearAlgebra
using LinearSolve
# using Plots
# using HomotopyContinuation
using StatsBase
using Quadmath

# Estimates polynomial coefficients for a specific dyadic interval transformation.
# Uses binomial expansion to compute coefficients that approximate the moment sequence
# over the transformed interval.
function estimate_ldp(i,cff,ldp_tol, last_interval)
  
  ## compute tilde_p_a(x;eps) on U_i
  eps = ldp_tol / 2 # approximating polynomial eps
  cff_1 = sum(abs.(cff)) # l1 norm of p(theta) ||cff||_1
  
	# degree of polynomial
	m = Int(ceil(exp(1+1/exp(1))*log(1/eps)))
	t = Int(floor(2^abs(i) * log(1/eps)))
	n = length(cff)-1 # largest degree 
	k_i = min(t,n)
	up = min(m-1,k_i)
	# right endpoint 
	a0 = !last_interval ? (1- 2.0^-abs(i)) : 1.0
	# interval length
	Ui_len = !last_interval ? (1-a0) : (2.0^(-abs(i)+1))
	
  # initialization for coefficient vector
	b = zeros(up+1)
	i_neg = (i < 0)
	
	# For better precision
	q_a0 = Quadmath.Float128(a0)
  q_Ui_len = Quadmath.Float128(Ui_len)

	for j in 0:up
	  # compute each b_j^{(i)} =  [sum_{k=j}^{k_i} sign(i)^k*c_k binom(k,j) a0^{k-j}] * (-1)^j* Ui_len^j
	  sum_val = zero(Quadmath.Float128)
	  
	  # recursive formula for binom(k, j) * a0^(k-j)
	  # [binom(k+1,j) * a0^(k+1-j)] = (k+1)/(k+1-j)*a0 * [binom(k, j) * a0^(k-j)]
	  # start with k = j : binom(j,j) * a0^0 = 1 
	  current_bterm = one(Quadmath.Float128)
	  
		for k in j:k_i
		  # inner term = sum_{k=j}^{k_i} sign(i)^k*c_k binom(k,j) a0^{k-j}
		  # sign(i)^k = -1 iff i < 0 and k is odd
		  current_sign = (i_neg && isodd(k)) ? Quadmath.Float128(-1.0) : Quadmath.Float128(1.0)
		  
		  sum_val += current_sign * Quadmath.Float128(cff[k+1]) * current_bterm
		  
		  # update current_bterm 
		  if k < k_i
		    current_bterm *= (k+1)/(k+1-j)*q_a0
		  end
			# b[j+1] += sign(i)^k * cff[k+1] * binomial(big(k), big(j)) * a0^(k-j) * Ui_len^j
		end
		val_j = sum_val * (-q_Ui_len)^j
		
		# Convert to standard float64
		val_j_f64 = Float64(val_j)
		if isnan(val_j_f64) || isinf(val_j_f64)
        error("estimate_ldp failed: Coefficient b[$(j+1)] became NaN. " * "Params: i=$i, eps=$eps, k_i=$k_i, a0=$a0")
    end
    b[j+1] =val_j_f64
	end
	
	## drop the largest high-degree tail whose abs-sum <= eps
	cs = cumsum(reverse(abs.(b)))
	jmax = findlast(<=(eps*cff_1), cs) # number of terms to drop
	if(jmax == nothing)
	  b_trunc = b # do not drop anything
  elseif (jmax == length(b))
    b_trunc = b[1:1] # keep constant term only
  else
    b_trunc = b[1:(length(b)-jmax)]
  end
  db_trunc = length(b_trunc)-1 # maxdegree of truncated polynomial
  
	return(b_trunc,a0,up,db_trunc)
end

# Maps a value x in [0,1] to a value in a specific interval determined by index i.
# The transformation uses dyadic intervals based on powers of 2.
function transform(x,i,last_interval)
  a0 = !last_interval ? (1- 2.0^-abs(i)) : 1
  Ui_len = !last_interval ? (1-a0) : (2.0^(-abs(i)+1))
  return(sign(i)*(a0 - x*Ui_len))
end

function get_transform_params(ind_val, last_interval)
  a0 = !last_interval ? (1- 2.0^-abs(ind_val)) : 1.0
  Ui_len = !last_interval ? (1.0-a0) : (2.0^(-abs(ind_val)+1))
  # coefficients for tx(x) = c0 + c1*x
  c0 = sign(ind_val) * a0
  c1 = -sign(ind_val) * Ui_len
  return c0,c1
end



function coeffs_dd_polys(ind_val, last_interval, coeffs_i, weight_option, V, supp, weight, h_power)
  # obtain theta(x)
  c0, c1 = get_transform_params(ind_val, last_interval)
  tx = Polynomial([c0,c1]) # theta(x) = c0 + c1*x
  
  # collect terms: 1-tx*supp[j]
  terms = [1 - tx*s for s in supp]
  s_len = length(supp)
  
  # denominator: product of (1 - theta * alpha_i)
  den = s_len > 0 ? prod(terms) : Polynomial([1.0])
  
  # num1: den * approx_poly
  num1 = den*Polynomial(coeffs_i)
  
  # num2: 
  # when V = Z, sum_{i=1}^s w_i* (1+ theta * alpha_i)* prod_{j!=i}(1-theta * alpha_j)
  # when V = N, sum_{i=1}^s w_i* prod_{j!=i}(1-theta * alpha_j)
  num2 = Polynomial([0.0])
  if s_len > 0
    prod1 = accumulate(*,terms) # prod_{j<=i} (1-theta * alpha_j)
    prod2 = accumulate(*,reverse(terms)) # prod_{j>=i} (1-theta * alpha_j)
    reverse!(prod2)
    
    # lo_i = prod_{j!=i}(1-theta * alpha_j)
    for i in 1:s_len
      if i==1
        lo_i = s_len > 1 ? prod2[2] : Polynomial([1.0]) # prod_{j>=2} (1-theta * alpha_j)
      elseif i == s_len
        lo_i = prod1[s_len-1] # prod_{j<=(s-1)} (1-theta * alpha_j)
      else
        lo_i = prod1[i-1] * prod2[i+1] # prod_{j<=i-1} (1-theta * alpha_j) * prod_{j>=i+1} (1-theta * alpha_j)
      end
      
      basis = (V == "Z") ? (1 + tx * supp[i]) : Polynomial([1.0])
      num2 += weight[i] * basis * lo_i
    end
  end # end of computing num2
  
  
  bracket_total = num1 - num2 
  
  
  if weight_option == "unweighted"
      num = -2 * bracket_total
  else
      # h(theta) = 1 - |theta|^h_power
      if ind_val > 0 # when theta > 0
          h_poly = 1 - tx^h_power
      else
          h_poly = 1 - (-tx)^h_power
      end
      num = -2 * h_poly * bracket_total
  end
  
  return (num_coeffs = Polynomials.coeffs(num), den_coeffs= Polynomials.coeffs(den))
end

# compute exact (weighted) directional derivative at (supp, weight) for comparison
# x is a vector where we want to evaluate directional derivatives 
function dd_poly_exact(x, ind_val, last_interval, cff, weight_option, V, supp, weight, h_power)
    
    # if x is a scalar, convert x to a vector of length 1
    x_vec = (x isa Number) ? [x] : x
    # output initialization
    results = zeros(Float64, length(x_vec))
    
    # get transformation 
    c0, c1 = get_transform_params(ind_val, last_interval)
    
    for (idx, x) in enumerate(x_vec)
      tx = c0 + c1 * x
      # tx = transform(x, ind_val, last_interval)
      
      # inner product
      # term1 = sum(cff[i] * tx^(i-1) for i in eachindex(cff))
      term1 = evalpoly(tx, cff)
      
      # term2 = sum_{i=1}^s w_i * (\sum_{j in V} (tx*supp_i)^{|j|})
      term2 = 0.0
      for i in 1:length(supp)
        ip_i = (V=="Z") ? (1+tx*supp[i])/(1-tx*supp[i]) : 1/(1-tx*supp[i])
        term2 += weight[i]*ip_i
      end
      
      dd_exact = -2.0*(term1 -term2)
      # if weighted, multiply dd_exact with (1-|tx|^h_power)
      if weight_option != "unweighted"
        dd_exact *= 1-abs(tx)^h_power
      end
      results[idx] = dd_exact
    end
    
    return(results)
end


function dd_poly_exact_original(theta, cff, weight_option, V, supp, weight, h_power)
    
    # if theta is a scalar, convert theta to a vector of length 1
    theta_vec = (theta isa Number) ? [theta] : theta
    # output initialization
    results = zeros(Float64, length(theta_vec))
    
    # get transformation 
    # c0, c1 = get_transform_params(ind_val, last_interval)
    
    for (idx, tx) in enumerate(theta_vec)
      # tx = transform(x, ind_val, last_interval)
      
      # inner product
      # term1 = sum(cff[i] * tx^(i-1) for i in eachindex(cff))
      term1 = evalpoly(tx, cff)
      
      # term2 = sum_{i=1}^s w_i * (\sum_{j in V} (tx*supp_i)^{|j|})
      term2 = 0.0
      for i in 1:length(supp)
        ip_i = (V=="Z") ? (1+tx*supp[i])/(1-tx*supp[i]) : 1/(1-tx*supp[i])
        term2 += weight[i]*ip_i
      end
      
      dd_exact = -2.0*(term1 -term2)
      # if weighted, multiply dd_exact with (1-|tx|^h_power)
      if weight_option != "unweighted"
        dd_exact *= 1-abs(tx)^h_power
      end
      results[idx] = dd_exact
    end
    
    return(results)
end

##### solvers for argmin_x num(x)/den(x) #####

function minimize_rf_SOS(ll, ul, num_coeffs, den_coeffs)
  
    solver = Clarabel.Optimizer
    # set up the optimizer
    DP.@polyvar x
		model = SOSModel(solver)
    set_string_names_on_creation(model, false)
    
    # create numerator and denominator polynomials
    # num = poly_from_coeffs(num_coeffs, x)
    # den = poly_from_coeffs(den_coeffs, x)
    num = sum(num_coeffs[i] * x^(i-1) for i in eachindex(num_coeffs))
    den = sum(den_coeffs[i] * x^(i-1) for i in eachindex(den_coeffs))
		#num, den = dd_polys(x, ind_val, p, weight_option, V, supp, weight)
		S = @set (ll <= x) && (x <= ul)
		
		# optimize
		@variable(model,s)
		@constraint(model,c, num >= s*den, domain = S)
		@objective(model, Max, s)
		set_silent(model)
		optimize!(model)
		
		# warn on "bad" termination
    if termination_status(model) != MOI.OPTIMAL
        @warn "SOS solve did not terminate OPTIMAL"
    end

     # values
     s_val = try value(s) catch; NaN end
     objval = try objective_value(model) catch; NaN end
     
    diff = abs(s_val - objval)
    if diff > 1e-6
        @warn "Consistency check failed: value(s) and objective_value(model) differ" s_val objval diff
        # Terminate and return 1000
        return Inf, NaN 
    end
    
    # otherwise form a moment matrix
    v = moment_matrix(model[:c])
    
    # Extract an atom assuming Rank 1
    pt = nothing
    pt = atomic_measure(v, FixedRank(1))
    
    if pt !== nothing
      x_atom = pt.atoms[1].center[1]
      # val_at_atom= evalpoly(x_atom, num_coeffs) / evalpoly(x_atom, den_coeffs)
      val_at_atom = num(x=>x_atom)/den(x=>x_atom)
      gap = val_at_atom - s_val
      
      # If Rank-1 is "good enough", return it immediately
        if abs(gap) <= 1e-4
            return (s_value = s_val, x_atom = x_atom, M_r = 1)
        end
      @warn "SOS optimization: Rank 1 validation failed" gap x_atom
    end
    
    # Only runs if Rank 1 failed validation or returned nothing
    # numerically determine the rank of the moment matrix V.Q    
    s_vals = svdvals(Matrix(v.Q)) 
    v_r = count(s -> s > 1e-6 * s_vals[1], s_vals)
    # extract all minimizers
    pt_r = atomic_measure(v, FixedRank(Int(v_r)))
    if pt_r !== nothing
      x_r = pt_r.atoms[1].center[1]
      # val_r = evalpoly(x_r, num_coeffs) / evalpoly(x_r, den_coeffs)
      val_r = num(x=>x_r)/den(x=>x_r)
      gap = val_r - s_val
      
      if abs(gap) <= 1e-4
            return (s_value = s_val, x_atom = x_r, M_r = Int(v_r))
      end
    end
    
    @error "SOS optimization: All attempts failed for interval [$ll, $ul]"
    return (s_value = Inf, x_atom = NaN, M_r = Int(v_r))
end

function minimize_rf_roots(ll, ul, num_coeffs, den_coeffs)
    
    
    # 1. Construct the polynomials f and g
    f = Polynomial(num_coeffs)
    g = Polynomial(den_coeffs)
    
    # 2. Compute the numerator of the derivative: f'g - fg'
    f_prime = Polynomials.derivative(f)
    g_prime = Polynomials.derivative(g)
    
    Q = f_prime * g - f * g_prime
    
    # 3. Extract coefficients of Q to find roots
    # q_coeffs = coefficients(Q, [x^i for i in 0:maxdegree(Q)])
    # all_roots = PolynomialRoots.roots(convert(Vector{Float64}, Q))
    c = collect(Float64,Polynomials.coeffs(Q))

    if abs(c[end]) < 1e-300 || abs(c[1] / c[end]) > 1e300 @warn "Scaling detected as dangerous." end
    all_roots = Polynomials.roots(Q)
    #all_roots = PolynomialRoots.roots(convert(Vector{Float64}, Polynomials.coeffs(Q)))
    
    # 4. Filter for real roots within the interval [ll, ul]
    # We allow a small imaginary tolerance for numerical noise
    candidate_xs = Float64[ll,ul]
    for r in all_roots
        #if abs(imag(r)) < 1e-9
            real_r = real(r)
            if ll <= real_r <= ul
                push!(candidate_xs, real_r)
            end
        #end
    end
    
    # 5. Evaluate the ratio at all candidates to find the global minimum
    best_x = NaN
    min_val = Inf  # Start at infinity so any real value is smaller

    for cx in candidate_xs
        # Evaluate ratio: num(x) / den(x)
        current_val = f(cx) / g(cx)
        
        if isless(current_val, min_val)
            min_val = current_val
            best_x = cx
        end
    end
    
    return (min_val = min_val, x_at_min = best_x)
end


function dd_optimize(p_dict, weight_option,V, supp, weight, delta, opt_type = "root", h_power = 1)
  
  J0 = (V == "Z") ? Int(length(p_dict) / 2) : length(p_dict)
  num_entries = (V == "Z") ? 2*J0 : J0
  gradients = zeros(num_entries)
  supports = zeros(num_entries)
	local_xs = zeros(num_entries)
	lower_bounds=zeros(num_entries)
	upper_bounds=zeros(num_entries)
	
	# set up sub-intervals
	intervals = (V == "Z") ? zip([1:1:2*J0;],append!([1:1:J0;], [-J0:1:-1;])) : zip([1:1:J0;],[1:1:J0;])
	
	
	for (ind_idx, ind_val) in intervals
    # println(ind_val)
    last_interval = (p_dict[ind_val][2] == 1) # if a0 = 1
    # create numerator and denominator polynomials
    num_coeffs, den_coeffs = coeffs_dd_polys(ind_val, last_interval, p_dict[ind_val][1], weight_option, V, supp, weight, h_power)

		# the search domain S
		if abs(ind_val) == J0 # end intervals
		  ll= !last_interval ? (2.0^(abs(J0)) * delta - 1.0) : (2.0^(abs(J0)-1)*delta)
		  ul = 1.0
		else
			ll = 0.0
			ul = 1.0
		end
		
		if(opt_type == "root")
		  res = minimize_rf_roots(ll,ul,num_coeffs, den_coeffs)
		  m_v, x_m = res.min_val, res.x_at_min
		else
		  res = minimize_rf_SOS(ll,ul,num_coeffs,den_coeffs)
		  m_v, x_m = res.s_value, res.x_atom
		end
		# println("ind_val: $ind_val, m_v: $m_v, x_m: $x_m")
		
      # store results		
		  gradients[ind_idx] = m_v
		  local_xs[ind_idx] = x_m
		  lower_bounds[ind_idx] = ll
		  upper_bounds[ind_idx] = ul
		  supports[ind_idx] = (x_m !== NaN) ? transform(x_m,ind_val,last_interval) : NaN
		# println("ind_val: $ind_val, grad: $(gradients[ind_idx]), supp: $(supports[ind_idx])")
	end # end of loopting through sub-intervals
	
	return(gradients=gradients, supports=supports, local_xs = local_xs, ll = lower_bounds, ul = upper_bounds)
end


  
  

# Performs a single Support Reduction (SR) step by solving for optimal weights given support points.
# Ensures non-negative weights by iteratively removing support points with negative weights.
# Parameters:
#   - supp: current support points
#   - weight: current weights
#   - r: input sequence
# Returns: (updated support, updated weights)
SR = function(V, supp, weight,r)
		validmeasure = false
		current = weight
	exponent = [0:1:length(r)-1;]
		while(!validmeasure)
		  # B[i,j]=<x_{alpha_i},x_{alpha_j}>
			B = zeros(length(supp), length(supp))
			for i in 1:length(supp)
				for j in 1:length(supp)
					B[i,j] = (V == "Z") ? (1+supp[i]*supp[j])/(1-supp[i]*supp[j]) : 1/(1-supp[i]*supp[j])
				end
			end
			
			# c[i] = <r, x_{alpha_i}>
			c = zeros(length(supp))
			for i in 1:length(supp)
				c[i] = (V == "Z") ? (2*sum((supp[i].^exponent) .* r) - r[1]) : sum( (supp[i].^exponent) .* r)
			end
			
			prob = LinearProblem(B,c)
			sol = LinearSolve.solve(prob)
			new = sol.u
			
			
			if(all( >=(0), new))
				validmeasure = true
				weight = new
			else
				t = - current ./(new .- current)
				pop!(t)
				t[t .< 0] .= Inf
				t[t .> 1] .= Inf
				bd = findmin(t)
				current = (1-bd[1]).*current + bd[1] .* new
				deleteat!(supp, bd[2])
				deleteat!(current,bd[2])
			end

		end
	return(supp, weight)
end




# Main gridless Support Reduction algorithm (SR1) for solving the moment problem.
# Iteratively adds support points by solving SOS optimization problems over multiple intervals
# until convergence or maximum iterations reached.
# Parameters:
#   - r: input sequence to match
#   - delta: precision parameter for interval bounds
#   - weight_option: "unweighted" or "weighted" (h(theta)=1-|theta|)
#   - V: "Z" for [-1,1] support, otherwise [0,1]
#   - maxit: maximum number of iterations
#   - init_supp: initial support points
#   - init_weight: initial weights
#   - tol: convergence tolerance
#   - opt_type: "root" for root-finding algorithm; otherwises it uses SOS optimization.
#   - verbose: print iteration details
# Returns: (final support points, final weights)
function SR1_gridless(r, delta, weight_option, V, maxit,init_supp, init_weight, tol, h_power = 1, opt_type = "root", verbose = false)
	
	 
  	# if(weight_option == "unweighted" && delta < 10e-12)
  	#	throw(ArgumentError("Delta is too small to use the unweighted SR algorithm"))
  	# end
  	
  	### low degree approximation over sub-intervals ###
  	# polynomial to be approximated: p(theta) = sum_{j=0}^n cff[j+1]*theta^j
  	if(V == "Z")
  	    cff = copy(r)
  	    cff[2:end] *= 2
    else
        cff = copy(r)
    end
    cff_1 = sum(abs.(cff)) # l1 norm of p(theta) ||cff||_1
    ldp_tol = tol / (6*cff_1)
    n = length(cff)-1 # highest degree in p(theta)
  	
  	# number of sub-intervals in [0,1]
  	J = Int(ceil(log2(n)))
  
    if(ceil(-log2(delta)) > J)
	    J0 = J+1
    else
  	  # 1-2^{-J0+1} <= 1-delta < 1 - 2^{-J0}
  	  # <-> log2(1/delta) < J0 <= 1+log2(1/delta) <-> J0 = ceil(log2(1/delta))
  	  J0 = Int(ceil(-log2(delta)))
    end
  	
  	
  	# create sub-intervals
    # obtain low-degree polynomial approximations
    if(V == "Z")
        id_int = append!([1:1:J0;], [-J0:1:-1;])
    else
        id_int = [1:1:J0;]
    end
    ldp_dict = Dict(id_int .=> [estimate_ldp(i,cff,ldp_tol,abs(i)==(J+1)) for i in id_int])
    
    if verbose
        println("$(length(id_int)) intervals with endpoints: ", [ldp_dict[i][2] for i in id_int])      
    end
    
    # initialization
    #exponents = [0:1:n-1;]
    supp = init_supp isa AbstractVector ? init_supp : Float64[init_supp]
    weight = init_weight isa AbstractVector ? init_weight : Float64[init_weight]
    conv = false
    #eftol = 6*tol*cff_1
    iter = 0
    min_grad_exact = Inf
    
    if verbose println("Starting SR1 Gridless: tol = $tol, maxit = $maxit") end
    
  
    while(iter < maxit && !conv)
        SRstep = SR(V, supp, weight,r) 
        supp = SRstep[1]
        weight = SRstep[2]
        
        if verbose
            println("--------------------------------")
            println("Iter $iter:")
            println("Current Support: ", supp)
            println("Current Weights: ", weight)
        end
        
        res = dd_optimize(ldp_dict, weight_option,V, supp, weight, delta, opt_type, h_power)
        
        # println("end of finding new support!!!")
        min_grad_approx, index = findmin(res.gradients)
        best_support_candidate = res.supports[index]
        
        # compute the exact directional derivative at best_support_candidate
        best_ind_val = id_int[index]
        x_m = res.local_xs[index]
        last_interval = (ldp_dict[best_ind_val][2] == 1)
        min_grad_exact = only(dd_poly_exact(x_m,best_ind_val,last_interval,cff,weight_option,V,supp,weight,h_power))
        
        if(min_grad_exact > -tol)
            conv = true
            if verbose println("min_grad: $min_grad_exact, -tol: $(-tol), conv: $conv") end
        else
          append!(supp, res.supports[index])
          append!(weight, 0)
        end
        
        
        if verbose
            println("End of finding a new support point")
            println("New point = $(res.supports[index]), Min Grad = $min_grad_exact, Min Grad (approx) = $min_grad_approx")
            println("New Support: ", supp)
            println("Current Weights: ", weight)
        end
        
        iter = iter + 1
    end
    # println(iter)
    if(!conv)
        println("failed to converge")
    end
  	
  	#pop!(supp)
  	#pop!(weight)
  	
  	return (
    	  support = supp, 
        weights = weight, 
        iterations = iter,
        min_grad = min_grad_exact,
        ldp_tol = ldp_tol,
        tol = tol,
        converged = conv
      )
end





function delta_h(r, delta, weight_option, V, supp, weight, ldp_tol, h_power=1, opt_type = "root", verbose = false)
  	
  	### low degree approximation over sub-intervals ###
  	# polynomial to be approximated: p(theta) = sum_{j=0}^n cff[j+1]*theta^j
  	if(V == "Z")
  	    cff = copy(r)
  	    cff[2:end] *= 2
    else
        cff = copy(r)
    end
    cff_1 = sum(abs.(cff)) # l1 norm of p(theta) ||cff||_1
    n = length(cff)-1 # highest degree in p(theta)
  	
  	# number of sub-intervals in [0,1]
  	J = Int(ceil(log2(n)))
  
    if(ceil(-log2(delta)) > J)
	    J0 = J+1
    else
  	  # 1-2^{-J0+1} <= 1-delta < 1 - 2^{-J0}
  	  # <-> log2(1/delta) < J0 <= 1+log2(1/delta) <-> J0 = ceil(log2(1/delta))
  	  J0 = Int(ceil(-log2(delta)))
    end
  	
  	
  	# create sub-intervals
    # obtain low-degree polynomial approximations
    if(V == "Z")
        id_int = append!([1:1:J0;], [-J0:1:-1;])
    else
        id_int = [1:1:J0;]
    end
    ldp_dict = Dict(id_int .=> [estimate_ldp(i,cff,ldp_tol,abs(i)==(J+1)) for i in id_int])
    
    if verbose
        println("$(length(id_int)) intervals with endpoints: ", [ldp_dict[i][2] for i in id_int])      
    end
    
    # initialization
    #exponents = [0:1:n-1;]
    supp = supp isa AbstractVector ? supp : Float64[supp]
    weight = weight isa AbstractVector ? weight : Float64[weight]
    res = dd_optimize(ldp_dict, weight_option,V, supp, weight, delta, opt_type, h_power)

    min_grad_approx, index = findmin(res.gradients)
    best_support_candidate = res.supports[index]

    # compute the exact directional derivative at best_support_candidate
    best_ind_val = id_int[index]
    x_m = res.local_xs[index]
    last_interval = (ldp_dict[best_ind_val][2] == 1)
    min_grad_exact = only(dd_poly_exact(x_m,best_ind_val,last_interval,cff,weight_option,V,supp,weight,h_power))

    epsilon_0 = 2*ldp_tol*cff_1
    
    # delta_h(mu) = max_theta -D_theta^h(mu) = -D_(theta*)^h(mu) 
    
    # output:
    # delta_h_candidate = -D_(theta_tilde)^h(mu)
    # delta_h_ldp = max_theta -D_theta^h(mu;ldp_tol) = -D_(theta_tilde)_^h(mu;ldp_tol) 
    delta_h_candidate = -min_grad_exact
    delta_h_ldp = -min_grad_approx

    # |D_(theta*)^h(mu) -D_theta*^h(mu;ldp_tol)|<=epsilon_0
    # -> -epsilon_0 + D_theta*^h(mu;ldp_tol)<=D_(theta*)^h(mu) 
    # -> -epsilon_0 + D_theta_tilde^h(mu;ldp_tol)<= D_(theta*)^h(mu) 
    # i.e., delta_h(mu)<= delta_h_ldp + epsilon_0
  	return (
      delta_h_candidate = delta_h_candidate,
      delta_h_ldp = delta_h_ldp,
      delta_h_upper = delta_h_ldp + epsilon_0,
      theta_tilde = best_support_candidate,
      ldp_tol = ldp_tol,
      epsilon_0 = epsilon_0
    )
end