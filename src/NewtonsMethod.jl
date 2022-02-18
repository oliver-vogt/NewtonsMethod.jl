module NewtonsMethod

using ForwardDiff, LinearAlgebra

function newtonroot(f, f′; x_0, tol=1E-7, maxiter=1000)
    # setup the algorithm
    x_old = x_0
    normdiff = Inf
    iter = 1
    while normdiff > tol && iter <= maxiter
        x_new = x_old - ( f(x_old) / f′(x_old) ) # use Newton's method
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    if iter <= maxiter
        return (value = x_old, normdiff=normdiff, iter=iter) # A named tuple
    else
        println("No unit root found with tolerance level $tol and Maximum Number of Iterations $maxiter")
        return nothing
    end
end

function newtonroot(f; x_0, tol=1E-7, maxiter=1000)
    # setup the algorithm
    x_old = x_0
    normdiff = Inf
    iter = 1
    D(f) = x -> ForwardDiff.derivative(f, x)
    f′ = D(f) # use auto diff
    while normdiff > tol && iter <= maxiter
        x_new = x_old - ( f(x_old) / f′(x_old) ) # use Newton's method
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    if iter <= maxiter
        return (value = x_old, normdiff=normdiff, iter=iter) # A named tuple
    else
        println("No unit root found with tolerance level $tol and Maximum Number of Iterations $maxiter")
        return nothing
    end
end

export newtonroot

end
