
#  USER GUIDE
#
# Brief descriptions of the functions in this file and their use.
#
# function julia11hypot(x, y)
# A direct copy of the original Julia 1.1 hypot code.
#
# function MyHypot1(x,y)
# Implements the Naive (Unfused) algorithm from the paper.
#
# function MyHypot2(x,y)
# Implements the Naive (Fused) algorithm from the paper.
#
# function MyHypot3(x,y)
# Implements the Corrected (Unfused) algorithm from the paper.
#
# function MyHypot4(x,y)
# Implements the Corrected (Fused) algorithm from the paper.
#
# function MyHypot5(x,y)
# Implements the Fully Corrected (Fused) algorithm from the paper.
#
# function testuniform(n, m)
# Runs n random trials of the first six algorithms with x,y~U(1,2) and computes
# percentage of incorrectly rounded results.
# The optional argument m rescales so that x~U(m,2*m).
# This code is used to generate the data in table 2.
#
# function multiscaletestuniform(n)
# Runs the testuniform code on a range of relative scales 2^k for k=0..11.
# Returns a matrix with the resulting error rates where each row corresponds
# to a relative scale of the operands and each column to one of the six
# methods.
# Input variable n is the number of trials to run at each scale.
# This code generates the data in table 2 by calling testuniform().
#
# function testnormal(n)
# Runs n random trials of the algorithms with x,y~N(0,1) and computes
# percentage of one ulp and two ulp errors.
# This code was used to create table 1 in the paper.

using Random
using Printf

function julia11hypot(x::T, y::T) where T<:Number #Julia 1.1 original hypot code.
    x = abs(x)
    y = abs(y)
    if x < y
        x, y = y, x
    end
    if iszero(x)
        r = y / oneunit(x)
    else
        r = y / x
    end
    rr = x * sqrt(1 + r * r)

    # Use type of rr to make sure that return type is the same for
    # all branches
    if isnan(r)
        isinf(x) && return oftype(rr, Inf)
        isinf(y) && return oftype(rr, Inf)
        return oftype(rr, r)
    else
        return rr
    end
end

function MyHypot1(x::T,y::T) where T<:AbstractFloat  # Naive (Unfused)
    # Take the absolute values of the inputs.
    absx, absy = abs(x), abs(y)
    # Return Inf if either or both inputs is Inf (Compliance with IEEE754)
    isinf(absx) && return absx
    isinf(absy) && return absy
    # Order the operands
    x = ifelse(absx < absy, absy, absx)
    y = ifelse(absx < absy, absx, absy)

    # Widely varying operands
    if y < x*sqrt(eps(T)/2)  #Note: This also gets y == 0
        return x
    end

    # Operands do not vary widely
    scale = eps(sqrt(floatmin(T)))  #Rescaling constant
    if x > sqrt(floatmax(T)/2)
        x = x*scale
        y = y*scale
        scale = inv(scale)
    elseif y < sqrt(floatmin(T))
        x = x/scale
        y = y/scale
    else
        scale = one(scale)
    end
    sqrt(x*x+y*y)*scale
end

function MyHypot2(x::T,y::T) where T<:AbstractFloat  # Naive (Fused)
    # Take the absolute values of the inputs.
    absx, absy = abs(x), abs(y)
    # Return Inf if either or both inputs is Inf (Compliance with IEEE754)
    isinf(absx) && return absx
    isinf(absy) && return absy
    # Order the operands
    x = ifelse(absx < absy, absy, absx)
    y = ifelse(absx < absy, absx, absy)

    # Widely varying operands
    if y <= x*sqrt(eps(T)/2)  #Note: This also gets y == 0
        return x
    end

    # Operands do not vary widely
    scale = eps(sqrt(floatmin(T)))  #Rescaling constant
    if x > sqrt(floatmax(T)/2)
        x = x*scale
        y = y*scale
        scale = inv(scale)
    elseif y < sqrt(floatmin(T))
        x = x/scale
        y = y/scale
    else
        scale = one(scale)
    end
    sqrt(fma(x,x,y*y))*scale
end

function MyHypot3(x::T,y::T) where T<:AbstractFloat  # Corrected (Unfused)
    # Take the absolute values of the inputs.
    absx, absy = abs(x), abs(y)
    # Return Inf if either or both inputs is Inf (Compliance with IEEE754)
    isinf(absx) && return absx
    isinf(absy) && return absy
    # Order the operands
    x = ifelse(absx < absy, absy, absx)
    y = ifelse(absx < absy, absx, absy)

    # Widely varying operands
    if y <= x*sqrt(eps(T)/2)  #Note: This also gets y == 0
        return x
    end

    # Operands do not vary widely
    scale = eps(sqrt(floatmin(T)))  #Rescaling constant
    if x > sqrt(floatmax(T)/2)
        x = x*scale
        y = y*scale
        scale = inv(scale)
    elseif y < sqrt(floatmin(T))
        x = x/scale
        y = y/scale
    else
        scale = one(scale)
    end
    h = sqrt(x*x+y*y)

    # This is a well balanced single branch code
    # delta = x-h
    #h += ( delta*(2*(x-y)) + (2*delta+y)*y - delta*delta)/(2*h)
    # End single branch

    # This is a well balanced twopart branch code
    if h <= 2*y
        delta = y-h
        h += (x*(2*delta + x)-(delta+2*(x-y))*delta)/(2*h)
    else
        delta = x-h
        h += ( 2*delta*(x-2*y) + (4*delta+y)*y - delta*delta)/(2*h)
    end
    # End two branch code
    h*scale
end


function MyHypot4(x::T,y::T) where T<:AbstractFloat  # Corrected (Fused)
    # Take the absolute values of the inputs.
    absx, absy = abs(x), abs(y)
    # Return Inf if either or both inputs is Inf (Compliance with IEEE754)
    isinf(absx) && return absx
    isinf(absy) && return absy
    # Order the operands
    x = ifelse(absx < absy, absy, absx)
    y = ifelse(absx < absy, absx, absy)

    # Widely varying operands
    if y <= x*sqrt(eps(T)/2)  #Note: This also gets y == 0
        return x
    end

    # Operands do not vary widely
    scale = eps(sqrt(floatmin(T)))  #Rescaling constant
    if x > sqrt(floatmax(T)/2)
        x = x*scale
        y = y*scale
        scale = inv(scale)
    elseif y < sqrt(floatmin(T))
        x = x/scale
        y = y/scale
    else
        scale = one(scale)
    end
    # Compute the first order terms
    x_sq = x*x
    y_sq = y*y
    sigma = x_sq+y_sq
    h = sqrt(sigma)

    # Compute tau using higher order terms and apply the correction
    sigma_e = y_sq - (sigma-x_sq)
    tau = fma(y,y,-y_sq) + fma(x,x,-x_sq) + sigma_e + fma(-h,h,sigma)
    h = fma(tau/h,.5,h)
    h*scale
end

function MyHypot5(x::T,y::T) where T<:AbstractFloat  # Fully Corrected (Fused)
    # Take the absolute values of the inputs.
    absx, absy = abs(x), abs(y)
    # Return Inf if either or both inputs is Inf (Compliance with IEEE754)
    isinf(absx) && return absx
    isinf(absy) && return absy
    # Order the operands
    x = ifelse(absx < absy, absy, absx)
    y = ifelse(absx < absy, absx, absy)

    # Widely varying operands
    if y <= x*sqrt(eps(T)/2)  #Note: This also gets y == 0
        return x
    end

    # Operands do not vary widely
    scale = eps(sqrt(floatmin(T)))  #Rescaling constant
    if x > sqrt(floatmax(T)/2)
        x = x*scale
        y = y*scale
        scale = inv(scale)
    elseif y < sqrt(floatmin(T))
        x = x/scale
        y = y/scale
    else
        scale = one(scale)
    end
    # Compute the first order terms
    x_sq = x*x
    y_sq = y*y
    sigma = x_sq+y_sq
    h = sqrt(sigma)

    # Compute the second order terms and apply the correction
    psum = fma(y,y,-y_sq) + fma(x,x,-x_sq) + (y_sq - (sigma-x_sq))
    h = fma((psum+fma(-h,h,sigma))/h,.5,h)

    # Phase one test for a correctly rounded result
    uh = eps(h)
    tau = psum+fma(-h,h,sigma)
    tauerror = 25*uh*uh
    if (tau >= 0)
        delta = uh
        correctly_rounded = (tau < (delta*h-tauerror))
    else
        # Since tau < 0 we use ulp_L(h) to define delta.
        delta = (significand(h) == 1) ? -uh/2 : -uh
        correctly_rounded = (tau > (delta*h+tauerror))
    end

    # If correct rounding cannot be verified at this point then we generate a
    # correctly rounded result with some additional work.
    if !correctly_rounded

        # Do the signofsum calculation for
        s = [-delta*h,fma(-h,h,sigma),y_sq-(sigma-x_sq),fma(x,x,-x_sq),fma(y,y,-y_sq),-delta*delta/4]
        while any(v->abs(v) > abs(s[6]/6),s[1:5])
            for i=1:5
                d = s[i]+s[i+1]
                s[i] = (abs(s[i]) > abs(s[i+1])) ? (s[i]-d)+s[i+1] : s[i] = (s[i+1]-d)+s[i]
                s[i+1] = d
            end
        end
        Phi = sign(s[6])

        # Apply any further needed corrections
        if Phi == sign(tau)
            h += delta
        elseif Phi == 0
        h += delta/2  # Ensure ties resolved by current mode.
        end
    end

    h*scale
end

function testuniform(n::Integer, m::Float64=1.0)
    # Runs n random trials of the algorithms with x,y~U(1,2) and computes
    # percentage of incorrectly rounded results.
    # The optional argument m rescales so that x~U(m,2*m)

    # Initialize
    Random.seed!(1234)
    anyulperrs = zeros(6,1)

    for i=1:n
        # Generate the random inputs
        a = ([1.0;1.0]+rand(2)).*[m;1]
        # Run each of the methods
        hnative = julia11hypot(a[1],a[2])
        hclib = ccall("hypot", Float64, (Float64, Float64),a[1],a[2])
        hmine1 = MyHypot1(a[1],a[2])
        hmine2 = MyHypot2(a[1],a[2])
        hmine3 = MyHypot3(a[1],a[2])
        hmine4 = MyHypot5(a[1],a[2])
        # Compute the corrrectly rounded result with MPFR
        abig, bbig = big(a[1]), big(a[2])
        hbig = convert(Float64,hypot(abig,bbig))
        # Find those that are incorrectly rounded and update the count.
        tmp=[abs(hnative-hbig);abs(hclib-hbig);abs(hmine1-hbig);
        abs(hmine2-hbig);abs(hmine3-hbig);abs(hmine4-hbig)]
        scale = 1/eps(hbig)
        tmp = tmp*scale
        anyulperrs += map(x->(x>0.0) ? 1.0 : 0.0 ,tmp)
    end

    # Rescale counts to a percentage and return.
    anyulperrs/n*100
end

function multiscaletestuniform(n::Integer=10000)
    # Runs the testuniform code on a range of relative scales 2^k for k=0..11.
    # Returns a matrix with the resulting error rates where each row corresponds
    # to a relative scale of the operands and each column to one of the six
    # methods.
    # Input variable n is the number of trials to run at each scale.
    A = zeros(6,12)
    for i=1:12
        A[:,i] = testuniform(n,2.0^(i-1))
    end
    A = A'
end

function testnormal(n::Integer)
    # Runs n random trials of the algorithms with x,y~N(0,1) and computes
    # percentage of one ulp and two ulp errors.

    # Initialize
    Random.seed!(1234)
    oneulperrs = zeros(6,1)
    twoulperrs = zeros(6,1)

    for i=1:n
        # Generate the random inputs
        a = randn(2)
        # Run each of the methods
        hnative = julia11hypot(a[1],a[2])
        hclib = ccall("hypot", Float64, (Float64, Float64),a[1],a[2])
        hmine1 = MyHypot1(a[1],a[2])
        hmine2 = MyHypot2(a[1],a[2])
        hmine3 = MyHypot3(a[1],a[2])
        hmine4 = MyHypot4(a[1],a[2])
        # Compute the corrrectly rounded result with MPFR
        abig, bbig = big(a[1]), big(a[2])
        hbig = convert(Float64,hypot(abig,bbig))
        # Find those that are incorrectly rounded and update the counts.
        tmp=[abs(hnative-hbig);abs(hclib-hbig);abs(hmine1-hbig);
        abs(hmine2-hbig);abs(hmine3-hbig);abs(hmine4-hbig)]
        scale = 1/eps(hbig)
        tmp = tmp*scale
        oneulperrs += map(x->(x==1.0) ? 1.0 : 0.0 ,tmp)
        twoulperrs += map(x->(x>1.0) ? 1.0 : 0.0 ,tmp)
    end

    # Rescale counts to percentages and print.
    M = round.(oneulperrs/n*100;digits=2)
    println("Percentage of 1 ulp deviations ",M)
    M = round.(twoulperrs/n*100;digits=2)
    println("Percentage of 2 ulp deviations ",M)
    println("Number of 2 ulp deviations ",twoulperrs)

end
