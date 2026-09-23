abstract type Role end
struct Loader  <: Role end
struct Compute <: Role end
struct Storer  <: Role end

# Kernel specialized by role
using KernelAbstractions

@kernel function staged_kernel(::Type{R}, A, B, C) where {R<:Role}
    i = @index(Global)
    @static if R <: Loader
        # e.g., prefetch from global to shared
        # rstore!(Shared, sharedA, rload(Global, globalA, i), i)
    elseif R <: Compute
        # math on shared tiles
    elseif R <: Storer
        # write back to global
    end
end

# Launch three specialized passes (or use cooperative groups in one kernel)
staged_kernel(Loader,  backend, ndrange)(A,B,C)
staged_kernel(Compute, backend, ndrange)(A,B,C)
staged_kernel(Storer,  backend, ndrange)(A,B,C)
