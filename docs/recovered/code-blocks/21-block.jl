using CUDA, KernelAbstractions

# Region wrappers
abstract type MemSpace end; struct Global<:MemSpace end; struct Shared<:MemSpace end
struct RegionView{S<:MemSpace,T,N}; data::CuDeviceArray{T,N}; end
GlobalView(A) = RegionView{Global,eltype(A),ndims(A)}(A)

# Role types
abstract type Role end; struct Loader<:Role end; struct Compute<:Role end

# Kernel that prefetches (Loader) then computes (Compute); simplified
@kernel function mykernel(::Type{R}, A_global, C_global) where {R<:Role}
    i = @index(Global)
    @static if R <: Loader
        # simulate load; in reality, use @cuStaticSharedMem and a SharedView
        @inbounds C_global[i] = A_global[i]
    elseif R <: Compute
        @inbounds C_global[i] = 2f0 * C_global[i] + 1f0
    end
end

# Launch pipeline
function run_pipeline(A::CuArray{Float32}, C::CuArray{Float32})
    backend = CUDADevice()
    n = length(A)
    ndr = n
    mykernel(Loader,  backend, ndr)(A, C)
    mykernel(Compute, backend, ndr)(A, C)
    synchronize(backend)
    return C
end
