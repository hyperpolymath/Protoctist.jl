abstract type MemSpace end
struct Global   <: MemSpace end
struct Shared   <: MemSpace end
struct Local    <: MemSpace end
struct Constant <: MemSpace end

# Region-typed view (phantom type S encodes intended region)
struct RegionView{S<:MemSpace,T,N}
    data::CuDeviceArray{T,N}   # CUDA device array underneath
end

# Constructors for clarity
GlobalView(A::CuDeviceArray)   = RegionView{Global,eltype(A),ndims(A)}(A)
SharedView(A::CuDeviceArray)   = RegionView{Shared,eltype(A),ndims(A)}(A)
ConstantView(A::CuDeviceArray) = RegionView{Constant,eltype(A),ndims(A)}(A)

# Region-aware loads/stores: specialize by S
@inline function rload(::Type{Global},  A::RegionView{Global,T,N}, I...) where {T,N}
    @inbounds return A.data[I...]
end
@inline function rstore!(::Type{Global}, A::RegionView{Global,T,N}, v, I...) where {T,N}
    @inbounds A.data[I...] = v
    return nothing
end

# Shared-memory helpers (CUDA.jl has @cuStaticSharedMem; we wrap it)
# For dynamic shared mem, prefer KernelAbstractions.@localmem; else statically size:
# buf = CUDA.@cuStaticSharedMem(T, blkdim)  # then wrap as RegionView{Shared}
