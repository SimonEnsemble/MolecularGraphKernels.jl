module Test_TypePiracy

using MolecularGraphKernels, Test

## https://discourse.julialang.org/t/pirate-hunter/20402 ##

function all_methods(mainmod::Module=Main)
    found_methods = Vector{Method}()

    # We begin with the assumption nothing in Core.Compiler should be indexed
    done_modules = Set{Module}()
    done_functions = Set()

    function index(mod::Module)
        mod ∈ done_modules && return
        push!(done_modules, mod)
        @info "Indexing" mod=mod

        for name in names(mod; all=true)
            val = try
                String(name)[1]=='#' && continue
                Core.eval(mod, name)
            catch err
                @warn "while working out what things was, had error" thing=name mod=mod exception=err
            end

            val isa Core.Compiler.BitArray{1} && continue  # This type is unhashable and unprintable

            index(val)
        end
    end

    function index_func(f)
        try
            f ∈ done_functions && return
            push!(done_functions, f)
            append!(found_methods, methods(f))
        catch err
            @warn exception=err
        end
    end

    function index(x)
        index_func(x) 
    end

    function index(x::Function)
        index_func(x)
        index_func(typeof(x))  # all functions have a type
    end

    # Invoke it
    index(mainmod)  # first arg is dummy

    @info "indexing complete" length(found_methods)
    return found_methods
end


union_types(U) = getproperty.(U, propertynames(U))

get_module(x) = get_module(typeof(x))
get_module(meth::Method) = meth.module
get_module(U::UnionAll) = get_module(U.body)
get_module(T::DataType) = T.name.module
function get_module(U::Union)
    modules = unique(get_module, union_types)
    if length(modules)==1
        return first(modules)
    else
        error("There is not a common module for all types in $U")
    end
end


hunt(mod::Module=Main; kwargs...) = hunt(all_methods(mod), mod; kwargs...)
function hunt(methods::Union{Base.MethodList, AbstractVector{Method}}, mod::Module; exclude::Vector{<: Type}=[])
    function is_pirate(meth::Method)
        method_module = get_module(meth)

        method_module ≠ mod && return false # not interested in reports of piracy from outside module mod

        method_pkg = Base.PkgId(method_module)

        meth.sig isa UnionAll && return false # TODO: Deal with this properly
        length(meth.sig.parameters) < 2 && return false  # can't pirate without args
        
        any(meth.sig.parameters[1] <: T for T in exclude) && return false # don't trigger on manually-excluded type methods

        function_type = meth.sig.parameters[1]
        local function_module
        try
            function_module = get_module(function_type)
            function_pkg = Base.PkgId(function_module)
            Base.PkgId(function_module) == method_pkg && return false  # can't pirate own function

            # Base and Core can't pirate each other
            ((function_pkg == Base.PkgId(Base) &&  method_pkg == Base.PkgId(Core)) ||
            function_pkg == Base.PkgId(Core) &&  method_pkg == Base.PkgId(Base)) && return false
        catch err
            # Not actually a big issue, as this can not lead to a different result
            @warn "Failed to find function's module" function_type exception=err
            function_module = "The module that defined $(function_type)"
        end

        function is_own(arg_type::DataType)
            try
                arg_module = get_module(arg_type)
                arg_pkg = Base.PkgId(arg_module)
                arg_pkg == method_pkg && return true
            catch err
                @warn "Failed to find arg types 's module" arg_type exception=err
            end

            return any(is_own, arg_type.parameters)
        end
        
        is_own(U::Union) = all(is_own, union_types(U))
        is_own(::UnionAll) = return true  # TODO: Deal with this properly
        is_own(x) = is_own(typeof(x)) # for Int, Symbol and other bit types that show up in type-params

        arg_types_type = meth.sig.parameters[2:end]
        any(is_own, arg_types_type) && return false
        @error "Pirate Found!" meth pirate=method_module victim=function_module
        return true  # has failed to prove that they are not a pirate.
    end

    return filter(is_pirate, methods)
end

@testset "Type Piracy Check" begin
    # define excluded first-parameter types
    exclude = [
        # seems that @kwdef causes an issue?
        Type{MolecularGraphKernels.VizGraphKwargs}
    ]
    
    # run test
    @test hunt(MolecularGraphKernels; exclude=exclude) == []
end

end
