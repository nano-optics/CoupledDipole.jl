function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end

function mymodel(;x = range(0, 10, length = 100), a = 1, b = 1, c = 2, fun = sin, kws...)
    DataFrame(x = x, y = a .* sin.(b .* x) .+ c, z = a .* fun.(b .* x) .+ c)
end

params = expand_grid(a=[0.1,0.2,0.3], b = [1,2,3], c = [0,0.5]);

all = map(eachrow(params)) do r
    mymodel(;fun = tanh, pairs(r)...)
end;


function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end

function mymodel(;x = range(0, 10, length = 100), a = 1, b = 1, c = 2, fun = sin, kws...)
    DataFrame(x = x, y = a .* sin.(b .* x) .+ c, z = a .* fun.(b .* x) .+ c)
end

function pmap_df(p, f, kws...; join=true)
   tmp = map(f, eachrow(p), kws...)
   all = reduce(vcat, tmp, source="id")
   if !join
    return all
   end
   p[!, :id] = 1:nrow(p)
   return DataFrames.leftjoin(p, all, on=:id)
end

function pmap_df2(p, f, kws...; join=true)
    all = mapreduce(f, (x,y) -> vcat(x,y, source="id", cols=:intersect), eachrow(p), kws...)
    if !join
     return all
    end
    p[!, :id] = 1:nrow(p)
    return DataFrames.leftjoin(p, all, on=:id)
 end

 
params = expand_grid(a=[0.1,0.2,0.3], b = [1,2,3], c = [0,0.5]);

pmap_df2(params, p -> mymodel(;p..., fun=tanh))





all = map(p -> mymodel(;p...), eachrow(params))
aa = reduce(vcat, all, source="id")

aha = mapreduce(p -> mymodel(;p...), vcat, eachrow(params))

aha




params[!, :id] = 1:nrow(params)

all_df = vcat(all...)
all_df[!, :id] = repeat(1:nrow(params), inner=100)

DataFrames.leftjoin(params, all_df, on=:id)


param1 = [1,2,3]; param2 = ["a", "b"];
df = rename!(DataFrame(Iterators.product(param1, param2)), ["param1", "param2"])
df[!, :result] .= 0.0; df

for r ∈ eachrow(df)
    r.result = mymodel(r.param1, r.param2)
end