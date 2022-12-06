Having recently come back to a Julia project after several months, I've just gone through a few too many missed guesses in a row, trying to initialise an array with an array of zeros with the same structure as an existing array. My thinking process went as follows,

```
# this is the array I'm starting from
A = SMatrix{3,5}(repeat([0.0+1im], 15))
# now I want to create a new array with same shape and type as A, but initialised at 0, i.e. equivalent to
B = 0*A
```

```
zeros(A) # a bit naive but nope, but I see there's a method zeros([T=Float64,] dims::Tuple)
dims(A) # nope
dim(A) # nope, that's in R
size(A) # aha, like Matlab
zeros(size(A))
```

but now I need the right type,

```
typeof(P) # hmm, doesn't look like it
# SMatrix{3, 5, ComplexF64, 15} (alias for SArray{Tuple{3, 5}, Complex{Float64}, 2, 15})
# ... after some time googling, it's called "eltype"

zeros(eltype(A), size(A)) # hurrah
```

I'm wondering if it might be helpful for beginners (or, in fact, for people who often switch between languages) to have a dummy package defining a few of those common missed guesses, such as:

```
dims() # alias for size(), or message suggesting to use size() if size() doesn't exist in the environment
length() # same for len()
type() # same for eltype()

```

