

using StaticArrays

A = SMatrix{2,2}(1.0, 2.0, 3.0, 4.0)

function test_alias(A)
    B = SMatrix{2,2}(1.0, 1.0, 1.0, 1.0)
    A += B * A
    return (A)
end


aa = test_alias(A)
