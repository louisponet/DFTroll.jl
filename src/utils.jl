cellcenter(crystal::Crystal) = Vec(ustrip.(sum(crystal.cell, dims=2) .|> a₀)...)
centered_index(n, mid) =  n - (n > mid) * 2mid

orthogonalize!(basis) = svd!(basis)
