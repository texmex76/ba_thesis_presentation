using QuantumOptics

k = 2π; ωr = 1
η = 10ωr; ωc = -10ωr; U0 = -1ωr

b_position = PositionBasis(0, 1, 32)
b_fock = FockBasis(16)
p = momentum(b_position)
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)

potential = x -> U0*cos(k*x)^2
H_int = (one(b_fock) ⊗ potentialoperator(b_position, potential))*ad*a
H_kin = (one(b_fock) ⊗ p^2) / k^2
H_cavity = -ωc*ad*a
H_pump = η*(a + ad)
H = H_kin + dense(H_int) + H_cavity + H_pump

E, ψ_states = eigenstates((H + dagger(H))/2, 3)
