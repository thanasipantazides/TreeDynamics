using TreeDynamics
using LinearAlgebra
using DifferentialEquations
using DiffEqCallbacks
using GLMakie

function main()
    n = 10
    C_Oi = zeros(3, 3, n)
    I_i = zeros(3, 3, n)
    w_i = zeros(3, n)
    l = zeros(n)
    k = zeros(n)
    c = zeros(n)
    a_i = zeros(3, n)
    a = [0; 0; 1]

    for i = 1:n
        # ang = pi / 2 * ((i - 1) / n)
        # ang = rand()*pi/3
        ang = 0
        C_Oi[:, :, i] = [
            1 0 0;
            0 cos(ang) -sin(ang);
            0 sin(ang) cos(ang)
        ]
        # C_Oi[:, :, i] = randr()

        I_i[:, :, i] = [
            1 0 0;
            0 1 0;
            0 0 4*i/n
        ] * 1e-3
        # w_i[:, i] = [1.0; 0.5; 0.1] * 10
        w_i[:,i] = rand(3)*10

        # k[i] = 2
        # k[i] = 2*(i - 1)/n
        k[i] = 1*(n - i)/n
        c[i] = 1e-2
        l[i] = 2*(n - i) / n
        a_i[:, i] = C_Oi[:, :, i] * a
    end

    # set some initial values
    p = TreeLinkage(
        n,
        l,
        k,
        c,
        I_i,
        a_i
    )

    x0 = [vec(w_i); vec(C_Oi)]

    # op = DiffEqArrayOperator(ones(length(x0), 1), update_func=dynamics!)

    # set a tspan
    tspan = (0, 10.0)

    cb = ManifoldProjection(proj_so3!, autodiff=AutoForwardDiff())

    # define a problem = ODEProblem(dynamics, initial_value, tspan)
    problem = ODEProblem(dynamics!, x0, tspan, p)
    # problem = ODEProblem(op, x0, tspan, p)

    println("Solving...")
    # sol = solve(problem, RK4(), callback=cb, dt=0.001, saveat=0.01)

    sol = solve(
        problem,
        Tsit5(),
        # callback=cb,
        saveat=0.01,
        reltol=1e-3,
        abstol=1e-3,
        maxiters=1e8,
        verbose=true
    )

    # sol = solve(
    #     problem,
    #     Rosenbrock23(autodiff = AutoFiniteDiff()),
    #     callback=cb,
    #     saveat=0.01,
    #     reltol = 1e-3,
    #     abstol = 1e-3,
    #     maxiters=1e8
    # )
    # sol = solve(
    #     problem,
    #     RK4(),
    #     callback=cb,
    #     saveat=0.01,
    #     reltol = 1e-3,
    #     abstol = 1e-3,
    #     maxiters=1e8
    # )

    # @show sol

    println("Reshaping...")
    # digest the result
    wh = sol[1:3*p.n, :]
    C_Oih = reshape(sol[3*p.n+1:end, :], 3, 3, p.n, length(sol.t))

    println(size(sol.t))

    println("Plotting...")
    # plot result
    GLMakie.activate!()
    fig = Figure()
    display(fig)

    time_slide(fig, sol, p)

end
