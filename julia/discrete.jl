
const ARRAY = Array
const VECTOR = Vector
const FILL = Base.fill
const RAND = Base.rand

@with_kw struct DiscreteState
    t_birth::VECTOR{Float32}
    t_death::VECTOR{Float32}
    
    # Number of infections in each host
    n_infections::VECTOR{UInt8}
    
    # Time of each infection
    # host X infection
    t_infection::ARRAY{Float32, 2}
    
    # Genes that make up each infection, stored directly as
    # sequences of sequences of allele IDs
    # host X infection X locus
    infection_genes::ARRAY{UInt16, 3}
    
    # Current expression index of each infection
    # host X infection
    expression_index::ARRAY{ExpressionIndex, 2}
    
    # Number of alleles at each locus
    n_alleles::Vector{AlleleId}
    
    # Immunity level for every allele
    # host X allele_id X locus
    immunity::ARRAY{ImmunityCount, 3}
end

function DiscreteState(p::Params)
    lifetime = VECTOR([draw_host_lifetime(p) for i in 1:p.n_hosts])
    t_birth = -RAND(Float32, p.n_hosts) .* lifetime
    t_death = t_birth + lifetime
    
    n_infections = FILL(UInt8(0), p.n_hosts)
    t_infection = FILL(NaN32, p.n_hosts, p.max_infection_count)
    infection_genes = FILL(AlleleId(0), p.n_hosts, p.max_infection_count, p.n_loci)
    expression_index = FILL(ExpressionIndex(0), p.n_hosts, p.max_infection_count)
    
    infection_hosts = sample(1:p.n_hosts, p.n_initial_infections, replace = false)
    n_infections[infection_hosts] .= 1
    t_infection[infection_hosts] .= 0.0
    infection_genes[infection_hosts, 1, :] = rand(
        AlleleId(1):AlleleId(p.n_alleles_per_locus_initial), p.n_initial_infections * p.n_loci
    )
    expression_index[infection_hosts, 1] .= 1
    
    DiscreteState(
        t_birth = t_birth,
        t_death = t_death,
        
        n_infections = n_infections,
        t_infection = t_infection,
        infection_genes = infection_genes,
        
        expression_index = expression_index,
        
        n_alleles = fill(AlleleId(p.n_alleles_per_locus_initial), p.n_loci),
        immunity = FILL(ImmunityCount(0), p.n_hosts, 2 * p.n_alleles_per_locus_initial, p.n_loci),
    )
end

function run_discrete(p::Params)
    s = DiscreteState(p)
    
    # TODO: use dt
    for t in 1:p.t_end
        println("stepping to t = $(t)")
        
        # Death and rebirth
        do_rebirth!(t, p, s)
        
        # Loss of immunity
        do_immunity_loss!(t, p, s) 
        
        # Expression transition
        do_transition!(t, p, s)
    end
end

function do_rebirth!(t, p, s)
    println("do_rebirth!()")
    
    dead_hosts = findall(s.t_death .< t)
    println("n dead: $(length(dead_hosts))")
    
    s.t_birth[dead_hosts] = s.t_death[dead_hosts]
    s.t_death[dead_hosts] = s.t_birth[dead_hosts] + [draw_host_lifetime(p) for i in 1:length(dead_hosts)]
    s.n_infections[dead_hosts] .= 0
    s.t_infection[dead_hosts, :, :] .= NaN32
    s.infection_genes[dead_hosts, :, :] .= 0
    s.expression_index[dead_hosts, :] .= 0
    s.immunity[dead_hosts, :, :] .= 0
end

function do_immunity_loss!(t, p, s)
    # TODO: adjust for dt
    p_loss = 1 - exp(-p.immunity_loss_rate)
    
    # Lose immunity one locus at a time
    for locus in 1:p.n_loci
        n_immunities = p.n_hosts * s.n_alleles[locus]
        n_loss = rand(Binomial(n_immunities, p_loss))
        println("n_loss($(locus)) = $(n_loss)")
        
        # Uniformly randomly sample indices in one-dimensional array
        indices = sample(1:n_immunities, n_loss, replace = false)
        
        # Create view of immunities for this locus as a one-dimensional array
        immunity_view = reshape(
            (@view(s.immunity[:, 1:s.n_alleles[locus], locus])),
            n_immunities
        )
        
        # Decrement immunity (leaving zeros at zero)
        immunity_view[indices] = max(zeros(ImmunityCount, n_loss), immunity_view[indices] .- ImmunityCount(1))
    end
end

function do_transition!(t, p, s)
    println("do_transition!()")
    # TODO: adjust for dt
    p_transition = 1 - exp(-p.transition_rate)
    
    n_infections = p.n_hosts * p.max_infection_count
    
    n_trans_raw = rand(Binomial(n_infections, p_transition))
    println("n_trans_raw: $(n_trans_raw)")
    
    # Uniformly randomly sample indices in one-dimensional array
    # Includes nonexistent infections
    indices_raw = sample(1:n_infections, n_trans_raw, replace = false)
    
    # Get one-dimensional views on infection
    expression_index_view = reshape(s.expression_index, n_infections)
    t_infection_view = reshape(s.t_infection, n_infections)
    
    # Filter indices down to infections that exist and are past the liver stage
    indices = indices_raw[
        (expression_index_view[indices_raw] .> 0) .& (t_infection_view[indices_raw] .+ p.t_liver_stage .< t)
    ]
    n_trans = length(indices)
    println("n_trans = $(n_trans)")
    
    # Compute host and infection indices from sampled indices
    host_indices = ((indices .- 1) .% p.n_hosts) .+ 1
    inf_indices = ((indices .- 1) ./ p.n_hosts) .+ 1
    
    # Construct matrix, with rows corresponding to infections, and columns
    # corresponding to genes in the expression sequence, so that
    # row i, column j is true if the host of infection i is *not immune* to gene j
    # *and* the expression index is greater than the current index.
end
