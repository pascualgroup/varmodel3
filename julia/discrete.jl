
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
    # host X locus X allele_id
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
        immunity = FILL(ImmunityCount(0), p.n_hosts, p.n_loci, 2 * p.n_alleles_per_locus_initial),
    )
end

function run_discrete(p::Params)
    s = DiscreteState(p)
    
    # TODO: use dt
    for t in 1:p.t_end
        println("stepping to t = $(t)")
        # Deaths
        do_rebirth(t, p, s)
        
    end
end

function do_rebirth(t, p, s)
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


