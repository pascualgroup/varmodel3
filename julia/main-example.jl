#!/usr/bin/env julia

include("parameters.jl")
include("output.jl")
include("model.jl")

function main()
    params = init_params()
    run(params)
end

function init_params()
    t_year = 360
    
    Params(
        implementation = DISCRETE_APPROXIMATION,
        dt = 1,
        
        output_db_filename = "output.sqlite",
        
        output_hosts = true,
        output_strains = false,
        output_genes = false,
    
        host_sampling_period = 30,
        host_sample_size = 100,
        expected_equilibrium = 10800,
    
        verification_on = true,
        verification_period = 30,
    
        rng_seed = nothing,
    
        t_year = t_year,
        t_end = 111 * t_year,
        t_burnin = 10 * t_year,
    
        n_hosts = 10000,
        n_initial_infections = 20,
        
        n_genes_initial = 9600,
        n_genes_per_strain = 60,
        
        n_loci = 2,
    
        n_alleles_per_locus_initial = 960,
    
        transmissibility = 0.5,
        coinfection_reduces_transmission = true,
    
        ectopic_recombination_rate = 1.8e-7,
        
        max_immunity_count = 100,
        immunity_loss_rate = 0.001,
    
        mutation_rate = 1.42e-8,
    
        t_liver_stage = 14.0,
        
        transition_rate = 1/6.0,
        
        moi_transmission_max = 9,
    
        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,
    
        immigration_on = true,
        immigration_rate_fraction = 0.0026,
        
        max_infection_count = 9,
        
        biting_rate_mean = 0.0005,
        daily_biting_rate_distribution = [
            727.4467846151298,
            727.4393612450948,
            727.4325889411169,
            727.4264099549741,
            727.4207697919934,
            727.415621378698,
            727.4109225065789,
            727.4066341838655,
            727.4027195424633,
            727.3991461557988,
            727.3958845491045,
            727.3929075701011,
            727.3901901476318,
            727.3877069206706,
            727.3854404373049,
            727.3833714190398,
            727.3814818628291,
            727.3797591633929,
            727.3781858290638,
            727.3767513027294,
            727.3754417791339,
            727.3742470060482,
            727.373156264827,
            727.3721608435912,
            727.3712522628946,
            727.3704224597545,
            727.3696654787001,
            727.3689742476798,
            727.3683430147117,
            727.3677674893947,
            727.3672417211801,
            727.3667619336727,
            727.3663244510045,
            727.3659248043872,
            727.3655603534133,
            727.3652275597987,
            727.3649234622708,
            727.3646461776384,
            727.3643934797981,
            727.3641622356696,
            727.3639511080177,
            727.3637589629986,
            727.363583454739,
            727.3634229259393,
            727.3632766441854,
            727.3631435260382,
            727.3630216824073,
            727.3629104575629,
            727.3628092897866,
            727.3627168895325,
            727.3626324156965,
            727.3625554462344,
            727.3624853029207,
            727.3624209385172,
            727.3623620610038,
            727.3623084222716,
            727.3622597741327,
            727.3622152974756,
            727.3621743424715,
            727.3621368742284,
            727.3621028578549,
            727.3620722015115,
            727.3620440410089,
            727.3620181121302,
            727.3619944370719,
            727.3619730380311,
            727.3619537775236,
            727.3619360329735,
            727.3619197415738,
            727.3619049060322,
            727.3618915290563,
            727.3618794111945,
            727.3618682568307,
            727.3618580478559,
            727.3618487682489,
            727.3618403946433,
            727.3618326855513,
            727.361825532844,
            727.3618189398356,
            727.3618129098402,
            727.3618074461718,
            727.3618025521449,
            727.3617982085909,
            727.3617941316172,
            727.3617902683498,
            727.3617866622487,
            727.3617833567746,
            727.3617803953877,
            727.3617778215483,
            727.3617756339031,
            727.3617735479852,
            727.3617715394117,
            727.3617696491916,
            727.3617679183337,
            727.361766387847,
            727.36176509874,
            727.3617640433582,
            727.3617630354819,
            727.3617620686814,
            727.3617611622192,
            727.3617603353581,
            727.3617596073606,
            727.3617589974891,
            727.3617584944656,
            727.3617580216186,
            727.3617575776527,
            727.361757166678,
            727.3617567928038,
            727.3617564601402,
            727.3617561727967,
            727.3617559236015,
            727.3617556925475,
            727.3617554785384,
            727.3617552809752,
            727.3617550992586,
            727.3617549327897,
            727.3617547809694,
            727.3617546431983,
            727.3617545188777,
            727.3617544074082,
            727.3617543081907,
            727.3617542206264,
            727.3617541449555,
            727.3617540835293,
            727.3617540346811,
            727.3617539962931,
            727.3617539662479,
            727.3617539424276,
            727.3617539227149,
            727.3617539049919,
            727.3617538871412,
            727.3617538670451,
            727.3617538425858,
            727.361753811646,
            727.3617537738965,
            727.3617537474546,
            727.3617537348399,
            727.3617537325969,
            727.3617537372701,
            727.3617537454039,
            727.3617537535431,
            727.3617537582319,
            727.361753756015,
            727.3617537434369,
            727.361753717042,
            727.3724037471579,
            727.5268077710775,
            728.0864145178148,
            729.3251572671344,
            731.4968073487536,
            734.8304504958904,
            739.530883006113,
            745.7798470548471,
            753.7374003255251,
            763.5432667797737,
            775.3181227037775,
            789.1648185537944,
            805.1694991609678,
            823.4026658605412,
            843.9201332226021,
            866.7638146467218,
            891.9625913159003,
            919.5328552491997,
            949.4791770721879,
            981.7947984292568,
            1016.462047063863,
            1053.4527878677054,
            1092.7287933154282,
            1134.2420804962405,
            1177.9352329712056,
            1223.7417686972058,
            1271.586488333063,
            1321.385823576554,
            1373.0481899965648,
            1426.474377689395,
            1481.557989614614,
            1538.1858350060556,
            1596.2383783398354,
            1655.5901954252295,
            1716.1104689905865,
            1777.6633600655653,
            1840.1086555595211,
            1903.302193301744,
            1967.0963367760687,
            2031.3404861338572,
            2095.881624850742,
            2160.5648137124863,
            2225.2336788899893,
            2289.730905359386,
            2353.898751476242,
            2417.579539346419,
            2480.6161210614,
            2542.852350880766,
            2604.133578830221,
            2664.307094181794,
            2723.222554156316,
            2780.732403841001,
            2836.692330896516,
            2890.961698608162,
            2943.403895355079,
            2993.8866995868498,
            3042.2826437464605,
            3088.469441815163,
            3132.3302955257896,
            3173.7541627318046,
            3212.636047611827,
            3248.877313010183,
            3282.386031099976,
            3313.0771415827603,
            3340.872649792165,
            3365.701830423097,
            3387.5014962180453,
            3406.2161935160534,
            3421.7982492523597,
            3434.207886830968,
            3443.4133447826616,
            3449.3910785803023,
            3452.125762536815,
            3451.6102642087735,
            3447.8456647148546,
            3440.8413064449724,
            3430.614866591087,
            3417.1921989429648,
            3400.607244367431,
            3380.901954494162,
            3358.126317298845,
            3332.3381011248593,
            3303.6027548742527,
            3271.9932376940906,
            3237.5898544258594,
            3200.4800352229054,
            3160.7581075251314,
            3118.525101007804,
            3073.888456226159,
            3026.9617486713532,
            2977.864468113825,
            2926.7216347446474,
            2873.6635404970575,
            2818.8254081437626,
            2762.3470533221125,
            2704.3725343964593,
            2645.049762708121,
            2584.530256205053,
            2522.9685548859484,
            2460.521950431697,
            2397.3501673424435,
            2333.614816724731,
            2269.47903106466,
            2205.107115354262,
            2140.664068947113,
            2076.3151717467226,
            2012.225597986068,
            1948.5599533289585,
            1885.481864506546,
            1823.1535527131925,
            1761.735398168373,
            1701.3855354260074,
            1642.2593872771913,
            1584.5092599210711,
            1528.2839005286014,
            1473.7280654538542,
            1420.9820941050098,
            1370.18148016111,
            1321.4564522395046,
            1274.93150987925,
            1230.725080533428,
            1188.9490371510747,
            1149.7083483954655,
            1113.1007404840648,
            1079.216405923785,
            1048.1252712769428,
            1019.7598018279631,
            993.9165197467398,
            970.3751943838186,
            948.9295386193717,
            929.391123488377,
            911.5883614896892,
            895.3653827807157,
            880.5805165213927,
            867.1049429031523,
            854.8216122001674,
            843.6239817952985,
            833.4152440450044,
            824.1073126248505,
            815.6200592988989,
            807.8805489013828,
            800.8224653018653,
            794.3853439908478,
            788.5142576870527,
            783.1591154189,
            778.2742934246289,
            773.8182842757668,
            769.753276847351,
            766.0447948383232,
            762.6614356572431,
            759.5745825082275,
            756.7581694336893,
            754.1884141329898,
            751.8436450213864,
            749.7041084506403,
            747.7517985465525,
            745.9702854162372,
            744.3446028339781,
            742.8610735238437,
            741.5072744901724,
            740.2718071194504,
            739.1443523017706,
            738.1154192275892,
            737.1763930198935,
            736.319424556843,
            735.5373046329322,
            734.8235126623326,
            734.1720784527512,
            733.5775363593688,
            733.0349177700982,
            732.5396908775579,
            732.0877100064706,
            731.6751869778658,
            731.2986745615367,
            730.9550375394026,
            730.64140500352,
            730.3551523764568,
            730.093887822083,
            729.8554303302728,
            729.6377873179922,
            729.4391415050799,
            729.2578346266531,
            729.0923533227427,
            728.9413161502168,
            728.8034620342182,
            728.6776402960493,
            728.5628002492863,
            728.4579825832777,
            728.3623129830245,
            728.2749920202244,
            728.1952916729266,
            728.1225478156992,
            728.0561530422575,
            727.9955529353682,
            727.9402412643857,
            727.8897587924921,
            727.8436829142603,
            727.8016300652412,
            727.7632454843513,
            727.7282102603887,
            727.696231842034,
            727.6670430612044,
            727.6404017866183,
            727.6160858564442,
            727.5938920325976,
            727.5736361061961,
            727.5551468559204,
            727.5382721502989,
            727.5228686172751,
            727.508809993699,
            727.495977180301,
            727.4842648045143,
            727.4735740153599,
            727.4638162084627,
            727.4549097979058,
        ],
    )
end

main()
