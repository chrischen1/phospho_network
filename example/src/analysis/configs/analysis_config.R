# Set parameters
result_path       = '~/Documents/workspace/phospho_network/example/script_files/pred_analysis'
algos             = 'en'  # en for elastic net, rf for random forest
alphas            = seq(0,1,0.05) #required for elastic net
i_penalty         = T     # required for elastic net use different penalty based on heat diffusion?
ncore             = 1     # number of cores used
outerfold         = 10
innerfold         = 5
scale_method      = "0-1"   # 0-1 or "scale" 
directional       = T       # Used except pred_choose is flat. Should only upstream nodes in the pathway be considered?
pred_choose       = 'all'   # method of choose different predictor : hf: by heat diffusion,connect:all connected nodes, direct: direct nodes, flat: all nodes in network
k                 = 0.001   # used if pred_choose is hf, a parameter to define the extend of predictors by keep nodes receive more heat than k*(heat_response_CNV)
max_level         = Inf     # used if pred_choose is up or connect, max level consdered for predictor selection

#LSF setting
cluster = T
queue   = 'short'
time    = '2:00'
mem     = '18000'

# Set inputs:
rna_filename          = '~/Documents/workspace/phospho_network/example/script_files/rna_processed.csv'
cnv_filename          = '~/Documents/workspace/phospho_network/example/script_files/cnv_processed.csv'
mut_filename          = '~/Documents/workspace/phospho_network/example/script_files/mutation_matrix_misc.csv'
mis_mut_filename      = '~/Documents/workspace/phospho_network/example/script_files/mutation_missense.csv'
heat_influence_file   = '~/Documents/workspace/temp_files/heat_influence.csv' #used if pred_choose is hf
network_file          = '~/Documents/workspace/temp_files/network.csv'

# target value input:
mdata_filename        = '~/Documents/workspace/phospho_network/example/script_files/rppa_processed.csv'
