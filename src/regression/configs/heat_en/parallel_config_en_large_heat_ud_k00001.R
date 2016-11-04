# Set parameters
result_path       = '~/Documents/workspace/phospho_network/script_files/output/en_large_heat_ud_k00001'
algos             = 'en'  # en for elastic net, rf for random forest
alphas            = seq(0,1,0.05) #required for elastic net
i_penalty         = T     # required for elastic net use different penalty based on heat diffusion?
ncore             = 1     # number of cores used
outerfold         = 10
innerfold         = 5
scale_method      = "0-1"   # 0-1 or "scale" 
directional       = F       # Used except pred_choose is flat. Should only upstream nodes in the pathway be considered?
pred_choose       = 'hf'    # method of choose different predictor : hf: by heat diffusion,connect:all connected nodes, direct: direct nodes, flat: all nodes in network
k                 = 0.00001   # used if pred_choose is hf, a parameter to define the extend of predictors by keep nodes receive more heat than k*(heat_response_CNV)
max_level         = Inf     # used if pred_choose is up or connect, max level consdered for predictor selection

#LSF setting
cluster = T
queue   = 'short'
time    = '2:00'
mem     = '38000'

# Set inputs:
rna_filename          = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/total_protein_processed.csv'
cnv_filename          = ''
mut_filename          = ''

mis_mut_filename      = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/mutation_matrix.csv'
heat_influence_file   = '~/Documents/workspace/temp_files/heat_matrix_ms.csv' #used if pred_choose is hf
network_file          = '~/Documents/workspace/temp_files/network_ms.csv'

# target value input:
mdata_filename        = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/msdata_processed.csv'
