#!/bin/bash
#SBATCH --job-name=RM2   # Job name
#SBATCH --mem=200gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=10              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/RM2_%j.log   # Standard output and error log
#SBATCH --error=logs/RM2_%j.err

cd ../../results/coordinate_analysis/repeat/mini_run

## For mini example
# nohup RepeatModeler -database ../../../../results/coordinate_analysis/repeat/references/saltmarsh
#     -threads 10 -LTRStruct >& saltmarsh_run.out &
nohup RepeatModeler -database ../../../../results/coordinate_analysis/repeat/references/flamingo
    -threads 10 -LTRStruct >& flamingo_run.out &


## For manuscript: amphibians
# nohup RepeatModeler -database references/tropical_clawed_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/African_clawed_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/frogs_and_toads
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Leishan_spiny_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Asiatic_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/African_clawed_frog2
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Yunnan_mustache_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Congo_dwarf_clawed_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Tungara_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Puerto_Rican_coqui
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Kenyan_clawed_frog
#       -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/plains_spadefoot_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Puerto_Rican_coqui
#       -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/hourglass_treefrog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/eastern_narrow-mouthed_toad
#       -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/wood_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Sardinian_treefrog
#       -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/plateau_brown_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/two-lined_caecilian
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/caecilians
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/caecilians2
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/common_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/common_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/axolotl
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Iberian_ribbed_newt
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/fire-bellied_toad
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/corroboree_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/mountain_yellow-legged_frog
#    -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/painted_frog
#    -threads 10 -LTRStruct >& run.out &


## For manuscript: non-amphibians
# nohup RepeatModeler -database references/diamondback_terrapin
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Aeolian_wall_lizard
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Far_Eastern_brook_lamprey
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/bony_fishes2
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/whale_shark
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Pacific_lamprey
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/blackcap
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Shaws_sea_snake
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/bony_fishes2
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Common_starling
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/gray_bichir
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/prairie_rattlesnake
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/reedfish
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Patagonian_moray_cod
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/smaller_spotted_catshark
#  -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/thorny_skate
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/lion
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/yellow-throated_sandgrouse
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/whitespotted_bambooshark
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/North_American_porcupine
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/sea_lamprey
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/epaulette_shark
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/West_African_lungfish
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Malagasy_flying_fox
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/greater_Indian_rhinoceros
# -threads 10 -LTRStruct >& run.out &
# nohup RepeatModeler -database references/Boesemans_rainbowfish
# -threads 10 -LTRStruct >& run.out &
