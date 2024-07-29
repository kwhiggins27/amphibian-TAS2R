#!/bin/bash
#SBATCH --job-name=RM   # Job name
#SBATCH --mem=100gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/RM_%j.log   # Standard output and error log
#SBATCH --error=logs/RM_%j.err


# BuildDatabase -name references/axolotl /lab/wengpj01/axolotl/axolotl_mygenome.fasta

#BuildDatabase -name references/bullfrog /lab/wengpj01/bullfrog/Lithobates_catesbeianus_bullfrog_2022.fasta
#
# BuildDatabase -name references/cane /lab/wengpj01/cane/canetoad.v2.2.fasta

#BuildDatabase -name references/tropical_clawed_frog	../../genomes/GCA_000004195.4_UCB_Xtro_10.0_genomic.fna
# BuildDatabase -name references/African_clawed_frog	../../genomes/GCA_001663975.1_Xenopus_laevis_v2_genomic.fna
# BuildDatabase -name references/frogs_and_toads	../../genomes/GCA_004786255.1_Pads_1.0_genomic.fna
# BuildDatabase -name references/Leishan_spiny_toad	../../genomes/GCA_009667805.1_ASM966780v1_genomic.fna
# BuildDatabase -name references/tropical_clawed_frog2	../../genomes/GCA_013368275.1_ASM1336827v1_genomic.fna
# BuildDatabase -name references/Asiatic_toad	../../genomes/GCA_014858855.1_ASM1485885v1_genomic.fna
# BuildDatabase -name references/African_clawed_frog2	../../genomes/GCA_017654675.1_Xenopus_laevis_v10.1_genomic.fna
# BuildDatabase -name references/Yunnan_mustache_toad	../../genomes/GCA_018994145.1_ASM1899414v1_genomic.fna
# BuildDatabase -name references/Congo_dwarf_clawed_frog	../../genomes/GCA_019447015.1_UCB_Hboe_1.0_genomic.fna
# BuildDatabase -name references/Tungara_frog	../../genomes/GCA_019512145.1_UCB_Epus_1.0_genomic.fna
# BuildDatabase -name references/Puerto_Rican_coqui	../../genomes/GCA_019857665.1_UCB_Ecoq_1.0_genomic.fna
# BuildDatabase -name references/Kenyan_clawed_frog	../../genomes/GCA_024363595.1_UCB_Xborealis_1_genomic.fna
# BuildDatabase -name references/plains_spadefoot_toad	../../genomes/GCA_027358695.2_aSpeBom1.2.pri_genomic.fna
# BuildDatabase -name references/painted_frog	../../genomes/GCA_027410445.1_aDisPic1.pri_genomic.fna
## BuildDatabase -name references/hourglass_treefrog	../../genomes/GCA_027789765.1_aDenEbr1.pat_genomic.fna
## BuildDatabase -name references/eastern_narrow-mouthed_toad	../../genomes/GCA_027917415.1_aGasCar1.hap2_genomic.fna
## BuildDatabase -name references/wood_frog	../../genomes/GCA_028564925.1_aRanSyl1.merge_genomic.fna
## BuildDatabase -name references/Sardinian_treefrog	../../genomes/GCA_029499605.1_aHylSar1.hap1_genomic.fna

# BuildDatabase -name references/plateau_brown_frog	../../genomes/GCA_029574335.1_ASM2957433v1_genomic.fna
# BuildDatabase -name references/two-lined_caecilian	../../genomes/GCA_901001135.2_aRhiBiv1.2_genomic.fna
# BuildDatabase -name references/caecilians	../../genomes/GCA_901765095.2_aMicUni1.2_genomic.fna
# BuildDatabase -name references/caecilians2	../../genomes/GCA_902459505.2_aGeoSer1.2_genomic.fna
# BuildDatabase -name references/common_toad	../../genomes/GCA_905171765.1_aBufBuf1.1_genomic.fna
# BuildDatabase -name references/common_frog	../../genomes/GCA_905171775.1_aRanTem1.1_genomic.fna


# #flag because big genome
# BuildDatabase -name references/axolotl	../../genomes/GCA_002915635.3_KH.fasta
# BuildDatabase -name references/Iberian_ribbed_newt	../../genomes/GCA_026652325.1_KH.fasta
# BuildDatabase -name references/fire-bellied_toad	../../genomes/GCA_027579735.1_KH.fasta
# BuildDatabase -name references/corroboree_frog	../../genomes/GCA_028390025.1_KH.fasta
# BuildDatabase -name references/mountain_yellow-legged_frog	../../genomes/GCA_029206835.1_KH.fasta

# BuildDatabase -name references/diamondback_terrapin ../../genomes/GCA_027887155.1_rMalTer1.hap1_genomic.fna
# BuildDatabase -name references/Aeolian_wall_lizard ../../genomes/GCA_027172205.1_rPodRaf1.pri_genomic.fna
# BuildDatabase -name references/Far_Eastern_brook_lamprey ../../genomes/GCA_015708825.1_ASM1570882v1_genomic.fna
#BuildDatabase -name references/bony_fishes3 ../../genomes/GCA_023029165.1_ASM2302916v1_genomic.fna
# BuildDatabase -name references/whale_shark ../../genomes/GCA_021869965.1_sRhiTyp1.1_genomic.fna
# BuildDatabase -name references/Pacific_lamprey ../../genomes/GCA_014621495.2_ETRf_v1_genomic.fna
# BuildDatabase -name references/blackcap ../../genomes/GCA_009819655.1_bSylAtr1.pri_genomic.fna
# BuildDatabase -name references/Shaws_sea_snake ../../genomes/GCA_019472885.1_HCur_v2_genomic.fna
#BuildDatabase -name references/bony_fishes2 ../../genomes/GCA_028583405.1_KU_S6_genomic.fna
# BuildDatabase -name references/Common_starling ../../genomes/GCA_023376015.1_Svulgaris_vAU_1.1_genomic.fna
# BuildDatabase -name references/gray_bichir ../../genomes/GCA_016835505.1_ASM1683550v1_genomic.fna
# BuildDatabase -name references/prairie_rattlesnake ../../genomes/GCA_003400415.2_UTA_CroVir_3.0_genomic.fna
# BuildDatabase -name references/reedfish ../../genomes/GCA_900747795.4_fErpCal1.3_genomic.fna
# BuildDatabase -name references/Patagonian_moray_cod ../../genomes/GCA_027704905.1_KU_S4_genomic.fna
# BuildDatabase -name references/thorny_skate ../../genomes/GCA_010909765.2_sAmbRad1.1.pri_genomic.fna
# BuildDatabase -name references/lion ../../genomes/GCA_018350215.1_P.leo_Ple1_pat1.1_genomic.fna
# BuildDatabase -name references/yellow-throated_sandgrouse ../../genomes/GCA_009769525.1_bPteGut1.pri_genomic.fna
# BuildDatabase -name references/whitespotted_bambooshark ../../genomes/GCA_004010195.1_ASM401019v1_genomic.fna
# BuildDatabase -name references/North_American_porcupine ../../genomes/GCA_028451465.1_mEreDor1.pri_genomic.fna
# BuildDatabase -name references/sea_lamprey ../../genomes/GCA_010993605.1_kPetMar1.pri_genomic.fna
# BuildDatabase -name references/epaulette_shark ../../genomes/GCA_020745735.1_sHemOce1.pat.X.cur._genomic.fna
# BuildDatabase -name references/West_African_lungfish ../../genomes/GCA_019279795.1_KH.fasta
# BuildDatabase -name references/Malagasy_flying_fox ../../genomes/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna
# BuildDatabase -name references/greater_Indian_rhinoceros ../../genomes/GCA_028646465.1_Rhinoceros_unicornis_HiC_genomic.fna
# BuildDatabase -name references/Boesemans_rainbowfish ../../genomes/GCA_017639745.1_fMelBoe1.pri_genomic.fna
# BuildDatabase -name references/rock_pigeon ../../genomes/GCA_001887795.1_colLiv2_genomic.fna


#
#BuildDatabase -name references/terribilis /lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/terribilis/P.terribilis.gapclosed.fasta
#
#BuildDatabase -name references/xenopus /lab/wengpj01/xenopus/Xenopus_Tropicalis_2022.fasta