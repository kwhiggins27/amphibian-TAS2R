path_in="/lab/solexa_public/Weng/230404_WIGTC-NOVASEQ1A_AHWMGTDSX5/FASTQ/L319_"

forward="_R1_001.fastq.gz"
reverse="_R2_001.fastq.gz"

A1_1="511_S10_L004"
A1_2="512_S11_L004"
A1_3="513_S12_L004"
A1_4="514_S13_L004"
A1_5="515_S14_L004"
A1_6="516_S15_L004"
A1_7="517_S16_L004"

A2_1="518_S17_L004"
A2_2="519_S18_L004"
A2_3="520_S19_L004"
A2_4="521_S20_L004"
A2_5="522_S21_L004"
A2_6="523_S22_L004"
A2_7="524_S23_L004"

A3_1="525_S24_L004"
A3_2="526_S25_L004"
A3_3="527_S26_L004"
A3_4="528_S27_L004"
A3_5="529_S28_L004"
A3_6="530_S29_L004"
A3_7="531_S30_L004"

PT1_1="409_S1_L001"
PT1_2="410_S2_L001"
PT1_3="411_S3_L001"
PT1_4="412_S4_L001"
PT1_5="413_S5_L001"
PT1_6="414_S6_L001"
PT1_7="415_S7_L001"

PT2_1="416_S8_L001"
PT2_2="417_S9_L001"
PT2_3="418_S10_L001"
PT2_4="419_S11_L001"
PT2_5="420_S12_L001"
PT2_6="421_S13_L001"
PT2_7="422_S14_L001"

PT3_1="423_S15_L001"
PT3_2="424_S16_L001"
PT3_3="425_S17_L001"
PT3_4="426_S18_L001"
PT3_5="427_S19_L001"
PT3_6="428_S20_L001"
PT3_7="429_S21_L001"

X1_1="430_S22_L002"
X1_2="431_S23_L002"
X1_3="432_S24_L002"
X1_4="433_S25_L002"
X1_5="434_S26_L002"
X1_6="435_S27_L002"
X1_7="436_S28_L002"

X2_1="437_S29_L002"
X2_2="438_S30_L002"
X2_3="439_S31_L002"
X2_4="440_S32_L002"
X2_5="441_S33_L002"
X2_6="442_S34_L002"
X2_7="443_S35_L002"

X3_1="444_S36_L002"
X3_2="445_S37_L002"
X3_3="446_S38_L002"
X3_4="447_S39_L002"
X3_5="448_S40_L002"
X3_6="449_S41_L002"
X3_7="450_S42_L002"

C1_1="451_S43_L003"
C1_2="452_S44_L003"
C1_3="453_S45_L003"
C1_4="454_S46_L003"
C1_5="455_S47_L003"
C1_6="456_S48_L003"
C1_7="457_S49_L003"

C2_1="458_S50_L003"
C2_2="459_S51_L003"
C2_3="460_S52_L003"
C2_4="461_S53_L003"
C2_5="462_S54_L003"
C2_6="463_S55_L003"
C2_7="464_S56_L003"

C3_1="465_S57_L003"
C3_2="466_S58_L003"
C3_3="467_S59_L003"
C3_4="468_S60_L003"
C3_5="469_S61_L003"
C3_6="470_S62_L003"
C3_7="471_S63_L003"

B1_1="472_S64_L004"
B1_2="473_S65_L004"
B1_3="474_S66_L004"
B1_4="475_S67_L004"
B1_5="476_S68_L004"
B1_6="477_S69_L004"
B1_7="478_S70_L004"

B2_1="479_S71_L004"
B2_2="480_S72_L004"
B2_3="481_S73_L004"
B2_4="482_S74_L004"
B2_5="483_S75_L004"
B2_6="484_S76_L004"
B2_7="485_S77_L004"

B3_1="486_S78_L004"
B3_2="487_S79_L004"
B3_3="488_S80_L004"
B3_4="489_S81_L004"
B3_5="490_S82_L004"
B3_6="491_S83_L004"
B3_7="492_S84_L004"



for var in "$A1_1" "$A1_2" "$A1_3" "$A1_4" "$A1_5" "$A1_6" "$A1_7" "$A2_1" "$A2_2" "$A2_3" "$A2_4" "$A2_5" "$A2_6" "$A2_7" "$A3_1" "$A3_2" "$A3_3" "$A3_4" "$A3_5" "$A3_6" "$A3_7"
do
  bsub -n 5 -R "rusage[mem=100000]" ./axolotl_trimmomatic_pt2.bat "$var"
done

###
#"$A1_1" "$A1_2" "$A1_3" "$A1_4" "$A1_5" "$A1_6" "$A1_7" "$A2_1" "$A2_2" "$A2_3" "$A2_4" "$A2_5" "$A2_6" "$A2_7" "$A3_1" "$A3_2" "$A3_3" "$A3_4" "$A3_5" "$A3_6" "$A3_7"
#"$PT1_1" "$PT1_2" "$PT1_3" "$PT1_4" "$PT1_5" "$PT1_6" "$PT1_7" "$PT2_1" "$PT2_2" "$PT2_3" "$PT2_4" "$PT2_5" "$PT2_6" "$PT2_7" "$PT3_1" "$PT3_2" "$PT3_3" "$PT3_4" "$PT3_5" "$PT3_6" "$PT3_7"
#"$X1_1" "$X1_2" "$X1_3" "$X1_4" "$X1_5" "$X1_6" "$X1_7" "$X2_1" "$X2_2" "$X2_3" "$X2_4" "$X2_5" "$X2_6" "$X2_7" "$X3_1" "$X3_2" "$X3_3" "$X3_4" "$X3_5" "$X3_6" "$X3_7"
#"$C1_1" "$C1_2" "$C1_3" "$C1_4" "$C1_5" "$C1_6" "$C1_7" "$C2_1" "$C2_2" "$C2_3" "$C2_4" "$C2_5" "$C2_6" "$C2_7" "$C3_1" "$C3_2" "$C3_3" "$C3_4" "$C3_5" "$C3_6" "$C3_7"
#"$B1_1" "$B1_2" "$B1_3" "$B1_4" "$B1_5" "$B1_6" "$B1_7" "$B2_1" "$B2_2" "$B2_3" "$B2_4" "$B2_5" "$B2_6" "$B2_7" "$B3_1" "$B3_2" "$B3_3" "$B3_4" "$B3_5" "$B3_6" "$B3_7"

#"$frog3_1" "$frog3_2" "$frog3_3" "$frog3_4" "$frog3_5" "$frog6_1" "$frog6_2" "$frog6_3" "$frog6_4" "$frog6_5" "$frog7_1" "$frog7_2" "$frog7_3" "$frog7_4" "$frog7_5"
