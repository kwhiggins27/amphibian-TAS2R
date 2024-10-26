#!shell
date
version

#specify data file, p-value threshold, # of threads to use, and log file
load -i birds_tas2r_numbers.fixedlabels.txt -t 12 -l birds.txt

#the phylogenetic tree structure with branch lengths
tree ((Dromaius-novaehollandiae:65.08832,Rhea-pennata:65.08832):42.5035,((((Tauraco-erythrolophus:74.96995,(((((Hemiprocne-comata:36.48226,Apus-apus:36.48226):16.06447,(Archilochus-colubris:4.19549,Calypte-anna:4.19549):48.35124):13.62058,Podargus-strigoides:66.16731):8.38664,Nyctibius-grandis:74.55395):1e-05,Caprimulgus-europaeus:74.55396):0.41599):1.296,(((((((Theristicus-caerulescens:64.34646,Ciconia-maguari:64.34646):1.06504,Phalacrocorax-aristotelis:65.4115):3.2747,Rhynochetos-jubatus:68.6862):1.21658,(Porphyrio-hochstetteri:60.84535,Grus-americana:60.84535):9.05743):2.05146,((Alca-torda:40.06087,(Sterna-hirundo:23.04832,Rissa-tridactyla:23.04832):17.01255):30.45868,Pluvialis-apricaria:70.51955):1.43469):0.61883,((Pterocles-gutturalis:63.49622,(Streptopelia-turtur:22.0472,Columba-livia:22.0472):41.44902):8.94862,Phoenicopterus-ruber:72.44484):0.12823):2.53742,((((((Merops-nubicus:62.73275,(((Colaptes-auratus:18.35185,Dryobates-pubescens:18.35185):17.3217,Indicator-indicator:35.67355):11.22986,Pogoniulus-pusillus:46.90341):15.82934):1.07822,Bucorvus-abyssinicus:63.81097):0.0598,Trogon-surrucura:63.87077):2.05072,Colius-striatus:65.92149):4.50672,((Gypaetus-barbatus:40.79859,((Haliaeetus-albicilla:27.95419,Accipiter-gentilis:27.95419):0.52198,Aquila-chrysaetos:28.47617):12.32242):19.8997,Gymnogyps-californianus:60.69829):9.72992):1e-05,(((((Melopsittacus-undulatus:32.36927,(Ara-ararauna:26.40085,(Amazona-aestiva:24.07135,Myiopsitta-monachus:24.07135):2.3295):5.96842):10.01763,Strigops-habroptila:42.3869):23.66972,((Rhegmatorhina-hoffmannsi:46.95689,((((((((Sylvia-borin:12.14381,Sylvia-atricapilla:12.14381):9.78774,Acrocephalus-scirpaceus:21.93155):0.45018,Prinia-subflava:22.38173):0.46108,Hirundo-rustica:22.84281):2.21723,(Parus-major:11.24516,Poecile-atricapillus:11.24516):13.81488):0.44567,((Certhia-americana:23.99299,(((Acridotheres-tristis:15.00999,Lamprotornis-superbus:15.00999):2.16848,Sturnus-vulgaris:17.17847):4.63848,((Ficedula-albicollis:14.20321,Oenanthe-melanoleuca:14.20321):1.7755,Erithacus-rubecula:15.97871):5.83824):2.17604):1.24098,((Diglossa-brunneiventris:20.4672,(((Corvus-moneduloides:6.45956,Corvus-cornix:6.45956):10.04272,((((Taeniopygia-guttata:12.48068,Lonchura-striata:12.48068):3.29017,(Vidua-chalybeata:5.47499,Vidua-macroura:5.47499):10.29586):0.27976,Corvus-hawaiiensis:16.05061):0.24588,Camarhynchus-parvulus:16.29649):0.20579):0.39624,Coloeus-monedula:16.89852):3.56868):3.74171,((Motacilla-alba:22.10999,((((Junco-hyemalis:4.54678,Zonotrichia-leucophrys:4.54678):4.28564,(Melospiza-georgiana:4.73478,Ammodramus-caudacutus:4.73478):4.09764):5.55445,((Geothlypis-trichas:6.40715,Setophaga-coronata:6.40715):7.56803,(Molothrus-ater:6.42673,Agelaius-phoeniceus:6.42673):7.54845):0.41169):5.82144,((Serinus-canaria:11.69604,Haemorhous-mexicanus:11.69604):6.1215,Fringilla-coelebs:17.81754):2.39077):1.90168):1.10266,Passer-domesticus:23.21265):0.99626):1.02506):0.27174):0.64574,Acanthisitta-chloris:26.15145):20.23534,(Lichenostomus-cassidix:30.26957,Malurus-cyaneus:30.26957):16.11722):0.5701):7.55693,Chiroxiphia-lanceolata:54.51382):11.5428):1e-05,((Falco-peregrinus:2.1397,((Falco-cherrug:0.29,Falco-rusticolus:0.29):0.92,Falco-biarmicus:1.21):0.9297):4.83106,Falco-naumanni:6.97076):59.08587):0.58922,Cariama-cristata:66.64585):3.78237):4.68227):1.15546):7.94685,(Otis-tarda:64.908,Cuculus-canorus:64.908):19.3048):6.63686,((((Anas-platyrhynchos:10.53724,(Aythya-marila:2.38746,Aythya-fuligula:2.38746):8.14978):23.8111,(((Anser-cygnoides:4.24954,Anser-indicus:4.24954):16.59046,(Cygnus-olor:0.01,Cygnus-atratus:0.01):20.83):13.21907,Oxyura-jamaicensis:34.05907):0.28927):1.13981,Cairina-moschata:35.48815):47.89678,((Coturnix-japonica:37.91099,((Lagopus-muta:19.96726,Meleagris-gallopavo:19.96726):15.8854,Gallus-gallus:35.85266):2.05833):8.43529,Numida-meleagris:46.34628):37.03865):7.46473):16.74216)

#search for one parameter model
#lambda -s

#search for 2 parameter model
lambdamu -s

