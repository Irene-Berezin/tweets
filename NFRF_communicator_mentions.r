# required libraries
library(tidyverse)
library(dplyr)

# importing csvs
adriandix <- read_csv("csvs/twitterpho/@adriandix_tweets.csv")
AlikaMD <- read_csv("csvs/twitterpho/@AlikaMD_tweets.csv")
AmyGreerKalisz <- read_csv("csvs/twitterpho/@AmyGreerKalisz_tweets.csv")
angie_rasmussen <- read_csv("csvs/twitterpho/@angie_rasmussen_tweets.csv")
AnnaBlakney <- read_csv("csvs/twitterpho/@AnnaBlakney_tweets.csv")
AntibioticDoc <- read_csv("csvs/twitterpho/@AntibioticDoc_tweets.csv")
asapscience <- read_csv("csvs/twitterpho/@asapscience_tweets.csv")
ASPphysician <- read_csv("csvs/twitterpho/@ASPphysician_tweets.csv")
atRachelGilmore <- read_csv("csvs/twitterpho/@atRachelGilmore_tweets.csv")
binhanv <- read_csv("csvs/twitterpho/@binhanv_tweets.csv")
BirinderNarang <- read_csv("csvs/twitterpho/@BirinderNarang_tweets.csv")
blackdocscanada <- read_csv("csvs/twitterpho/@blackdocscanada_tweets.csv")
BogochIsaac <- read_csv("csvs/twitterpho/@BogochIsaac_tweets.csv")
bornk <- read_csv("csvs/twitterpho/@bornk_tweets.csv")
carlyweeks <- read_csv("csvs/twitterpho/@carlyweeks_tweets.csv")
CaulfieldTim <- read_csv("csvs/twitterpho/@CaulfieldTim_tweets.csv")
CDCofBC <- read_csv("csvs/twitterpho/@CDCofBC_tweets.csv")
cdube_sante <- read_csv("csvs/twitterpho/@cdube_sante_tweets.csv")
cfpcceo <- read_csv("csvs/twitterpho/@cfpcceo_tweets.csv")
ChiefSciCan <- read_csv("csvs/twitterpho/@ChiefSciCan_tweets.csv")
cmcovidtf <- read_csv("csvs/twitterpho/@cmcovidtf_tweets.csv")
CMOH_Alberta <- read_csv("csvs/twitterpho/@CMOH_Alberta_tweets.csv")
CMOH_NL <- read_csv("csvs/twitterpho/@CMOH_NL_tweets.csv")
conquercovid19 <- read_csv("csvs/twitterpho/@conquercovid19_tweets.csv")
COVID_19_Canada <- read_csv("csvs/twitterpho/@COVID_19_Canada_tweets.csv")
COVIDSciOntario <- read_csv("csvs/twitterpho/@COVIDSciOntario_tweets.csv")
CPHO_Canada <- read_csv("csvs/twitterpho/@CPHO_Canada_tweets.csv")
ctouzin <- read_csv("csvs/twitterpho/@ctouzin_tweets.csv")
CTV_AvisFavaro <- read_csv("csvs/twitterpho/@CTV_AvisFavaro_tweets.csv")
DeNovo_Fatima <- read_csv("csvs/twitterpho/@DeNovo_Fatima_tweets.csv")
deonandan <- read_csv("csvs/twitterpho/@deonandan_tweets.csv")
drfisman <- read_csv("csvs/twitterpho/@drfisman_tweets.csv")
Dr_ChrisSimpson <- read_csv("csvs/twitterpho/@Dr_ChrisSimpson_tweets.csv")
drgigiosler <- read_csv("csvs/twitterpho/@drgigiosler_tweets.csv")
DrKaliBarrett <- read_csv("csvs/twitterpho/@DrKaliBarrett_tweets.csv")
drmwarner <- read_csv("csvs/twitterpho/@drmwarner_tweets.csv")
drsusanshaw <- read_csv("csvs/twitterpho/@drsusanshaw_tweets.csv")
DrVivianS <- read_csv("csvs/twitterpho/@DrVivianS_tweets.csv")
egpayne <- read_csv("csvs/twitterpho/@egpayne_tweets.csv")
epdevilla <- read_csv("csvs/twitterpho/@epdevilla_tweets.csv")
ErnieHudsonPEI <- read_csv("csvs/twitterpho/@ErnieHudsonPEI_tweets.csv")
everetthindley <- read_csv("csvs/twitterpho/@everetthindley_tweets.csv")
First10EM <- read_csv("csvs/twitterpho/@First10EM_tweets.csv")
GermHunterMD <- read_csv("csvs/twitterpho/@GermHunterMD_tweets.csv")
glenpyle <- read_csv("csvs/twitterpho/@glenpyle_tweets.csv")
heysciencesam <- read_csv("csvs/twitterpho/@heysciencesam_tweets.csv")
hgagneTVA <- read_csv("csvs/twitterpho/@hgagneTVA_tweets.csv")
IDEpiPhD <- read_csv("csvs/twitterpho/@IDEpiPhD_tweets.csv")
imgrund <- read_csv("csvs/twitterpho/@imgrund_tweets.csv")
iPreetBrar <- read_csv("csvs/twitterpho/@iPreetBrar_tweets.csv")
IrfanDhalla <- read_csv("csvs/twitterpho/@IrfanDhalla_tweets.csv")
j_mcelroy <- read_csv("csvs/twitterpho/@j_mcelroy_tweets.csv")
jasonfherring <- read_csv("csvs/twitterpho/@jasonfherring_tweets.csv")
jfrketich <- read_csv("csvs/twitterpho/@jfrketich_tweets.csv")
jkwan_md <- read_csv("csvs/twitterpho/@jkwan_md_tweets.csv")
Johnrockdoc <- read_csv("csvs/twitterpho/@Johnrockdoc_tweets.csv")
JuliaWongCBC <- read_csv("csvs/twitterpho/@JuliaWongCBC_tweets.csv")
juliegreenMLA <- read_csv("csvs/twitterpho/@juliegreenMLA_tweets.csv")
Justin_Ling <- read_csv("csvs/twitterpho/@Justin_Ling_tweets.csv")
jyangstar <- read_csv("csvs/twitterpho/@jyangstar_tweets.csv")
KashPrime <- read_csv("csvs/twitterpho/@KashPrime_tweets.csv")
KatharineSmart <- read_csv("csvs/twitterpho/@KatharineSmart_tweets.csv")
Kevin_Parent <- read_csv("csvs/twitterpho/@Kevin__Parent_tweets.csv")
KindrachuckJason <- read_csv("csvs/twitterpho/@KindrachukJason_tweets.csv")
KrishanaSankar <- read_csv("csvs/twitterpho/@KrishanaSankar_tweets.csv")
kwadwo777 <- read_csv("csvs/twitterpho/@kwadwo777_tweets.csv")
LaurenPelley <- read_csv("csvs/twitterpho/@LaurenPelley_tweets.csv")
LisaBarrettID <- read_csv("csvs/twitterpho/@LisaBarrettID_tweets.csv")
McGillOSS <- read_csv("csvs/twitterpho/@McGillOSS_tweets.csv")
MerrimanPaul <- read_csv("csvs/twitterpho/@MerrimanPaul_tweets.csv")
MichaelSchwandt <- read_csv("csvs/twitterpho/@MichaelSchwandt_tweets.csv")
MLAStefanson <- read_csv("csvs/twitterpho/@MLAStefanson_tweets.csv")
moirawyton <- read_csv("csvs/twitterpho/@moirawyton_tweets.csv")
moriartylabs <- read_csv("csvs/twitterpho/@moriartylabs_tweets.csv")
MPaiMD <- read_csv("csvs/twitterpho/@MPaiMD_tweets.csv")
NaheedD <- read_csv("csvs/twitterpho/@NaheedD_tweets.csv")
NathanStall <- read_csv("csvs/twitterpho/@NathanStall_tweets.csv")
NightShiftMD <- read_csv("csvs/twitterpho/@NightShiftMD_tweets.csv")
NoLore <- read_csv("csvs/twitterpho/@NoLore_tweets.csv")
OttawaHealth <- read_csv("csvs/twitterpho/@OttawaHealth_tweets.csv")
paimadhu <- read_csv("csvs/twitterpho/@paimadhu_tweets.csv")
PattyHajdu <- read_csv("csvs/twitterpho/@PattyHajdu_tweets.csv")
picardonhealth <- read_csv("csvs/twitterpho/@picardonhealth_tweets.csv")
RicharLisa <- read_csv("csvs/twitterpho/@RicharLisa_tweets.csv")
roussin_brent <- read_csv("csvs/twitterpho/@roussin_brent_tweets.csv")
sabaeitizaz <- read_csv("csvs/twitterpho/@sabaeitizaz_tweets.csv")
sabiVM <- read_csv("csvs/twitterpho/@sabiVM_tweets.csv")
SammyG_MD <- read_csv("csvs/twitterpho/@SammyG_MD_tweets.csv")
sarperotto <- read_csv("csvs/twitterpho/@sarperotto_tweets.csv")
SciChefCan <- read_csv("csvs/twitterpho/@SciChefCan_tweets.csv")
sciencemonkeyca <- read_csv("csvs/twitterpho/@sciencemonkeyca_tweets.csv")
ScienceUpFirst <- read_csv("csvs/twitterpho/@ScienceUpFirst_tweets.csv")
sdbaral <- read_csv("csvs/twitterpho/@sdbaral_tweets.csv")
shandro <- read_csv("csvs/twitterpho/@shandro_tweets.csv")
SharkawyMD <- read_csv("csvs/twitterpho/@SharkawyMD_tweets.csv")
shazmamithani <- read_csv("csvs/twitterpho/@shazmamithani_tweets.csv")
ShephardDorothy <- read_csv("csvs/twitterpho/@ShephardDorothy_tweets.csv")
srinmurthy99 <- read_csv("csvs/twitterpho/@srinmurthy99_tweets.csv")
SteiniBrown <- read_csv("csvs/twitterpho/@SteiniBrown_tweets.csv")
theresaboyle <- read_csv("csvs/twitterpho/@theresaboyle_tweets.csv")
thisisourshotca <- read_csv("csvs/twitterpho/@thisisourshotca_tweets.csv")
TorontoIDDOC <- read_csv("csvs/twitterpho/@TorontoIDDOC_tweets.csv")
UbakaOgbogu <- read_csv("csvs/twitterpho/@UbakaOgbogu_tweets.csv")
VaxHuntersCan <- read_csv("csvs/twitterpho/@VaxHuntersCan_tweets.csv")
VeraEtches <- read_csv("csvs/twitterpho/@VeraEtches_tweets.csv")
VikCBC <- read_csv("csvs/twitterpho/@VikCBC_tweets.csv")
wickdchiq <- read_csv("csvs/twitterpho/@wickdchiq_tweets.csv")
zachchurchhill <- read_csv("csvs/twitterpho/@zachchurchill_tweets.csv")
zchangla <- read_csv("csvs/twitterpho/@zchagla_tweets.csv")
DrKathleenRoss1 <- read_csv("csvs/twitterpho/DrKathleenRoss1_tweets.csv")

# Your original list of datasets with their names
df_list <- list(
  adriandix, AlikaMD, AmyGreerKalisz, angie_rasmussen, AnnaBlakney, AntibioticDoc,
  asapscience, ASPphysician, atRachelGilmore, binhanv, BirinderNarang, blackdocscanada,
  BogochIsaac, bornk, carlyweeks, CaulfieldTim, CDCofBC, cdube_sante, cfpcceo, ChiefSciCan,
  cmcovidtf, CMOH_Alberta, CMOH_NL, conquercovid19, COVID_19_Canada, COVIDSciOntario,
  CPHO_Canada, ctouzin, CTV_AvisFavaro, DeNovo_Fatima, deonandan, drfisman, Dr_ChrisSimpson,
  drgigiosler, DrKaliBarrett, drmwarner, drsusanshaw, DrVivianS, egpayne, epdevilla, ErnieHudsonPEI,
  everetthindley, First10EM, GermHunterMD, glenpyle, heysciencesam, hgagneTVA, IDEpiPhD, imgrund,
  iPreetBrar, IrfanDhalla, j_mcelroy, jasonfherring, jfrketich, jkwan_md, Johnrockdoc, JuliaWongCBC,
  juliegreenMLA, Justin_Ling, jyangstar, KashPrime, KatharineSmart, Kevin_Parent, KindrachuckJason,
  KrishanaSankar, kwadwo777, LaurenPelley, LisaBarrettID, McGillOSS, MerrimanPaul, MichaelSchwandt,
  MLAStefanson, moirawyton, moriartylabs, MPaiMD, NaheedD, NathanStall, NightShiftMD, NoLore, OttawaHealth,
  paimadhu, PattyHajdu, picardonhealth, RicharLisa, roussin_brent, sabaeitizaz, sabiVM, SammyG_MD,
  sarperotto, SciChefCan, sciencemonkeyca, ScienceUpFirst, sdbaral, shandro, SharkawyMD, shazmamithani,
  ShephardDorothy, srinmurthy99, SteiniBrown, theresaboyle, thisisourshotca, TorontoIDDOC, UbakaOgbogu,
  VaxHuntersCan, VeraEtches, VikCBC, wickdchiq, zachchurchhill, zchangla, DrKathleenRoss1
)

remove_retweets <- function(df) {
  df %>%
    filter(!startsWith(text, "RT")) %>%
    select(text) 
}

# Applying function
transformed_list <- lapply(df_list, transform_data)

# Assign the names of the datasets to the transformed list
names(transformed_list) <- c("adriandix", "AlikaMD", "AmyGreerKalisz", "angie_rasmussen", "AnnaBlakney", 
                             "AntibioticDoc", "asapscience", "ASPphysician", "atRachelGilmore", "binhanv", 
                             "BirinderNarang", "blackdocscanada", "BogochIsaac", "bornk", "carlyweeks", 
                             "CaulfieldTim", "CDCofBC", "cdube_sante", "cfpcceo", "ChiefSciCan", "cmcovidtf", 
                             "CMOH_Alberta", "CMOH_NL", "conquercovid19", "COVID_19_Canada", "COVIDSciOntario", 
                             "CPHO_Canada", "ctouzin", "CTV_AvisFavaro", "DeNovo_Fatima", "deonandan", "drfisman", 
                             "Dr_ChrisSimpson", "drgigiosler", "DrKaliBarrett", "drmwarner", "drsusanshaw", 
                             "DrVivianS", "egpayne", "epdevilla", "ErnieHudsonPEI", "everetthindley", "First10EM", 
                             "GermHunterMD", "glenpyle", "heysciencesam", "hgagneTVA", "IDEpiPhD", "imgrund", 
                             "iPreetBrar", "IrfanDhalla", "j_mcelroy", "jasonfherring", "jfrketich", "jkwan_md", 
                             "Johnrockdoc", "JuliaWongCBC", "juliegreenMLA", "Justin_Ling", "jyangstar", "KashPrime", 
                             "KatharineSmart", "Kevin_Parent", "KindrachuckJason", "KrishanaSankar", "kwadwo777", 
                             "LaurenPelley", "LisaBarrettID", "McGillOSS", "MerrimanPaul", "MichaelSchwandt", 
                             "MLAStefanson", "moirawyton", "moriartylabs", "MPaiMD", "NaheedD", "NathanStall", 
                             "NightShiftMD", "NoLore", "OttawaHealth", "paimadhu", "PattyHajdu", "picardonhealth", 
                             "RicharLisa", "roussin_brent", "sabaeitizaz", "sabiVM", "SammyG_MD", "sarperotto", 
                             "SciChefCan", "sciencemonkeyca", "ScienceUpFirst", "sdbaral", "shandro", "SharkawyMD", 
                             "shazmamithani", "ShephardDorothy", "srinmurthy99", "SteiniBrown", "theresaboyle", 
                             "thisisourshotca", "TorontoIDDOC", "UbakaOgbogu", "VaxHuntersCan", "VeraEtches", 
                             "VikCBC", "wickdchiq", "zachchurchhill", "zchangla", "DrKathleenRoss1")

# Use list2env to assign each element of transformed_list to a variable in the global environment
list2env(transformed_list, envir = .GlobalEnv)


# Create a new dataframe which contains the # of instances each user is mentioned
tweet_counts <- data.frame(
    User = c(
        "adriandix", "AlikaMD", "AmyGreerKalisz", "angie_rasmussen", "AnnaBlakney", "AntibioticDoc",
        "asapscience", "ASPphysician", "atRachelGilmore", "binhanv", "BirinderNarang",
        "blackdocscanada", "BogochIsaac", "bornk", "carlyweeks", "CaulfieldTim", "CDCofBC",
        "cdube_sante", "cfpcceo", "ChiefSciCan", "cmcovidtf", "CMOH_Alberta", "CMOH_NL",
        "conquercovid19", "COVID_19_Canada", "COVIDSciOntario", "CPHO_Canada", "ctouzin",
        "CTV_AvisFavaro", "DeNovo_Fatima", "deonandan", "drfisman", "Dr_ChrisSimpson",
        "drgigiosler", "DrKaliBarrett", "drmwarner", "drsusanshaw", "DrVivianS", "egpayne",
        "epdevilla", "ErnieHudsonPEI", "everetthindley", "First10EM", "GermHunterMD", "glenpyle",
        "heysciencesam", "hgagneTVA", "IDEpiPhD", "imgrund", "iPreetBrar", "IrfanDhalla",
        "j_mcelroy", "jasonfherring", "jfrketich", "jkwan_md", "Johnrockdoc", "JuliaWongCBC",
        "juliegreenMLA", "Justin_Ling", "jyangstar", "KashPrime", "KatharineSmart", "Kevin_Parent",
        "KindrachuckJason", "KrishanaSankar", "kwadwo777", "LaurenPelley", "LisaBarrettID",
        "McGillOSS", "MerrimanPaul", "MichaelSchwandt", "MLAStefanson", "moirawyton",
        "moriartylabs", "MPaiMD", "NaheedD", "NathanStall", "NightShiftMD", "NoLore",
        "OttawaHealth", "paimadhu", "PattyHajdu", "picardonhealth", "RicharLisa", "roussin_brent",
        "sabaeitizaz", "sabiVM", "SammyG_MD", "sarperotto", "SciChefCan", "sciencemonkeyca",
        "ScienceUpFirst", "sdbaral", "shandro", "SharkawyMD", "shazmamithani", "ShephardDorothy",
        "srinmurthy99", "SteiniBrown", "theresaboyle", "thisisourshotca", "TorontoIDDOC",
        "UbakaOgbogu", "VaxHuntersCan", "VeraEtches", "VikCBC", "wickdchiq", "zachchurchhill",
        "zchangla", "DrKathleenRoss1"
    ),
    Count = c(
        nrow(adriandix), nrow(AlikaMD), nrow(AmyGreerKalisz), nrow(angie_rasmussen),
        nrow(AnnaBlakney), nrow(AntibioticDoc), nrow(asapscience), nrow(ASPphysician),
        nrow(atRachelGilmore), nrow(binhanv), nrow(BirinderNarang), nrow(blackdocscanada),
        nrow(BogochIsaac), nrow(bornk), nrow(carlyweeks), nrow(CaulfieldTim), nrow(CDCofBC),
        nrow(cdube_sante), nrow(cfpcceo), nrow(ChiefSciCan), nrow(cmcovidtf), nrow(CMOH_Alberta),
        nrow(CMOH_NL), nrow(conquercovid19), nrow(COVID_19_Canada), nrow(COVIDSciOntario),
        nrow(CPHO_Canada), nrow(ctouzin), nrow(CTV_AvisFavaro), nrow(DeNovo_Fatima),
        nrow(deonandan), nrow(drfisman), nrow(Dr_ChrisSimpson), nrow(drgigiosler),
        nrow(DrKaliBarrett), nrow(drmwarner), nrow(drsusanshaw), nrow(DrVivianS), nrow(egpayne),
        nrow(epdevilla), nrow(ErnieHudsonPEI), nrow(everetthindley), nrow(First10EM),
        nrow(GermHunterMD), nrow(glenpyle), nrow(heysciencesam), nrow(hgagneTVA), nrow(IDEpiPhD),
        nrow(imgrund), nrow(iPreetBrar), nrow(IrfanDhalla), nrow(j_mcelroy), nrow(jasonfherring),
        nrow(jfrketich), nrow(jkwan_md), nrow(Johnrockdoc), nrow(JuliaWongCBC), nrow(juliegreenMLA),
        nrow(Justin_Ling), nrow(jyangstar), nrow(KashPrime), nrow(KatharineSmart), nrow(Kevin_Parent),
        nrow(KindrachuckJason), nrow(KrishanaSankar), nrow(kwadwo777), nrow(LaurenPelley),
        nrow(LisaBarrettID), nrow(McGillOSS), nrow(MerrimanPaul), nrow(MichaelSchwandt),
        nrow(MLAStefanson), nrow(moirawyton), nrow(moriartylabs), nrow(MPaiMD), nrow(NaheedD),
        nrow(NathanStall), nrow(NightShiftMD), nrow(NoLore), nrow(OttawaHealth), nrow(paimadhu),
        nrow(PattyHajdu), nrow(picardonhealth), nrow(RicharLisa), nrow(roussin_brent),
        nrow(sabaeitizaz), nrow(sabiVM), nrow(SammyG_MD), nrow(sarperotto), nrow(SciChefCan),
        nrow(sciencemonkeyca), nrow(ScienceUpFirst), nrow(sdbaral), nrow(shandro), nrow(SharkawyMD),
        nrow(shazmamithani), nrow(ShephardDorothy), nrow(srinmurthy99), nrow(SteiniBrown),
        nrow(theresaboyle), nrow(thisisourshotca), nrow(TorontoIDDOC), nrow(UbakaOgbogu),
        nrow(VaxHuntersCan), nrow(VeraEtches), nrow(VikCBC), nrow(wickdchiq), nrow(zachchurchhill),
        nrow(zchangla), nrow(DrKathleenRoss1)
    )
)

# Calculating total number of tweets across all datasets
total_tweets <- sum(tweet_counts$Count)

# Adding a column for the percentage
tweet_counts <- tweet_counts %>%
    mutate(percentage = Count / total_tweets * 100)

#run to see the df
tweet_counts 
