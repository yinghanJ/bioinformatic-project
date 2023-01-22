library("plyr")
library("dplyr")
library(data.table)
library(tidyverse)
c2_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c2_snps_nafld.txt")
c5_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c5_snps_nafld.txt")
c8_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c8_snps_nafld.txt")
c19_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c19_snps_nafld.txt")
c20_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c20_snps_nafld.txt")
c22_snps <- fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/SNPs_NAFLD/c22_snps_nafld.txt")
data_final<-fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/data_ukbb_cubn_pnpla3_07042021_final.txt")
dia_data_final<-fread(file="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/data_dm_ukbb_no_egfr_outliers_07042021_final.txt")

merge_c2 <-inner_join(c2_snps, data_final, by ='eid')
merge_c19 <-inner_join(c19_snps, data_final, by ='eid')
merge_c20 <-inner_join(c20_snps, data_final, by ='eid')
merge_c22 <-inner_join(c22_snps, data_final, by ='eid')
merge_c5 <-inner_join(c5_snps, data_final, by ='eid')
merge_c8 <-inner_join(c8_snps, data_final, by ='eid')


c2_snp_data <- merge_c2 %>% select(eid,sex,age,ln_gfr,rs1260326_T,rs780094_T)
c19_snp_data <- merge_c19 %>% select(eid,sex,age,ln_gfr,rs58542926_T, rs10401969_C,rs17751061_T, rs73001065_C,rs4808199_A,rs12977524_G,rs73002956_G,rs77215230_A,rs641738_T,rs12975366_C)
c20_snp_data <- merge_c20 %>% select(eid,sex,age,ln_gfr,rs2295490_G)
c22_snp_data <- merge_c22 %>% select(eid,sex,age,ln_gfr,rs738409_G.x,rs738408_T,rs3747207_A,rs78569621_T)
c5_snp_data <- merge_c5 %>% select(eid,sex,age,ln_gfr,rs17056583_C,rs10515775_T,rs17719770_C,rs12188266_G)
c8_snp_data <- merge_c8 %>% select(eid,sex,age,ln_gfr,rs17321515_G,rs2954029_T,rs2980859_C,rs28601761_G)

#model1 the snps and the ln gfr
model1_c2_rs1260326_T <- lm(ln_gfr ~ rs1260326_T,data=c2_snp_data%>% select(eid,rs1260326_T,ln_gfr))
model1_c2_rs780094_T <- lm(ln_gfr ~ rs780094_T,data=c2_snp_data%>% select(eid,rs780094_T,ln_gfr))

model1_c19_rs58542926_T <- lm(ln_gfr ~ rs58542926_T,data=c19_snp_data%>% select(eid,rs58542926_T,ln_gfr))
model1_c19_rs10401969_C <- lm(ln_gfr ~ rs10401969_C,data=c19_snp_data%>% select(eid,rs10401969_C,ln_gfr))
model1_c19_rs17751061_T <- lm(ln_gfr ~ rs17751061_T,data=c19_snp_data%>% select(eid,rs17751061_T,ln_gfr))
model1_c19_rs73001065_C <- lm(ln_gfr ~ rs73001065_C,data=c19_snp_data%>% select(eid,rs73001065_C,ln_gfr))
model1_c19_rs4808199_A <- lm(ln_gfr ~ rs4808199_A,data=c19_snp_data%>% select(eid,rs4808199_A,ln_gfr))
model1_c19_rs12977524_G <- lm(ln_gfr ~rs12977524_G,data=c19_snp_data%>% select(eid,rs12977524_G,ln_gfr))
model1_c19_rs73002956_G <- lm(ln_gfr ~ rs73002956_G,data=c19_snp_data%>% select(eid,rs73002956_G,ln_gfr))
model1_c19_rs77215230_A <- lm(ln_gfr ~ rs77215230_A,data=c19_snp_data%>% select(eid,rs77215230_A,ln_gfr))
model1_c19_rs641738_T <- lm(ln_gfr ~ rs641738_T,data=c19_snp_data%>% select(eid,rs641738_T,ln_gfr))
model1_c19_rs12975366_C <- lm(ln_gfr ~ rs12975366_C,data=c19_snp_data%>% select(eid,rs12975366_C,ln_gfr))

model1_c20_rs2295490_G <- lm(ln_gfr ~ rs2295490_G,data=c20_snp_data%>% select(eid,rs2295490_G,ln_gfr))

model1_c5_rs17056583_C <- lm(ln_gfr ~ rs17056583_C,data=c5_snp_data%>% select(eid,rs17056583_C,ln_gfr))
model1_c5_rs10515775_T <- lm(ln_gfr ~ rs10515775_T,data=c5_snp_data%>% select(eid,rs10515775_T,ln_gfr))
model1_c5_rs17719770_C <- lm(ln_gfr ~ rs17719770_C,data=c5_snp_data%>% select(eid,rs17719770_C,ln_gfr))
model1_c5_rs12188266_G <- lm(ln_gfr ~ rs12188266_G,data=c5_snp_data%>% select(eid,rs12188266_G,ln_gfr))

model1_c8_rs17321515_G <- lm(ln_gfr ~ rs17321515_G,data=c8_snp_data%>% select(eid,rs17321515_G,ln_gfr))
model1_c8_rs2954029_T <- lm(ln_gfr ~ rs2954029_T,data=c8_snp_data%>% select(eid,rs2954029_T,ln_gfr))
model1_c8_rs2980859_C <- lm(ln_gfr ~ rs2980859_C,data=c8_snp_data%>% select(eid,rs2980859_C,ln_gfr))
model1_c8_rs28601761_G <- lm(ln_gfr ~ rs28601761_G,data=c8_snp_data%>% select(eid,rs28601761_G,ln_gfr))

model1_c22_rs738409_G.x <- lm(ln_gfr ~ rs738409_G.x,data=c22_snp_data%>% select(eid,rs738409_G.x,ln_gfr))
model1_c22_rs738408_T <- lm(ln_gfr ~ rs738408_T,data=c22_snp_data%>% select(eid,rs738408_T,ln_gfr))
model1_c22_rs3747207_A <- lm(ln_gfr ~ rs3747207_A,data=c22_snp_data%>% select(eid,rs3747207_A,ln_gfr))
model1_c22_rs78569621_T <- lm(ln_gfr ~ rs78569621_T,data=c22_snp_data%>% select(eid,rs78569621_T,ln_gfr))

#model2 sex age ln gfr in final data set
model2_c2_rs1260326_T <- lm(ln_gfr ~ rs1260326_T + sex + age,data=c2_snp_data%>% select(eid,rs1260326_T,sex,age,ln_gfr))
model2_c2_rs780094_T <- lm(ln_gfr ~ rs780094_T+ sex + age,data=c2_snp_data%>% select(eid,rs780094_T,sex,age,ln_gfr))

model2_c19_rs58542926_T <- lm(ln_gfr ~ rs58542926_T+ sex + age,data=c19_snp_data%>% select(eid,rs58542926_T,sex,age,ln_gfr))
model2_c19_rs10401969_C <- lm(ln_gfr ~ rs10401969_C+ sex + age,data=c19_snp_data%>% select(eid,rs10401969_C,sex,age,ln_gfr))
model2_c19_rs17751061_T <- lm(ln_gfr ~ rs17751061_T+ sex + age,data=c19_snp_data%>% select(eid,rs17751061_T,sex,age,ln_gfr))
model2_c19_rs73001065_C <- lm(ln_gfr ~ rs73001065_C+ sex + age,data=c19_snp_data%>% select(eid,rs73001065_C,sex,age,ln_gfr))
model2_c19_rs4808199_A <- lm(ln_gfr ~ rs4808199_A+ sex + age,data=c19_snp_data%>% select(eid,rs4808199_A,sex,age,ln_gfr))
model2_c19_rs12977524_G <- lm(ln_gfr ~rs12977524_G+ sex + age,data=c19_snp_data%>% select(eid,rs12977524_G,sex,age,ln_gfr))
model2_c19_rs73002956_G <- lm(ln_gfr ~ rs73002956_G+ sex + age,data=c19_snp_data%>% select(eid,rs73002956_G,sex,age,ln_gfr))
model2_c19_rs77215230_A <- lm(ln_gfr ~ rs77215230_A+ sex + age,data=c19_snp_data%>% select(eid,rs77215230_A,sex,age,ln_gfr))
model2_c19_rs641738_T <- lm(ln_gfr ~ rs641738_T+ sex + age,data=c19_snp_data%>% select(eid,rs641738_T,sex,age,ln_gfr))
model2_c19_rs12975366_C <- lm(ln_gfr ~ rs12975366_C+ sex + age,data=c19_snp_data%>% select(eid,rs12975366_C,sex,age,ln_gfr))

model2_c20_rs2295490_G <- lm(ln_gfr ~ rs2295490_G+ sex + age,data=c20_snp_data%>% select(eid,rs2295490_G,sex,age,ln_gfr))

model2_c5_rs17056583_C <- lm(ln_gfr ~ rs17056583_C+ sex + age,data=c5_snp_data%>% select(eid,rs17056583_C,sex,age,ln_gfr))
model2_c5_rs10515775_T <- lm(ln_gfr ~ rs10515775_T+ sex + age,data=c5_snp_data%>% select(eid,rs10515775_T,sex,age,ln_gfr))
model2_c5_rs17719770_C <- lm(ln_gfr ~ rs17719770_C+ sex + age,data=c5_snp_data%>% select(eid,rs17719770_C,sex,age,ln_gfr))
model2_c5_rs12188266_G <- lm(ln_gfr ~ rs12188266_G+ sex + age,data=c5_snp_data%>% select(eid,rs12188266_G,sex,age,ln_gfr))

model2_c8_rs17321515_G <- lm(ln_gfr ~ rs17321515_G+ sex + age,data=c8_snp_data%>% select(eid,rs17321515_G,sex,age,ln_gfr))
model2_c8_rs2954029_T <- lm(ln_gfr ~ rs2954029_T+ sex + age,data=c8_snp_data%>% select(eid,rs2954029_T,sex,age,ln_gfr))
model2_c8_rs2980859_C <- lm(ln_gfr ~ rs2980859_C+ sex + age,data=c8_snp_data%>% select(eid,rs2980859_C,sex,age,ln_gfr))
model2_c8_rs28601761_G <- lm(ln_gfr ~ rs28601761_G+ sex + age,data=c8_snp_data%>% select(eid,rs28601761_G,sex,age,ln_gfr))

model2_c22_rs738409_G.x <- lm(ln_gfr ~ rs738409_G.x+ sex + age,data=c22_snp_data%>% select(eid,rs738409_G.x,sex,age,ln_gfr))
model2_c22_rs738408_T <- lm(ln_gfr ~ rs738408_T+ sex + age,data=c22_snp_data%>% select(eid,rs738408_T,sex,age,ln_gfr))
model2_c22_rs3747207_A <- lm(ln_gfr ~ rs3747207_A+ sex + age,data=c22_snp_data%>% select(eid,rs3747207_A,sex,age,ln_gfr))
model2_c22_rs78569621_T <- lm(ln_gfr ~ rs78569621_T+ sex + age,data=c22_snp_data%>% select(eid,rs78569621_T,sex,age,ln_gfr))

library(broom)
tidy(model2_c2_rs1260326_T)
augment(model2_c2_rs1260326_T,age)
library(stargazer)
#stargazer(model1_c2_rs1260326_T,model1_c2_rs780094_T,title="123",out="L:/LovbeskyttetMapper/Masters Thesis Nov22/02_Yinghan_UKBB/modeltest.txt")
#the result is not it i need

model_loop_test <- c2_snp_data %>% 
group_by(eid) %>%
  do(
    lm(egfr ~ fu_time_yrs, data = .) %>% tidy(conf.int = T) %>% slice(2)
  )
