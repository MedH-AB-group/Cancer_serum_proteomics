## ----Libraries, message=FALSE, warning=FALSE, include=FALSE----------------------------------
library("readr")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ChAMP)
library(pastecs)
library(car)
library(ggpubr)
library(summarytools)
library(DT)


## ----read_data-------------------------------------------------------------------------------
###Reading the needed data
data <- read.csv2("data/protein_screening.csv")
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        mutate_at("sample_id", as.character) %>%
        arrange(group) 
## Reading and Cleaning metadata
metadata <- read.csv2("data/MetaDATA_Ghada.csv",  header=TRUE, check.names = FALSE) %>%
     dplyr::rename(resection_type = `Resect type GBC`) %>%
     mutate(Death = gsub('alive', '0', Death)) %>%
     mutate(Death = gsub('dead', '1', Death)) %>%
     mutate(Alb = na_if(Alb, "m")) %>%
     mutate(CRP0 = na_if(CRP0, "m")) %>%
     mutate_if(is.character, as.numeric) %>%
     mutate_at("resection_type", as.character)
###Table prepared for statistics on bmi, sex etc (table 1)
stat <- read.csv2("results/Paper_illustratios/stat_table.csv") 


## ----stat_summary_controls, message=FALSE, warning=FALSE, paged.print=FALSE, include=FALSE----
### statistical description of the metadata: age, bmi and proteins (rootinely measured)
datatable(stat.desc(metadata %>% 
                      filter(DiagPostop == "0") %>%
                      dplyr::select(AgeatDiagn, BMI, CRP0, Alb, Brb0, CA199)))
metadata %>% filter(DiagPostop == "0" & CA199 == 34)
dfSummary(metadata %>%
            filter(DiagPostop == "0") %>%
            dplyr::select(AgeatDiagn, BMI, CRP0, Alb, Brb0, CA199))


## ----stats-----------------------------------------------------------------------------------
## Statistics table
datatable(stat,
          extensions = 'Buttons', 
    options = list(dom = 'Bfrtip', 
                   buttons = c('excel', "csv")))




## ----sex, fig.cap="Proportions of females (F) and males (M) in the cohort"-------------------
###Defining bmi and age categories then plotting them
metadata %>% 
   mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(sex = as.character(Gender)) %>%
    mutate(sex = gsub('0', 'M', sex)) %>%
    mutate(sex = gsub('1', 'F', sex)) %>%
  group_by(DiagPostop, sex) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, sex) %>%
  distinct(count, DiagPostop, sex) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = sex)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_fill_brewer(palette = "Paired")


## ----age_cat, fig.cap="Age categories in the cohort"-----------------------------------------
###Visualizing age categories
metadata %>% 
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(age_cat = cut(AgeatDiagn, breaks = c(0, 20, 30, 40, 50, 60, 70, 80, 100), labels = c("<20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", ">80"))) %>%
  group_by(DiagPostop, age_cat) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, age_cat) %>%
  distinct(count, DiagPostop, age_cat) %>%
  ggplot(aes(fill = age_cat, y = count, x = DiagPostop)) +
  geom_bar(position = "stack", stat = "identity")


## ----bmi_category, fig.cap="BMI categories in the cohort"------------------------------------
##Visualize bmi categories
metadata %>% 
  mutate_at("BMI", as.integer) %>%
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
   mutate(bmi_cat = cut(BMI, breaks = c(0, 18, 25, 30, 35, 90), labels = c("<18", "18-25", "25-30", "30-35", ">35"))) %>%
  mutate(bmi_cat_description = cut(BMI, breaks = c(0, 18, 25, 30, 35, 90), labels = c("Underweight", "Normal", "Overweight", "Light_obesity", "Extreme_obesity"))) %>%
  group_by(DiagPostop, bmi_cat) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, bmi_cat) %>%
  distinct(count, DiagPostop, bmi_cat) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = bmi_cat)) +
  geom_bar(position = "stack", stat = "identity")


## ----survival, fig.cap="Overall survival after surgery in both groups"-----------------------
##survival after surgery
metadata %>% 
  group_by(DiagPostop, Death) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(Death = gsub('0', 'Survival', Death)) %>%
  mutate(Death = gsub('1', 'Death', Death)) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, Death) %>%
  distinct(count, DiagPostop, Death) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = Death)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))


## ----Resection_type, fig.cap="Type of the resection in both groups"--------------------------
##Type of resection in the two groups, visualization
metadata %>% 
  group_by(DiagPostop, resection_type) %>%
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(resection_type = gsub('0', 'Explorative laparotomy', resection_type)) %>%
  mutate(resection_type = gsub('1', 'Cholecystectomy', resection_type)) %>%
  mutate(resection_type = gsub('2', 'Cholecystectomy+BD', resection_type)) %>%
  mutate(resection_type = gsub('3', 'Bisegmentectomy', resection_type)) %>%
  mutate(resection_type = gsub('4', 'Bisegmentectomy+BD', resection_type)) %>%
  mutate(resection_type = gsub('5', 'Extended Cholecystectomy', resection_type)) %>%
  mutate(resection_type = gsub('6', 'Major resection', resection_type)) %>%
  mutate(resection_type = gsub('7', 'Other', resection_type)) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, resection_type) %>%
  distinct(count, DiagPostop, resection_type) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = resection_type)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_fill_brewer(palette = "Dark2")


## ----Cholangiography, fig.cap="Perfmed cholangiography preoperatively"-----------------------
### Cholangiography visualization
metadata %>% 
  group_by(DiagPostop, Cholangio) %>%
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(Cholangio = gsub('0', 'None', Cholangio)) %>%
  mutate(Cholangio = gsub('1', 'PTC', Cholangio)) %>%
  mutate(Cholangio = gsub('3', 'ERCP', Cholangio)) %>%
  mutate(Cholangio = gsub('4', 'both', Cholangio)) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, Cholangio) %>%
  distinct(count, DiagPostop, Cholangio) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = Cholangio)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_fill_brewer(palette = "Dark2")


## ----Radiology, fig.cap="Performed radiology preoperatively"---------------------------------
###Radiology visualisation
metadata %>% 
  group_by(DiagPostop, Radiology) %>%
  mutate(DiagPostop = gsub('0', 'Control', DiagPostop)) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  mutate(Radiology = gsub('1', 'CT', Radiology)) %>%
  mutate(Radiology = gsub('2', 'MRI', Radiology)) %>%
  mutate(Radiology = gsub('3', 'both', Radiology)) %>%
  mutate(Radiology = gsub('0', 'None', Radiology)) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, Radiology) %>%
  distinct(count, DiagPostop, Radiology) %>%
  ggplot(aes(y = count, x = DiagPostop, fill = Radiology)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_color_manual(values=c("red", "blue", "yellow", "cyan"))


## ----tumor_extension, fig.cap="Proportion of different types of tumor extensions in the GBC group"----
### Tumor extensions in GBC groups; visualisation
metadata %>% 
  group_by(DiagPostop, pT) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  filter(DiagPostop == "Cancer") %>%
  mutate(pT = gsub('0', 'in situ', pT)) %>%
  mutate(pT = gsub('2', '2a', pT)) %>%
  mutate(pT = gsub('3', '2b', pT)) %>%
  mutate(pT = gsub('4', '3', pT)) %>%
  mutate(pT = gsub('5', '4', pT)) %>%
  mutate(count = n()) %>%
  select(count, DiagPostop, pT) %>%
  distinct(count, DiagPostop, pT) %>%
  rename(pT = "tumor_extension") %>%
  ggplot(aes(y = count, x = DiagPostop, fill = tumor_extension)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_fill_brewer(palette = "Dark2")


## ----resection_margin, fig.cap="Resection margin of the performed resection surgery on GBC patients"----
### Visualising the resection margin in GBC
metadata %>% 
  group_by(DiagPostop, pR) %>%
  mutate(DiagPostop = gsub('4', 'Cancer', DiagPostop)) %>%
  filter(DiagPostop == "Cancer") %>%
  mutate(pR = gsub('0', 'R0', pR)) %>%
  mutate(pR = gsub('1', 'R1', pR)) %>%
  
  mutate(count = n()) %>%
  select(count, DiagPostop, pR) %>%
  distinct(count, DiagPostop, pR) %>%
  rename(pR = "Microscopically_tumor_positive_resection_margin") %>%
  ggplot(aes(y = count, x = DiagPostop, fill = Microscopically_tumor_positive_resection_margin)) +
  geom_bar(position = "stack", stat = "identity") +
 scale_fill_brewer(palette = "Dark2")


## ----SVD, message=FALSE, warning=FALSE, paged.print=FALSE, fig.keep='last', echo=FALSE, results='hide'----
###Singlar value decomposition analysis
champ.SVD(beta = as.data.frame(t(data %>%
                                   # filter(group == 0) %>%
                                   arrange(sample_id) %>%
                                   select(-sample_id))), 
## pd is the metadata (group and other clinical information as a data-frame)
          pd = as.data.frame(metadata %>% 
                              mutate_at("BMI", as.integer) %>%
                              mutate(bmi_cat = cut(BMI, breaks = c(0, 18, 25, 30, 35, 90), labels = c("<18", "18-25", "25-30", "30-35", ">35"))) %>%
                              mutate(age_cat = cut(AgeatDiagn, breaks = c(0, 20, 30, 40, 50, 60, 70, 80, 100), labels = c("<20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", ">80"))) %>%
                              dplyr::rename("sample_id" = `biobank number`) %>%
                              arrange(sample_id) %>%
                              inner_join(as.data.frame(data) %>%
                                          mutate_at("sample_id", as.integer) %>%
                                          select(sample_id, group),
                                        by = "sample_id") %>%
                             # filter(group == 0) %>%
                             mutate_at("Death", as.factor) %>%
                             # select(group, Gender, age_cat, bmi_cat, Diabetes, PSC, Cirrhosis)))
                             select(group, CRP0, Alb, pT, CA199, Brb0, `pN#pos`, Death)))



