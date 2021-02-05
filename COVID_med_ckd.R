rm(list = ls())
library('openxlsx')
library('tidyr', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
library('dplyr', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
library('logistf', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
options(stringsAsFactors = FALSE)


PATIENT_DATE = '08312020'
DATA_CREATION_DATE = '09012020'
ANALYSIS_DATE = format(Sys.time(), "%m%d%Y")


data_dir = 'J:\\OHIA_Investigators\\YiDing\\covid19\\Data'
result_dir = 'J:\\OHIA_Investigators\\YiDing\\covid19\\results'


time_string = 'before'
data_dir = file.path(data_dir, PATIENT_DATE)
result_dir = file.path(result_dir, PATIENT_DATE)
data_file = file.path(data_dir, sprintf('patient_data_%s_%s.csv', time_string, DATA_CREATION_DATE))
df = read.csv(data_file, row.names = 1)


enrichment_test <- function(risk_factors, outcome_string, 
                            cohort_string, 
                            df, 
                            time_string, 
                            membership_string,
                            adjust_string = 'adj_SA', race_ethnicity = 'All'){
    # subset data frame 
    condition =  ((df[outcome_string ] == 1) | (df[cohort_string] == 1))
    if (membership_string == 'NotNewToUCLA'){condition = (df['NewToUCLA'] == 0) & condition}
    
    if (race_ethnicity == 'Hispanic'){
        condition = condition & (df['SIRE_Ethnicity'] == 'Hispanic or Latino')
    }else if (race_ethnicity == 'NonHispanic'){
        condition = condition & (df['SIRE_Ethnicity'] == 'Not Hispanic or Latino')
    }else if (race_ethnicity == 'NonHispanicWhite'){
        condition = condition & ((df['SIRE_Ethnicity'] == 'Not Hispanic or Latino') & (df[,'SIRE_Race'] == 'White or Caucasian'))
    }else if (race_ethnicity == 'All'){
        condition = condition
    }else{
        stop('No this Ethnicity')
    }       

    analysis_df = df[condition,]
    # get covariates
    if (adjust_string == 'crude'){
        covariates = c()
    }else if (adjust_string == 'adj_SA'){ 
        covariates = c('age','age2', 'Sex')
    }else if (adjust_string == 'adj_SAM'){
        covariates = c('age','age2','Sex', 'NewToUCLA')
    }else if(adjust_string == 'adj_SA_krf'){
        covariates = c( 'age','age2','Sex', known_risk_factors)
    }else if(adjust_string == 'adj_SAM_krf'){
        covariates = c( 'age','age2','Sex', 'NewToUCLA', known_risk_factors)
    }else{
        stop('adjustment type not recognized')
    }
    
    if (cohort_string == 'DDRSAMPLE'){
        covariates = c('AgeGroup', covariates)
        adjust_string = paste0(adjust_string, '_', 'MC')
    }
        
    # initiate progress bar 
    total = length(risk_factors)
    pb = txtProgressBar(min = 0, max = total, style = 3)
    # hold results for each risk factor
    results_df = data.frame('Membership' = character(), 'Outcome' = character(), 'Cohort' = character(), 'Adjust'= character(), 'ckd/medicine'= character(), 
                'Odds Ratio' = numeric(), '2.5%'= numeric(), '97.5%'= numeric(), 'ORCI' = character(),'Pvalue (Firth)'= numeric(), 
                             '# ckd/medicine in Outcome+'=integer(), '# Outcome+'=integer(),
                              '# ckd/medicine in Outcome-'=integer(), '# Outcome-'=integer(),
                             'Freq ckd/medicine in Outcome+'= numeric(), 'Freq ckd/medicine in Outcome-'= numeric(),
                             'Freq ckd/medicine in SEVERE'= numeric(),'Freq ckd/medicine in INPATIENT'= numeric(),
                            'Freq ckd/medicine in COVID'= numeric(),'Freq ckd/medicine in TESTED'= numeric(),
                             '# ckd/medicine in SEVERE'=integer(), '# SEVERE' = integer(),
                             '# ckd/medicine in INPATIENT'=integer(), '# INPATIENT'=integer(),
                             '# ckd/medicine in COVID'=integer(),'# COVID'=integer(),
                             '# ckd/medicine in TESTED'=integer(), '# TESTED'=integer(),
                             'age coef' = numeric(), 'age 2.5%' = numeric(),'age 97.5%' = numeric(), 'age pvalue' = numeric(),
                             'age2 coef' = numeric(), 'age2 2.5%' = numeric(),'age2 97.5%' = numeric(), 'age2 pvalue' = numeric(),
                             'SexFemale' = numeric(), 'SexFemale 2.5%' = numeric(),'SexFemale 97.5%' = numeric(), 'SexFemale pvalue' = numeric())
    
    for (i in seq(total)){
        # update progress bar 
        setTxtProgressBar(pb, i)
        # run logistic regression
        rf = risk_factors[i]
        model_covariates = c(covariates, rf)
        f = as.formula( paste0(outcome_string, '~' ,paste(model_covariates, collapse = ' + '), sep = ' '))
        model = logistf(formula = f, data = analysis_df)
        # get odds ratio, ci and pvalue 
        p = model$prob[rf]
        oddsratio = round(exp(coef(model)[rf]),2)
        ci_lower = round(exp(confint(model)[rf,1]),2)
        ci_upper = round(exp(confint(model)[rf,2]),2)
        orci_string = sprintf('%s [%s, %s]', oddsratio, ci_lower, ci_upper)

        # age 
        age_coef = signif(exp(coef(model)['age']),3)
        age_lower = signif(exp(confint(model)['age',1]),3)
        age_upper = signif(exp(confint(model)['age',2]),3)
        age_p = model$prob['age']

        # age2
        age2_coef = signif(exp(1000*coef(model)['age2']),3)
        age2_lower = signif(exp(1000*confint(model)['age2',1]),3)
        age2_upper = signif(exp(1000*confint(model)['age2',2]),3)
        age2_p = model$prob['age2']

        # SexFemale
        SexFemale_coef = signif(exp(- coef(model)['SexMale']),3)
        SexFemale_lower = signif(exp(- confint(model)['SexMale',2]),3)
        SexFemale_upper = signif(exp(- confint(model)['SexMale',1]),3)
        SexFemale_p = model$prob['SexMale']
        
        # frequency in positive and negative outcome
        n_pos_rf = sum(analysis_df[rf] * analysis_df[outcome_string])
        n_neg_rf = sum(analysis_df[rf] * (1-analysis_df[outcome_string]))
        n_pos = sum(analysis_df[outcome_string])
        n_neg = sum((1-analysis_df[outcome_string]))
        freq_pos_rf = round(n_pos_rf/n_pos*100, 2)
        freq_neg_rf = round(n_neg_rf/n_neg*100, 2)
        
        # frequency in tested
        n_TESTED = sum(df['TESTED'])
        n_TESTED_rf = sum(df['TESTED']*df[rf])
        freq_TESTED_rf = round(n_TESTED_rf/n_TESTED*100,2)
        # frequency in covid
        n_COVID = sum(df['COVID'])
        n_COVID_rf = sum(df['COVID']*df[rf])
        freq_COVID_rf = round(n_COVID_rf/n_COVID*100,2)
        # frequency in inpatient
        n_INPATIENT = sum(df['INPATIENT'])
        n_INPATIENT_rf = sum(df['INPATIENT']*df[rf])
        freq_INPATIENT_rf = round(n_INPATIENT_rf/n_INPATIENT*100,2)
        # frequency in severe
        n_SEVERE = sum(df['SEVERE'])
        n_SEVERE_rf = sum(df['SEVERE']*df[rf])
        freq_SEVERE_rf = round(n_SEVERE_rf/n_SEVERE*100,2)
        
        result_row= list(membership_string, outcome_string, cohort_string, adjust_string, rf, 
                      oddsratio, ci_lower, ci_upper,orci_string, p,
                      n_pos_rf, n_pos, 
                      n_neg_rf, n_neg,
                      freq_pos_rf, freq_neg_rf, 
                      freq_SEVERE_rf, freq_INPATIENT_rf, freq_COVID_rf, freq_TESTED_rf,
                      n_SEVERE_rf, n_SEVERE,
                      n_INPATIENT_rf, n_INPATIENT,
                      n_COVID_rf, n_COVID,
                      n_TESTED_rf, n_TESTED,
                      age_coef, age_lower, age_upper, age_p, 
                      age2_coef, age2_lower, age2_upper, age2_p,
                      SexFemale_coef, SexFemale_lower, SexFemale_upper, SexFemale_p)
        results_df[i,] = result_row
    } 
  
    results_df = data.frame(results_df)
    colnames(results_df) = c('Membership', 'Outcome', 'Cohort', 'Adjust', 'ckd/medicine', 
                             'Odds Ratio', '2.5%', '97.5%', 'ORCI','Pvalue (Firth)', 
                             '# ckd/medicine in Outcome+', '# Outcome+',
                              '# ckd/medicine in Outcome-', '# Outcome-',
                             'Freq ckd/medicine in Outcome+', 'Freq ckd/medicine in Outcome-',
                             'Freq ckd/medicine in SEVERE','Freq ckd/medicine in INPATIENT','Freq ckd/medicine in COVID','Freq ckd/medicine in TESTED',
                             '# ckd/medicine in SEVERE', '# SEVERE',
                             '# ckd/medicine in INPATIENT', '# INPATIENT',
                             '# ckd/medicine in COVID','# COVID',
                             '# ckd/medicine in TESTED', '# TESTED',
                              'age OR' , 'age 2.5%' ,'age 97.5%' , 'age pvalue' ,
                             'age2 OR' , 'age2 2.5%' ,'age2 97.5%' , 'age2 pvalue' ,
                             'SexFemale OR' , 'SexFemale 2.5%' ,'SexFemale 97.5%' , 'SexFemale pvalue' )

    return(results_df)
}

membership_string = 'NotNewToUCLA'
known_risk_factors= c('chf', 'diabetes', 'hyperlipidemia', 'hypertension', 'obesity', 'ckd', 'copd', 'chd')
compare_list = list( c('COVID', 'TESTED'), c('INPATIENT', 'COVID'), c('SEVERE', 'INPATIENT'))
risk_factors = c('ImmunoSuppressants', 'steroids', 'ACEInhibitors', 'ARBs', 'AntiCoagulants', 'NSAID', 'SSRI')
sheet_list = list()
race_ethnicity_group = c('Hispanic', 'NonHispanicWhite', 'All')
for(compare in compare_list){
    outcome_string = compare[1]
    cohort_string = compare[2]
    count = 1
    results_df = data.frame('Membership' = character(), 'Outcome' = character(), 'Cohort' = character(), 'Adjust'= character(), 'ckd/medicine'= character(), 
                'Odds Ratio' = numeric(), '2.5%'= numeric(), '97.5%'= numeric(), 'ORCI' = character(),'Pvalue (Firth)'= numeric(), 
                             '# ckd/medicine in Outcome+'=integer(), '# Outcome+'=integer(),
                              '# ckd/medicine in Outcome-'=integer(), '# Outcome-'=integer(),
                             'Freq ckd/medicine in Outcome+'= numeric(), 'Freq ckd/medicine in Outcome-'= numeric(),
                             'Freq ckd/medicine in SEVERE'= numeric(),'Freq ckd/medicine in INPATIENT'= numeric(),
                            'Freq ckd/medicine in COVID'= numeric(),'Freq ckd/medicine in TESTED'= numeric(),
                             '# ckd/medicine in SEVERE'=integer(), '# SEVERE' = integer(),
                             '# ckd/medicine in INPATIENT'=integer(), '# INPATIENT'=integer(),
                             '# ckd/medicine in COVID'=integer(),'# COVID'=integer(),
                             '# ckd/medicine in TESTED'=integer(), '# TESTED'=integer(),
                             'age OR' = numeric(), 'age 2.5%' = numeric(),'age 97.5%' = numeric(), 'age pvalue' = numeric(),
                             'age2 OR' = numeric(), 'age2 2.5%' = numeric(),'age2 97.5%' = numeric(), 'age2 pvalue' = numeric(),
                             'SexFemale OR' = numeric(), 'SexFemale 2.5%' = numeric(),'SexFemale 97.5%' = numeric(), 'SexFemale pvalue' = numeric())
    
    for (race_ethnicity in race_ethnicity_group){
        results_df[count,1] = race_ethnicity
        sub_results_df = enrichment_test(risk_factors, outcome_string, 
                            cohort_string, 
                            df, 
                            time_string, 
                            membership_string,
                            adjust_string = 'adj_SA_krf', race_ethnicity)
        results_df[(count+1):(count + nrow(sub_results_df)),] = sub_results_df
        count = count + nrow(sub_results_df) + 1
        
        
    }
    colnames(results_df) = c('Membership', 'Outcome', 'Cohort', 'Adjust', 'ckd/medicine', 
                             'Odds Ratio', '2.5%', '97.5%', 'ORCI','Pvalue (Firth)', 
                             '# ckd/medicine in Outcome+', '# Outcome+',
                              '# ckd/medicine in Outcome-', '# Outcome-',
                             'Freq ckd/medicine in Outcome+', 'Freq ckd/medicine in Outcome-',
                             'Freq ckd/medicine in SEVERE','Freq ckd/medicine in INPATIENT','Freq ckd/medicine in COVID','Freq ckd/medicine in TESTED',
                             '# ckd/medicine in SEVERE', '# SEVERE',
                             '# ckd/medicine in INPATIENT', '# INPATIENT',
                             '# ckd/medicine in COVID','# COVID',
                             '# ckd/medicine in TESTED', '# TESTED',
                              'age OR' , 'age 2.5%' ,'age 97.5%' , 'age pvalue' ,
                             'age2 OR' , 'age2 2.5%' ,'age2 97.5%' , 'age2 pvalue' ,
                             'SexFemale OR' , 'SexFemale 2.5%' ,'SexFemale 97.5%' , 'SexFemale pvalue' )
    sheet_list[[paste0(outcome_string, 'in', cohort_string)]] = results_df
}

output_file = file.path(result_dir, sprintf('medication_krf_%s.xlsx', ANALYSIS_DATE))

write.xlsx(sheet_list, file = output_file) 

membership_string = 'NotNewToUCLA'
known_risk_factors= c('chf', 'diabetes', 'hyperlipidemia', 'hypertension', 'obesity', 'ckd', 'copd', 'chd')
compare_list = list( c('COVID', 'TESTED'), c('INPATIENT', 'COVID'), c('SEVERE', 'INPATIENT'))
risk_factors = c('ckd_any', 'ckd_1', 'ckd_2', 'ckd_3', 'ckd_4', 'ckd_5')
sheet_list = list()
race_ethnicity_group = c('Hispanic', 'NonHispanicWhite', 'All')
for(compare in compare_list){
    outcome_string = compare[1]
    cohort_string = compare[2]
    count = 1
    results_df = data.frame('Membership' = character(), 'Outcome' = character(), 'Cohort' = character(), 'Adjust'= character(), 'ckd/medicine'= character(), 
                'Odds Ratio' = numeric(), '2.5%'= numeric(), '97.5%'= numeric(), 'ORCI' = character(),'Pvalue (Firth)'= numeric(), 
                             '# ckd/medicine in Outcome+'=integer(), '# Outcome+'=integer(),
                              '# ckd/medicine in Outcome-'=integer(), '# Outcome-'=integer(),
                             'Freq ckd/medicine in Outcome+'= numeric(), 'Freq ckd/medicine in Outcome-'= numeric(),
                             'Freq ckd/medicine in SEVERE'= numeric(),'Freq ckd/medicine in INPATIENT'= numeric(),
                            'Freq ckd/medicine in COVID'= numeric(),'Freq ckd/medicine in TESTED'= numeric(),
                             '# ckd/medicine in SEVERE'=integer(), '# SEVERE' = integer(),
                             '# ckd/medicine in INPATIENT'=integer(), '# INPATIENT'=integer(),
                             '# ckd/medicine in COVID'=integer(),'# COVID'=integer(),
                             '# ckd/medicine in TESTED'=integer(), '# TESTED'=integer(),
                             'age OR' = numeric(), 'age 2.5%' = numeric(),'age 97.5%' = numeric(), 'age pvalue' = numeric(),
                             'age2 OR' = numeric(), 'age2 2.5%' = numeric(),'age2 97.5%' = numeric(), 'age2 pvalue' = numeric(),
                             'SexFemale OR' = numeric(), 'SexFemale 2.5%' = numeric(),'SexFemale 97.5%' = numeric(), 'SexFemale pvalue' = numeric())
    
    
    for (race_ethnicity in race_ethnicity_group){
        results_df[count,1] = race_ethnicity
        sub_results_df = enrichment_test(risk_factors, outcome_string, 
                            cohort_string, 
                            df, 
                            time_string, 
                            membership_string,
                            adjust_string = 'adj_SA_krf', race_ethnicity)
        results_df[(count+1):(count + nrow(sub_results_df)),] = sub_results_df
        count = count + nrow(sub_results_df) + 1
        
        
    }
    colnames(results_df) = c('Membership', 'Outcome', 'Cohort', 'Adjust', 'ckd/medicine', 
                             'Odds Ratio', '2.5%', '97.5%', 'ORCI','Pvalue (Firth)', 
                             '# ckd/medicine in Outcome+', '# Outcome+',
                              '# ckd/medicine in Outcome-', '# Outcome-',
                             'Freq ckd/medicine in Outcome+', 'Freq ckd/medicine in Outcome-',
                             'Freq ckd/medicine in SEVERE','Freq ckd/medicine in INPATIENT','Freq ckd/medicine in COVID','Freq ckd/medicine in TESTED',
                             '# ckd/medicine in SEVERE', '# SEVERE',
                             '# ckd/medicine in INPATIENT', '# INPATIENT',
                             '# ckd/medicine in COVID','# COVID',
                             '# ckd/medicine in TESTED', '# TESTED',
                              'age OR' , 'age 2.5%' ,'age 97.5%' , 'age pvalue' ,
                             'age2 OR' , 'age2 2.5%' ,'age2 97.5%' , 'age2 pvalue' ,
                             'SexFemale OR' , 'SexFemale 2.5%' ,'SexFemale 97.5%' , 'SexFemale pvalue' )
    sheet_list[[paste0(outcome_string, 'in', cohort_string)]] = results_df
}

output_file = file.path(result_dir, sprintf('ckd_krf_%s.xlsx', ANALYSIS_DATE))

write.xlsx(sheet_list, file = output_file) 
