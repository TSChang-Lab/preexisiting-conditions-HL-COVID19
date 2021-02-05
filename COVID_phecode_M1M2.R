rm(list = ls())
library('openxlsx')
library('tidyr', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
library('dplyr', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
library('logistf', lib = '\\\\netshares.ccnet.ucla.edu/ctx_profiles/_EEG_XDR/_Redirect/YiDing/Documents/R/win-library/3.6')
options(stringsAsFactors = FALSE)

PATIENT_DATE = '08312020'
DATA_CREATION_DATE = '09012020'
ANALYSIS_DATE  = format(Sys.time(), "%m%d%Y")


data_dir = 'J:\\OHIA_Investigators\\YiDing\\covid19\\Data'
result_dir = 'J:\\OHIA_Investigators\\YiDing\\covid19\\results'
dir.create(file.path(result_dir, PATIENT_DATE))


load_data <- function(data_dir, time_string, PATIENT_DATE, DATA_CREATION_DATE){
    data_file = file.path(data_dir, PATIENT_DATE, sprintf('patient_data_%s_%s.csv', time_string, DATA_CREATION_DATE))
    print(data_file)
    df = read.csv(data_file, row.names = 1)
    
    
    col_names = colnames(df)
    for (i in seq(length(col_names))){
        var = col_names[i]
        if (substring(var,1,1) == 'X'){
            var = sprintf( 'X%.2f', as.numeric(substring(var,2)))
            col_names[i] = var
        }
    }
    col_names
    colnames(df) = col_names

    return(df)
}

time_string = 'before'
membership_string = 'NotNewToUCLA'
known_risk_factors= c('chf', 'diabetes', 'hyperlipidemia', 'hypertension', 'obesity', 'ckd', 'copd', 'chd')

df = load_data(data_dir, time_string, PATIENT_DATE, DATA_CREATION_DATE)

enrichment_test <- function(phecode_list, known_risk_factors, df, 
                            outcome_string, cohort_string, time_string, 
                            membership_string,adjust_string = 'adj_SA', race_ethnicity){
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
        stop('no this!')
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
    
    if (cohort_string %in% c('DDRSAMPLE', 'ALL', 'ALLNEG')){
        covariates = c('AgeGroup', covariates)
        adjust_string = paste0(adjust_string, '_', 'MC')
    }

        
    # initiate progress bar 
    total = length(phecode_list)
    pb = txtProgressBar(min = 0, max = total, style = 3)
    # hold results for each risk factor
    results_df = data.frame('Membership' = character(), 'Outcome' = character(), 'Cohort' = character(), 'Adjust Type'= character(), 
                            'Phecode'= character(), 'Phecode Name'= character(), 'Phecode Category'= character(), 
                            'Odds Ratio' = numeric(), '2.5%'= numeric(), '97.5%'= numeric(), 'ORCI' = character(),'Pvalue'= numeric(), 
                            '# Phecode in Outcome+'=integer(), '# Outcome+'=integer(),
                            '# Phecode in Outcome-'=integer(), '# Outcome-'=integer(),
                            'Freq Phecode in Outcome+'= numeric(), 'Freq Phecode in Outcome-'= numeric(),
                            'Freq Phecode in SEVERE'= numeric(),'Freq Phecode in INPATIENT'= numeric(),
                            'Freq Phecode in COVID'= numeric(),'Freq Phecode in TESTED'= numeric(),
                            '# Phecode in SEVERE'=integer(), '# SEVERE' = integer(),
                            '# Phecode in INPATIENT'=integer(), '# INPATIENT'=integer(),
                            '# Phecode in COVID'=integer(),'# COVID'=integer(),
                            '# Phecode in TESTED'=integer(), '# TESTED'=integer(),
                             'age OR' = numeric(), 'age 2.5%' = numeric(),'age 97.5%' = numeric(), 'age pvalue' = numeric(),
                             'age2 OR' = numeric(), 'age2 2.5%' = numeric(),'age2 97.5%' = numeric(), 'age2 pvalue' = numeric(),
                             'SexFemale OR' = numeric(), 'SexFemale 2.5%' = numeric(),'SexFemale 97.5%' = numeric(), 'SexFemale pvalue' = numeric())
    
    count = 0
    for (i in seq(total)){
        # update progress bar 
        setTxtProgressBar(pb, i)
        # run logistic regression
        target_phecode = phecode_list[i]
        target_phecode_name = phecodes_map[target_phecode, 'Phecode.Name']
        target_phecode_category = phecodes_map[target_phecode, 'Phecode.Category']
        if (sum(analysis_df[,target_phecode]) ==0) next
        
        model_covariates = c(covariates, target_phecode)
        f = as.formula( paste0(outcome_string, '~' ,paste(model_covariates, collapse = ' + '), sep = ' '))
        model = logistf(formula = f, data = analysis_df)
        # get odds ratio, ci and pvalue 
        p = model$prob[target_phecode]
        oddsratio = round(exp(coef(model)[target_phecode]),2)
        ci_lower = round(exp(confint(model)[target_phecode,1]),2)
        ci_upper = round(exp(confint(model)[target_phecode,2]),2)
        orci_string = sprintf('%s [%s, %s]', oddsratio, ci_lower, ci_upper)
        count = count + 1

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
        n_pos_phecode = sum(analysis_df[target_phecode] * analysis_df[outcome_string])
        n_neg_phecode = sum(analysis_df[target_phecode] * (1-analysis_df[outcome_string]))
        n_pos = sum(analysis_df[outcome_string])
        n_neg = sum((1-analysis_df[outcome_string]))
        freq_pos_phecode = round(n_pos_phecode/n_pos*100, 2)
        freq_neg_phecode = round(n_neg_phecode/n_neg*100, 2)
        
        # frequency in tested
        n_TESTED = sum(df['TESTED'])
        n_TESTED_phecode = sum(df['TESTED']*df[target_phecode])
        freq_TESTED_phecode = round(n_TESTED_phecode/n_TESTED*100,2)
        # frequency in covid
        n_COVID = sum(df['COVID'])
        n_COVID_phecode = sum(df['COVID']*df[target_phecode])
        freq_COVID_phecode = round(n_COVID_phecode/n_COVID*100,2)
        # frequency in inpatient
        n_INPATIENT = sum(df['INPATIENT'])
        n_INPATIENT_phecode = sum(df['INPATIENT']*df[target_phecode])
        freq_INPATIENT_phecode = round(n_INPATIENT_phecode/n_INPATIENT*100,2)
        # frequency in severe
        n_SEVERE = sum(df['SEVERE'])
        n_SEVERE_phecode = sum(df['SEVERE']*df[target_phecode])
        freq_SEVERE_phecode = round(n_SEVERE_phecode/n_SEVERE*100,2)
#         # frequency in DDRSAMPLE
#         n_DDRSAMPLE = sum(df['DDRSAMPLE'])
#         n_DDRSAMPLE_phecode = sum(df['DDRSAMPLE']*df[target_phecode])
#         freq_DDRSAMPLE_phecode = round(n_DDRSAMPLE_phecode/n_DDRSAMPLE*100,2)
        
        result_row= list(membership_string, outcome_string, cohort_string, adjust_string, 
                         target_phecode, target_phecode_name, target_phecode_category,  
                      oddsratio, ci_lower, ci_upper, orci_string, p,
                      n_pos_phecode, n_pos, 
                      n_neg_phecode, n_neg,
                      freq_pos_phecode, freq_neg_phecode, 
                      freq_SEVERE_phecode, freq_INPATIENT_phecode, freq_COVID_phecode, freq_TESTED_phecode, #freq_DDRSAMPLE_phecode,
                      n_SEVERE_phecode, n_SEVERE,
                      n_INPATIENT_phecode, n_INPATIENT,
                      n_COVID_phecode, n_COVID,
                      n_TESTED_phecode, n_TESTED,
                      age_coef, age_lower, age_upper, age_p, 
                      age2_coef, age2_lower, age2_upper, age2_p,
                      SexFemale_coef, SexFemale_lower, SexFemale_upper, SexFemale_p)
#                         n_DDRSAMPLE_phecode, n_DDRSAMPLE)
        results_df[count,] = result_row
    } 
  
    results_df = data.frame(results_df)
    results_df = results_df[order(results_df['Pvalue']),]
    results_df[ 'Nominal Sig.'] = (results_df['Pvalue'] <0.05)*1
    results_df[ 'Bonferr Sig.'] = (results_df['Pvalue'] <0.05/nrow(results_df))*1
    colnames(results_df) = c('Membership' , 'Outcome', 'Cohort', 'Adjust Type', 
            'Phecode', 'Phecode Name', 'Phecode Category', 
            'Odds Ratio' , '2.5%', '97.5%','ORCI','Pvalue (Firth)', 
            '# Phecode in Outcome+', '# Outcome+',
            '# Phecode in Outcome-', '# Outcome-',
            'Freq Phecode in Outcome+', 'Freq Phecode in Outcome-',
            'Freq Phecode in SEVERE','Freq Phecode in INPATIENT',
            'Freq Phecode in COVID','Freq Phecode in TESTED',
            '# Phecode in SEVERE', '# SEVERE',
            '# Phecode in INPATIENT', '# INPATIENT',
            '# Phecode in COVID','# COVID',
            '# Phecode in TESTED', '# TESTED', 
             'age OR' , 'age 2.5%' ,'age 97.5%' , 'age pvalue' ,
            'age2 OR' , 'age2 2.5%' ,'age2 97.5%' , 'age2 pvalue' ,
            'SexFemale OR' , 'SexFemale 2.5%' ,'SexFemale 97.5%' , 'SexFemale pvalue',
            'Nominal Sig.','Bonferr Sig.' )
    return(results_df)

}


phecodes_map = read.csv(file.path(data_dir, 'phecodes_map.csv'))
phecodes_map = data.frame(phecodes_map)
rownames(phecodes_map) = sapply(phecodes_map['Phecode'], function(x) sprintf('X%.2f', x))      
file.path(data_dir, 'phecodes_map.csv')

time_string = 'before'
membership_string = 'NotNewToUCLA'
known_risk_factors= c('chf', 'diabetes', 'hyperlipidemia', 'hypertension', 'obesity', 'ckd', 'copd', 'chd')

df = load_data(data_dir, time_string, PATIENT_DATE, DATA_CREATION_DATE)

compare_list = list(c('SEVERE', 'INPATIENT'),c('INPATIENT', 'COVID'),c('COVID', 'TESTED'))
sheet_list = list()

race_ethnicity_group = c('Hispanic','NonHispanicWhite')

for (race_ethnicity in race_ethnicity_group){
     writeLines(sprintf('\n--------%s--------\n', race_ethnicity))
    for (compare in compare_list){
        outcome_string = compare[1]
        cohort_string = compare[2]
        
        file_name = sprintf('covid_phecode_firth_logistic_%sin%s_%s_M1M2_%s_%s.xlsx', 
                            outcome_string, cohort_string, time_string, ANALYSIS_DATE, race_ethnicity)

        output_file = file.path(result_dir, PATIENT_DATE, file_name)
        
        
        sheet_list = list()

        # run M1
        writeLines(sprintf('\n--------M1: %sin%s--------\n', outcome_string, cohort_string))
        M1_phecodes = colnames(df) [substring(colnames(df), 1,1) == 'X']
        M1_results = enrichment_test(M1_phecodes, known_risk_factors, df, 
                                    outcome_string, cohort_string, time_string, 
                                    membership_string,adjust_string = 'adj_SA', race_ethnicity)
        sheet_list[['M1']] = M1_results

        # save results
        write.xlsx(sheet_list, file = output_file) 

        # run M2
        writeLines(sprintf('\n--------M2: %sin%s--------\n', outcome_string, cohort_string))
        M2_phecodes = M1_results[M1_results['Bonferr Sig.'] == 1, 'Phecode']
        if(length(M2_phecodes)>0){
            M2_results = enrichment_test(M2_phecodes, known_risk_factors, df, 
                                    outcome_string, cohort_string, time_string, 
                                    membership_string,adjust_string = 'adj_SA_krf', race_ethnicity)

        }else{
            M2_results = data.frame(matrix(0,2,2))
        }
        sheet_list[['M2']] = M2_results

        # save results
        write.xlsx(sheet_list, file = output_file) 
    }
}
                    


