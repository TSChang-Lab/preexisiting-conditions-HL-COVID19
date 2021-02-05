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

known_risk_factors= c('chf', 'diabetes', 'hyperlipidemia', 'hypertension', 'obesity', 'ckd', 'copd', 'chd')

load_data <- function(data_dir, time_string, PATIENT_DATE, DATA_CREATION_DATE){
    data_file = file.path(data_dir, PATIENT_DATE, sprintf('patient_data_%s_%s.csv', time_string, DATA_CREATION_DATE))
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
df = load_data(data_dir, time_string, PATIENT_DATE, DATA_CREATION_DATE)

# get hispanic and nonhispanic white only
condition = (df['SIRE_Ethnicity'] == 'Hispanic or Latino')
condition = condition | ((df['SIRE_Ethnicity'] == 'Not Hispanic or Latino') & (df[,'SIRE_Race'] == 'White or Caucasian'))
df = df[condition,]


# covariates_list = list(baseline = c('HispanicORLatino', 'Sex', 'age', 'age2'), 
#                        diabetes = c('HispanicORLatino', 'Sex', 'age', 'age2', 'diabetes'),
#                        ckd =      c('HispanicORLatino', 'Sex', 'age', 'age2', 'ckd'),
#                        both =     c('HispanicORLatino', 'Sex', 'age', 'age2', 'ckd', 'diabetes'))

covariates_list = list(baseline = c('HispanicORLatino', 'Sex', 'age', 'age2'), 
                       krf = c('HispanicORLatino', 'Sex', 'age', 'age2', known_risk_factors))

compare_list = list(c('COVID', 'TESTED'), c('INPATIENT', 'COVID'), c('SEVERE', 'INPATIENT'))

membership = 'NotNewToUCLA'
results_df = data.frame( 'Membership' = character(), 'Outcome' = character(), 'Cohort' = character(), 'Model' = character(), 
                         'His. Odds Ratio' =numeric(), '2.5%' = numeric(), '97.5%' = numeric(), 'ORCI' = character(), 'Pvalue (Firth)' = numeric(), 
                         'Freq His. in Outcome+' = numeric(), 'Freq His. in Outcome-' = numeric(),
                         '# His. in Outcome+' = integer(), '# NHW in Outcome+'=integer(), '# Outcome+' = integer(),
                         '# His. in Outcome-' = integer(), '# NHW in Outcome-'=integer(),'# Outcome-' = integer(),
                             'age OR' = numeric(), 'age 2.5%' = numeric(),'age 97.5%' = numeric(), 'age pvalue' = numeric(),
                             'age2 OR' = numeric(), 'age2 2.5%' = numeric(),'age2 97.5%' = numeric(), 'age2 pvalue' = numeric(),
                             'SexFemale OR' = numeric(), 'SexFemale 2.5%' = numeric(),'SexFemale 97.5%' = numeric(), 'SexFemale pvalue' = numeric())
    

count = 0
for(compare in compare_list){
    outcome_string = compare[1]
    cohort_string = compare[2]
    
    condition =  (df['NewToUCLA'] == 0) & ((df[outcome_string ] == 1) | (df[cohort_string] == 1))
    analysis_df = df[condition,]
    
    for (model_name in names(covariates_list)){
        covariates = covariates_list[[model_name]]
        if (cohort_string == 'DDRSAMPLE'){
            covariates = c('AgeGroup', covariates)
        }
        f = as.formula( paste0(outcome_string, '~' ,paste(covariates, collapse = ' + '), sep = ' '))
        model = logistf(formula = f, data = analysis_df)
        
        p = model$prob['HispanicORLatino']
        oddsratio = round(exp(coef(model)['HispanicORLatino']),2)
        ci_lower = round(exp(confint(model)['HispanicORLatino',1]),2)
        ci_upper = round(exp(confint(model)['HispanicORLatino',2]),2)
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
        

        
        n_pos_his = sum(analysis_df['HispanicORLatino'] * analysis_df[outcome_string])
        n_neg_his = sum(analysis_df['HispanicORLatino'] * (1-analysis_df[outcome_string]))
        
        n_pos = sum(analysis_df[outcome_string])
        n_neg = sum((1-analysis_df[outcome_string]))
        
        freq_pos_his = round(n_pos_his/n_pos, 2)
        freq_neg_his = round(n_neg_his/n_neg, 2)

        n_pos_nhw = sum( (1-analysis_df['HispanicORLatino']) * analysis_df[outcome_string])
        n_neg_nhw = sum( (1-analysis_df['HispanicORLatino']) * (1-analysis_df[outcome_string]))
        result_row= list(membership, outcome_string, cohort_string, model_name, oddsratio, ci_lower, ci_upper, orci_string, p,
                      freq_pos_his, freq_neg_his, n_pos_his, n_pos_nhw, n_pos, n_neg_his,n_neg_nhw, n_neg,
                      age_coef, age_lower, age_upper, age_p, 
                      age2_coef, age2_lower, age2_upper, age2_p,
                      SexFemale_coef, SexFemale_lower, SexFemale_upper, SexFemale_p)
        count = count + 1
        results_df[count,] = result_row
    } 
    count = count + 1
    results_df[count,] = NA
}

colnames(results_df) = c('Membership', 'Outcome', 'Cohort', 'Model', 
                         'His. Odds Ratio', '2.5%', '97.5%', 'ORCI', 'Pvalue (Firth)', 
                         'Freq His. in Outcome+', 'Freq His. in Outcome-',
                         '# His. in Outcome+', '# NHW in Outcome+', '# Outcome+',
                         '# His. in Outcome-', '# NHW in Outcome-', '# Outcome-',
            'age OR' , 'age 2.5%' ,'age 97.5%' , 'age pvalue' ,
            'age2 OR' , 'age2 2.5%' ,'age2 97.5%' , 'age2 pvalue' ,
            'SexFemale OR' , 'SexFemale 2.5%' ,'SexFemale 97.5%' , 'SexFemale pvalue')

    # save results
    file_name = sprintf('covid_Hispanic_%s.xlsx', ANALYSIS_DATE)

    output_file = file.path(result_dir, PATIENT_DATE, file_name)
    write.xlsx(results_df, file = output_file) 




