import pandas as pd
import csv
import re
import matplotlib.pyplot as plt


# This function opens a .csv input_IPDA_file that has been exported from the BioRad ddPCR analyzer
# and exports an excel table that states which samples passed/failed the data quality control tests
def IPDA_quality_control(input_IPDA_file, minimum_required_droplets):

    # create empty lists for each column of your QC table that you need to use for analysis
    wells = []
    samples = []
    targets = []
    concentrations = []
    double_pos = []
    Ch1_pos = []
    Ch2_pos = []
    double_neg = []
    AcceptedDroplets = []
    Accepted_droplet_test = []
    below_thirty_percent_pos_test_passed = []

    # open the input_IPDA_file and add data of interest to lists
    file_handle = open(input_IPDA_file, 'r')
    csv_reader = csv.reader(file_handle, delimiter=',')
    
    # to skip first row which are just titles the print statement prevents a linter error
    header_row = next(csv_reader)
    header_row

    for row in csv_reader:
        
        # this is based on the regular output from the ddPCR analyzer
        wells.append(row[0])
        samples.append(row[3])
        targets.append(row[5])
        concentrations.append(row[7])

        # the droplet channel results and droplet counts will be needed multiple times
        double_positives = int(row[16])
        single_channel1 = int(row[17])
        single_channel2 = int(row[18])
        double_pos.append(double_positives)
        negatives = int(row[19])
        num_droplets = int(row[21])
        positive_droplets = single_channel1 + single_channel2 + double_positives

        Ch1_pos.append(single_channel1)
        Ch2_pos.append(single_channel2)
        double_neg.append(negatives)
        AcceptedDroplets.append(num_droplets)

        if (num_droplets >= minimum_required_droplets):
            Accepted_droplet_test.append('passed')
        else:
            Accepted_droplet_test.append('failed')

        if positive_droplets <= 0.3 * (positive_droplets * negatives):
            below_thirty_percent_pos_test_passed.append('passed')
        else:
            below_thirty_percent_pos_test_passed.append('failed')

    file_handle.close()

    # create data frame that will store all data from the quality control
    quality_control = pd.DataFrame({
        'Well': wells,
        'Sample': samples,
        'Target': targets,
        'Concentration': concentrations,
        'Ch1+Ch2+': double_pos,
        'Ch1+Ch2-': Ch1_pos,
        'Ch1-Ch2+': Ch2_pos,
        'Ch1-Ch2-': double_neg,
        'Number of Droplets': AcceptedDroplets,
        'More than 10,000 droplets?': Accepted_droplet_test,
        'Below 30% positives?': below_thirty_percent_pos_test_passed})
    
    # for better visibility and easier continuation of analysis, sort data
    # so that the same samples and Targets will stick together
    quality_control = quality_control.sort_values(by=['Sample', 'Target'], axis=0)

    # sometimes, people export not just the single well data
    # but select "Both" single and merged data. So in order to make our program compatible with either file input
    # we will exclude any merged data if present. Such data has a 'Well' starting with 'M'
    quality_control = quality_control.loc[~quality_control['Well'].str.contains('M')]

    # export to new Excel file
    quality_control.to_excel('output_files/Analyzed_IPDA_data.xlsx', sheet_name='Quality Control', index=False)

    return quality_control


# this is where the second function starts. Write an exclude_outliers function that you can call over and over.


# This function accesses a data frame that contains the quality control of an input_IPDA_file
# and copies only the samples that passed the data quality control tests
# it then performs the complete IPDA analysis
# by correcting the IPDA RPP30 concentration for the actual concentration that the HIV reaction was run at
# and calculating the % of unsheared DNA
# using this information, Gag and Env HIV reactions will be normalized
# Finally Env and Gag data will be combined to calculate the % intact HIV
# and the HIV copies per million cells
def IPDA_normalized_to_housekeeping_gene(input_IPDA_file, DNA_concentrations_used_input_file, minimum_required_droplets):

    # use the cleaned and quality controlled dataframe and select only data that passed QC tests, then remove the QC test data
    analyzed_IPDA_data = IPDA_quality_control(input_IPDA_file, minimum_required_droplets)
    analyzed_IPDA_data = analyzed_IPDA_data[analyzed_IPDA_data['More than 10,000 droplets?'] == 'passed']
    # analyzed_IPDA_data = analyzed_IPDA_data[~analyzed_IPDA_data['More than 10,000 droplets?']]
    analyzed_IPDA_data = analyzed_IPDA_data[analyzed_IPDA_data['Below 30% positives?'] == 'passed']
    # analyzed_IPDA_data = analyzed_IPDA_data[~analyzed_IPDA_data['Below 30% positives?']]
    # analyzed_IPDA_data = analyzed_IPDA_data[~analyzed_IPDA_data['Below 30% positives?']]
    analyzed_IPDA_data = analyzed_IPDA_data.drop(['Number of Droplets', 'More than 10,000 droplets?', 'Below 30% positives?'], axis=1)
    # somehow, 'Concentration' is an object instead of a number, so we need to convert it to float
    # sometimes the ddPCR Machine outputs the text "No Call" if the concentration is too low.
    # "No Call" cannot be converted to float, that's why we turn this into a 0
    analyzed_IPDA_data['Concentration'].replace('No Call', '0', inplace=True)
    analyzed_IPDA_data['Concentration'] = analyzed_IPDA_data['Concentration'].astype(float)
    # split dataframe into 4 new dataframes, one for each Target used
    rpp30_fam_df = analyzed_IPDA_data.loc[analyzed_IPDA_data['Target'].str.contains('Rpp30 Fam|Rpp30 Shear', case=False, flags=re.IGNORECASE, regex=True)]
    rpp30_vic_df = analyzed_IPDA_data.loc[analyzed_IPDA_data['Target'].str.contains('vic', case=False, flags=re.IGNORECASE, regex=True)]
    gag_df = analyzed_IPDA_data.loc[analyzed_IPDA_data['Target'].str.contains('Psi|Gag', case=False, flags=re.IGNORECASE, regex=True)]
    env_df = analyzed_IPDA_data.loc[analyzed_IPDA_data['Target'].str.contains('env', case=False, flags=re.IGNORECASE, regex=True)]
    # correct RPP30 concentration for HIV concentration
    DNA_excel_file = pd.ExcelFile(DNA_concentrations_used_input_file)
    DNA_concentration_dataframe = DNA_excel_file.parse('Sheet1')

    # we merge the dataframe of the DNAconcentration_input_file
    rpp30_fam_df = pd.merge(rpp30_fam_df, DNA_concentration_dataframe, on='Sample')

    # Now we can correct the concentration values from the old rpp30_fam_df
    # using the two new columns that we got from the DNA_concentration_dataframe
    actual_fam_conc = (rpp30_fam_df['Concentration'] / rpp30_fam_df['DNA conc I used [ng/µL] for RPP30'])
    dilution_factor = actual_fam_conc * rpp30_fam_df['DNA conc I used [ng/µL] for HIV Gag Env reactions']
    rpp30_fam_df['Corrected concentration [ng/µL]'] = dilution_factor
    
    # perform the same DNA correction steps on the rpp30_vic_df
    rpp30_vic_df = pd.merge(rpp30_vic_df, DNA_concentration_dataframe, on='Sample')
    actual_vic_conc = (rpp30_vic_df['Concentration'] / rpp30_vic_df['DNA conc I used [ng/µL] for RPP30'])
    rpp30_vic_df['Corrected concentration [ng/µL]'] = actual_vic_conc * rpp30_vic_df['DNA conc I used [ng/µL] for HIV Gag Env reactions']
    
    # drop the columns that you don't need anymore
    rpp30_fam_df = rpp30_fam_df.drop(['DNA conc I used [ng/µL] for RPP30', 'DNA conc I used [ng/µL] for HIV Gag Env reactions'], axis=1)
    rpp30_vic_df = rpp30_vic_df.drop(['DNA conc I used [ng/µL] for RPP30', 'DNA conc I used [ng/µL] for HIV Gag Env reactions'], axis=1)

    # now using the corrected concentrations, we will combine the data of the rpp30_fam_df and  rpp30_vic_df
    # to find the average concentrations, as well as the %unsheared DNA
    # which is the fraction of double positive droplets
    # for this, we first need to take the mean of the corrected concentrations of each replicate
    # to combine the data of all replicates for one given sample
    rpp30_combined_df = pd.merge(rpp30_fam_df, rpp30_vic_df, on=['Well', 'Sample'], suffixes=['_fam', '_vic'])
    dataframe_for_rpp30_means = rpp30_combined_df.groupby(['Sample'], as_index=False).mean(numeric_only=True)
    dataframe_for_rpp30_means = dataframe_for_rpp30_means[['Sample', 'Corrected concentration [ng/µL]_fam', 'Corrected concentration [ng/µL]_vic']]
    dataframe_for_rpp30_stdevs = rpp30_combined_df.groupby(['Sample'], as_index=False).std(numeric_only=True)
    dataframe_for_rpp30_stdevs = dataframe_for_rpp30_stdevs[['Sample', 'Corrected concentration [ng/µL]_fam', 'Corrected concentration [ng/µL]_vic']]

    # now that we have the average and stdev for each corrected concentration, let's rename the columns
    # and then merge all three dataframes so that we have single replicate data, means, and stdevs in one dataframe
    dataframe_for_rpp30_means = dataframe_for_rpp30_means.rename(columns={'Corrected concentration [ng/µL]_fam': 'fam_mean'})
    dataframe_for_rpp30_means = dataframe_for_rpp30_means.rename(columns={'Corrected concentration [ng/µL]_vic': 'vic_mean'})
    dataframe_for_rpp30_stdevs = dataframe_for_rpp30_stdevs.rename(columns={'Corrected concentration [ng/µL]_fam': 'fam_stdev'})
    dataframe_for_rpp30_stdevs = dataframe_for_rpp30_stdevs.rename(columns={'Corrected concentration [ng/µL]_vic': 'vic_stdev'})
    rpp30_combined_df = pd.merge(rpp30_combined_df, dataframe_for_rpp30_means, on=['Sample'])
    rpp30_combined_df = pd.merge(rpp30_combined_df, dataframe_for_rpp30_stdevs, on=['Sample'])

    # For each Sample, check for outliers where the 'Corrected concentration [ng/µL]' is outside of mean±(2*stDev) and exclude those
    # test for each row if it's outside the mean±2*StDev range for that sample
    fam_lower = (rpp30_combined_df['fam_mean'] - (rpp30_combined_df['fam_stdev']) * 2)
    fam_upper = (rpp30_combined_df['fam_mean'] + (rpp30_combined_df['fam_stdev']) * 2)
    fam_test1 = rpp30_combined_df['Corrected concentration [ng/µL]_fam'] >= fam_lower
    fam_test2 = rpp30_combined_df['Corrected concentration [ng/µL]_fam'] <= fam_upper
    vic_lower = (rpp30_combined_df['vic_mean'] - (rpp30_combined_df['vic_stdev']) * 2)
    vic_upper = (rpp30_combined_df['vic_mean'] + (rpp30_combined_df['vic_stdev']) * 2)
    vic_test1 = rpp30_combined_df['Corrected concentration [ng/µL]_vic'] >= vic_lower
    vic_test2 = rpp30_combined_df['Corrected concentration [ng/µL]_vic'] <= vic_upper
    rpp30_combined_df['Outlier_test passed'] = fam_test1 & fam_test2 & vic_test1 & vic_test2
    rpp30_combined_df['Outlier_test passed'] = rpp30_combined_df['Outlier_test passed'].map({True: 'passed', False: 'failed'})
    rpp30_combined_df = rpp30_combined_df.drop(['Ch1+Ch2+_vic', 'Ch1+Ch2-_vic', 'Ch1-Ch2+_vic', 'Ch1-Ch2-_vic'], axis=1)
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Ch1+Ch2+_fam': 'Ch1+Ch2+'})
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Ch1+Ch2-_fam': 'Ch1+Ch2-'})
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Ch1-Ch2+_fam': 'Ch1-Ch2+'})
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Ch1-Ch2-_fam': 'Ch1-Ch2-'})
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Corrected concentration [ng/µL]_fam': 'Corr conc fam'})
    rpp30_combined_df = rpp30_combined_df.rename(columns={'Corrected concentration [ng/µL]_vic': 'Corr conc vic'})

    # now that you have excluded any outliers, continue to create the means table again
    no_outliers = rpp30_combined_df[rpp30_combined_df['Outlier_test passed'] == 'passed']

    no_outliers = no_outliers[['Sample', 'Ch1+Ch2+', 'Ch1+Ch2-', 'Ch1-Ch2+', 'Corr conc fam', 'Corr conc vic']]
    dataframe_for_rpp30_means = no_outliers.groupby('Sample', as_index=False).mean(numeric_only=True)
    
    average_shear_vic = (dataframe_for_rpp30_means['Corr conc fam'] + dataframe_for_rpp30_means['Corr conc vic']) / 2
    dataframe_for_rpp30_means['Average RPP30'] = average_shear_vic
    d_pos_vic = dataframe_for_rpp30_means['Ch1+Ch2+']

    # note, this step differs slightly from how we usually manually calculate the '%Unsheared averaged'
    # usually, we do this for each sample individually rather than on the average droplet counts across channels.
    # the values are very similar, however, and do not affect the final result
    unsheared = d_pos_vic / (d_pos_vic + ((dataframe_for_rpp30_means['Ch1+Ch2-'] + dataframe_for_rpp30_means['Ch1-Ch2+']) / 2))
    dataframe_for_rpp30_means['%Unsheared averaged'] = unsheared

    # now that we have done all the calculations that required us to work on each sample individually
    # we can add the two new columns to our original rpp30_combined_df that shows all technical replicates
    # also add the concentration from the env table to gag for later calculations
    env_df = env_df[['Well', 'Sample', 'Concentration']]
    env_df = env_df.rename(columns={'Concentration': 'Env conc'})
    gag_df = gag_df.rename(columns={'Concentration': 'Gag conc'})
    gag_df = pd.merge(gag_df, env_df, on=['Well', 'Sample'])
    dataframe_for_rpp30_means = dataframe_for_rpp30_means[['Sample', 'Average RPP30', '%Unsheared averaged']]
    gag_df = pd.merge(gag_df, dataframe_for_rpp30_means, on=['Sample'])
    
    # we calculate the 3'Deleted copies per million CD4+ T cells
    # by multiplying the concentration by 1 million
    # and correcting for the actual number of copies we expect per cell compared to copies of RPP30
    # which is half of the Average RPP30 (since human genes have 2 alleles on paternal and maternal chromosomes per cell)
    # and HIV integrates into only one chromosome per cell.
    totl_gag_conc = gag_df['Gag conc']
    totl_env_conc = gag_df['Env conc']
    expected_copies = (gag_df['Average RPP30'] / 2)
    gag_df["3'Deleted/hypermutated/M (Gag/M)"] = (totl_gag_conc * 1000000) / expected_copies
    gag_df["5'Deleted/M (Env/M)"] = (totl_env_conc * 1000000) / expected_copies
    
    # at this stage, we are going to perform an outlier exclusion step
    # to exclude any rows in which gag_df["3'Deleted/hypermutated/M (Gag/M)"] or gag_df["5'Deleted/M (Env/M)"]
    # are more than 2 standard deviations away from the mean
    dataframe_HIV_means = gag_df.groupby(['Sample'], as_index=False).mean(numeric_only=True)
    dataframe_HIV_means = dataframe_HIV_means[['Sample', "3'Deleted/hypermutated/M (Gag/M)", "5'Deleted/M (Env/M)"]]
    dataframe_HIV_means = dataframe_HIV_means.rename(columns={"3'Deleted/hypermutated/M (Gag/M)": 'Gag/M_mean'})
    dataframe_HIV_means = dataframe_HIV_means.rename(columns={"5'Deleted/M (Env/M)": 'Env/M_mean'})
    dataframe_HIV_stdev = gag_df.groupby(['Sample'], as_index=False).std(numeric_only=True)
    dataframe_HIV_stdev = dataframe_HIV_stdev[['Sample', "3'Deleted/hypermutated/M (Gag/M)", "5'Deleted/M (Env/M)"]]
    dataframe_HIV_stdev = dataframe_HIV_stdev.rename(columns={"3'Deleted/hypermutated/M (Gag/M)": 'Gag/M_stdev'})
    dataframe_HIV_stdev = dataframe_HIV_stdev.rename(columns={"5'Deleted/M (Env/M)": 'Env/M_stdev'})

    gag_df = pd.merge(gag_df, dataframe_HIV_means, on=['Sample'])
    gag_df = pd.merge(gag_df, dataframe_HIV_stdev, on=['Sample'])

    # For each Sample, check for outliers where the '3'Deleted/hypermutated/M (Gag/M)' or the  "5'Deleted/M (Env/M)"
    # are outside of mean±(2*stDev) and exclude those
    gag_lower = (gag_df['Gag/M_mean'] - (gag_df['Gag/M_stdev']) * 2)
    gag_upper = (gag_df['Gag/M_mean'] + (gag_df['Gag/M_stdev']) * 2)
    gag_test1 = gag_df["3'Deleted/hypermutated/M (Gag/M)"] >= gag_lower
    gag_test2 = gag_df["3'Deleted/hypermutated/M (Gag/M)"] <= gag_upper
    env_lower = (gag_df['Env/M_mean'] - (gag_df['Env/M_stdev']) * 2)
    env_upper = (gag_df['Env/M_mean'] + (gag_df['Env/M_stdev']) * 2)
    env_test1 = gag_df["5'Deleted/M (Env/M)"] >= env_lower
    env_test2 = gag_df["5'Deleted/M (Env/M)"] <= env_upper
    gag_df['Outlier_test passed'] = gag_test1 & gag_test2 & env_test1 & env_test2
    gag_df['Outlier_test passed'] = gag_df['Outlier_test passed'].map({True: 'passed', False: 'failed'})
    no_HIV_outliers = gag_df[gag_df['Outlier_test passed'] == 'passed']

    # to later calculate the % and number of intact HIV, per million cells, we will need to know
    # the concentration copies of intact Gag (currently, 'Concentration' is the sum of the concentrations of intact as well as Gag single positives)
    d_pos_HIV = no_HIV_outliers['Ch1+Ch2+']
    s_pos_gag = no_HIV_outliers['Ch1+Ch2-']
    s_pos_env = no_HIV_outliers['Ch1-Ch2+']
    intact_gag_conc = d_pos_HIV * totl_gag_conc / (d_pos_HIV + s_pos_gag)
    intact_env_conc = d_pos_HIV * totl_env_conc / (d_pos_HIV + s_pos_env)

    # this is the concentration of intact provirus as determined using the totl_gag_conc
    no_HIV_outliers['D+ accord. to Gag rxn'] = intact_gag_conc
    no_HIV_outliers['D+ accord. to Env rxn'] = intact_env_conc
    intact_conc = (intact_gag_conc + intact_env_conc) / 2
    no_HIV_outliers['intact concentration'] = intact_conc
    no_HIV_outliers['intact concentration/M'] = (intact_conc * 1000000) / expected_copies

    # correct the intact concentration for the % of DNA shearing
    no_HIV_outliers['Intact/M'] = no_HIV_outliers['intact concentration/M'] / no_HIV_outliers['%Unsheared averaged']
    
    # calculate % intact
    total_HIV = no_HIV_outliers['Intact/M'] + no_HIV_outliers["3'Deleted/hypermutated/M (Gag/M)"] + no_HIV_outliers["5'Deleted/M (Env/M)"]
    intact_per_total_copies = no_HIV_outliers['Intact/M'] / total_HIV
    no_HIV_outliers['Intact [%]'] = intact_per_total_copies * 100

    # rename unclear column names and .drop columns that no longer make sense before you export the table
    no_HIV_outliers = no_HIV_outliers.rename(columns={'intact concentration': 'non-corr. intact conc.'})
    no_HIV_outliers = no_HIV_outliers.drop(['Well', 'Target'], axis=1)

    # save new dataframes to Excel file
    with pd.ExcelWriter('output_files/Analyzed_IPDA_data.xlsx', engine='openpyxl', mode='a') as writer:
        rpp30_combined_df.to_excel(writer, sheet_name='RPP30 Analysis', index=False)
        gag_df.to_excel(writer, sheet_name='HIV Analysis', index=False)

    return no_HIV_outliers.copy()


# since we usually only care about the actual results of the analysis,
# we'll export those results as a separate sheet, ready to be copy-pasted into GraphPad Prism
# and we also plot the data as .png files here
def export_analyzed_IPDA_as_Excel_and_png(input_IPDA_file, DNA_concentrations_used_input_file, minimum_required_droplets):

    # use analyzed dataframe only which is what the IPDA_normalized_to_housekeeping_gene function returns
    summary_data_to_be_exported = IPDA_normalized_to_housekeeping_gene(input_IPDA_file, DNA_concentrations_used_input_file, minimum_required_droplets)

    # select only the 'Sample', "3'Deleted/hypermutated/M (Gag/M)", "5'Deleted/M (Env/M), 'Corrected Intact Concentration/M'
    summary_data_to_be_exported = summary_data_to_be_exported[['Sample', "3'Deleted/hypermutated/M (Gag/M)", "5'Deleted/M (Env/M)", 'Intact/M', 'Intact [%]']]

    # save to Excel sheet two next to Quality Control
    with pd.ExcelWriter('output_files/Analyzed_IPDA_data.xlsx', engine='openpyxl',
                        mode='a') as writer:
        summary_data_to_be_exported.to_excel(writer, sheet_name='Summary', index=False)
    
    grouped_df = summary_data_to_be_exported.groupby('Sample', as_index=False).sum()
    no_percent_intact = grouped_df.drop(['Intact [%]'], axis=1)
    percent_intact = grouped_df[['Sample', 'Intact [%]']]

    # create bar plots
    no_percent_intact.plot.bar(x='Sample', rot=45, title='Copies per million CD4+ T cells')
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.4, top=0.9)
    plt.xlabel('Sample', labelpad=15)
    plt.savefig('output_files/copies_per_million.png')
    percent_intact.plot.bar(x='Sample', rot=45, title='% Intact HIV')
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.4, top=0.9)
    plt.xlabel('Sample', labelpad=15)
    plt.savefig('output_files/percent_intact.png')
    
    return summary_data_to_be_exported
