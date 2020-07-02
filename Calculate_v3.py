#!/usr/bin/env python

import numpy as np
import argparse, sys

"""
This program will take in log.txt files as outputted by the plotter
and calculate the Scale factor, efficiency, and purity, along with
propagating systematic errors.
"""

def sum_of_squares(list):
    """Performs a simple sum of squares calculation, used for
    error propagation.
    """
    total=0
    sum=0
    for elem in list:
        square=np.square(elem)
        sum=sum+square
    total=np.sqrt(sum)
    return total

def open_log_file(file_name, cut_of_interest):
    """Takes in a log.txt file, returns the process_names and
    number of events after the final applied cut in list form.
    """
    final_cut_values_list = []
    with open(file_name, "r") as input_file:
        data_raw=input_file.readlines()
        for elem in data_raw:
            if "Process" in elem: #This looks for the line with the process names and saves them.
                temp_var1 = elem.split("\\")
                process_names = temp_var1[0]
                process_names_list = process_names.split(" & ")
            if str(cut_of_interest) in elem:
                temp_var2 = elem.split("\\\\")#The +3 accounts for the first 3 lines in the log file before the "Process" line.
                final_cut_values = temp_var2[0]
                final_cut_values_list = final_cut_values.split(" & ")
    return  process_names_list, final_cut_values_list
   
def check_for_signal(input_file_name, signal_name):
    """Checks if there is a signal. If there is, it prints the names
    and returns a list where the value refers to the index of signal
    in the MC list.
    """
    process_name_list, _ = open_log_file(input_file_name, 0)
#    process_name_list = process_name_list.split(" & ")
    del process_name_list[0] #Removes "process"
    del process_name_list[0] # Removes "data"
    if signal_name == None:
        print("Please rerun the program, specifying the signal using '-s'. Available options are:", process_name_list)
        sys.exit()
    process_present = False
    signal_process_index_values = []
    signal_process_names = []
    for i,elem in enumerate(process_name_list):
        for name in signal_name:
            if str(elem) == str(name):
                signal_process_names.append(elem)
                signal_process_index_values.append(i) #The -1 acconuts for us removing the "process" from the list of process names
                process_present = True
    if process_present == True:
        print("signal(s):", signal_process_names)
    if process_present == False:
        print("signal not found in processes. Available options are:", process_name_list)
        sys.exit()
    return signal_process_index_values

def get_for_cut_names(file_name):
    """Takes in a log.txt file, returns a list of cuts applied.
    """
    list_of_cuts = []
    cuts_reached = False
    with open(file_name, "r") as input_file:
        data_raw=input_file.readlines()
        for elem in data_raw:
            if "Process" in elem:
                cuts_reached = True
            if "tabular" in elem:
                cuts_reached = False
            if cuts_reached == True:
                cut_name = elem.split(" & ")
                list_of_cuts.append(cut_name[0])
    del list_of_cuts[0] #Gets rid of first cut name which is "process"
    return list_of_cuts
                

def convert_opened_file(process_names_list, cut_values_list, signal_process_index_values):
    """Takes the open_log_file output and returns the values of data, data error,
    signal, signal error, background, and background error in float form.
    """
    background_list, background_err_list, signal_list, signal_err_list = [], [], [], []
    del process_names_list[0] #gets rid of first value, which is just "process"
    del cut_values_list[0] #gets rid of first value, which is just the name of the cut applied
    data = float(cut_values_list[0])
    data_err = np.sqrt(data)
    del cut_values_list[0] #remove data since we have saved it to a variable, leaving only MC
    for i, elem in enumerate(cut_values_list):
        if i in signal_process_index_values:
            temp = elem.split(" $\\pm$ ")
            signal_list.append(float(temp[0]))
            signal_err_list.append(float(temp[1]))
        else:
            temp = elem.split(" $\\pm$ ")
            background_list.append(float(temp[0]))
            background_err_list.append(float(temp[1]))
    signal = sum(signal_list)
    signal_err = sum_of_squares(signal_err_list)
    background = sum(background_list)
    background_err = sum_of_squares(background_err_list)
    return data, data_err, signal, signal_err, background, background_err
   
def scale_factor(data, data_err, signal, signal_err, background, background_err):
    """Takes in data, data error, signal, signal error, background, and
    background error and returns the correctly calculated scale factor and error.
    """
    SF = (data-background)/signal
    num_err=sum_of_squares([data_err,background_err])
    SF_err=SF*np.sqrt(np.square(num_err/(data-background))+np.square(signal_err/signal))
    return [SF, SF_err]
    
def efficiency(before_cut_data, before_cut_data_err, before_cut_background, before_cut_background_err, after_cut_data, after_cut_data_err, after_cut_background, after_cut_background_err):
    """Takes in data, data error, signal, signal error, background, and
    background error band returns the correctly calculated efficiency and error.
    """
    efficiency_value =(after_cut_data-after_cut_background)/(before_cut_data-before_cut_background)
    num_err = sum_of_squares([after_cut_data_err,after_cut_background_err])
    denom_err = sum_of_squares([before_cut_data_err,before_cut_background_err])
    efficiency_error_value = efficiency_value*np.sqrt(np.square(num_err/(after_cut_data-after_cut_background))+np.square(denom_err/(before_cut_data-before_cut_background)))
    return ['{:.2E}'.format(efficiency_value), '{:.2E}'.format(efficiency_error_value)]
    
def purity(signal, signal_err, background, background_err):
    """
    """
    purity_value=signal/(signal+background)
    return purity_value
    
def total_event_statistics(data, data_err, signal, signal_err, background, background_err):
    data_to_MC = data/(signal+background)
    data_to_MC_err = data_to_MC*np.sqrt(np.square(data_err/data)+np.square(sum_of_squares([signal_err,background_err])/(signal+background)))
    print("Total data: %.1f +/- %.1f" % (data, data_err))
    print("Total signal: %.1f +/- %.1f" % (signal, signal_err))
    print("Total background: %.1f +/- %.1f" % (background, background_err))
    print("Global data/MC ratio: %f +/- %f" % (data_to_MC, data_to_MC_err))
    
def main():

    #This defines the parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", nargs="+", help="The input file(s)")
    parser.add_argument("-s", "--signal", dest="signal", nargs="+", help=" Specify the signal sample here")
    parser.add_argument("-cs", "--central_selection", dest="central_selection", help="Central selection sample")
    parser.add_argument("-fs", "--final_selection", dest="final_selection", help="Final selection sample")
    parser.add_argument("-su", "--systematic_uncertainty", dest="systematic_uncertainty,", help="Systematic uncertainty to be propagated")
    args = parser.parse_args()

    #Makes sure there is an appropriate input file and signal specified.
    if args.input == None:
        print("Please provide input log.txt files, exiting program.")
        sys.exit()
    else:
        print("input file(s) is:", args.input)
    
    #Make sure there is a signal process
    signal_index_values = check_for_signal(args.input[0], args.signal)
    
    #Decides what the final cut is
    cuts_in_file = get_for_cut_names(args.input[0])
    final_cut_selection = str(cuts_in_file[-1])
    if args.final_selection is not None:
        final_cut_selection = str(args.final_selection)
        if final_cut_selection not in cuts_in_file:
            print("Final selection cut not found, avaiable options are: ", cuts_in_file)
            sys.exit()
    print("Final selection cut is ", final_cut_selection)
    
    #Scale factor is always divided by central scale factor so if we dont have one, set it to 1.
    if args.central_selection is None:
        cs_scale_factor = 1
        print("No central selection specified, efficiency & central selection scale factor therefore not calculated")
    else:
        print("central selection cut is:", args.central_selection)
        if cuts_in_file.index(final_cut_selection) <= cuts_in_file.index(str(args.central_selection)):#Makes sure central selection cut comes before final selection cut
            print("Make sure your central selection cut comes after your final selection cut.")
            print("cut options are: ", cuts_in_file)
            sys.exit()
    
    #Main for loop of the program:
    for i,files in enumerate(args.input):
        if (args.central_selection is not None):
            if str(args.central_selection) not in cuts_in_file:
                print("central selection cut not found, avaiable options are: ", cuts_in_file)
                sys.exit()
            cs_process_names_list, cs_final_cut_values_list = open_log_file(files, args.central_selection)
            cs_data, cs_data_err, cs_signal, cs_signal_err, cs_background, cs_background_err = convert_opened_file(cs_process_names_list, cs_final_cut_values_list, signal_index_values)
            cs_scale_factor, cs_scale_factor_err = scale_factor(cs_data, cs_data_err, cs_signal, cs_signal_err, cs_background, cs_background_err)
            cs_purity = purity(cs_signal, cs_signal_err, cs_background, cs_background_err)
            print("Central selection Scale Factor is: %f +/- %f" % (cs_scale_factor, cs_scale_factor_err))
            print("Central Selection purity is: %f" % cs_purity)
            total_event_statistics(cs_data, cs_data_err, cs_signal, cs_signal_err, cs_background, cs_background_err)
            print("***************************")
        input_process_names_list, input_final_cut_values_list = open_log_file(files, final_cut_selection)
        input_data, input_data_err, input_signal, input_signal_err, input_background, input_background_err= convert_opened_file(input_process_names_list, input_final_cut_values_list, signal_index_values)
        input_scale_factor = scale_factor(input_data, input_data_err, input_signal, input_signal_err, input_background, input_background_err)
        input_purity = purity(input_signal, input_signal_err, input_background, input_background_err)
        print("%s scale factor is: %f +/- %f (central selection SF is accounted for)" % (files, (input_scale_factor[0]/cs_scale_factor), input_scale_factor[1]))
        print("%s purity is: %f" % (files, input_purity))
        if args.central_selection is not None:
            input_efficiency = efficiency(cs_data, cs_data_err, cs_background, cs_background_err, input_data, input_data_err, input_background, input_background_err)
            print("%s efficiency is : %s +/- %s" % (files, input_efficiency[0], input_efficiency[1]))
        total_event_statistics(input_data, input_data_err, input_signal, input_signal_err, input_background, input_background_err)
        print("***************************")
    
if __name__ == '__main__':
    main()

