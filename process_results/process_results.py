import subprocess
import sys
import os.path
import argparse
from collections import OrderedDict

#constants
EPSILON = 0.0001
NUM_PARAMETERS = 7
NUM_METRICS = 3
MAX_LF,AVG_EPK,MAX_EPK,K,B,T,U = 0,1,2,3,4,5,6
OPTIONS = ["--cuda-max-lf","--cuda-avg-epk","--cuda-max-epk","-K","-B","-t","--ultra-thresh"]
#what percentage difference in performance is considered significant
SIGNIFICANCE_LEVEL = 5
#significance level but for comparing midpoint to neighbours
MIDPOINT_SIGNIFICANCE_LEVEL = 5
#options that use only one dataset
ONE_DATASET = [-3,-1,0,1,3]
#options that take the midpoint of two datasets
TWO_DATASETS = [-4,-2,2]

"""
COMPUTATIONAL HELPER FUNCTIONS
============================================================================
Functions that help to compute certain results needed for the main functions
"""

#compare two floats
def f_equals(a,b):
    return abs(a-b) < EPSILON

#see if two dictionaries have the same parameter values
def same_parameters(dict1,dict2):
    keys = list(dict1.keys())
    all_same = True
    for i in range(NUM_PARAMETERS):
        all_same = all_same and f_equals(dict1[keys[i]],dict2[keys[i]])
    return all_same

#get the average of all values in the given list
def average(ls):
    if(len(ls) == 0):
        print("EMPTY LIST")
        return 0.0
    ave = 0.0
    for i in ls:
        ave += i
    return (ave/len(ls))

#modify the value of p1 based on the instructions given.
#increase: whether or not we want to increase or decrease the parameter
def binary_modify(p1,p2,increase):
    multiplier = 0
    result = p1
    if(increase and p2 is None):
        #double the parameter value
        multiplier = 2.0
    elif(increase and not p2 is None):
        if(p2 - p1 > EPSILON):
            #take the midpoint of current and old
            result = float(p1 + p2)/2.0
        else:
            #not restricted by old value
            multiplier = 2.0
    elif(not increase and p2 is None):
        #halve the parameter value
        multiplier = 0.5
    elif(not increase and not p2 is None):
        if(p2 - p1 < -1 * EPSILON):
            result = float(p1 + p2)/2.0
        else:
            multiplier = 0.5
    #multiply by 2 or 0.5 if no bound found
    if(multiplier != 0):
        result = p1 * multiplier
    return result

#either halve, double or take midpoint of params_1 and params_2 for the parameter specified by this warning
#1st list is the current values, 2nd list is the old values,
#3rd list tells whether a parameter has been modified, 4th list is a list of tokens from a warning
def heed_warning(params_1,params_2,parameter_modified,msg):
    msg_length = len(msg)
    updated = False
    for i in range(msg_length):
        token = msg[i].strip(".")
        #encounter suggestion to increase or decrease a variable
        if(token in ["increasing","decreasing"] and i != msg_length -1):
            parameter = msg[i+1].strip(".")
            #check whether the parameter is a supported option
            supported = False
            for j in range(len(OPTIONS)):
                if(parameter == OPTIONS[j]):
                    supported = True
                    break
            if((not supported) or parameter_modified[j]):
                break
            #proceed to modify parameter value
            updated = True
            multiplier = 0
            increase = False
            param_1,param_2 = params_1[j],None
            if(not params_2 is None):
                param_2 = params_2[j]
            if(token == "increasing"):
                increase = True
            params_1[j] = binary_modify(param_1,param_2,increase)
            parameter_modified[j] = True
    return params_1,parameter_modified,updated

#see whether the provided midpoint is better. has to be better than both
def is_midpoint_better(averages,params,midpoint,test_1,test_2):
    midpoint_test = [float(i) for i in averages[midpoint]]
    midpoint_params = [float(i) for i in params[midpoint]]
    #calculate performance comparison
    test_1_diffs = [(midpoint_test[i] - test_1[i]) for i in range(NUM_METRICS)]
    test_2_diffs = [(midpoint_test[i] - test_2[i]) for i in range(NUM_METRICS)]
    test_1_percent_diffs = [100*(test_1_diffs[i]/test_1[i]) for i in range(NUM_METRICS)]
    test_2_percent_diffs = [100*(test_2_diffs[i]/test_2[i]) for i in range(NUM_METRICS)]
    midpoint_performance = 0
    #only focus on alignment time for now
    for i in range(NUM_METRICS - 2):
        if(test_1_percent_diffs[i] < -1 * MIDPOINT_SIGNIFICANCE_LEVEL and test_2_percent_diffs[i] < -1 * MIDPOINT_SIGNIFICANCE_LEVEL):
            #significant improvement in performance over both
            midpoint_performance = 1
        elif(test_1_percent_diffs[i] > MIDPOINT_SIGNIFICANCE_LEVEL or test_2_percent_diffs[i] > MIDPOINT_SIGNIFICANCE_LEVEL):
            #significant worsening in performance over either test
            midpoint_performance = -1
    return midpoint_performance

"""
handle_direction
==============================================================================
handle the value of the direction variable produced by evaluate_adjustments
direction: the action to take (see evaluate_adjustments)
param_index: the index of the parameter we want to modify
main_params: line number in averages.txt file of the main test we are looking at
second_params: optional second test in case we are taking a midpoint
"""

def handle_direction(direction,param_index,main_params,second_params):
    #sanity check
    if((not second_params is None) and (direction in ONE_DATASET)):
        print("second_params should be None.")
        return
    if((second_params is None) and (direction in TWO_DATASETS)):
        print("second_params shouldn't be None.")
        return

    #always initialise
    main_value = main_params[param_index]
    final_value = 0.0
    stop_test = False

    #do different things for different values of direction
    if(direction == 4):
        print("Error occurred deciding which direction to go in.")
    elif(abs(direction) == 3):
        #halve or double
        multiplier = 2 if direction == 3 else 0.5
        final_value = main_value * multiplier
    elif(abs(direction) == 2 or direction == -4):
        #take midpoint of the two tests
        second_value = second_params[param_index]
        final_value = float(main_value + second_value)/2.0
    elif(abs(direction) == 1):
        #take the value of test_1 or test_2
        final_value = main_value if direction == 1 else second_params[param_index]
        stop_test = True
    elif(direction == 0):
        #take the midpoint parameter
        final_value = main_value
        stop_test = True
    else:
        print(f"Invalid direction: {direction}")

    #edit the value in the list and return it
    final_params = [i for i in main_params]
    final_params[param_index] = final_value

    return final_value,final_params,stop_test

"""
FILE I/O HELPER FUNCTIONS
==============================================================
Functions that perform routine file I/O for the main functions
"""

#try to open file for read
def open_file(filename,mode):
    try:
        FILE = open(filename,mode)
    except FileNotFoundError:
        print(f"File {filename} not found")
        sys.exit()
    return FILE

#extract parameter values from a result file such as {useful/raw}.txt
def read_result_file(filename):
    parameter_values = []
    #use test output
    TEST_FILE = open_file(filename,"r")
    test_data = TEST_FILE.readlines()
    if(len(test_data) < 4):
        print(f"{filename} does not have enough data or profile variable should be true")
        return
    for i in range(len(test_data)-1,len(test_data)-5,-1):
        if( "max-lf:" in test_data[i].split(" ")):
            test_parameters = test_data[i].split(", ")
            break
    for i in range(NUM_PARAMETERS):
        test_data = test_parameters[i].split(": ")
        try:
            parameter_values.append(float(test_data[1]))
        except ValueError:
            print(f"Could not convert parameter data from {filename}.")
            return
    return parameter_values

#extract parameter values from a .profile file
def read_profile_file(filename):
    parameter_values = []
    #only read if it is a .profile file
    if(filename.split(".")[1] != "profile"):
        print("Please enter the name of a file ending in .profile for this option.")
        return
    PROFILE = open_file(filename,"r")
    profile_data = PROFILE.readline().strip()
    PROFILE.close()
    profile_values = profile_data.split(" ")
    try:
        for i in profile_values:
            if( i != "" ):
                parameter_values.append(float(i))
    except ValueError:
        print(f"Could not convert data from {filename} to float.")
        return
    return parameter_values

#write values in the parameter_values list to filename
def write_param_values(filename,parameter_values):
    PARAM_FILE = open_file(filename,"w")
    for i in range(len(parameter_values)):
        if(i > 2):
            parameter_values[i] = int(parameter_values[i])
        PARAM_FILE.write(str(parameter_values[i])+" ")
    PARAM_FILE.close()

#read average performance results from averages.txt file
def read_average_file(filename):
    averages = []
    params = []
    AVE_FILE = open_file(filename,"r")
    line = AVE_FILE.readline()
    if(line == ""):
        print(f"{filename} is empty!")
        return averages
    #split the line by colons
    line_parts = line.strip().split(": ")
    line_averages = line_parts[1]
    line_params = line_parts[2]
    while (line != ""):
        averages.append(line_averages.split(" "))
        params.append(line_params.split(" "))
        line = AVE_FILE.readline()
        if(line != ""):
            line_parts = line.strip().split(": ")
            line_averages = line_parts[1]
            line_params = line_parts[2]
    AVE_FILE.close()
    return averages,params

def write_average_values(filename,average_values,group_params):
    if(len(average_values) != len(group_params)):
        print(f"Length of average_values ({len(average_values)}) and group_params ({len(group_params)}) don't match!")
        return
    AVE_FILE = open_file(filename,"w")
    i = 0
    for j in range(len(average_values)):
        result = average_values[i]
        AVE_FILE.write(f"{i}: ")
        AVE_FILE.write(" ".join(["{0:.2f}".format(i) for i in result]))
        AVE_FILE.write(": ")
        AVE_FILE.write(" ".join(["{0:.2f}".format(i) for i in group_params[i]]))
        AVE_FILE.write("\n")
        i += 1
    AVE_FILE.close()

#write a single line to the file
def write_line(filename,line):
    STOP_FILE = open_file(filename,"w")
    STOP_FILE.write(str(line))
    STOP_FILE.close()

"""
MAIN FUNCTIONS
==========================================
Functions to be called in the main program
"""

"""
calculate_averages
============================================================================================================
calculate average runtime for different performance metrics for each group of tests with the same parameters
"""

def calculate_averages():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',help='Full path to the parameters.txt file that contains parameters and results for different tests')
    args = parser.parse_args(sys.argv[2:])
    filename = args.filename
    #remove number at beginning of each line
    FILE = open_file(filename,"r")
    data = FILE.readlines()
    FILE.close()
    for i in range(len(data)):
        result = data[i]
        for j in range(len(result)):
            if(result[j] == ":"):
                data[i] = result[(j+1):]
                break

    #remove spaces and newline characters
    data_split = [i.replace(" ","").rstrip().split(',') for i in data]

    #transform into dictionary
    data_dict = []
    for result in data_split:
        result_dict = {}
        for i in range(len(result)):
            result_split = result[i].split(":")
            key = result_split[0]
            val = result_split[1]
            result_dict[key] = float(val)
        data_dict.append(result_dict)

    #compute average over datasets with identical parameters. assumes that identical datasets are adjacent
    if(len(data_dict) == 0):
        print("Data Dictionary empty")
        return
    else:
        keys = list(data_dict[0].keys())

    #the parameters used for each group of runs
    group_params = []
    #temporary list containing lists of different test metric values. used for computing averages
    metric_values = [[] for i in range(NUM_METRICS)]
    #list that contains lists of average metric values for each group of tests
    all_averages = []

    #load data from the first entry
    #counter for which metric we are looking at
    k=0
    for j in range(NUM_PARAMETERS,NUM_PARAMETERS+NUM_METRICS):
        metric_values[k].append(data_dict[0][keys[j]])
        k += 1
    #counter for which group of tests we are looking at
    test_num=-1

    #calculate and store the average value for each metric for each test group
    for i in range(len(data_dict)):
        #if the next test is a new test group, calculate average for the current group and save it
        if(i == len(data_dict) - 1 or not same_parameters(data_dict[i],data_dict[i+1])):
            #temporary variable that stores averages for each metric for one group of tests.
            group_averages = []
            #calculate averages for each metric and store them in group_averages
            for j in range(NUM_METRICS):
                group_averages.append(average(metric_values[j]))
            all_averages.append(group_averages)
            #add the test's parameters to the group_params list
            test_num += 1
            group_params.append([])
            for k in range(NUM_PARAMETERS):
                group_params[test_num].append(data_dict[i][keys[k]])
            #reset temporary variable
            metric_values = [[] for l in range(NUM_PARAMETERS)]
        #index for positions in the metric_values list
        k=0
        for j in range(NUM_PARAMETERS,NUM_PARAMETERS+NUM_METRICS):
            metric_values[k].append(data_dict[i][keys[j]])
            k += 1
    #take the path to the parameters file and replace parameters with averages
    output_file = filename.split("/")[:-1]
    output_file.append("averages.txt")
    output_filename = "/".join(output_file)
    write_average_values(output_filename,all_averages,group_params)

"""
evaluate_adjustments
=================================================================================================================================
Compare performance metrics before and after adjusting parameters. Tells whether to stop, keep doubling/halving, or take midpoint
Should be used when only 1 parameter is adjusted at a time (running average without -o flag)
Test 1 is the newer test in which a parameter from test 2 was adjusted
"""

def evaluate_adjustments():
    #parse arguments for the function
    parser = argparse.ArgumentParser()
    parser.add_argument('average',help='Path to the averages.txt file that contains average test metrics')
    parser.add_argument('test_1',type=int,help='Line number of the more recent test in the averages.txt file to compare with a previous test')
    parser.add_argument('test_2',type=int,help='Line number of a previous test in the averages.txt to be compared against')
    parser.add_argument('streak',type=int,help='How many doubles/halves we have done in a row')
    parser.add_argument('num_worse',type=int,help='How many times have our changes led to worse performance')
    parser.add_argument('--w_limit','-w',type=int,default=2,help='Maximum value of num_worse before we stop searching')
    parser.add_argument('--s_limit','-s',type=int,default=4,help='Maximum value of streak before we stop searching')
    parser.add_argument('--midpoint','-m',type=int,default=-1,help='Line number of the test in the averages.txt file which a midpoint was taken between test_1 and test_2')
    args = parser.parse_args(sys.argv[2:])

    #save the args in local variables
    streak = args.streak
    num_worse = args.num_worse
    w_limit = args.w_limit
    s_limit = args.s_limit
    midpoint = args.midpoint

    #list that contains all the argument values to this function for the next iteration. contents are written to a evaluate.params file
    evaluate_params = []

    #read data from files
    averages,params = read_average_file(args.average)

    #extract the relevant lines to work with
    index_1 = args.test_1
    index_2 = args.test_2
    try:
        #convert string to float
        test_1 = [float(i) for i in averages[index_1]]
        test_2 = [float(i) for i in averages[index_2]]
        params_1 = [float(i) for i in params[index_1]]
        params_2 = [float(i) for i in params[index_2]]
    except IndexError:
        print("Indices out of range.")
        return

    #calculate performance comparisons
    diffs = [(test_1[i] - test_2[i]) for i in range(NUM_METRICS)]
    percent_diffs = [100*(diffs[i]/test_2[i]) for i in range(NUM_METRICS)]

    #1 means test 1 has better performance, 0 means no significant difference, -1 means worse performance
    performance = 0
    #only focus on alignment time for now
    for i in range(NUM_METRICS - 2):
        if(percent_diffs[i] > SIGNIFICANCE_LEVEL):
            #significant worsening in performance
            performance = -1
        elif(percent_diffs[i] < -1 * SIGNIFICANCE_LEVEL):
            #significant improvement in performance
            performance = 1

    #find the index of the parameter that was updated
    i = 0
    while i < len(params_1):
        if(not f_equals(params_1[i],params_2[i])):
            break
        i += 1

    #arguments that we will pass in to the handle_direction function
    direction = 4
    midpoint_params = [float(i) for i in params[midpoint]]
    main_params = None
    second_params = None

    """
    -4: midpoint of two tests, -3: halve newer value, -2: midpoint towards test_2 (older one), -1: take test_2's value,
    0: take the midpoint value, 1: take test_1's value, 2: midpoint towards test_1, 3: double, 4: error
    """
    #based on arguments of the function and performance decide the course of action to take
    if(num_worse == w_limit or streak == s_limit):
        #already reached limit so see which one to take
        if(midpoint == -1):
            #no midpoint provided so take the best of the new or old
            direction = 1 if percent_diffs[0] < 0 else -1
            main_params = params_1 if direction == 1 else params_2
        else:
            #midpoint provided so take the best of the three
            midpoint_performance = is_midpoint_better(averages,params,midpoint,test_1,test_2)
            if(midpoint_performance == 1):
                #midpoint is the best so choose it
                direction = 0
                main_params = midpoint_params
            else:
                #midpoint isn't the best so choose between test_1 and test_2
                direction = 1 if percent_diffs[0] < 0 else -1
                main_params = params_1 if direction == 1 else params_2
    elif(midpoint == -1):
        #No midpoint provided, so simply compare the two tests
        if(performance == 1):
            #performance has improved (test 1 is better than test 2)
            if(streak > 0):
                #on a streak so keep doubling/halving
                direction = 3 if params_1[i] > params_2[i] else -3
                #always double/halve the newer test parameters
                main_params = params_1
                streak += 1
            else:
                print("Why are you here if you're not on a streak?")
                return
        elif(performance == -1 or performance == 0):
            #check the midpoint of new and old
            direction = -4
            main_params = params_1
            second_params = params_2
            num_worse += 1
            streak = 0
    else:
        #Midpoint provided so find out whether it is better
        midpoint_performance = is_midpoint_better(averages,params,midpoint,test_1,test_2)
        if(midpoint_performance != 1):
            #check the midpoint of new and old
            num_worse += 1
        #take the midpoint no matter what
        direction = 2 if percent_diffs[0] < 0 else  -2
        main_params = midpoint_params
        second_params = params_1 if direction == 2 else params_2

    #sanity check
    if(main_params is None):
        print("main_params is None!")
        return

    #make adjustments to parameter values and return them
    final_value,final_params,stop_test = handle_direction(direction,i,main_params,second_params)

    #take the path to the averages file and replace it with a different filename
    output_profile = output_stop_test = output_evaluate_params = args.average.split("/")[:-1]
    output_profile.append("new.profile")
    output_stop_test.append("stop.test")
    output_evaluate_params.append("evaluate.parameters")
    #write values to the evaluate_params list and
    #data stored in format: test_1, test_2, streak, num_worse, w_limit, s_limit, midpoint
    if direction in ONE_DATASET:
        #no midpoint but next test is the new test_1
        evaluate_params = [index_1+1,index_1,streak,num_worse,w_limit,s_limit]
    elif direction in TWO_DATASETS:
        #midpoint is the next test
        evaluate_params = [index_1,index_2,streak,num_worse,w_limit,s_limit,index_1+1]
    else:
        print("Hmmm something wrong with direction")
    output_filename = "/".join(output_profile)
    output_stop_filename = "/".join(output_stop_test)
    output_evaluate_filename = "/".join(output_evaluate_params)

    #write new parameters to the file
    write_param_values(output_filename,final_params)
    write_line(output_stop_filename,stop_test)
    write_line(output_evaluate_filename,", ".join(evaluate_params))


    print("\n====================================================================================================")
    print(f"Test 1 averages: {test_1} Parameters used: {params_1}")
    print(f"Test 2 averages: {test_2} Parameters used: {params_2}")
    print(f"Absolute differences: {diffs}")
    print(f"Percentage differences: {percent_diffs}")
    print(f"Changed parameter '{OPTIONS[i]}' to: {final_value}. test 1: {params_1[i]}, test 2: {params_2[i]}")
    print(f"New parameters: {final_params}")
    print("====================================================================================================\n")

"""
adjust_warnings
===============================================================================
adjusts parameters by either halving or doubling based on error messages found
"""

def adjust_warnings():
    #parse arguments for the function
    parser = argparse.ArgumentParser()
    parser.add_argument('warnings',help='Path to the warning_{test_number}.txt file that contains warnings for a test')
    parser.add_argument('parameters',help='Path to the .profile/stderr output file that contains the parameter values for a test')
    parser.add_argument('-x',action='store_false',help='Whether or not the parameter file is a .profile file or not. Default is True')
    parser.add_argument('-o',action='store_false',help='Whether or not to update only one or several parameters at once. Default is True')
    args = parser.parse_args(sys.argv[2:])
    warning_file = args.warnings
    param_file = args.parameters
    profile = args.x
    one_update = args.o
    #read and parse data
    FILE = open_file(warning_file,"r")
    data = FILE.readlines()
    #if there were no warnings:
    if(len(data) > 1 and data[0] == "No warnings."):
        print("No warnings found.")
        return None
    FILE.close()
    tokenized_data = [line.split(" ")[1:] for line in data]
    #whether or not the parameter in this position has been modified already
    parameter_modified = [False for i in range(NUM_PARAMETERS)]

    #parameter values to be used for next test. initialize based on either given profile or output of the test.
    parameter_values = []
    old_parameter_values = None
    if(profile == False):
        parameter_values = read_result_file(param_file)
    else:
        #use profile values
        parameter_values = read_profile_file(param_file)
        if(os.path.isfile(param_file + ".old")):
            #read values from the old parameter file as well so that we have upper/lower bounds for binary search
            old_parameter_values = read_profile_file(param_file + ".old")
    if(parameter_values is None):
        print("Failed to parse parameter file. Have a nice day.")
        return

    if(len(parameter_values) < NUM_PARAMETERS):
        print(f"Something went wrong with reading {param_file}. Not enough parameters read.")
        return

    #write old parameters to new file. allows upper and lower bounds for binary search
    if(param_file.split(".")[1] != "profile"):
        #we have read data from the result file so create a profile and store it there
        param_file = "/".join(param_file.split("/")[:-2]) +  "/test_" + param_file.split("/")[0] + ".profile"

    #save old values before continuing
    old_param_file = param_file + ".old"
    write_param_values(old_param_file,parameter_values)

    #automatically adjust -B option to specified limit
    if(len(tokenized_data) > 0 and ("accommodate" in tokenized_data[0])):
        for i in range(len(tokenized_data[0])):
            if(tokenized_data[0][i] == "bases."):
                #extract the value and remove the trailing 'M'
                try:
                    value_mil = float(tokenized_data[0][i-1][:-1])
                except ValueError:
                    print("Could not convert to float.")
                    break;
                value = int(value_mil * 1000000)
                parameter_values[B] = value
                parameter_modified[B] = True

    #loop through warnings and adjust parameters if there is a suggestion to
    for msg in tokenized_data:
        parameter_values,parameter_modified,updated = heed_warning(parameter_values,old_parameter_values,parameter_modified,msg)
        #if an update has been made and we only want one update
        if(one_update and updated):
            break

    #print helpful message notifying if no parameters were modified
    modified = False
    for i in parameter_modified:
        modified = modified or i
    if(not modified):
        print("No parameters modified.")
        return(parameter_values)
    #write new parameter values to the same parameter file
    write_param_values(param_file,parameter_values)
    print("New values: ",parameter_values,"stored in ",param_file)
    return(parameter_values)

"""
modify_parameter
=============================================================================================================
Directly modify the value of a parameter from a given .profile file by specifying the multiplier or new value
Used for halving or doubling the value of the given parameter
"""

def modify_parameter():
    #parse arguments for the function
    parser = argparse.ArgumentParser()
    parser.add_argument('profile',help='Path to the .profile file that contains the parameter values to be changed')
    parser.add_argument('--old','-o',default=None,help='Optional path to the name for the file that will store the current parameter values specified by [profile file]')
    parser.add_argument('parameter_number',type=int,help='The number of the parameter to be modified. MAX_LF,AVG_EPK,MAX_EPK,K,B,T,ULTRA_THRESH are represented by 0-6 respectively')
    parser.add_argument('--scale_factor','-s',type=float,default=0.0,help='The factor to multiply with the current value of the specified parameter')
    parser.add_argument('--value','-v',type=float,default=0.0,help='The value to assign to the specified parameter. Overrides effect of --scale-factor')
    args = parser.parse_args(sys.argv[2:])
    #read data from files
    current_params = read_profile_file(args.profile)
    #sanity checks
    if (args.parameter_number < 0 or args.parameter_number > NUM_PARAMETERS):
        print(f"Parameter number {args.parameter_number} out of range.")
        return None
    if args.scale_factor < 0:
        print(f"Parameter number {args.scale_factor} out of range.")
        return None
    if args.value < 0:
        print(f"Parameter number {args.value} out of range.")
        return None
    #default behaviour if no output file for the old values is specified
    if args.old is None:
        write_param_values(args.profile + ".old",current_params)
    else:
        write_param_values(args.old,current_params)
    #retrieve the desired value and modify it
    old_value = current_params[args.parameter_number]
    print("Old: ",old_value)
    #modify parameter value based on args provided
    if(args.value > 0):
        new_value = args.value
    elif(args.scale_factor > 0):
        new_value = old_value * args.scale_factor
    else:
        print("Please set one of the -s or -v flags.")
        return None
    current_params[args.parameter_number] = new_value
    write_param_values(args.profile,current_params)
    print("New: ",new_value)
    return new_value

"""
midpoint()
============================================================================================================================================
Take the midpoint of two parameter values or go in the best direction and create a new .profile file that uses the midpoint as the new value
"""

def midpoint():
    #parse arguments for the function
    parser = argparse.ArgumentParser()
    parser.add_argument('parameters',help='Path to the parameters.txt file that contains the parameter values for each test')
    parser.add_argument('test_1',type=int,help='Line number of the more recent test in the averages.txt file to compare with a previous test')
    parser.add_argument('test_2',type=int,help='Line number of a previous test in the averages.txt to be compared against')
    parser.add_argument('--midpoint','-m',type=int,help='Line number of the test in which a midpoint was taken between test_1 and test_2')
    args = parser.parse_args(sys.argv[2:])
    print("I don't do anything yet :(")
    return

"""
main
===================================================================================
main program which parses command line arguments and executes the correct functions
"""

parser = argparse.ArgumentParser(
    description='Processes results from running f5c',
    usage=
    '''python3 /path/to/process_results.py [mode] ([other arguments])

    modes:
            'average' :            (for calculating average performance metrics for tests)
            'warning' :            (for adjusting algorithm parameters based on warning messages from f5c)
            'evaluate':            (for evaluating whether a change in parameter led to better or worse performance)
            'modify'  :            (for directly modifying the value of a given parameter)
            'midpoint':            (for taking the midpoint between a given set of parameters)

    usage:
            'average [result file]'
                [result file] :        (path to parameters.txt file generated by param_test.sh)

            'evaluate [average file] [test 1] [test 2] [streak] [num_worse] (-w [w_limit] -s [s_limit] -m [midpoint])'
                [average file] :       (evaluate mode only. path to averages.txt file to compare output from)
                [test 1] :             (line number of more recent test in the average file to compare with a previous test)
                [test 2] :             (line number of previous test in the average file to be compared against)
                [streak] :             (how many doubles/halves we have done in a row)
                [num_worse] :          (how many times have our changes led to worse performance)
                '-w [w_limit]' :       (maximum value of num_worse before we stop searching)
                '-s [s_limit]' :       (maximum value of streak before we stop searching)
                '-m [midpoint]' :      (Line number of the test in which a midpoint was taken between test_1 and test_2)

            'warning [warning file] [parameter file] ([other flags])'
                [warning file] :       (path to warning_{test_number}.txt file generated by param_test.sh)
                [parameter file] :     (path to .profile (or f5c stderr output file in warning mode) that lists parameters used in the test)
                '-x' :                 (optional boolean flag to specify whether using a .profile file or not. default is True)
                '-o' :                 (optional boolean flag to specify whether to only update one parameter or not. default is True)

            'modify [profile file] (-o [old profile]) [parameter number] [direction] (-s [scale factor]/-v [value]):
                [profile file] :       (path to the .profile file that contains the parameter values to be changed)
                '-o [old profile]' :   (optional path to the name for the file that will store the current parameter values specified by [profile file])
                [parameter number] :   (the number of the parameter to be modified)
                '-s [scale factor]' :  (optional factor to multiply with the current value of the specified parameter)
                '-v [value]' :         (optional value to assign to the specified parameter. overrides effect of --scale-factor)

            'midpoint [parameters] [test_1] [test_2] (-m [midpoint])'
                [parameters]
                [test_1]
                [test_2]
                '-m [midpoint]'
            '''
)

parser.add_argument('mode', help='Functionality you would like to use. Options are: averages,warnings')
if(len(sys.argv) < 2):
    parser.print_usage()
    sys.exit()
args = parser.parse_args(sys.argv[1:2])
if(args.mode == 'average'):
    calculate_averages()
elif(args.mode == 'evaluate'):
    evaluate_adjustments()
elif(args.mode == 'warning'):
    adjust_warnings()
elif(args.mode == 'modify'):
    modify_parameter()
elif(args.mode == 'midpoint'):
    midpoint()
else:
    parser.print_usage()
