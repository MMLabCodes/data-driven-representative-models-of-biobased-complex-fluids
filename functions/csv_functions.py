# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:15:14 2022

@author: 983045
"""
import requests 
from bs4 import BeautifulSoup
import csv
import sys
'''
This file contains functions to webscrape information from the pubchem database.

Functions contained in this file are (indented functions are used by the one above):
    
    get_csv_info(csv_file, True)
        get_pubchem_info(names, info, error_file_path="error_log.txt")
        check_error_file(error_file_path)
        
    is_csv_full(csv_file_path)

Further descriptions of each and their usagae can be found after each function.
'''
def get_pubchem_info(names, info, error_file_path="data_clean_error_log.txt"):
    urls = []
    info_list = []
    number_of_errors = 0

    with open(error_file_path, 'w') as error_file:
        for i in range(len(names)):
            molecule = names[i]
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + molecule + "/xml"
            urls.append(url)

        for i in range(len(urls)):
            url = urls[i]
            name = url.split('/')[7]

            data = requests.get(url)
            html = BeautifulSoup(data.content, "xml")
            mw_tag = html.find(name="PC-Urn_label", string=info)

            try:
                mw_parents = mw_tag.find_parent("PC-InfoData")
            except AttributeError:
                print("Logging error for", name)
                error_file.write(f"Error for {name}\n")
                number_of_errors += 1
                info_list.append("N/A")
            else:
                SMILES = mw_parents.find('PC-InfoData_value_sval').string
                print(SMILES)
                info_list.append(SMILES)

        print("Total", info, "searched for", len(urls))
        print("Total errors", number_of_errors)
        print("Success rate is", ((len(urls) - number_of_errors) / len(urls)) * 100, "%")

    return info_list
'''
USAGE: a = get_pubchem_info(names, info)
    names = list of molecules that data is wanted for
    info = this is the info you want to receive
        different info that can be searched for is:
            - "SMILES"
            - "Molecular Weight"
            - "Molecular Formula"
            - Other things can be searched for, but those have not been tested.
        To search for this info, the second argument must be written as above including quotation marks.
        
Returns:
    Returns a list of the desired information for the list of molecules passed to the function.
    
Note:
    This function is used within the the "get_csv_info" to fill in the molecular weights, molecular formula and SMILES section of the csv.
    See the next function to see how the data obtained with this function is written to file.
'''
def get_csv_info(csv_file, fill_csv, error_file_path="data_clean_error_log.txt"):
    if fill_csv == False:
        print("No information required for csv")
        return()
    key = [0, 0, 0, 0, 0] # molecules, mw, area, smi, formula 0 means write it, 1 means data already not
    information = []
    with open(csv_file, 'r') as file:
      reader = csv.reader(file)
      information = []
      for item in reader:
          information.append(item)
    file.close()
    for i in range(len(information)):
        if i == 0:
            pass
        else:
            if information[i][0] != '': # if name col not empty
                key[0] = 1
            if information[i][1] != '': # if mw col not empty
                key[1] = 1
            if information[i][2] != '': # if area col not empty
                key[2] = 1
            if information[i][3] != '': # if smi col not empty
                key[3] = 1
            if information[i][4] != '': # if form col not empty
                key[4] = 1      
    if 0 in key: # this little bit passes through if no info to be added
        pass
    else: # if no info to be added - function exits
        check_error_file(error_file_path)
    names = []
    formulas = []
    mw = []
    smiles = []
    for i in range(len(information)):
        if i == 0:
            pass
        else:
            name = information[i][0]
            names.append(name)
    if key[3] == 0:        
        smiles = get_pubchem_info(names, "SMILES")
    if key[4] == 0:
        formulas = get_pubchem_info(names, "Molecular Formula")
    if key[1] == 0:
        mws = get_pubchem_info(names, "Molecular Weight")
    for i in range(len(information)):
        if i == 0:
            pass
        else:
            if key[4] == 0:
                information[i][4] = formulas[i-1]
            if key[1] == 0:
                information[i][1] = mws[i-1]
            if key[3] == 0:
                information[i][3] = smiles[i-1]
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = information[0]
        # write the header
        writer.writerow(header)
        for i in range(len(information)):
            if i == 0:
                pass
            else:
                data = information[i]
                writer.writerow(data)
    f.close()
    return("CSV file:", csv_file, " is now updated with Molecular weights, Molecular formulas and Smiles")
'''
USAGE: get_csv_info(csv_file_path)
    csv_file = This is the csv file coming from GC-MS data, which has been treated to get peak areas and molecule names
    
    csv_file must have following headers (in this order):
        molecule_name, molecular weight, area, smiles, molecular formula
    
Returns:
    The function returns an updated csv file with mws, formulas and smiles
    
Notes:
    For the bio-oil model generator to work, the peak areas must included in the initial csv.
'''
def check_error_file(error_file_path):
    with open(error_file_path, 'r') as error_file:
        errors = error_file.read().strip()
        if errors:
            print("Error file is not empty. Stopping the code.")
            print("Blankspaces in the csv must be filled")
            print("")
            sys.exit(1)  # Exit with a non-zero status code to indicate an error
        else:
            print("Error file is empty. Continuning.")
'''
USAGE: check_error_file(error_file_path)

Returns:
    A system exit and notification if more data is required
    Nothing is no more dats is required in the csv file
'''
def is_csv_full(csv_file_path):
    with open(csv_file_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        header = next(reader)  # Assuming the first row is the header
        for row in reader:
            for value in row:
                if not value.strip():
                    return(False)  # If any value (after stripping whitespace) is empty, return False
    return(True)  # All values are non-empty
'''
USAGE: is_csv_full(csv_file_path)

Returns:
    True: if csv is full
    False: is data is missing
'''
