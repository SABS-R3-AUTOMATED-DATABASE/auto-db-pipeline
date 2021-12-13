"""
This file contains the code to retrieve potential/possible PDB IDs from a url (e.g. the url of a paper) 
and return those potential PDB IDs as a list. The final (fifth) function get_pdbs utilizes
the first four functions. Please read the docstrings in the functions for more information. 

Note that there will likely be false positives (potential PDB IDs that are not actual PDBs).
"""


from bs4 import BeautifulSoup
import requests
import re


def get_txt(url):
    """
    Gets all the text from the html of a url of a paper. 
    Returns the prettified text as a string.
    
    Requires the following libraries:
    * requests
    * BeautifulSoup from bs4
    """
    
    # Get the response to the request to the url
    r = requests.get(url, allow_redirects=True)
    
    # Get the beautiful soup 
    soup = BeautifulSoup(r.content, 'html.parser')
    
    # Convert the html in the beautiful soup into text (type string)
    soup_txt = soup.prettify()
    
    return soup_txt



def get_alpha_numeric(string):
    """
    Takes a string and returns as a list only those that are alpha-numeric 
    using regular expressions.
    
    Requires the following libraries:
    * re
    
    Misc note:
    There are two equivalent ways of doing this:
    finding all alpha-numeric strings, and,
    splitting on all non-alpha-numeric strings.
    We arbitrarily use the first version. 
    Hint: we can test that these give equivalent results.
    """
    
    # Find all alphabetical and numeric strings
    alpha_numeric = re.findall('[a-zA-Z0-9]+', string)

    # Split on all non alphabetical and numeric strings
    # alpha_numeric = re.split('[^a-zA-Z0-9]+', string)
    
    return alpha_numeric
    
    

def filter_pdb(string_list):
    """
    Takes a list of strings and returns a sub list of strings 
    of only those strings that follow the rules of PDB IDs, given below:
    
    1. All characters are alphabetical or numeric
    1. Length of four.
    2. First charater is numeric. 
    3. First character is in the range of 1-9 inclusive (greater than 0). 

    These rules are taken from [protopedia](https://proteopedia.org/wiki/index.php/PDB_code).
    
    Returns a list.
    
    Misc note:
    Checking the first rule is actually redundant in our pipeline because we filter
    only the alpha-numeric strings in the "get_alpha_numeric" function.

    """
    output_list = list()
    
    for string in string_list:
        
        # check if string is alpha-numeric
        if string.isalnum():

            # check length is 4
            if len(string) == 4:  

                # check first character is numeric
                if string[0].isnumeric():  

                    # check the numeric character is > 0.
                    if int(string[0]) > 0:  

                        # We have met all the PDB conditions
                        output_list.append(string)

    return output_list



def standardize(pdb_list):
    """
    Takes a list of possible PDB IDs and does the following:
    1. Converts them to upper case.
    2. Removes duplicates.
    3. Sorts in alphabetical-numerical order.
    
    Returns a list.
    
    Misc note:
    We can remove any of these steps and this function will not throw an error.
    """
    
    # 1. Convert to upper case
    pdb_list = [pdb.upper() for pdb in pdb_list]
    
    # 2. Remove duplicates
    pdb_list = list(set(pdb_list))
    
    # 3. Sort alphabetically / numerically
    pdb_list.sort()
    
    return pdb_list



def get_poss_pdbs(url):
    """
    Gets the possible PDBs from the url of a paper.
    
    Returns a list.
    """
    
    # Get html of the url as a string
    string = get_txt(url)
    
    # Filter out the alpha-numeric strings
    string_list = get_alpha_numeric(string)
    
    # Filter out those that follow the PDB ID rules
    pdb_list = filter_pdb(string_list)
    
    # Standardize the list of PDB IDs
    pdb_list = standardize(pdb_list)
    
    return pdb_list

