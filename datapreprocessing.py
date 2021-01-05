import pandas as pd 
import numpy as np
from rdkit import Chem  
from rdkit.Chem import rdMolDescriptors,Draw,AllChem 
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display
import urllib.request
import requests
import os
import csv
from pprint import pprint
import re
import itertools
#import util 
#sys.setrecursionlimit(10000)
from type_conversion import smarts_to_smiles, circular_fingerprints, similarity_check, running_through_all_reaction_templates
from rdkit import DataStructs
from rdkit import RDLogger
from type_conversion import moltosmiles
import numpy as np
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.Chem.EState as EState
import rdkit.Chem.rdPartialCharges as rdPartialCharges
import rdkit.Chem.rdChemReactions as rdRxns
import rdkit.Chem.rdChemReactions as ReactionFromSmarts
from rdkit.Chem import rdChemReactions
from chemical_module import Chemical
from reaction import Reaction
from rdkit.Chem import rdMolDescriptors,Draw,AllChem 
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display




att_dtype = np.float32
RDLogger.DisableLog('rdApp.*')

def kegg_data():
    kegg_reactions = pd.read_csv("ATLAS_Curated_KEGG_Reactions.csv")

def reaction_rules():
    return pd.read_csv("retrorules_rr02_rp2_flat_retro.csv")

def SMILE_rxn_rules():
    return pd.read_csv("SMILE_reaction_rules.csv")

def breaking(string):
    """We have long kegg ID strings but since we know that they can be only of length 6 we 
    are splitting them into parts of 6."""
    return [string[start:start+6] for start in range(0, len(string), 6)]

def conversion(data):
    """
    Kegg reactions are not directly writing molecular compounds but they reference each compound by 
    using something called kegg ID, which is unique for each compound. 
    """
    temp= []
    for i in range(0,data.shape[0]):
        temp.append((data.REACTION[i]).split(" "))
    su = ""
    for i in range(0,len(data.REACTION)):
        su +=str(data.REACTION[i])
    import re 
    res = re.sub('\(.*?\)','', su).split(" ")
    temp = []
    for i in range(0,len(res)):
        if "C" in res[i]:
            temp.append(res[i])
    stage2 = ""
    for i in range(0,len(temp)):
        stage2 +=str(temp[i])
    stage3 = []
    stage3 = stage2.split("+")
    stage4 = []
    for i in range(0,len(stage3)):
        if len(stage3[i])==6:
            stage4.append(stage3[i])
        if len(stage3[i])> 6:
            stage4.extend(breaking(stage3[i]))
    return list(set(stage4))

def collecting_reactions_kegg(data):
    """
    Objective:- This function makes a call to the kegg api makes a list of lists with index 1 returning the Kegg ID and second index 
    being the kegg compund formula. 
    Input:- I have made a function call within this function which already does that, 
    NOT NEEDED(please input a list containing the kegg id which you wish to query, these need to input one at a time 
    without any seperator or plus, stoichometric coefficeints).   
    Output:- [kegg ID, corresponding SMILES string]
    Return Type : Lists of lists 
    """
    def ses(url_t):
        with urllib.request.urlopen(url_t) as kegg:
            data = kegg.read()

        try:
            textfile = open('demo20.txt', 'wb')
            textfile.write(data)
            textfile.close()
            demo = Chem.MolFromMolFile('demo20.txt')
            return Chem.MolToSmiles(demo)

        except:
            os.remove('demo20.txt')

    list_url = []
    stage4 = conversion(data)
    unprocessed = []
    print(stage4[0:10])
    for i in range(0,len(stage4)):
        url_t = 'http://rest.kegg.jp/get/%s/mol'%(stage4[i])
        if i != 17:
            try:
                list_url.append([stage4[i],ses(url_t)])
                print(len(list_url))
            except :  # This is the correct syntax
                unprocessed.append(stage4[i])
                pass 
    return list_url, unprocessed

print("this is working")


def split_into_reactions_products(string):
    """
    Input:-  would be a string which will check the position of <=> special character and based on that split the 
    data into reactents and products. Another issue which we did not want to run into was dealing with spacing 
    stoichometric coefficients, so this function would also be removing them as well. 
    
    Output:- 2 strings, first one would be the reactent and other one would be product. 
    """
    reactent,product = [],[]
    #string = string.split(" ")
    index = string.index("<=>")
    reactent.append(string[:index -1])
    product.append(string[index+3:])
    reactent, product =  re.sub(r'\([^)]*\)', '',reactent[0]), re.sub(r'\([^)]*\)', '',product[0])
    return " ".join(reactent.split()), " ".join(product.split())

def smiles(kegg_id,dictionary):
    """
    Input : 
        1. Kegg ID 
        2. Dictionary with the kegg-smiles conversion already defined
    Output:- 
        Return type : SMILE string 
    """
    return dictionary[kegg_id]

def kegg_dictionary():
    """
    Input : The location of this smiles file with only 2 columns, 
    Column 1: Name of the kegg Compound 
    Column 2: Corresponding SMILE String 
    
    Output : Format Dictionary 
        1. keys are the names of compounds, 
        2. values are the SMILE strings
    
    """
    
    class ListDict(dict):
        """ Dictionary who's values are lists. """
        def __missing__(self, key):
            value = self[key] = []
            return value

    filename = 'keggsmiles.csv'

    lstdct = ListDict()
    with open(filename, "rt") as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            key, value = row[:2]
            lstdct[key].append(value)
    return lstdct 

def dictionary_mapping(lis,dictionary):
    """
    Input: String Format
    
    The input must either be the the whole of reactents and whole of products, which must be in string format 
    seperated by + sign without having stoichometric coefficients"""
    #my_string = ""
    my_string = []
    plus = lis.split("+")
    for i in range(0,len(plus)):
        #my_string += str(smiles(plus[i],dictionary))
        my_string.append(smiles(plus[i],dictionary))
    return my_string

def original_kegg_reactions(kegg_reactions,converting_dictionary): 
    """
    MAIN function:- Kegg Datset
    This function should be able to generate the smiles of the origial kegg 
    Return type: List 1: Reactents, List 2: Products 
    """
    #oro = unconverted_kegg_reactions()
    lstdct = converting_dictionary()
    REA, PRO = [],[]
    for i in range(0,kegg_reactions.REACTION.shape[0]):
        #print(i)
        r, p = split_into_reactions_products(kegg_reactions.REACTION[i])
        #print(p)
        REA.append(dictionary_mapping(r,lstdct))
        PRO.append(dictionary_mapping(p,lstdct))
    return REA, PRO 

def reducing_reaction_rules(data,cutoff = 50):
    """
    Input: Extracted reaction rules which should have 2 columns one with smarts format and more 
    importantly a smile format with column name "SMILES".,
    optional argument: Cutoff value for how many times should that rule at least occur in the reactions 
    database.
    Output: reduced set of reaction rules
    Example :-
    SMILE_rxn_rules = pd.read_csv("/Users/aasimwani/Downloads/MS research/extractedData/SMILE_reaction_rules.csv")
    reducing_reaction_rules(SMILE_rxn_rules,cutoff = 100)
    """
    #SMILE_rxn_rules = pd.read_csv("/Users/aasimwani/Downloads/MS research/extractedData/SMILE_reaction_rules.csv")
    number_list = np.array(data["REACTION_RULES"])
    (unique, counts) = np.unique(number_list, return_counts=True)
    frequencies = np.asarray((unique, counts)).T
    bot = pd.DataFrame(frequencies,columns = ["strings","counts"])
    counted = bot.sort_values(by ='counts',ascending = False )
    total = np.sum(counted["counts"])
    return counted.loc[counted["counts"]> cutoff]

def classififes_negative_reaction():
    """
    This self contained function should be converting the kegg reactions smiles and using more then one product to create;
    this function would sort the list of products based on the length of the smile string. The first one would be classified as the 
    more important string classified as positive strings and others which are smaller in length strings would be classifies as negative 
    strings. There exist multiple null reactions in the datacase wich will also need to be removed.
    """
    REA, PRO = original_kegg_reactions(kegg_reactions(),converting_dictionary=kegg_dictionary())
    temp = []
    for i in range(0,len(PRO)):
        current_length = len(PRO[i])
        if current_length >1:
            sorted_product = list(reversed(sorted(PRO[i], key=len)))
            for j in range(0,current_length):
                if sorted_product[j]:
                    if j == 0 :
                        temp.append([REA[i],sorted_product[j],"Positive"])
                    else:
                        temp.append([REA[i],sorted_product[j],"Negative"])
                else:
                    pass 
        else:
            temp.append([REA[i],PRO[i],"Positive"])
    df = pd.DataFrame(temp)
    df.columns = ["Reactent","Product","Type"]
    return df

def generating_reaction_fingerprints():
    def circular_fingerprint(smile_string):
        mol = Chem.MolFromSmiles(smile_string)
        bi = {}
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, bitInfo=bi)
        return fp.ToBitString()

    def encoding(current_element):
        return [int(item) for item in circular_fingerprint(current_element)]

    def making_product_reactents_morgan_fingerprints(a_a):
        merging_element_df = []
        temporary_merging_element = np.zeros((1,2048))
        for i in range(0,len(a_a)):
            for j in range(0,len(a_a[i])):
                item = a_a[i][j]
                if item != []:
                    temporary_merging_element+=encoding(item[0])
                #print(len(merging_element_df))
                temporary_merging_element = np.zeros((1,2048))
            merging_element_df.append(temporary_merging_element)
        return merging_element_df

    REA, PRO = original_kegg_reactions(kegg_reactions,converting_dictionary=kegg_dictionary)
    a  = making_product_reactents_morgan_fingerprints(REA)
    reactent_fingerprints = pd.DataFrame([item[0] for item in a])
    b  = making_product_reactents_morgan_fingerprints(PRO)
    product_fingerprints = pd.DataFrame([item[0] for item in b])
    reaction_fingerprints = np.add(reactent_fingerprints,product_fingerprints)
    return reaction_fingerprints


def generating_products_from_reaction_rules(molefile,cutoff):
    """
    Input - Molefile of the 
    Number_of_rules - number of reaction rules which you want to generate reactents 
    """
    #morphine = Chem.MolFromMolFile("/Users/aasimwani/Downloads/MolData/morphine.mol")
    ruleDF = pd.read_csv("retrorules_rr02_rp2_flat_retro.csv")
    SMILE_rxn_rules = pd.read_csv("SMILE_reaction_rules.csv")
    extracted_reaction_rules = reducing_reaction_rules(SMILE_rxn_rules,cutoff)
    rules = pd.DataFrame(pd.merge(extracted_reaction_rules.reset_index(),ruleDF.reset_index(),how = "inner",on = "index")["Rule"])
    item = []
    negative_reactions = []
    for i in range(0,rules.shape[0]):
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts(rules["Rule"][i])
        reacts = Chem.AddHs(Chem.MolFromSmiles(molefile))
        products = rxn.RunReactants((reacts,))
        item.append(products)
        if products == ():
            negative_reactions.append([molefile,rdChemReactions.ReactionFromSmarts(rules['Rule'][i])])
    formatted = list(itertools.chain(*item))
    let = []
    for i in range(0,len(formatted)):
        let.append(Chem.MolToSmiles(formatted[i][0]))
    a = pd.DataFrame(let)
    a.columns = ["ar1"]
    only_uniques = a["ar1"].unique()
    return only_uniques, negative_reactions

def please_work(initial,target,cutoff):
    def reaction_rules_final(cutoff):
        data = pd.read_csv("retrorules_rr02_rp3_nohs/retrorules_rr02_flat_all.tsv",sep='\t')
        number_list = np.array(data["Rule_SMARTS"])
        (unique, counts) = np.unique(number_list, return_counts=True)
        frequencies = np.asarray((unique, counts)).T
        bot = pd.DataFrame(frequencies,columns = ["strings","counts"])
        counted = bot.sort_values(by ='counts',ascending = False )
        unique_strings = np.array(counted.loc[counted["counts"]> cutoff]["strings"])
        return unique_strings

    def state_space_generation(molefile,unique_strings):
        item = []
        from rdkit.Chem import rdChemReactions
        molefile = morphine()
        for i in range(0,len(unique_strings)):
            rxn = rdChemReactions.ReactionFromSmarts(unique_strings[i])
            reacts = Chem.MolFromSmiles(molefile)
            products = rxn.RunReactants((reacts,))
            if products != (): 
                item.append(products)
        formatted = list(itertools.chain(*item))
        print("Number of Possible poducts generated are", len(formatted))
        string_for_morphine = [Chem.MolToSmiles(item[0]) for item in formatted]
        return string_for_morphine

    def searching_state_space_for_desired_product(possible_reactants,pattern):
        totalParts = []
        for item in possible_reactants:
            subparts = (item.split("."))
            totalParts.append(subparts)
        found_it = []
        positions = []
        for i in range(0,len(totalParts)):
            m = Chem.MolFromSmiles(pattern)
            patt = Chem.MolFromSmiles(totalParts[i][0])
            if patt.HasSubstructMatch(m) and m.HasSubstructMatch(patt):
            #if patt.HasSubstructMatch(m):
                found_it.append(Chem.MolFromSmiles(totalParts[i][0]))
                positions.append(i)
        return found_it, positions
    
    r1 = reaction_rules_final(49)
    possible_products = state_space_generation(initial,r1)
    r3,r4 = searching_state_space_for_desired_product(possible_products,target)
    return ("Index position where target molecule is found in state space:",r4 )

def get_hydrogen_atoms(my_mol):
    my_mol = Chem.MolFromSmiles(moltosmiles(my_mol))
    my_mol_with_explicit_h = Chem.AddHs(my_mol)
    return my_mol_with_explicit_h.GetNumAtoms() - my_mol_with_explicit_h.GetNumHeavyAtoms()


def reversing_smarts():
    data = pd.read_csv("retrorules_rr02_rp3_nohs/retrorules_rr02_flat_all.tsv",sep='\t')
    number_list = np.array(data["Rule_SMARTS"])
    reversed_strings = []
    for i in range(0,len(number_list)):
        reaction_smarts_retro = number_list[i].split(">>")[1] + ">>" + number_list[i].split(">>")[0]
        reversed_strings.append(str(reaction_smarts_retro))    
    return reversed_strings
"""
    all_values = []
    current_max = []
    for i in range(0,len(REA)):
        for item in REA[i]:
            if item != []:
                temp = get_hydrogen_atoms(item[0])
            current_max.append(temp)
        all_values.append(max(current_max))

    Y = pd.DataFrame(all_values)
    salutaridine =  "COC1=C[C@]23CCN(C)[C@H](Cc4ccc(OC)c(O)c24)C3=CC1=O"
r_reticulene = "CN1CCC2=CC(=C(C=C2[C@@H]1CC3=CC(=C(C=C3)OC)O)O)OC"
s_reticulene = "COc1ccc(C[C@@H]2N(C)CCc3cc(OC)c(O)cc23)cc1O"
s_n_methyl_coculaurine = "COc1cc2CCN(C)[C@@H](Cc3ccc(O)cc3)c2cc1O"
#s_coculaurine = "COc1cc2CCN[C@@H](Cc3ccc(O)cc3)c2cc1O"
s_coculaurine= "COc1cc2c(cc1O)[C@H](Cc1ccc(O)cc1)NCC2"
#COC1=C(C=C2C(NCCC2=C1)CC3=CC=C(C=C3)O)O
#norococulaurine = "Oc1ccc(C[C@@H]2NCCc3cc(O)c(O)cc23)cc1"
norococulaurine = "Oc1ccc(CC2NCCc3cc(O)c(O)cc32)cc1"
                   
s3_hydroxy_methylcoclaurine = "COc1cc2CCN(C)[C@@H](Cc3ccc(O)c(O)c3)c2cc1O"
salutaridnal = "COC1=C[C@]23CCN(C)[C@H](Cc4ccc(OC)c(O)c42)C3=C[C@@H]1O"
                
thebaine = "COC1=CC=C2[C@H]3Cc4ccc(OC)c5O[C@@H]1[C@]2(CCN3C)c45"

glucose = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
glucose_6p = "O[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O"

    """

def pushing_through_reaction_for_reaction_centers(status,reactent=None,product = None):
    def extracting(reactant,product):
        R2 = "*"
        temp_reactant, temp_product = [],[]
        for reactant_index in range(0,len(reactant)):
            try:
                if reactant[reactant_index] != []:
                    m = Chem.MolFromSmiles(reactant[reactant_index][0][0])
                    if m is not None:
                        reactant_smarts = Chem.MolToSmarts(m,isomericSmiles=True)
                        temp_reactant.append(reactant_smarts.format(R2))
                    else:
                        pass 
            except:
                IndexError
        part1 = ".".join(temp_reactant)
        print(index)
        for product_index in range(0,len(product)):
            try:
                if product[product_index] != []:
                    m = Chem.MolFromSmiles(product[product_index][0][0])
                    if m is not None:
                        product_smarts = Chem.MolToSmarts(m,isomericSmiles=True)
                        temp_product.append(product_smarts.format(R2))
                    if m is None:
                        pass
            except:
                IndexError
        part2 = ".".join(temp_product)
        #print(part2)
        return '>>'.join([part2,part1])

    if status == "training":
        REA, PRO = original_kegg_reactions(kegg_reactions,converting_dictionary=kegg_dictionary)
        extracting_reaction_core_smarts = []
        for index in range(0,len(REA)):
            r_part = ".".join([item[0] for item in REA[index] if item != []])
            p_part = ".".join([item[0] for item in PRO[index] if item != []])
            r_part = str(r_part).replace('*', "C")
            p_part = str(p_part).replace('*', "C")
            extracting_reaction_core_smarts.append(extracting(r_part,p_part))
    if status == "testing":
        extracting_reaction_core_smarts = []
        r_part = ".".join([item[0] for item in reactent if item != []])
        p_part = ".".join([item[0] for item in product if item != []])
        r_part = str(r_part).replace('*', "C")
        p_part = str(p_part).replace('*', "C")
        extracting_reaction_core_smarts.append(extracting(r_part,p_part))
    
    return extracting_reaction_core_smarts

def getting_reaction_cores():
    REA, PRO = original_kegg_reactions(kegg_reactions,converting_dictionary=kegg_dictionary)
    extracting_reaction_core_smarts = pushing_through_reaction_for_reaction_centers(REA,PRO)
    import os
    from chemical_module import Chemical
    from reaction import Reaction
    from rdkit.Chem import rdChemReactions
    list_of_cores = []
    for i in range(0,len(extracting_reaction_core_smarts)):
        rxn = rdChemReactions.ReactionFromSmarts(extracting_reaction_core_smarts[i])
        rxn_smiles = AllChem.ReactionToSmiles(rxn)
        rxn = Reaction(rxn_smiles, verbose=True)
        core = rxn.find_core()
        list_of_cores.append(core)
    return list_of_cores


def mol_level_descriptors(mol):
	'''
	Given an RDKit mol, returns a list of molecule-level descriptors 
	and their names
	returns: (labels, attributes)
	'''
	
	labels = [label for (label, f) in Descriptors._descList]
	attributes = [f(mol) for (label, f) in Descriptors._descList]
	
	return (labels, attributes)

def atom_level_descriptors(mol, include = ['functional'], asOneHot = False, ORIGINAL_VERSION = False):
	'''
	Given an RDKit mol, returns an N_atom-long list of lists,
	each of which contains atom-level descriptors and their names
	returns: (label, attributes)
	'''

	attributes = [[] for i in mol.GetAtoms()]
	labels = []
	if 'functional' in include:

		[attributes[i].append(x[0]) \
			for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
		labels.append('Crippen contribution to logp')

		[attributes[i].append(x[1]) \
			for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
		labels.append('Crippen contribution to mr')

		[attributes[i].append(x) \
			for (i, x) in enumerate(rdMolDescriptors._CalcTPSAContribs(mol))]
		labels.append('TPSA contribution')

		[attributes[i].append(x) \
			for (i, x) in enumerate(rdMolDescriptors._CalcLabuteASAContribs(mol)[0])]
		labels.append('Labute ASA contribution')

		[attributes[i].append(x) \
			for (i, x) in enumerate(EState.EStateIndices(mol))]
		labels.append('EState Index')

		rdPartialCharges.ComputeGasteigerCharges(mol)
		[attributes[i].append(float(a.GetProp('_GasteigerCharge'))) \
			for (i, a) in enumerate(mol.GetAtoms())]
		labels.append('Gasteiger partial charge')

		# Gasteiger partial charges sometimes gives NaN
		for i in range(len(attributes)):
			if np.isnan(attributes[i][-1]):
				attributes[i][-1] = 0.0

		[attributes[i].append(float(a.GetProp('_GasteigerHCharge'))) \
			for (i, a) in enumerate(mol.GetAtoms())]
		labels.append('Gasteiger hydrogen partial charge')

		# Gasteiger partial charges sometimes gives NaN
		for i in range(len(attributes)):
			if np.isnan(attributes[i][-1]):
				attributes[i][-1] = 0.0
	
	if 'structural' in include:
		[attributes[i].extend(atom_structural(mol.GetAtomWithIdx(i), asOneHot = asOneHot, ORIGINAL_VERSION = ORIGINAL_VERSION)) \
			for i in range(len(attributes))]
		labels.append('--many structural--')

	return (labels, attributes)

def pushing_through_reaction_for_reaction_centers(REA,PRO):
    def extracting(reactant,product):
        R2 = "*"
        temp_reactant, temp_product = [],[]
        for reactant_index in range(0,len(reactant)):
            try:
                if reactant[reactant_index] != []:
                    m = Chem.MolFromSmiles(reactant[reactant_index][0][0])
                    if m is not None:
                        reactant_smarts = Chem.MolToSmarts(m,isomericSmiles=True)
                        temp_reactant.append(reactant_smarts.format(R2))
                    else:
                        pass 
            except:
                IndexError
        part1 = ".".join(temp_reactant)
        print(index)
        for product_index in range(0,len(product)):
            try:
                if product[product_index] != []:
                    m = Chem.MolFromSmiles(product[product_index][0][0])
                    if m is not None:
                        product_smarts = Chem.MolToSmarts(m,isomericSmiles=True)
                        temp_product.append(product_smarts.format(R2))
                    if m is None:
                        pass
            except:
                IndexError
        part2 = ".".join(temp_product)
        #print(part2)
        return '>>'.join([part2,part1])
    extracting_reaction_core_smarts = []
    for index in range(0,len(REA)):
        extracting_reaction_core_smarts.append(extracting(REA[index],PRO[index]))
    return extracting_reaction_core_smarts

def training_dataset_features():
    def getting_reaction_cores():
        REA, PRO = original_kegg_reactions(kegg_reactions,converting_dictionary=kegg_dictionary)
        extracting_reaction_core_smarts = pushing_through_reaction_for_reaction_centers(REA,PRO)
        import os
        from chemical_module import Chemical
        from reaction import Reaction
        from rdkit.Chem import rdChemReactions
        list_of_cores = []
        for i in range(0,len(extracting_reaction_core_smarts)):
            rxn = rdChemReactions.ReactionFromSmarts(extracting_reaction_core_smarts[i])
            rxn_smiles = AllChem.ReactionToSmiles(rxn)
            rxn = Reaction(rxn_smiles, verbose=True)
            core = rxn.find_core()
            list_of_cores.append(core)
        return list_of_cores

    def generating_product_reactant_features():
        list_of_cores = getting_reaction_cores()
        product_array = []
        reactant_array = []
        for iteration in range(0,len(list_of_cores)):
            pro_nju = atom_level_descriptors(Chem.MolFromSmiles(list_of_cores[iteration].split(">>")[0]))
            rea_nju = atom_level_descriptors(Chem.MolFromSmiles(list_of_cores[iteration].split(">>")[1]))
            index = 0 
            product_addition = np.zeros((1,7))
            reactant_addition = np.zeros((1,7))
            for i in range(0,len(pro_nju[1])):
                product_addition = np.add(product_addition,pro_nju[1][i])
            product_array.append(product_addition)
            for j in range(0,len(rea_nju[1])):
                reactant_addition = np.add(reactant_addition,rea_nju[1][j])
            reactant_array.append(reactant_addition)
        return reactant_array,product_array

    reactent_features, product_features = generating_product_reactant_features() 
    reactant_array = [item[0] for item in reactent_features]
    product_array = [item[0] for item in product_features]
    reactant_features = pd.DataFrame(reactant_array)
    product_features = pd.DataFrame(product_array)
    reactant_columns = ['reactant_Crippen contribution to logp',
                        'reactant_Crippen contribution to mr',
                        'reactant_TPSA contribution',
                        'Labute ASA contribution',
                        'reactant_EState Index',
                        'reactant_Gasteiger partial charge',
                        'reactant_Gasteiger hydrogen partial charge']
    product_columns = ['product_Crippen contribution to logp', 
                        'product_Crippen contribution to mr', 
                        'product_TPSA contribution', 
                        'product_Labute ASA contribution', 
                        'product_EState Index', 
                        'product_Gasteiger partial charge', 
                        'product_Gasteiger hydrogen partial charge']
    reactant_features.columns = reactant_columns 
    product_features.columns = product_columns
    merged_features = pd.concat([reactant_features,product_features],axis =1)
    #merged_features.drop(columns = ["reactant_EState Index","product_EState Index"],inplace = True)
    return merged_features

"""
path = [morphine, morphinone,oripavine, thebaine, salutaridnol,salutaridine,
        r_reticuline,s_reticulene,s3_hydroxy_methylcoclaurine,s_n_methyl_coculaurine,
       s_coculaurine,norococulaurine, dopamine, tyramine, l_tyrosine]
"""

def generating_product_rxn_templates(rxn,second_product):
    """
    This function will take in a reaction rule and add a product to the existing product 
    Example 
    
    rxn rule : 
        '([#8&v2&H1:1]-[#8&v2&H1:2])>>([#8&v2&H0:1]=[#8&v2&H0:2])'
    
    smile form of rxn rule:
        '[OH:1][OH:2]>>[O:1]=[O:2]'
    
    second_product : morphine molecule "CN1CC[C@]23C4=C5C=CC(O)"
    
    generating_product_rxn_templates(rxn_rule , morphine)
    
    Output: 
    '[OH:1][OH:2].CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5>>[O:1]=[O:2]'
    
    """
    if type(second_product) != str:
        return "Type error, the second_product needs to be in smiles format"
    else:
        #Since the reaction is initially in smarts format it needs to be read into reaction format 
        rxn_smiles = rdChemReactions.ReactionFromSmarts(rxn)
        #The reaction format is not directly useful we will convert into smiles fromat 
        reaction_smiles = AllChem.ReactionToSmiles(rxn_smiles)
        #We want to split smiles strings into reaction/product becuase we wanted to add a new product to 
        # existing product 
        initial_product, reactant = reaction_smiles.split(">>")
        product = ".".join([initial_product,second_product])
        return ">>".join([product,reactant])

def setting_reaction_rules():
    data = pd.read_csv("retrorules_rr02_rp3_nohs/retrorules_rr02_flat_all.tsv",sep='\t')
    r_50 = reaction_rules_final(data,50)
    r_49 = reaction_rules_final(data,49)
    r_1 = reaction_rules_final(data,1)
    r_0 = reaction_rules_final(data,0)
    r_30 = reaction_rules_final(data,30)
    r_20 = reaction_rules_final(data,20)
    r_10 = reaction_rules_final(data,10)
    return r_50, r_49, r_30, r_20, r_10, r_1, r_0 


def one_generate_features(current_rule):
    rxn = rdChemReactions.ReactionFromSmarts(current_rule)
    rxn_smiles = AllChem.ReactionToSmiles(rxn)
    rxn = Reaction(rxn_smiles, verbose=True)
    core = rxn.find_core()
    return core


def model_development():
    def get_hydrogen_atoms(my_mol):
        my_mol = Chem.MolFromSmiles(moltosmiles(my_mol))
        my_mol_with_explicit_h = Chem.AddHs(my_mol)
        return my_mol_with_explicit_h.GetNumAtoms() - my_mol_with_explicit_h.GetNumHeavyAtoms()

    all_values = []
    current_max = []
    for i in range(0,len(REA)):
        for item in REA[i]:
            if item != []:
                temp = get_hydrogen_atoms(item[0])
            current_max.append(temp)
        all_values.append(max(current_max))
        
    merged_features = training_dataset_features()
    X = pd.concat([reactent_fingerprints,merged_features],axis =1)
    Y = pd.DataFrame(all_values)
    from keras.models import Sequential
    from keras.layers import Dense
    model = Sequential()
    model.add(Dense(1, input_dim=2060, activation='sigmoid'))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    hist = model.fit(X.values,Y.values, epochs=2, batch_size=20,verbose=1, 
                    validation_data=(X.values,Y.values), shuffle=True)

def generating_atomic_features(extracting_reaction_core_smarts):
    import os
    from chemical_module import Chemical
    from reaction import Reaction
    from rdkit.Chem import rdChemReactions
    list_of_cores = []
    rxn = rdChemReactions.ReactionFromSmarts(extracting_reaction_core_smarts)
    rxn_smiles = AllChem.ReactionToSmiles(rxn)
    rxn = Reaction(rxn_smiles, verbose=True)
    core = rxn.find_core()
    list_of_cores.append(core)
    return list_of_cores 

def testing_atomic_features(smile_string_list):
    temp = []
    current_sum = np.zeros((1,7))
    for i in range(0,len(smile_string_list)):
        for j in range(0,len(smile_string_list[i])):
            atom_level = pd.DataFrame(np.array(atom_level_descriptors(smile_string_list[i]),asOneHot=True)[1])
            sum_values =  atom_level.sum()
            current_sum = np.sum(current_sum,sum_values)
    return current_sum

def getting_reaction_cores(reaction):
    extracting_reaction_core_smarts = pushing_through_reaction_for_reaction_centers(reaction.split(">>")[1],reaction.split(">>")[0])
    import os
    from chemical_module import Chemical
    from reaction import Reaction
    from rdkit.Chem import rdChemReactions
    list_of_cores = []
    for i in range(0,len(extracting_reaction_core_smarts)):
        rxn = rdChemReactions.ReactionFromSmarts(extracting_reaction_core_smarts[i])
        rxn_smiles = AllChem.ReactionToSmiles(rxn)
        rxn = Reaction(rxn_smiles, verbose=True)
        core = rxn.find_core()
        list_of_cores.append(core)
    return list_of_cores
def one_generate_features(current_rule):
    rxn = rdChemReactions.ReactionFromSmarts(current_rule)
    rxn_smiles = AllChem.ReactionToSmiles(rxn)
    rxn = Reaction(rxn_smiles, verbose=True)
    core = rxn.find_core()
    return core

def testing_atomic_features(cores):
    product_addition = np.zeros((1,7))
    reactant_addition = np.zeros((1,7))
    pro_nju = atom_level_descriptors(Chem.MolFromSmiles(cores.split(">>")[0]))
    rea_nju = atom_level_descriptors(Chem.MolFromSmiles(cores.split(">>")[1]))
    for i in range(0,len(pro_nju[1])):
        product_addition = np.add(product_addition,pro_nju[1][i])
    for j in range(0,len(rea_nju[1])):
        reactant_addition = np.add(reactant_addition,rea_nju[1][j])
    P2_removing_e_state_index = np.concatenate([product_addition[0][0:3],product_addition[0][4:]])
    P1_removing_e_state_index = np.concatenate([reactant_addition[0][0:3],reactant_addition[0][4:]])
    return np.concatenate([P1_removing_e_state_index,P2_removing_e_state_index])

def mol_level_descriptors(mol):
    '''
    Given an RDKit mol, returns a list of molecule-level descriptors 
    and their names
    returns: (labels, attributes)
    '''
    
    labels = [label for (label, f) in Descriptors._descList]
    attributes = [f(mol) for (label, f) in Descriptors._descList]
    
    return (labels, attributes)

def atom_level_descriptors(mol, include = ['functional'], asOneHot = False, ORIGINAL_VERSION = False):
    '''
    Given an RDKit mol, returns an N_atom-long list of lists,
    each of which contains atom-level descriptors and their names
    returns: (label, attributes)
    '''

    attributes = [[] for i in mol.GetAtoms()]
    labels = []
    if 'functional' in include:

        [attributes[i].append(x[0]) \
            for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
        labels.append('Crippen contribution to logp')

        [attributes[i].append(x[1]) \
            for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
        labels.append('Crippen contribution to mr')

        [attributes[i].append(x) \
            for (i, x) in enumerate(rdMolDescriptors._CalcTPSAContribs(mol))]
        labels.append('TPSA contribution')

        [attributes[i].append(x) \
            for (i, x) in enumerate(rdMolDescriptors._CalcLabuteASAContribs(mol)[0])]
        labels.append('Labute ASA contribution')

        [attributes[i].append(x) \
            for (i, x) in enumerate(EState.EStateIndices(mol))]
        labels.append('EState Index')

        rdPartialCharges.ComputeGasteigerCharges(mol)
        [attributes[i].append(float(a.GetProp('_GasteigerCharge'))) \
            for (i, a) in enumerate(mol.GetAtoms())]
        labels.append('Gasteiger partial charge')

        # Gasteiger partial charges sometimes gives NaN
        for i in range(len(attributes)):
            if np.isnan(attributes[i][-1]):
                attributes[i][-1] = 0.0

        [attributes[i].append(float(a.GetProp('_GasteigerHCharge'))) \
            for (i, a) in enumerate(mol.GetAtoms())]
        labels.append('Gasteiger hydrogen partial charge')

        # Gasteiger partial charges sometimes gives NaN
        for i in range(len(attributes)):
            if np.isnan(attributes[i][-1]):
                attributes[i][-1] = 0.0
    
    if 'structural' in include:
        [attributes[i].extend(atom_structural(mol.GetAtomWithIdx(i), asOneHot = asOneHot, ORIGINAL_VERSION = ORIGINAL_VERSION)) \
            for i in range(len(attributes))]
        labels.append('--many structural--')

    return (labels, attributes)

def circular_fingerprint(smile_string):
    mol = Chem.MolFromSmiles(smile_string)
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, bitInfo=bi)
    return fp.ToBitString()

def encoding(current_element):
    return [int(item) for item in circular_fingerprint(current_element)]